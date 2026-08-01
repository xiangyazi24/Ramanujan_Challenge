#!/usr/bin/env python3
"""Verify: (a) collision pairs (r,r+d) with b_r=b_{r+d}=c!=0 force y_r = c*sigma_d(r) mod p;
(b) sigma_d pairwise distinct as rational functions (symbolic, d,d' <= 10)."""
import sympy as sp

# --- (b) symbolic over Q ---
r = sp.symbols('r')
a = sp.expand((2*r+1)*(17*r**2+17*r+5))/ (r+1)**3
beta = -r**3/(r+1)**3
# A_d, B_d: b_{r+d} = A_d b_r + B_d b_{r-1}; recursion in d with shifted coefficients:
# A_{d+1}(r) = a(r+d) A_d(r) + beta(r+d) A_{d-1}(r), likewise B.
A = {0: sp.Integer(1), 1: a}
B = {0: sp.Integer(0), 1: beta}
DMAX = 10
for d in range(1, DMAX):
    ashift = a.subs(r, r+d); bshift = beta.subs(r, r+d)
    A[d+1] = sp.cancel(ashift*A[d] + bshift*A[d-1])
    B[d+1] = sp.cancel(ashift*B[d] + bshift*B[d-1])
print("symbolic sigma_d distinctness (Psi = (1-A_d)B_d' - (1-A_d')B_d):")
bad = []
for d in range(1, DMAX+1):
    for dp in range(d+1, DMAX+1):
        Psi = sp.cancel((1-A[d])*B[dp] - (1-A[dp])*B[d])
        if Psi == 0:
            bad.append((d,dp))
print("  degenerate pairs (d,d') with Psi==0:", bad if bad else "NONE (all distinct, d,d'<=10)")

# --- (a) numeric mod p=101 ---
p = 101
b = [1,5]
for n in range(1,p+2):
    b.append(((2*n+1)*(17*n*n+17*n+5)*b[n] - n**3*b[n-1]) * pow((n+1)**3, -1, 10**60) % 10**60 if False else 0)
# redo exactly with integers
b = [1,5]
for n in range(1, p+2):
    num = (2*n+1)*(17*n*n+17*n+5)*b[n] - n**3*b[n-1]
    q, rem = divmod(num, (n+1)**3); assert rem == 0
    b.append(q)
bm = [x % p for x in b[:p]]
# transfer A_d,B_d mod p evaluated at base point rr: iterate
def AB(d, rr, p):
    Ad_1, Ad = 1 % p, 0  # A_0=1; track pair (A_{k}, A_{k-1})
    # use vectors: (A_k, B_k) with recursion A_{k+1} = a(rr+k)A_k + beta(rr+k)A_{k-1}
    Aprev, Acur = 1, None; Bprev, Bcur = 0, None
    # k=1: a(rr), beta(rr)
    def coefs(x):
        den = pow((x+1)**3 % p, -1, p)
        return ((2*x+1)*(17*x*x+17*x+5) % p)*den % p, (-(x**3) % p)*den % p
    a1, b1 = coefs(rr % p)
    Acur, Bcur = a1, b1
    for k in range(1, d):
        ak, bk = coefs((rr+k) % p)
        Acur, Aprev = (ak*Acur + bk*Aprev) % p, Acur
        Bcur, Bprev = (ak*Bcur + bk*Bprev) % p, Bcur
    return Acur, Bcur
ok, tested, skipped = 0, 0, 0
for rr in range(1, p-1):
    for d in range(1, p-1-rr):
        if bm[rr] != 0 and bm[rr] == bm[rr+d]:
            c = bm[rr]
            Ad, Bd = AB(d, rr, p)
            # sanity: b_{r+d} = A_d b_r + B_d b_{r-1}
            assert (Ad*bm[rr] + Bd*bm[rr-1]) % p == bm[rr+d] % p, (rr,d)
            if Bd % p == 0:
                skipped += 1; continue
            y_pred = c * (1 - Ad) % p * pow(Bd, -1, p) % p
            tested += 1
            if y_pred == bm[rr-1] % p: ok += 1
print(f"numeric p={p}: collision pairs tested {tested}, y=c*sigma_d(r) verified {ok}, B_d(r)=0 skipped {skipped}")
