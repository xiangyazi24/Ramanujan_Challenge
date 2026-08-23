#!/usr/bin/env python3
"""Exact symbolic and finite-field audit for Q2362.

This uses the definitions in the prompt literally:
  A_m(r) = N_{m+1}(r),
  B_m(r) = R_m(r),
  W_m = A_m B_{m-1} - A_{m-1} B_m.

It also prints D_m=A_m-B_m because the Q2355 factorization N_3-R_2
belongs to D_2, not to W_1 under the displayed definition.
"""
from __future__ import annotations

from math import gcd
import sympy as sp

r = sp.symbols("r")


def P(t):
    return sp.expand(34*t**3 + 51*t**2 + 27*t + 5)


def Q(m):
    return sp.expand((2*r + 2*m + 1) * (
        3*r**2 + 2*r*m + m**2 + 3*r + m + 1
    ))


def primitive_factor_string(poly):
    poly = sp.Poly(sp.expand(poly), r, domain=sp.ZZ)
    content, facs = sp.factor_list(poly.as_expr(), r)
    pieces = []
    if content != 1:
        pieces.append(str(content))
    for f, e in facs:
        fs = str(sp.expand(f))
        pieces.append(f"({fs})" if e == 1 else f"({fs})^{e}")
    return content, facs, " * ".join(pieces) if pieces else "1"


def factor_integer(n):
    n = int(n)
    sign = -1 if n < 0 else 1
    fac = sp.factorint(abs(n))
    body = " * ".join(f"{p}^{e}" if e != 1 else str(p) for p, e in sorted(fac.items()))
    if not body:
        body = "1"
    return ("-" if sign < 0 else "") + body


# Apéry gap continuants N_h.
N = {0: sp.Integer(0), 1: sp.Integer(1)}
for h in range(1, 8):
    N[h+1] = sp.expand(P(r+h) * N[h] - (r+h)**6 * N[h-1])

# Racah gap continuants R_m.
R = {-1: sp.Integer(0), 0: sp.Integer(1)}
for m in range(1, 7):
    R[m] = sp.expand(Q(m) * R[m-1] - (r+m)**6 * R[m-2])

A = {m: N[m+1] for m in range(0, 7)}
B = {m: R[m] for m in range(0, 7)}
W = {}
D = {}

print("Q2362_SYMBOLIC_BEGIN")
print("INDEXING_CHECK")
print("A0=", A[0])
print("B0=", B[0])
print("W1_expected=P(r+1)-Q1")
print("Q2355_difference_is_D2=A2-B2=N3-R2")

for m in range(1, 6):
    W[m] = sp.expand(A[m]*B[m-1] - A[m-1]*B[m])
    D[m] = sp.expand(A[m] - B[m])
    print(f"\nM={m} H={m+1}")
    print(f"DEG_A={sp.degree(A[m], r)} DEG_B={sp.degree(B[m], r)}")
    print(f"DEG_W={sp.degree(W[m], r)}")
    c, facs, fs = primitive_factor_string(W[m])
    print(f"CONTENT_W={c}")
    print(f"W_FACTOR={fs}")
    print("W_FACTOR_LIST=" + repr([(str(sp.expand(f)), int(e), int(sp.degree(f, r))) for f,e in facs]))
    print(f"DEG_D={sp.degree(D[m], r)}")
    cd, facsd, fsd = primitive_factor_string(D[m])
    print(f"CONTENT_D={cd}")
    print(f"D_FACTOR={fsd}")
    print("D_FACTOR_LIST=" + repr([(str(sp.expand(f)), int(e), int(sp.degree(f, r))) for f,e in facsd]))
    print(f"GCD_A_W={sp.factor(sp.gcd(sp.Poly(A[m],r,domain=sp.ZZ),sp.Poly(W[m],r,domain=sp.ZZ)).as_expr())}")
    print(f"GCD_A_D={sp.factor(sp.gcd(sp.Poly(A[m],r,domain=sp.ZZ),sp.Poly(D[m],r,domain=sp.ZZ)).as_expr())}")

# Verify the exact mixed-Wronskian recurrence.
print("\nWRONSKIAN_RECURRENCE_CHECKS")
for m in range(2, 6):
    delta = sp.expand(P(r+m) - Q(m))
    rhs = sp.expand(delta*A[m-1]*B[m-1] + (r+m)**6*W[m-1])
    print(f"m={m} delta_factor={sp.factor(delta)} ok={sp.expand(W[m]-rhs)==0}")

# Resultants against every irreducible factor of W_m and D_m.  We print
# exact values and complete factorizations only when factorint terminates;
# m<=5 is small enough for the reduced factors found here.
print("\nREDUCED_FACTOR_RESULTANTS")
for family_name, family in (("W", W), ("D", D)):
    for m in range(1, 6):
        Am = sp.Poly(A[m], r, domain=sp.ZZ)
        _, facs = sp.factor_list(sp.Poly(family[m], r, domain=sp.ZZ).as_expr(), r)
        for idx, (f,e) in enumerate(facs, start=1):
            fp = sp.Poly(f, r, domain=sp.ZZ)
            res = sp.resultant(Am.as_expr(), fp.as_expr(), r)
            print(f"{family_name}{m}_FACTOR{idx}_DEG={fp.degree()} EXP={e}")
            print(f"{family_name}{m}_FACTOR{idx}={sp.expand(f)}")
            print(f"RES_A{m}_{family_name}{m}F{idx}={res}")
            # Complete factorization for manageable values; factorint itself
            # is exact and will fail loudly rather than label a composite prime.
            print(f"RESFAC_A{m}_{family_name}{m}F{idx}={factor_integer(res)}")


def primes_upto(n):
    return list(sp.primerange(2, n+1))


def apery_row_mod_p(p):
    b = [0]*p
    b[0] = 1
    if p > 1:
        b[1] = 5 % p
    for n in range(1, p-1):
        den = pow(n+1, 3, p)
        assert den != 0
        pn = (34*n**3 + 51*n**2 + 27*n + 5) % p
        b[n+1] = ((pn*b[n] - pow(n,3,p)*b[n-1]) * pow(den, -1, p)) % p
    return b


def Pmod(t,p):
    return (34*t**3 + 51*t**2 + 27*t + 5) % p


def Qmod(r0,m,p):
    return ((2*r0+2*m+1) * (3*r0*r0 + 2*r0*m + m*m + 3*r0 + m + 1)) % p


def gap_data_mod(r0,h,p):
    # N_0,N_1 through N_h.
    nprev, ncur = 0, 1
    Ns = [nprev,ncur]
    for j in range(1,h):
        nnext = (Pmod(r0+j,p)*ncur - pow(r0+j,6,p)*nprev) % p
        Ns.append(nnext)
        nprev,ncur = ncur,nnext
    # R_-1,R_0 through R_{h-1}.
    rm1, rcur = 0,1
    Rs = [rcur]
    for j in range(1,h):
        rnext = (Qmod(r0,j,p)*rcur - pow(r0+j,6,p)*rm1) % p
        Rs.append(rnext)
        rm1,rcur = rcur,rnext
    m=h-1
    A_m = Ns[h]
    A_prev = Ns[h-1]
    B_m = Rs[m]
    B_prev = Rs[m-1] if m>=1 else 0
    Wm = (A_m*B_prev - A_prev*B_m) % p if m>=1 else None
    return A_m,A_prev,B_m,B_prev,Wm


def Kmod(r0,s0,p):
    total=0
    for k in range(0,min(r0,s0)+1):
        # D(n,k)=binom(n,k)binom(n+k,k), all integer then reduce.
        total += sp.binomial(r0,k)*sp.binomial(r0+k,k)*sp.binomial(s0,k)*sp.binomial(s0+k,k)
    return int(total % p)

print("\nPHYSICAL_SCAN_P_LE_500")
all_pairs=[]
zero_w=[]
zero_w_nonmate=[]
short_pairs=[]
for p in primes_upto(500):
    if p < 5:
        continue
    b=apery_row_mod_p(p)
    zeros=[i for i in range(1,p-1) if b[i]==0]
    for ii,r0 in enumerate(zeros):
        for s0 in zeros[ii+1:]:
            h=s0-r0
            if h<2:
                continue
            Am,Aprev,Bm,Bprev,Wm=gap_data_mod(r0,h,p)
            assert Am==0
            kval=Kmod(r0,s0,p)
            # Under b_r=0, B_m=0 iff K=0 in the nonwrapping range;
            # verify directly in every scanned case.
            assert (Bm==0)==(kval==0)
            mate=(r0+s0+1==p)
            rec=(p,r0,s0,h,Wm,Bm,kval,mate)
            all_pairs.append(rec)
            if h<=6:
                short_pairs.append(rec)
            if Wm==0:
                zero_w.append(rec)
                if not mate:
                    zero_w_nonmate.append(rec)

print(f"TOTAL_PHYSICAL_PAIRS={len(all_pairs)}")
print(f"TOTAL_SHORT_H_LE_6={len(short_pairs)}")
print(f"W_ZERO_COUNT={len(zero_w)}")
print(f"W_ZERO_NONMATE_COUNT={len(zero_w_nonmate)}")
print("W_ZERO_RECORDS="+repr(zero_w))
print("W_ZERO_NONMATE_RECORDS="+repr(zero_w_nonmate))
print("SHORT_PAIR_RECORDS="+repr(short_pairs))
print("Q2362_SYMBOLIC_END")
