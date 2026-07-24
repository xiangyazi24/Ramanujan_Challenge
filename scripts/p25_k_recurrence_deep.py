#!/usr/bin/env python3
"""
P2.5 deep holonomy test for the Delannoy-basis coefficients f(k).

Part 1: Is the inverse triangle B^{-1} proper hypergeometric?
        B(N,k) = 2^k C(2k,k) C(N,k) C(N+k,k);  f = B^{-1} Qhat.
        If B^{-1}(k,N) is hypergeometric in both k and N, then f(k) is a
        definite sum (proper hypergeometric)*(holonomic) => holonomic by
        Zeilberger/WZ theory, and the empty searches just mean the minimal
        recurrence exceeds (r,d)=(8,24).

Part 2: One-shot maximal guessing at (r,d) = (12,22), (16,16), (8,34)
        mod p with data to k=KMAX. A recurrence of any order r0<=r and
        degree d0<=d embeds (zero-padded) in the null space, so nullity=0
        rules out the whole box.
"""
from fractions import Fraction as F
from math import comb
import sys, time

KMAX = 320

def M_entries(n):
    m11 = (-2*n-5)*(n+3)**2*(136*n**4+1424*n**3+5548*n**2+9551*n+6141)
    m12 = 384*n**6+6384*n**5+44168*n**4+162698*n**3+336377*n**2+369933*n+169011
    m13 = -(480*n**4+4980*n**3+19210*n**2+32690*n+20730)
    m21 = (n+2)**2*(n+3)**2*(4*n+10)*(48*n**3+386*n**2+1017*n+879)
    m22 = (n+2)**2*(-272*n**5-3848*n**4-21732*n**3-61184*n**2-85761*n-47808)
    m23 = (n+2)**2*(320*n**3+2540*n**2+6610*n+5640)
    m31 = (-4*n-10)*(n+2)**2*(n+3)**2*(32*n**4+302*n**3+1037*n**2+1530*n+813)
    m32 = (n+2)**2*(192*n**6+2984*n**5+19116*n**4+64452*n**3+120256*n**2+117279*n+46476)
    m33 = (n+2)**2*(-16*n**5-408*n**4-2912*n**3-8884*n**2-12254*n-6240)
    return [[m11,m12,m13],[m21,m22,m23],[m31,m32,m33]]

def delta_H(n):
    return -2*(n+2)**2*(n+3)**2*(2*n+5)*(2*n+7)**2

def B(N, k):
    if k < 0 or k > N:
        return 0
    return 2**k * comb(2*k, k) * comb(N, k) * comb(N+k, k)

# ----------------------------------------------------------------------
# Part 1: inverse triangle structure
# ----------------------------------------------------------------------
print("Part 1: inverse triangle B^{-1}")
KT = 12
# invert lower-triangular integer matrix over Q
Binv = [[F(0)]*(KT+1) for _ in range(KT+1)]
for k in range(KT+1):
    Binv[k][k] = F(1, B(k, k))
    for N in range(k-1, -1, -1):
        # sum_{m=N..k} Binv[k][m] B(m,N) = 0 for N<k
        s = sum(Binv[k][m] * B(m, N) for m in range(N+1, k+1))
        Binv[k][N] = -s / B(N, N)

# candidate closed form: Binv(k,N) = (-1)^{k+N} (2N+1) C(k+N,k-N) C(2N,N) ... ?
# test hypergeometricity: ratio in N should be rational; try to match
# the Legendre inverse: for T(N,k)=C(N,k)C(N+k,k),
#   T^{-1}(k,N) = (-1)^{k+N} (2N+1) * (k+N)!k!(k-N)!^{-1} ... check numerically.
print("  Binv(k,N) for k<=6 (exact):")
for k in range(7):
    print("   ", [str(Binv[k][N]) for N in range(k+1)])

# guess: Binv(k,N) = (-1)^{k+N} (2N+1) * C(k,N) C(k+N,N) / (2^k C(2k,k) * (k+N+1))?
def guess1(k, N):
    return F((-1)**(k+N) * (2*N+1) * comb(k, N) * comb(k+N, N),
             2**k * comb(2*k, k) * (k+N+1))
ok1 = all(Binv[k][N] == guess1(k, N) for k in range(KT+1) for N in range(k+1))
print(f"  guess (-1)^(k+N)(2N+1)C(k,N)C(k+N,N) / (2^k C(2k,k)(k+N+1)): "
      f"{'MATCH' if ok1 else 'no'}")
if not ok1:
    # print ratios to see structure
    print("  ratio Binv(k,N)/guess1(k,N):")
    for k in range(6):
        print("   ", [str(Binv[k][N]/guess1(k,N)) for N in range(k+1)])
sys.stdout.flush()

# ----------------------------------------------------------------------
# Part 2: data to KMAX and one-shot maximal guesses
# ----------------------------------------------------------------------
t0 = time.time()
print(f"\nPart 2: computing f(k) for k=0..{KMAX} ...")
sys.stdout.flush()
r1 = [F(1), F(0), F(0)]
Qe1 = [F(1)]
for N in range(KMAX):
    M = M_entries(N)
    d = F(delta_H(N))
    MH = [[F(M[i][j]) / d for j in range(3)] for i in range(3)]
    r1 = [sum(r1[i]*MH[i][k] for i in range(3)) for k in range(3)]
    Qe1.append(r1[0])
print(f"  trajectory done ({time.time()-t0:.1f}s)")
sys.stdout.flush()

P0 = (1 << 61) - 1

def frac_modp(x, p):
    return (x.numerator % p) * pow(x.denominator % p, p - 2, p) % p

# triangular inversion directly mod p (much faster than exact for large KMAX)
Q_mod = [frac_modp(x, P0) for x in Qe1]
Bkk_inv = [pow(B(k, k) % P0, P0-2, P0) for k in range(KMAX+1)]
f_mod = []
t1 = time.time()
for K in range(KMAX+1):
    rhs = Q_mod[K]
    for k in range(K):
        rhs = (rhs - f_mod[k] * (B(K, k) % P0)) % P0
    f_mod.append(rhs * Bkk_inv[K] % P0)
print(f"  mod-p inversion done ({time.time()-t1:.1f}s)")
sys.stdout.flush()

def nullity_modp(A, ncols, p):
    m = len(A)
    prow = 0
    for col in range(ncols):
        piv = next((r for r in range(prow, m) if A[r][col] % p != 0), None)
        if piv is None:
            continue
        A[prow], A[piv] = A[piv], A[prow]
        inv = pow(A[prow][col], p - 2, p)
        A[prow] = [(x * inv) % p for x in A[prow]]
        for r in range(m):
            if r != prow and A[r][col] % p != 0:
                fac = A[r][col]
                Ar = A[r]; Ap = A[prow]
                A[r] = [(Ar[j] - fac * Ap[j]) % p for j in range(ncols)]
        prow += 1
        if prow == m:
            break
    return ncols - prow

def one_shot(r, d, kmin=5, extra=8):
    nunk = (r + 1) * (d + 1)
    kmax_row = KMAX - r
    nrows = min(kmax_row - kmin + 1, nunk + extra)
    if nrows < nunk + extra:
        print(f"  ({r},{d}): not enough data")
        return
    t = time.time()
    A = []
    for k in range(kmin, kmin + nrows):
        row = []
        for j in range(r + 1):
            fk = f_mod[k + j]
            kp = 1
            for _ in range(d + 1):
                row.append(fk * kp % P0)
                kp *= k
        A.append(row)
    nl = nullity_modp(A, nunk, P0)
    print(f"  box r<={r}, d<={d}: nullity = {nl}  "
          f"({'RECURRENCE EXISTS' if nl else 'ruled out'}, {time.time()-t:.1f}s)")
    sys.stdout.flush()
    return nl

print("\n  One-shot maximal boxes (any recurrence with r0<=r, d0<=d embeds):")
one_shot(12, 22)
one_shot(16, 16)
one_shot(8, 34)
one_shot(20, 12)
print(f"\nTotal {time.time()-t0:.1f}s")
