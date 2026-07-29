#!/usr/bin/env python3
"""
P2.5 deep holonomy test, stage 2.

Verified: B^{-1}(k,N) = (-1)^{k+N}(2N+1) C(k,N)C(k+N,N)
                        / (2^k C(2k,k)(k+N+1) C(2k-N,k-N)^2)
is proper hypergeometric, so f(k) = sum_N B^{-1}(k,N) Qhat_N is holonomic
by Zeilberger theory. Boxes (16,16), (20,12) are ruled out, so the minimal
recurrence is big. Push data to k=1200 (mod p) and test big boxes:
(6,160), (10,100), (14,70), (24,40).
"""
from fractions import Fraction as F
from math import comb
import sys, time

KMAX = 1200
P0 = (1 << 61) - 1

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

t0 = time.time()

# --- verify the B^{-1} closed form exactly up to k=20 ---
def B(N, k):
    if k < 0 or k > N:
        return 0
    return 2**k * comb(2*k, k) * comb(N, k) * comb(N+k, k)

def Binv_closed(k, N):
    return F((-1)**(k+N) * (2*N+1) * comb(k, N),
             2**k * comb(2*k, k) * (k+N+1) * comb(k+N, N))

KT = 20
ok = True
for k in range(KT+1):
    for N in range(k+1):
        s = sum(Binv_closed(k, m) * B(m, N) for m in range(N, k+1))
        want = F(1) if N == k else F(0)
        if s != want:
            ok = False
            print(f"  closed form FAILS at k={k}, N={N}")
print(f"B^-1 closed form verified as two-sided inverse up to k={KT}: "
      f"{'YES' if ok else 'NO'}  ({time.time()-t0:.1f}s)")
sys.stdout.flush()

# --- trajectory mod p ---
print(f"computing Qhat^e1 mod p for N=0..{KMAX} ...")
sys.stdout.flush()
t1 = time.time()

def inv(x):
    return pow(x % P0, P0 - 2, P0)

r1 = [1, 0, 0]
Q_mod = [1]
for N in range(KMAX):
    M = M_entries(N)
    dinv = inv(delta_H(N))
    MH = [[M[i][j] * dinv % P0 for j in range(3)] for i in range(3)]
    r1 = [sum(r1[i]*MH[i][k] for i in range(3)) % P0 for k in range(3)]
    Q_mod.append(r1[0])
print(f"  trajectory mod p done ({time.time()-t1:.1f}s)")
sys.stdout.flush()

# --- triangular inversion mod p, using iterative row update for B(K,k) ---
t1 = time.time()
f_mod = []
Brow = [1]  # B(K,k) mod p for current K, k=0..K
for K in range(KMAX + 1):
    if K > 0:
        # update B(K,k) from B(K-1,k): C(N,k)C(N+k,k) with N: K-1->K
        # B(K,k)/B(K-1,k) = [C(K,k)/C(K-1,k)] * [C(K+k,k)/C(K-1+k,k)]
        #                 = [K/(K-k)] * [(K+k)/K] = (K+k)/(K-k)
        Brow = [Brow[k] * ((K + k) % P0) % P0 * inv(K - k) % P0
                for k in range(K)]
        # new diagonal entry B(K,K) = 2^K C(2K,K)^2
        Brow.append(2**K % P0 * pow(comb(2*K, K) % P0, 2, P0) % P0)
    rhs = Q_mod[K]
    for k in range(K):
        rhs = (rhs - f_mod[k] * Brow[k]) % P0
    f_mod.append(rhs * inv(Brow[K]) % P0)
    if K % 200 == 0:
        print(f"  inversion K={K}  ({time.time()-t1:.1f}s)")
        sys.stdout.flush()
# sanity: f(1)=5749/3136, f(2)=16811771/4572288
assert f_mod[1] == 5749 * inv(3136) % P0
assert f_mod[2] == 16811771 * inv(4572288) % P0
print(f"  mod-p inversion done, sanity OK ({time.time()-t1:.1f}s)")
sys.stdout.flush()

# --- big-box nullity tests ---
def nullity_modp(A, ncols):
    m = len(A)
    prow = 0
    for col in range(ncols):
        piv = next((r for r in range(prow, m) if A[r][col] != 0), None)
        if piv is None:
            continue
        A[prow], A[piv] = A[piv], A[prow]
        iv = pow(A[prow][col], P0 - 2, P0)
        A[prow] = [x * iv % P0 for x in A[prow]]
        Ap = A[prow]
        for r in range(m):
            if r != prow and A[r][col] != 0:
                fac = A[r][col]
                Ar = A[r]
                A[r] = [(Ar[j] - fac * Ap[j]) % P0 for j in range(ncols)]
        prow += 1
        if prow == m:
            break
    return ncols - prow

def one_shot(r, d, kmin=5, extra=10):
    nunk = (r + 1) * (d + 1)
    kmax_row = KMAX - r
    if kmax_row - kmin + 1 < nunk + extra:
        print(f"  box ({r},{d}): not enough data (need {nunk+extra} rows)")
        return None
    t = time.time()
    A = []
    for k in range(kmin, kmin + nunk + extra):
        row = []
        for j in range(r + 1):
            fk = f_mod[k + j]
            kp = 1
            for _ in range(d + 1):
                row.append(fk * kp % P0)
                kp *= k
                kp %= P0
        A.append(row)
    nl = nullity_modp(A, nunk)
    print(f"  box r<={r}, d<={d}: nullity = {nl}  "
          f"({'RECURRENCE EXISTS' if nl else 'ruled out'}, {time.time()-t:.1f}s)")
    sys.stdout.flush()
    return nl

print("\nBig-box tests on f(k) (any r0<=r, d0<=d recurrence embeds):")
for (r, d) in [(6, 160), (10, 100), (14, 70), (24, 40)]:
    one_shot(r, d)

print(f"\nTotal {time.time()-t0:.1f}s")
