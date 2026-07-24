#!/usr/bin/env python3
"""
P2.5: Search for hypergeometric summand of Q̂_N.

If Q̂_N = Σ_k a(N,k), we need a(N,k) hypergeometric in N,k.

Strategy:
1. Compute exact rational Q̂_N for N=0..30
2. Test ansatz: Q̂_N = Σ_{k=0}^N f(k) × 2^k C(2k,k)C(N,k)C(N+k,k)
   where f(k) is a polynomial or simple rational function.
3. If that fails, try more complex ansätze.

The Delannoy model: D_N² = Σ_k 2^k C(2k,k)C(N,k)C(N+k,k)
"""
from fractions import Fraction as F
from math import comb

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

NMAX = 30

# Compute all three columns of the normalized CMF trajectory from e₁
# Q̂_N^{(j)} = e_j · Π M_H · e₁ for j=1,2,3
print("Computing Q̂_N values (exact fractions)...")
rows = [[F(1), F(0), F(0)],   # e₁
        [F(0), F(1), F(0)],   # e₂
        [F(0), F(0), F(1)]]   # e₃

qhat = {j: [rows[j][0]] for j in range(3)}

for N in range(NMAX):
    M = M_entries(N)
    d = F(delta_H(N))
    MH = [[F(M[i][j]) / d for j in range(3)] for i in range(3)]

    for j in range(3):
        r = rows[j]
        new_r = [sum(r[i]*MH[i][k] for i in range(3)) for k in range(3)]
        rows[j] = new_r
        qhat[j].append(new_r[0])

# Also compute from q-row and p-row
q = [F(33750), F(-36000), F(9000)]
p = [F(30921), F(-32972), F(8240)]

Qhat_q = []  # Q̂_N from q-row
Phat_p = []  # P̂_N from p-row

for N in range(NMAX + 1):
    Qhat_q.append(sum(q[j] * qhat[j][N] for j in range(3)))
    Phat_p.append(sum(p[j] * qhat[j][N] for j in range(3)))

print(f"Q̂_0 (e₁) = {qhat[0][0]}")
print(f"Q̂_1 (e₁) = {qhat[0][1]}")
print(f"Q̂_2 (e₁) = {qhat[0][2]}")
print(f"Q̂_3 (e₁) = {qhat[0][3]}")
print()

# Now test: Q̂_N^{e₁} = Σ_k f(k) · base(N,k)?
# where base(N,k) = 2^k C(2k,k) C(N,k) C(N+k,k)

def delannoy_summand(N, k):
    """F(N,k) = 2^k C(2k,k) C(N,k) C(N+k,k)"""
    if k < 0 or k > N:
        return F(0)
    return F(2**k * comb(2*k, k) * comb(N, k) * comb(N+k, k))

# D_N² = Σ_k base(N,k)
for N in range(10):
    dn2 = sum(delannoy_summand(N, k) for k in range(N+1))
    print(f"D_{N}² = {dn2}", end="")
    # Also print Q̂_N^{e₁}
    print(f"   Q̂_{N} = {qhat[0][N]}")

# If Q̂_N = Σ_k f(k) · base(N,k), then for each N:
# Q̂_N - Σ_{k=0}^N f(k) base(N,k) = 0
# This is a system for f(0), f(1), ..., f(N)
# Using N=0: f(0) · base(0,0) = Q̂_0 → f(0) = Q̂_0 / base(0,0) = 1/1 = 1
# Using N=1: f(0)·base(1,0) + f(1)·base(1,1) = Q̂_1
# etc.

print("\n=== Solving for f(k) with fixed polynomial guess ===")
f_vals = []
for K in range(NMAX + 1):
    # Use N = K equation
    rhs = qhat[0][K]
    for k in range(K):
        rhs -= f_vals[k] * delannoy_summand(K, k)
    b_KK = delannoy_summand(K, K)
    if b_KK == 0:
        print(f"  base({K},{K}) = 0, can't determine f({K})")
        break
    f_vals.append(rhs / b_KK)

print("f(k) values:")
for k in range(min(15, len(f_vals))):
    print(f"  f({k}) = {f_vals[k]}")
    # Check if f is a polynomial in k
    if k >= 1:
        ratio = f_vals[k] / f_vals[k-1] if f_vals[k-1] != 0 else "inf"
        if isinstance(ratio, F):
            print(f"    f({k})/f({k-1}) = {ratio} = {float(ratio):.6f}")

# Verify: these f(k) should reproduce Q̂_N for ALL N, not just N ≥ k
print("\n=== Verification ===")
for N in range(NMAX + 1):
    val = sum(f_vals[k] * delannoy_summand(N, k) for k in range(min(N+1, len(f_vals))))
    err = val - qhat[0][N]
    if err != 0:
        print(f"  N={N}: FAIL, residual = {err}")
        break
    else:
        if N < 10 or N % 5 == 0:
            print(f"  N={N}: OK")

# Check if f(k) is a polynomial
print("\n=== Checking if f(k) is polynomial ===")
# Compute finite differences
diffs = [list(f_vals)]
for order in range(1, min(8, len(f_vals))):
    new_diffs = []
    for i in range(len(diffs[-1]) - 1):
        new_diffs.append(diffs[-1][i+1] - diffs[-1][i])
    diffs.append(new_diffs)
    if all(d == 0 for d in new_diffs):
        print(f"  f(k) is polynomial of degree {order - 1}")
        break
    else:
        print(f"  Δ^{order} f: first few = {[float(d) for d in new_diffs[:5]]}")

# Check if f(k) has a hypergeometric pattern: f(k+1)/f(k) rational in k
print("\n=== Checking if f(k) is hypergeometric ===")
for k in range(min(12, len(f_vals) - 1)):
    if f_vals[k] != 0:
        r = f_vals[k+1] / f_vals[k]
        print(f"  f({k+1})/f({k}) = {float(r):.15f}")

# Also check for the q-row trajectory
print("\n\n=== Same analysis for Q̂_N from q-row ===")
f_vals_q = []
for K in range(NMAX + 1):
    rhs = Qhat_q[K]
    for k in range(K):
        rhs -= f_vals_q[k] * delannoy_summand(K, k)
    b_KK = delannoy_summand(K, K)
    if b_KK == 0:
        break
    f_vals_q.append(rhs / b_KK)

print("f_q(k) values:")
for k in range(min(10, len(f_vals_q))):
    print(f"  f_q({k}) = {f_vals_q[k]} = {float(f_vals_q[k]):.6f}")

print("\nDone.")
