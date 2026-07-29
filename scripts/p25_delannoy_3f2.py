#!/usr/bin/env python3
"""
P2.5: Test whether Φ(w) = Σ_k c_k w^k is a ₃F₂ at 32w.

Key insight: F(z) = Σ u₁(N) z^N, and with w = z/(1-z)²:
F(z) = 1/(1-z) · Φ(w)  where Φ(w) = Σ_k α₁(k) · 2^k C(2k,k)² w^k

Both z=ρ and z=1/ρ map to w=1/32, z=1→∞.
So Φ has only 3 singular points: 0, 1/32, ∞ → could be ₃F₂.

If Φ(w) = ₃F₂(a,b,c; d,e; 32w), then:
c_{k+1}/c_k = 32(a+k)(b+k)(c+k)/((d+k)(e+k)(k+1))

Test: compute c_k = α₁(k) · 2^k C(2k,k)², then check if ratio
is a rational function of k of the form above.
"""
from fractions import Fraction as Fr
from math import comb
from mpmath import mp, mpf, identify

mp.dps = 50

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
        return Fr(0)
    return Fr(2**k * comb(2*k, k) * comb(N, k) * comb(N+k, k))

KMAX = 60
print(f"Computing CMF trajectories for N=0..{KMAX}...", flush=True)

rows = [[Fr(1), Fr(0), Fr(0)],
        [Fr(0), Fr(1), Fr(0)],
        [Fr(0), Fr(0), Fr(1)]]
u = {j: [rows[j][0]] for j in range(3)}

for N in range(KMAX):
    M = M_entries(N)
    d = Fr(delta_H(N))
    MH = [[Fr(M[i][j]) / d for j in range(3)] for i in range(3)]
    for j in range(3):
        r = rows[j]
        new_r = [sum(r[i]*MH[i][k] for i in range(3)) for k in range(3)]
        rows[j] = new_r
        u[j].append(new_r[0])
    if (N+1) % 20 == 0:
        print(f"  N={N+1} done", flush=True)

# Triangular inversion: α₁(k)
print(f"\nTriangular inversion for α₁(k)...", flush=True)
alpha1 = []
for K in range(KMAX + 1):
    rhs = u[0][K]
    for k in range(K):
        rhs -= alpha1[k] * B(K, k)
    alpha1.append(rhs / B(K, K))

# c_k = α₁(k) · 2^k · C(2k,k)²
print(f"\nComputing c_k = α₁(k) · 2^k · C(2k,k)²...", flush=True)
c = []
for k in range(KMAX + 1):
    ck = alpha1[k] * Fr(2**k * comb(2*k, k)**2)
    c.append(ck)

# Compute ratios c_{k+1}/c_k
print(f"\nRatios c_{{k+1}}/c_k:")
ratios = []
for k in range(KMAX):
    r = c[k+1] / c[k]
    ratios.append(r)
    if k < 15 or k % 10 == 0:
        print(f"  k={k}: c_{{k+1}}/c_k = {float(r):.15f}", flush=True)

# If ₃F₂: ratio = 32(a+k)(b+k)(c+k)/((d+k)(e+k)(k+1))
# So ratio/32 = (a+k)(b+k)(c+k)/((d+k)(e+k)(k+1))
# At large k: ratio → 32.
# Define R(k) = ratio(k)/32 - 1 → 0 as k → ∞
# R(k) = (a+b+c-d-e-1)/k + O(1/k²)

print(f"\nChecking if ratio/32 is rational in k:")
print(f"  ratio/32 for large k:")
for k in range(KMAX - 5, KMAX):
    print(f"  k={k}: ratio/32 = {float(ratios[k])/32:.15f}")

# More sophisticated: compute (ratio/32)(k)(k+1) = (a+k)(b+k)(c+k)/((d+k)(e+k))
# Let P(k) = (ratio/32)(k) · (k+1) · (d+k)(e+k) = (a+k)(b+k)(c+k)
# First compute R(k) = ratio(k) · (k+1) / 32
print(f"\n'Numerator polynomial' R(k) = ratio(k)·(k+1)/32:")
R = []
for k in range(KMAX):
    rk = ratios[k] * Fr(k+1) / Fr(32)
    R.append(rk)
    if k < 10 or k % 10 == 0:
        print(f"  k={k}: R(k) = {float(rk):.15f}", flush=True)

# If ₃F₂: R(k) = (a+k)(b+k)(c+k)/((d+k)(e+k))
# This is a degree-3-over-degree-2 rational function.
# At large k: R(k) ~ k + (a+b+c-d-e) + O(1/k)

# Compute R(k)/k for large k → 1 + const/k
print(f"\nR(k)/k for large k (should approach 1 if ₃F₂):")
for k in range(max(1, KMAX-10), KMAX):
    print(f"  k={k}: R/k = {float(R[k])/k:.15f}")

# Compute R(k) - k for large k → a+b+c-d-e + O(1/k)
print(f"\nR(k) - k (should approach constant if ₃F₂):")
for k in range(max(1, KMAX-10), KMAX):
    print(f"  k={k}: R-k = {float(R[k]) - k:.15f}")

# Compute successive differences to check if R(k) is exactly rational in k
# If R(k) = p(k)/q(k) with p degree 3, q degree 2,
# then R is determined by 6 parameters.
# Interpolate from 7 values.
print(f"\n--- Fitting rational function ---")
print(f"R(k) = (a+k)(b+k)(c+k)/((d+k)(e+k))")
print(f"= (k³ + s₂k² + s₁k + s₀) / (k² + t₁k + t₀)")
print(f"where s₂=a+b+c, s₁=ab+ac+bc, s₀=abc, t₁=d+e, t₀=de")

# From R(k) = P(k)/Q(k), P(k) = R(k)·Q(k)
# P(k) = k³ + s₂k² + s₁k + s₀
# Q(k) = k² + t₁k + t₀
# So: R(k)·(k² + t₁k + t₀) = k³ + s₂k² + s₁k + s₀
# Rearranging: R(k)·k² + R(k)·t₁·k + R(k)·t₀ = k³ + s₂k² + s₁k + s₀
# We have 5 unknowns: s₂, s₁, s₀, t₁, t₀
# Use least squares from many k values

# Method: R(k)(k² + t₁k + t₀) - k³ - s₂k² - s₁k - s₀ = 0
# → R(k)k² + t₁·R(k)·k + t₀·R(k) - k³ - s₂·k² - s₁·k - s₀ = 0
# → t₁·R(k)·k + t₀·R(k) - s₂·k² - s₁·k - s₀ = k³ - R(k)·k²
# Unknowns: t₁, t₀, s₂, s₁, s₀

# Linear system: for each k
# [R(k)·k,  R(k),  -k²,  -k,  -1] · [t₁, t₀, s₂, s₁, s₀]^T = k³ - R(k)·k²

n_unknowns = 5
# Use k = 1,...,n_unknowns+3 for overdetermined system
rows_data = []
rhs_data = []
for k in range(1, 20):
    rk = R[k]
    row = [rk * k, rk, -Fr(k**2), -Fr(k), Fr(-1)]
    rows_data.append(row)
    rhs_data.append(Fr(k**3) - rk * Fr(k**2))

# Gaussian elimination on exact fractions
n = n_unknowns
# Use first n equations
aug = [rows_data[i][:] + [rhs_data[i]] for i in range(n)]
for col in range(n):
    pivot_row = None
    for r in range(col, n):
        if aug[r][col] != 0:
            pivot_row = r
            break
    if pivot_row is None:
        print(f"Singular at col {col}!")
        break
    aug[col], aug[pivot_row] = aug[pivot_row], aug[col]
    for r in range(n):
        if r != col and aug[r][col] != 0:
            factor = aug[r][col] / aug[col][col]
            for c2 in range(n + 1):
                aug[r][c2] -= factor * aug[col][c2]

sol = [aug[i][n] / aug[i][i] for i in range(n)]
t1, t0, s2, s1, s0 = sol

print(f"\nSolution (if ₃F₂):")
print(f"  t₁ = d+e = {float(t1):.15f} (exact: {t1})")
print(f"  t₀ = de  = {float(t0):.15f} (exact: {t0})")
print(f"  s₂ = a+b+c = {float(s2):.15f} (exact: {s2})")
print(f"  s₁ = ab+ac+bc = {float(s1):.15f} (exact: {s1})")
print(f"  s₀ = abc = {float(s0):.15f} (exact: {s0})")

# Verify with additional data
print(f"\nVerification:")
max_err = Fr(0)
for k in range(1, KMAX):
    predicted = (Fr(k**3) + s2*Fr(k**2) + s1*Fr(k) + s0) / (Fr(k**2) + t1*Fr(k) + t0)
    actual = R[k]
    err = abs(predicted - actual)
    if err > max_err:
        max_err = err
    if k < 5 or k == KMAX - 1:
        print(f"  k={k}: predicted={float(predicted):.15f}, actual={float(actual):.15f}, err={float(err):.2e}")

print(f"  max error over k=1..{KMAX-1}: {float(max_err):.2e}")

if max_err < Fr(1, 10**10):
    print(f"\n*** ₃F₂ HYPOTHESIS CONFIRMED ***")
    # Find a, b, c from s₂, s₁, s₀ (roots of t³ - s₂t² + s₁t - s₀ = 0)
    # Find d, e from t₁, t₀ (roots of t² - t₁t + t₀ = 0)
    import numpy as np
    abc_roots = np.roots([1, -float(s2), float(s1), -float(s0)])
    de_roots = np.roots([1, -float(t1), float(t0)])
    print(f"  a, b, c = {abc_roots}")
    print(f"  d, e = {de_roots}")

    # Also try to identify as simple fractions
    for name, val in [("a+b+c", s2), ("ab+ac+bc", s1), ("abc", s0), ("d+e", t1), ("de", t0)]:
        print(f"  {name} = {val} = {float(val):.15f}")
        # Check if it's a simple fraction p/q with small q
        for denom in range(1, 200):
            numer = round(float(val) * denom)
            if abs(val - Fr(numer, denom)) < Fr(1, 10**10):
                print(f"    → {numer}/{denom}")
                break
else:
    print(f"\n  ₃F₂ hypothesis REJECTED (max error too large)")
    # Try higher degree: maybe it's ₄F₃ or has polynomial corrections
    print(f"\n--- Trying ₄F₃ ansatz ---")
    # R(k) = (a+k)(b+k)(c+k)(f+k)/((d+k)(e+k)(g+k)(k+1)) would give
    # ratio = 32 · P4(k)/P4(k) — but then the ratio already has (k+1) in denominator
    # Actually for ₄F₃: ratio = 32(a+k)(b+k)(c+k)(f+k)/((d+k)(e+k)(g+k)(k+1))
    # So R(k) = ratio·(k+1)/32 = (a+k)(b+k)(c+k)(f+k)/((d+k)(e+k)(g+k))
    # This is degree 4 over degree 3.

    # Let's check: R(k) ~ k + const for large k already...
    # For ₄F₃: R(k) ~ k + const + O(1/k).
    # For ₃F₂: R(k) ~ k + const + O(1/k).
    # Same leading behavior!

    # The difference: ₃F₂ has R = deg3/deg2, ₄F₃ has R = deg4/deg3.

    # Try R(k) = (k⁴+...)/( k³+...)
    # 7 unknowns: s₃,s₂,s₁,s₀, t₂,t₁,t₀
    print("Fitting degree-4/degree-3 rational function...")
    n_unknowns2 = 7
    rows_data2 = []
    rhs_data2 = []
    for k in range(1, 30):
        rk = R[k]
        row = [rk * Fr(k**3), rk * Fr(k**2), rk * Fr(k), rk,
               -Fr(k**3), -Fr(k**2), -Fr(k)]
        rows_data2.append(row)
        rhs_data2.append(Fr(k**4) - rk * Fr(k**3) + Fr(0))
        # Wait, need to be more careful

    # R(k) · (k³ + t₂k² + t₁k + t₀) = k⁴ + s₃k³ + s₂k² + s₁k + s₀
    # R(k)·k³ + t₂·R(k)·k² + t₁·R(k)·k + t₀·R(k) = k⁴ + s₃·k³ + s₂·k² + s₁·k + s₀
    # → t₂·R(k)·k² + t₁·R(k)·k + t₀·R(k) - s₃·k³ - s₂·k² - s₁·k - s₀ = k⁴ - R(k)·k³
    rows_data2 = []
    rhs_data2 = []
    for k in range(1, 30):
        rk = R[k]
        row = [rk*Fr(k**2), rk*Fr(k), rk, -Fr(k**3), -Fr(k**2), -Fr(k), Fr(-1)]
        rows_data2.append(row)
        rhs_data2.append(Fr(k**4) - rk * Fr(k**3))

    aug2 = [rows_data2[i][:] + [rhs_data2[i]] for i in range(n_unknowns2)]
    for col in range(n_unknowns2):
        pivot_row = None
        for r in range(col, n_unknowns2):
            if aug2[r][col] != 0:
                pivot_row = r
                break
        if pivot_row is None:
            print(f"Singular at col {col}!")
            break
        aug2[col], aug2[pivot_row] = aug2[pivot_row], aug2[col]
        for r in range(n_unknowns2):
            if r != col and aug2[r][col] != 0:
                factor = aug2[r][col] / aug2[col][col]
                for c2 in range(n_unknowns2 + 1):
                    aug2[r][c2] -= factor * aug2[col][c2]

    sol2 = [aug2[i][n_unknowns2] / aug2[i][i] for i in range(n_unknowns2)]
    t2_4, t1_4, t0_4, s3_4, s2_4, s1_4, s0_4 = sol2

    print(f"\n₄F₃ Solution:")
    print(f"  t₂ = {float(t2_4):.15f}")
    print(f"  t₁ = {float(t1_4):.15f}")
    print(f"  t₀ = {float(t0_4):.15f}")
    print(f"  s₃ = {float(s3_4):.15f}")
    print(f"  s₂ = {float(s2_4):.15f}")
    print(f"  s₁ = {float(s1_4):.15f}")
    print(f"  s₀ = {float(s0_4):.15f}")

    # Verify
    max_err2 = Fr(0)
    for k in range(1, KMAX):
        numer = Fr(k**4) + s3_4*Fr(k**3) + s2_4*Fr(k**2) + s1_4*Fr(k) + s0_4
        denom = Fr(k**3) + t2_4*Fr(k**2) + t1_4*Fr(k) + t0_4
        predicted = numer / denom
        actual = R[k]
        err = abs(predicted - actual)
        if err > max_err2:
            max_err2 = err
    print(f"  max error (₄F₃): {float(max_err2):.2e}")

    if max_err2 < Fr(1, 10**10):
        print(f"\n*** ₄F₃ HYPOTHESIS CONFIRMED ***")
        import numpy as np
        abcf_roots = np.roots([1, -float(s3_4), float(s2_4), -float(s1_4), float(s0_4)])
        deg_roots = np.roots([1, -float(t2_4), float(t1_4), -float(t0_4)])
        print(f"  a, b, c, f = {abcf_roots}")
        print(f"  d, e, g = {deg_roots}")

        for name, val in [("s₃", s3_4), ("s₂", s2_4), ("s₁", s1_4), ("s₀", s0_4),
                          ("t₂", t2_4), ("t₁", t1_4), ("t₀", t0_4)]:
            print(f"  {name} = {val}")
            for denom in range(1, 200):
                numer = round(float(val) * denom)
                if abs(val - Fr(numer, denom)) < Fr(1, 10**10):
                    print(f"    → {numer}/{denom}")
                    break

print("\nDone.")
