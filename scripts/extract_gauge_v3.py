#!/usr/bin/env python3
"""Extract gauge r(n) with enough precision.

Key realization: r(n) = -16 * (n+a1)...(n+a7) is a polynomial × (-16),
since the Poincaré root c=-16 with denominator = 1.

Strategy:
1. Compute matrix product to N=30 with mp.dps=1500 (entries are ~10^225,
   giving ~1275 digits of relative precision)
2. Cancel dominant mode using Richardson extrapolation
3. Extract r(n) values for n=5..25 (~20 values, overdetermined for 7 params)
4. Use the 7 values of r(n) at specific n to find the roots a_i
"""
from mpmath import mp, mpf, nstr, matrix, log10, fabs
mp.dps = 1500

def M_mat(n):
    n = mpf(n)
    m11 = (-2*n-5)*(n+3)**2 * (136*n**4 + 1424*n**3 + 5548*n**2 + 9551*n + 6141)
    m12 = 384*n**6 + 6384*n**5 + 44168*n**4 + 162698*n**3 + 336377*n**2 + 369933*n + 169011
    m13 = -480*n**4 - 4980*n**3 - 19210*n**2 - 32690*n - 20730
    m21 = (n+2)**2*(n+3)**2*(4*n+10)*(48*n**3 + 386*n**2 + 1017*n + 879)
    m22 = (n+2)**2*(-272*n**5 - 3848*n**4 - 21732*n**3 - 61184*n**2 - 85761*n - 47808)
    m23 = (n+2)**2*(320*n**3 + 2540*n**2 + 6610*n + 5640)
    m31 = (-4*n-10)*(n+2)**2*(n+3)**2*(32*n**4 + 302*n**3 + 1037*n**2 + 1530*n + 813)
    m32 = (n+2)**2*(192*n**6 + 2984*n**5 + 19116*n**4 + 64452*n**3 + 120256*n**2 + 117279*n + 46476)
    m33 = (n+2)**2*(-16*n**5 - 408*n**4 - 2912*n**3 - 8884*n**2 - 12254*n - 6240)
    return matrix([[m11, m12, m13], [m21, m22, m23], [m31, m32, m33]])

N_max = 35
T = matrix([[1,0,0],[0,1,0],[0,0,1]])
s = {j: [] for j in range(3)}

for N in range(N_max + 1):
    for j in range(3):
        s[j].append(T[0,j])
    if N < N_max:
        T = T * M_mat(N)

# Check precision: how many digits do we have at N=30?
val30 = s[0][30]
if val30 != 0:
    bits = int(float(log10(fabs(val30))))
    print(f"log10(|s0(30)|) ≈ {bits}")
    print(f"mp.dps = {mp.dps}, estimated relative precision ≈ {mp.dps - bits} digits")

# Compute eigenvector ratios f_j(N) = s_j(N)/s_0(N)
# These converge to A_j^dom / A_0^dom as N → ∞
# Use Aitken Δ² extrapolation

f1 = {}; f2 = {}
for N in range(5, N_max+1):
    if s[0][N] != 0:
        f1[N] = s[1][N] / s[0][N]
        f2[N] = s[2][N] / s[0][N]

# Check: how fast does f1(N) converge?
print("\n=== s1/s0 convergence ===")
for N in [10, 15, 20, 25, 30]:
    if N in f1 and N+1 in f1:
        diff = fabs(f1[N+1] - f1[N])
        rel = diff / fabs(f1[N]) if f1[N] != 0 else diff
        digits = -float(log10(rel)) if rel > 0 else 999
        print(f"  N={N}: f1 = {nstr(f1[N], 30)}, stability: {digits:.1f} digits")

# Richardson extrapolation: f(N) ≈ R + α ρ^N where ρ = c_gauge/c_dom ≈ 16/543.5 ≈ 0.029
# Use THREE-point extrapolation at N, N+1, N+2:
# ρ = (f(N+2)-f(N+1)) / (f(N+1)-f(N))
# R = f(N) - (f(N+1)-f(N))^2 / (f(N+2)-2f(N+1)+f(N))

# Then do iterated Aitken (Aitken on Aitken) for even more precision
def aitken(vals):
    """Apply Aitken Δ² to a sequence of values."""
    result = {}
    keys = sorted(vals.keys())
    for i in range(len(keys)-2):
        N = keys[i]
        d1 = vals[keys[i+1]] - vals[keys[i]]
        d2 = vals[keys[i+2]] - 2*vals[keys[i+1]] + vals[keys[i]]
        if d2 != 0:
            result[N] = vals[N] - d1**2 / d2
    return result

R1_aitken = aitken(f1)
R2_aitken = aitken(f2)

# Check stability of Aitken-level-1
print("\n=== Aitken level 1: s1/s0 ===")
keys1 = sorted(R1_aitken.keys())
for N in [10, 15, 20, 25]:
    if N in R1_aitken and N+1 in R1_aitken:
        diff = fabs(R1_aitken[N+1] - R1_aitken[N])
        rel = diff / fabs(R1_aitken[N]) if R1_aitken[N] != 0 else diff
        digits = -float(log10(rel)) if rel > 0 else 999
        print(f"  N={N}: R1 = {nstr(R1_aitken[N], 40)}, stability: {digits:.1f} digits")

# Apply iterated Aitken (level 2)
R1_aitken2 = aitken(R1_aitken)
R2_aitken2 = aitken(R2_aitken)

print("\n=== Aitken level 2: s1/s0 ===")
keys2 = sorted(R1_aitken2.keys())
for N in keys2[-5:]:
    if N+1 in R1_aitken2:
        diff = fabs(R1_aitken2[N+1] - R1_aitken2[N])
        rel = diff / fabs(R1_aitken2[N]) if R1_aitken2[N] != 0 else diff
        digits = -float(log10(rel)) if rel > 0 else 999
        print(f"  N={N}: R1'' = {nstr(R1_aitken2[N], 40)}, stability: {digits:.1f} digits")

# Use best estimate from Aitken level 2
if keys2:
    R1_best = R1_aitken2[keys2[-1]]
else:
    R1_best = R1_aitken[sorted(R1_aitken.keys())[-1]]

R2_aitken2_keys = sorted(R2_aitken2.keys())
if R2_aitken2_keys:
    R2_best = R2_aitken2[R2_aitken2_keys[-1]]
else:
    R2_best = R2_aitken[sorted(R2_aitken.keys())[-1]]

print(f"\n  Best R1 = {nstr(R1_best, 50)}")
print(f"  Best R2 = {nstr(R2_best, 50)}")

# Cancel dominant mode
g = [s[1][N] - R1_best * s[0][N] for N in range(N_max+1)]

# Extract gauge ratios
print("\n=== Gauge-mode ratios r(N) = g(N+1)/g(N) ===")
r_vals = {}
for N in range(3, 30):
    if g[N] != 0 and g[N+1] != 0:
        r = g[N+1] / g[N]
        r_vals[N] = r
        rn7 = r / (mpf(-16) * mpf(N)**7)
        if N <= 25:
            print(f"  N={N:2d}: r(N) = {nstr(r, 20)},  r/(-16*N^7) = {nstr(rn7, 15)}")

# Check: is r(N)/(-16) a degree-7 polynomial?
# If so, 8th finite difference should be 0
print("\n=== Finite differences of r(N)/(-16) ===")
vals = [r_vals[N] / mpf(-16) for N in range(3, 25) if N in r_vals]
diffs = [list(vals)]
for d in range(1, 11):
    prev = diffs[-1]
    new_diff = [prev[i+1] - prev[i] for i in range(len(prev)-1)]
    diffs.append(new_diff)
    if len(new_diff) > 0:
        max_val = max(fabs(x) for x in new_diff)
        ref = max(fabs(x) for x in diffs[0])
        ratio = float(max_val / ref) if ref > 0 else 0
        print(f"  Diff order {d}: max |Δ^{d}| / |r/16| = {ratio:.3e}",
              "← ZERO!" if ratio < 1e-100 else "")

# Also try: compute r(N) from SECOND cancellation (using s_2)
g2 = [s[2][N] - R2_best * s[0][N] for N in range(N_max+1)]
print("\n=== Cross-check: ratios from s2-cancelled ===")
for N in range(5, 20):
    if g2[N] != 0 and g2[N+1] != 0:
        r2 = g2[N+1] / g2[N]
        if N in r_vals and r_vals[N] != 0:
            rel_diff = fabs(r2 - r_vals[N]) / fabs(r_vals[N])
            digits = -float(log10(rel_diff)) if rel_diff > 0 else 999
            print(f"  N={N}: agreement = {digits:.0f} digits")
