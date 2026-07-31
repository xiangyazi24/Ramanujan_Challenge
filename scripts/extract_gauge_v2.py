#!/usr/bin/env python3
"""Extract the c=-16 gauge mode using Richardson-extrapolated dominant cancellation.

Key insight: s_j(N)/s_0(N) = R_j + alpha_j * rho^N + O(rho'^N)
where rho = c_1/c_3 ≈ 0.029 and rho' = c_2/c_3 ≈ 0.00087.
Aitken Delta^2 extrapolation removes the rho^N term, giving R_j to ~154 digits.
Then g(N) = s_j(N) - R_j * s_0(N) is dominated by the gauge mode.
"""
from mpmath import mp, mpf, nstr, matrix, power, fac, log10, inf
mp.dps = 300  # Extra precision for the extrapolation

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

# Compute matrix product and extract column solutions
N_max = 80
T = matrix([[1,0,0],[0,1,0],[0,0,1]])
# Store all three row-0 entries for each N
s = {j: [] for j in range(3)}

for N in range(N_max + 1):
    for j in range(3):
        s[j].append(T[0,j])
    if N < N_max:
        T = T * M_mat(N)

# Compute f_j(N) = s_j(N) / s_0(N) at several N values around N=60-70
print("=== Computing eigenvector ratios via Richardson extrapolation ===")

# Use N=65,66,67 for Aitken extrapolation
for jj in [1, 2]:
    f = {}
    for N in range(55, 76):
        f[N] = s[jj][N] / s[0][N]
    
    # Aitken Delta^2: R = f(N) - (f(N+1)-f(N))^2 / (f(N+2)-2f(N+1)+f(N))
    R_vals = {}
    for N in range(55, 73):
        d1 = f[N+1] - f[N]
        d2 = f[N+2] - 2*f[N+1] + f[N]
        if abs(d2) > 0:
            R_vals[N] = f[N] - d1**2 / d2
    
    # Check stability of Richardson estimates
    Ns = sorted(R_vals.keys())
    print(f"\n  s_{jj}/s_0 Richardson estimates:")
    for N in [55, 60, 65, 70]:
        if N in R_vals:
            print(f"    N={N}: {nstr(R_vals[N], 40)}")
    
    # Check convergence by looking at differences between successive Richardson estimates
    if len(Ns) >= 2:
        for i in range(max(0, len(Ns)-3), len(Ns)-1):
            N1, N2 = Ns[i], Ns[i+1]
            diff = abs(R_vals[N2] - R_vals[N1])
            if diff > 0:
                digits = -float(log10(diff / abs(R_vals[N2])))
            else:
                digits = 300
            print(f"    Stability N={N1}→{N2}: {digits:.0f} agreeing digits")

# Use the best Richardson estimate (N=70) for cancellation
R1 = None
R2 = None
for jj in [1, 2]:
    f = {}
    for N in range(55, 76):
        f[N] = s[jj][N] / s[0][N]
    d1 = f[71] - f[70]
    d2 = f[72] - 2*f[71] + f[70]
    R = f[70] - d1**2 / d2
    if jj == 1:
        R1 = R
    else:
        R2 = R

print(f"\n  R1 (s1/s0 limit) = {nstr(R1, 50)}")
print(f"  R2 (s2/s0 limit) = {nstr(R2, 50)}")

# Now cancel the dominant mode using the precise R values
g = [s[1][N] - R1 * s[0][N] for N in range(N_max+1)]

# The cancelled solution g(N) is dominated by the gauge (c=-16) mode.
# Compute r(N) = g(N+1) / g(N)
print("\n=== Gauge-mode ratios r(N) = g(N+1)/g(N) ===")
print("  Expected: r(N) ~ -16 * N^7 for large N\n")

for N in range(3, 55):
    if g[N] != 0 and g[N+1] != 0:
        r = g[N+1] / g[N]
        rn7 = r / (mpf(-16) * mpf(N)**7)
        if N <= 20 or N % 5 == 0:
            print(f"  N={N:2d}: r(N)/(-16*N^7) = {nstr(rn7, 30)}")

# The ratio r(N)/(-16*N^7) should converge to 1 as N → ∞.
# But more precisely, r(N) = -16 * (N+a1)(N+a2)...(N+a7) / ((N+b1)(N+b2)...(N+b7))
# or r(N) = -16 * P(N) where P is a degree-7 polynomial if the denominator is constant.

# Let's check: is r(N)/(-16) a POLYNOMIAL in N?
# Compute the "polynomial remainder" by taking successive differences.
# For a degree-d polynomial, d+1 successive differences give 0.
print("\n=== Testing if r(N)/(-16) is a polynomial ===")
vals = []
for N in range(3, 25):
    if g[N] != 0 and g[N+1] != 0:
        r = g[N+1] / g[N]
        vals.append(r / mpf(-16))

# Apply 8 successive differences (degree 7 → 8th diff = 0)
diffs = [list(vals)]
for d in range(1, 10):
    prev = diffs[-1]
    new_diff = [prev[i+1] - prev[i] for i in range(len(prev)-1)]
    diffs.append(new_diff)
    if len(new_diff) > 2:
        # Check if this level is approximately 0
        max_val = max(abs(x) for x in new_diff)
        ref = max(abs(x) for x in diffs[0])
        ratio = float(max_val / ref) if ref != 0 else float(max_val)
        print(f"  Diff order {d}: max |Δ^{d}| / |r/16| = {ratio:.3e}", 
              "← ZERO!" if ratio < 1e-100 else "")

