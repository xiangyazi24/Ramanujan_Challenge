#!/usr/bin/env python3
"""Extract the c=-16 gauge mode by cancelling the dominant mode.

Strategy:
1. Compute 3 independent forward solutions from matrix columns
2. Cancel the dominant (c≈-543) mode using large-N ratios
3. The remaining signal is the gauge (c=-16) mode
4. Compute r(n) = gauge(n+1)/gauge(n) → rational function
"""
from mpmath import mp, mpf, nstr, matrix
mp.dps = 200

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

# Compute matrix product T_N = M(0)·M(1)·...·M(N-1)
N_max = 60
T = matrix([[1,0,0],[0,1,0],[0,0,1]])
T_vals = [None] * (N_max + 1)
T_vals[0] = matrix([[1,0,0],[0,1,0],[0,0,1]])

for N in range(N_max):
    T = T * M_mat(N)
    T_vals[N+1] = T.copy()

# Three forward solutions: s_j(N) = T_N[0,j] (first row, columns 0,1,2)
s = [[T_vals[N][0,j] for N in range(N_max+1)] for j in range(3)]

# At large N, s_j(N) ≈ A_j * h_dominant(N) + B_j * h_gauge(N) + C_j * h_recessive(N)
# Dominant mode ratio: s_j(N) / s_0(N) → A_j / A_0
N_ref = 50  # Reference point for dominant mode ratio
ratio_1 = s[1][N_ref] / s[0][N_ref]  # should approach A_1/A_0
ratio_2 = s[2][N_ref] / s[0][N_ref]  # should approach A_2/A_0

# Also check stability of these ratios
ratio_1_49 = s[1][49] / s[0][49]
ratio_2_49 = s[2][49] / s[0][49]
print("=== Dominant mode eigenvector ratios ===")
print(f"  s1/s0 at N=49: {nstr(ratio_1_49, 20)}")
print(f"  s1/s0 at N=50: {nstr(ratio_1, 20)}")
print(f"  s2/s0 at N=49: {nstr(ratio_2_49, 20)}")
print(f"  s2/s0 at N=50: {nstr(ratio_2, 20)}")
print(f"  Stability: {nstr(abs(ratio_1 - ratio_1_49)/abs(ratio_1), 5)}")

# Cancel dominant mode:
# g_j(N) = s_j(N) - (A_j/A_0) * s_0(N)  for j=1,2
# These two solutions have the dominant mode cancelled.
g1 = [s[1][N] - ratio_1 * s[0][N] for N in range(N_max+1)]
g2 = [s[2][N] - ratio_2 * s[0][N] for N in range(N_max+1)]

# g1 and g2 are dominated by the gauge (c=-16) mode for large N.
# Check: ratio g1(N+1)/g1(N) should approach r(N) ~ -16 N^7
print("\n=== Gauge-mode ratios from g1 ===")
for N in [5, 8, 10, 15, 20, 25, 30, 35, 40, 45]:
    if g1[N] != 0:
        r = g1[N+1] / g1[N]
        rn7 = r / (mpf(N)**7)
        print(f"  g1({N+1})/g1({N}) = {nstr(r, 15)},  / N^7 = {nstr(rn7, 10)}")

# Also try with g2
print("\n=== Gauge-mode ratios from g2 ===")
for N in [5, 8, 10, 15, 20, 25, 30, 35, 40, 45]:
    if g2[N] != 0:
        r = g2[N+1] / g2[N]
        rn7 = r / (mpf(N)**7)
        print(f"  g2({N+1})/g2({N}) = {nstr(r, 15)},  / N^7 = {nstr(rn7, 10)}")

# If the dominant mode cancellation worked, g1 is gauge + recessive.
# For N not too large, the gauge dominates over the recessive.
# The ratio g1(N+1)/g1(N) should be close to r(N) = h(N+1)/h(N).

# Now identify r(N) as a rational function.
# r(N) / (-16) should be a ratio of degree-7 polynomials.
# Compute r(N) / (-16) for several N values:
print("\n=== r(N)/(-16) from g1 ===")
r_values = {}
for N in range(3, 40):
    if g1[N] != 0 and g1[N+1] != 0:
        r = g1[N+1] / g1[N]
        r_values[N] = r / mpf(-16)
        if N <= 20:
            print(f"  N={N:2d}: r/(-16) = {nstr(r_values[N], 25)}")

# Try to match r(N)/(-16) = P(N)/Q(N) where P,Q are degree-7 polynomials.
# Use polynomial interpolation: if r/(-16) = P/Q, then r/(-16) * Q(N) - P(N) = 0.
# This gives a system of linear equations.
# With 16 unknowns (8 for P + 8 for Q, minus 1 normalization = 15),
# we need at least 15 sample points.

# Alternatively, try to identify the FACTORED form:
# r(N) = -16 * (N+a1)(N+a2)...(N+a7) / ((N+b1)(N+b2)...(N+b7))
# The ai and bi should be half-integers or integers matching the recurrence structure.

# Check: is r(N)/(-16) close to a simple product of shifted factorials?
# From the recurrence, we expect parameters like (1/2, 1, 3/2, 2, 5/2, 3, 7/2)
# or similar _7F_6 parameters.

# Let's try to find zeros of r(N) (where h(N+1) = 0 → r(N) = 0 → N = -ai):
# r(N) = 0 means N + ai = 0 for some i, i.e., ai is a positive integer or half-integer.
# Check: r(-k) = 0 for which k?
# We can't evaluate r at negative N directly, but we can check:
# If ai = k (integer), then r(k-ai) = r(0) involves N=0, etc.

# Actually, let me just try fitting. Using 15 values of r/(-16) at N=3,...,17,
# fit P(N) and Q(N).

# Set Q(N) = N^7 + b6 N^6 + ... + b0 (monic degree 7)
# P(N) = N^7 + a6 N^6 + ... + a0 (monic degree 7, since r/(-16) ~ N^7/N^7 = 1 leading)

# Wait, r ~ -16 N^7, so r/(-16) ~ N^7. If P and Q are both monic degree 7,
# then P/Q ~ 1 for large N. But r/(-16) = N^7 * (1 + ...). Hmm.

# Actually, r(N) = -16 * (N^7 + ...) / (1 + ...) would give r/(-16) ~ N^7.
# So P is degree 7 (monic or LC=1) and Q is degree 0 (constant)?
# No, that would make r(N) polynomial.

# Let me check: is r(N)/(-16) close to a polynomial (i.e., Q = constant)?
# For N=10: r/(-16) ≈ 1618043 and N^7 = 10^7 = 10000000. So r/(-16)/N^7 ≈ 0.16.
# Not approaching 1. So r/(-16) grows SLOWER than N^7.

# Actually from the ratios: r(N)/(-16N^7) is decreasing. At N=45:
if 45 in r_values:
    print(f"\n  N=45: r/(-16*N^7) = {nstr(r_values[45]/(mpf(45)**7), 10)}")

# Hmm, maybe r(N) is NOT -16 N^7 + lower. Let me reconsider.
# The Poincaré root c=-16 means the growth of the gauge solution is
# proportional to (-16)^N · (correction involving N).
# But the "N^7" part is absorbed into the Pochhammer factors.

# Actually, for a ₇F₆ hypergeometric:
# h(n+1)/h(n) = (-16) * (n+a1)...(n+a7) / ((n+b1)...(n+b7))
# For large n: (n+ai)/(n+bi) → 1, so r(n) → -16.
# NOT -16*n^7!

# Wait, that can't be right either. The Poincaré root is c=-16 with
# the scaling q_n ~ c^n * n^{7n} * .... The n^{7n} part comes from
# ∏ n^7 = (n!)^7. So:
# h(n) = (-16)^n * (n!)^7 * R(n)
# where R(n) is a slowly-varying rational function.
# h(n+1)/h(n) = (-16) * (n+1)^7 * R(n+1)/R(n) ≈ -16 * n^7 * (1+7/n+...).

# But if h(n) = (-16)^n * ∏_{k=1}^n (k+a1)...(k+a7) / ∏(k+b1)...(k+b7),
# then h(n+1)/h(n) = (-16) * (n+1+a1)...(n+1+a7)/((n+1+b1)...(n+1+b7))
# For large n, this is ≈ -16 * n^0 = -16 (ratio of monic degree-7 polys → 1).

# This contradicts r(n) ~ -16 n^7!

# The discrepancy is that in the standard hypergeometric framework, the
# FIRST parameter in the ratio r(n) = (n+a1)...(n+a7)/((n+b1)...(n+b7))
# gives r(n) → 1 for large n. The Pochhammer growth is:
# ∏_{k=1}^n (k+a) / (k+b) ~ n^{a-b} * const.
# So h(n) = (-16)^n * ∏_{k=1}^n [(k+a1)...(k+a7) / ((k+b1)...(k+b7))]
# grows as (-16)^n * n^{Σ(ai-bi)}.

# For the Poincaré root c=-16 with degree pattern (28,21,14,7):
# q_n ~ h(n) grows as c^n * n^{cn} (superexponential) where the n^{cn}
# comes from the degree gap.

# Hmm, I'm confusing myself. Let me go back to basics.
# The recurrence has degree pattern (28,21,14,7). This means:
# c₃(N) ~ K₃ N^7, c₂(N) ~ K₂ N^14, c₁(N) ~ K₁ N^21, c₀(N) ~ K₀ N^28.

# For a solution q_N ~ μ^N * N^(7N):
# The ratio q_{N+1}/q_N ~ μ * (N+1)^7 ~ μ * N^7.
# So the ratio DOES scale as N^7, and the "Poincaré root" μ = -16 is the
# coefficient in front of N^7.

# So r(N) = h(N+1)/h(N) ~ -16 * N^7 * (1 + O(1/N)).
# The gauge h(N) = ∏_{k=0}^{N-1} r(k).

# For r(N) to be rational and ~ -16 N^7:
# r(N) = -16 * P(N) / Q(N) where P has degree 7 and Q has degree 0?
# Then r(N) ~ -16 * N^7 / constant, which works.
# But Q being degree 0 means r(N) is a polynomial times -16.
# r(N) = -16 * polynomial(N).

# Alternatively: P has degree d_P and Q has degree d_Q with d_P - d_Q = 7.
# Could be d_P = 14, d_Q = 7, for instance.

# From the actual data: r(N)/(-16) is growing but not as N^7.
# At N=10: r/(-16) ≈ 1.618e6, while N^7 = 1e7. So r/(-16)/N^7 ≈ 0.16.
# At N=20: r/(-16) ≈ 1.495e9, while N^7 = 1.28e9. So r/(-16)/N^7 ≈ 1.17.

# Hmm, it's approaching 1? Let me check more carefully.
# Actually wait, the data above shows r/(-16N^7) DECREASING from 2.7 at N=3 to 0.058 at N=27.
# But that was from the backward recursion, which selected the WRONG mode.
# The g1 data (dominant-mode-cancelled) should be different.

# Let me just look at the g1 ratios more carefully.
