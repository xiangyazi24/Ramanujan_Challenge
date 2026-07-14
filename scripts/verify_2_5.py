#!/usr/bin/env python3
"""Problem 2.5: Catalan's constant G from 3×3 CMF.

Compute the matrix product, extract scalar recurrence, find the summation-lift
factorization, and identify the order-2 kernel.
"""
from mpmath import mp, mpf, catalan, nstr, log, pi, matrix
from fractions import Fraction as F

mp.dps = 200

# ---- Matrix entries m_{ij}(n) from the challenge ----
def M(n):
    """3x3 matrix M(n) with exact rational entries."""
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

    return matrix([[m11, m12, m13],
                   [m21, m22, m23],
                   [m31, m32, m33]])

# Initial matrix A (2×3)
A = matrix([[mpf(30921), mpf(-32972), mpf(8240)],
            [mpf(33750), mpf(-36000), mpf(9000)]])

# ---- Compute matrix products and verify convergence ----
print("=== Problem 2.5: Catalan's Constant G ===")
print(f"G = {nstr(catalan, 50)}")

# Compute A * M(0) * M(1) * ... * M(N-1)
N_max = 80
prod = M(0)
for n in range(1, N_max):
    prod = prod * M(n)

result = A * prod
print(f"\nA * M(0)...M({N_max-1}):")
for j in range(3):
    Pj = result[0, j]
    Qj = result[1, j]
    ratio = Pj / Qj
    err = float(abs(ratio - catalan))
    print(f"  j={j+1}: P/Q = {nstr(ratio, 30)}, err = {err:.3e}")

# ---- Extract scalar recurrence via direct computation ----
# Compute the sequence q_N = Q_{N,1} (first column of Q row)
print("\n=== Scalar sequence extraction ===")

# We need to track the 3 columns independently
# T_N = M(0) * M(1) * ... * M(N-1), then A * T_N gives P_{N,j} and Q_{N,j}
# Q_{N,j} = row 2 of A * T_N, column j

# For the scalar recurrence, we need Q_{N,1} for many N values.
# Build incrementally.

q_vals = []  # q[N] = Q_{N,1}
T = matrix([[1,0,0],[0,1,0],[0,0,1]])  # identity
for N in range(60):
    AT = A * T
    q_val = AT[1, 0]  # Q_{N,1}
    q_vals.append(q_val)
    T = T * M(N)

print("First few q_N values:")
for i in range(min(10, len(q_vals))):
    print(f"  q[{i}] = {nstr(q_vals[i], 20)}")

# Try to find the recurrence order and coefficients
# For order 3: alpha_3 q_{N} + alpha_2 q_{N-1} + alpha_1 q_{N-2} + alpha_0 q_{N-3} = 0
# We need q_N for N = 0, 1, 2, ..., and solve for polynomial coefficients alpha_i(N)

# Method: use sufficient q values to determine the recurrence
# With order 3 and degree d coefficients, we need 3+d+1 equations per degree determination

# First verify the order is 3 by checking determinant
print("\n=== Checking recurrence order ===")
# For order k: the Casorati determinant of k+1 consecutive values should be zero
for order in range(2, 6):
    pass

# Simpler approach: fit the recurrence at multiple N values
# For order 3: q_{N+3} + c_2(N) q_{N+2} + c_1(N) q_{N+1} + c_0(N) q_N = 0
# If c_i are polynomials of degree d, we have 3(d+1) unknowns

# Try Poincaré analysis: ratio q_{N+1}/q_N for large N
print("\n=== Poincaré analysis ===")
for N in [10, 20, 30, 40, 50]:
    if N+1 < len(q_vals) and q_vals[N] != 0:
        ratio = q_vals[N+1] / q_vals[N]
        print(f"  q[{N+1}]/q[{N}] = {nstr(ratio, 15)}")

# The q_N values grow very fast. Let me look at the growth rate.
print("\n=== Growth rate ===")
for N in range(1, min(30, len(q_vals))):
    if q_vals[N] != 0 and q_vals[N-1] != 0:
        ratio = q_vals[N] / q_vals[N-1]
        print(f"  q[{N}]/q[{N-1}] = {nstr(ratio, 12)}")

# Also compute the difference sequence v_N = q_{N+1} - q_N
# If L = L_1 * (S-1), then v_N satisfies L_1
print("\n=== Difference sequence v_N = q_{N+1} - q_N ===")
v_vals = [q_vals[i+1] - q_vals[i] for i in range(len(q_vals)-1)]
print("First few v_N values:")
for i in range(min(10, len(v_vals))):
    print(f"  v[{i}] = {nstr(v_vals[i], 20)}")

print("\nRatio v[N+1]/v[N]:")
for N in range(1, min(20, len(v_vals)-1)):
    if v_vals[N] != 0:
        ratio = v_vals[N+1] / v_vals[N]
        print(f"  v[{N+1}]/v[{N}] = {nstr(ratio, 15)}")

# If v_N satisfies a 2-term recurrence, successive ratios should approach a limit
# The limit should be one of the Poincaré roots of L_1, i.e., eigenvalues of C_2

# Check: do the ratios approach (3+2√2)² = 17+12√2 ≈ 33.97?
from mpmath import sqrt
silver_sq = (3 + 2*sqrt(mpf(2)))**2
print(f"\n(3+2√2)² = 17+12√2 = {nstr(silver_sq, 15)}")
print(f"Compare to large-N ratio v[N+1]/v[N]")

# Or maybe the Poincaré root is -16(17+12√2) = -(1+√2)⁴ * 16?
# From the proof.tex: χ(t) = (t-1)(t²-34t+1), roots 1, 17±12√2
# After dividing by (S-1), the Poincaré roots of L_1 are 17±12√2
# These are (1+√2)⁴ and (1-√2)⁴ = (√2-1)⁴
pp_root = 17 + 12*sqrt(mpf(2))
print(f"17+12√2 = {nstr(pp_root, 15)}")
print(f"(1+√2)⁴ = {nstr((1+sqrt(mpf(2)))**4, 15)}")

# But wait — the Poincaré roots of the ORDER-3 recurrence for q_N are
# c_1/(-16), c_2/(-16), c_3/(-16) normalized. The raw roots are -16, -16(17±12√2).
# After dividing by (S-1) (which removes the root 1 from the normalized poly),
# the remaining roots of L_1 are the roots of t²-34t+1 = 0, i.e. 17±12√2.
# The dominant root is 17+12√2 ≈ 33.97.

# So v[N+1]/v[N] should → 17+12√2 ≈ 33.97... Let's check
# Wait, we need to be more careful. The original Poincaré roots are for q, not v.
# If q satisfies L q = 0 with char poly (t+16)(t²+544t+256),
# then q_N ~ c₁·(-16)^N + c₂·(-16(17-12√2))^N + c₃·(-16(17+12√2))^N
# The dominant root is -16(17+12√2) ≈ -543.5
# So q_{N+1}/q_N → -16(17+12√2) for the dominant mode.

dominant = -16*(17+12*sqrt(mpf(2)))
print(f"\nDominant Poincaré root: -16(17+12√2) = {nstr(dominant, 15)}")

# For v_N = q_{N+1} - q_N, under the dominant mode:
# v_N ~ c₃·(-16(17+12√2))^N · [(-16(17+12√2)) - 1]
# So v_{N+1}/v_N → -16(17+12√2) as well. Same dominant root.

# The (S-1) removes the NEUTRAL mode (Poincaré root 1 in the normalized poly,
# which is root -16 in the raw poly). So in v_N, the -16 mode is absent.
# v_N has only the two silver-ratio modes.

# Numerically verify: p_N/q_N should approach G
print("\n=== Verifying p_N/q_N → G ===")
T = matrix([[1,0,0],[0,1,0],[0,0,1]])
for N in range(70):
    AT = A * T
    p_val = AT[0, 0]
    q_val = AT[1, 0]
    if N >= 50 and q_val != 0:
        ratio = p_val / q_val
        err = float(abs(ratio - catalan))
        print(f"  N={N}: P/Q - G = {err:.3e}")
    T = T * M(N)

# ---- OEIS check: Q_N denominators ----
print("\n=== OEIS check ===")
# q[0] should be the first element of the Q sequence
# From A * identity = A itself, Q_{0,1} = A[1,0] = 33750
# Actually let me print all q values again
print("Raw q_N values (integer check):")
for i in range(8):
    print(f"  q[{i}] = {nstr(q_vals[i], 30)}")
