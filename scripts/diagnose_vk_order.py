#!/usr/bin/env python3
"""Diagnose the recurrence order of v_k = q_{k+1} - q_k for Problem 2.5.

Use Casorati determinant ratios: if v_k satisfies an order-r recurrence with
polynomial coefficients, then the ratio of consecutive order-r Casorati determinants
is a rational function of k.

Also try: find the order-3 recurrence for q_N directly using the matrix structure.
"""
from mpmath import mp, mpf, nstr, matrix, catalan, log

mp.dps = 200

def M(n):
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

A_mat = matrix([[mpf(30921), mpf(-32972), mpf(8240)],
                [mpf(33750), mpf(-36000), mpf(9000)]])

N_max = 60
q_vals = []
T = matrix([[1,0,0],[0,1,0],[0,0,1]])
for N in range(N_max + 1):
    AT = A_mat * T
    q_vals.append(AT[1, 0])
    T = T * M(N)

v_vals = [q_vals[k+1] - q_vals[k] for k in range(N_max)]

# === Casorati determinant ratio diagnostic ===
print("=== Order-2 Casorati diagnostic ===")
print("If v_k has order 2: C2(k+1)/C2(k) should be a rational function of k\n")

def casorati2(k):
    return v_vals[k] * v_vals[k+2] - v_vals[k+1]**2

for k in range(5, 40):
    c2k = casorati2(k)
    c2k1 = casorati2(k+1)
    if c2k != 0:
        ratio = c2k1 / c2k
        print(f"  k={k:2d}: C2(k+1)/C2(k) = {nstr(ratio, 20)}, log10|C2| = {float(log(abs(c2k), 10)):.1f}")

print("\n=== Order-3 Casorati diagnostic (for v_k) ===")
print("If v_k has order 3: C3(k+1)/C3(k) should be a rational function of k\n")

def casorati3(k):
    return matrix([
        [v_vals[k], v_vals[k+1], v_vals[k+2]],
        [v_vals[k+1], v_vals[k+2], v_vals[k+3]],
        [v_vals[k+2], v_vals[k+3], v_vals[k+4]]
    ])

for k in range(5, 30):
    d3k = mp.det(casorati3(k))
    d3k1 = mp.det(casorati3(k+1))
    if d3k != 0:
        ratio = d3k1 / d3k
        print(f"  k={k:2d}: C3(k+1)/C3(k) = {nstr(ratio, 20)}, log10|C3| = {float(log(abs(d3k), 10)):.1f}")

# === Direct approach: compute the ORDER-3 recurrence for q_N ===
# q_N satisfies: α₃(N) q_{N+3} + α₂(N) q_{N+2} + α₁(N) q_{N+1} + α₀(N) q_N = 0
# The coefficients α_i(N) come from the matrix M(N).

# Method: for each N, we have the relation
# r_{N+1} = r_N · M(N)
# where r_N is a 1×3 row vector. The scalar recurrence for the first component
# q_N = r_{N,1} can be derived by computing:
#
# | q_N    q_{N+1}  q_{N+2}  q_{N+3} |   | 1    M11(N)  P11(N,N+1)  P11(N,N+1,N+2) |
# | q'_N   q'_{N+1} q'_{N+2} q'_{N+3}| = | ...  ...     ...          ...             |
# | q''_N  ...      ...      ...      |   | ...  ...     ...          ...             |
#
# The coefficient of q_{N+3} is det(Casorati_3(N)), etc.

# Simpler: compute the recurrence numerically.
# For order 3: we need 4 coefficients α₀,...,α₃ as polynomials in N.
# Key insight: normalize by dividing by q_N to make the system better conditioned.

print("\n\n=== DIRECT order-3 recurrence for q_N ===")
print("Testing: α₃ q_{N+3} + α₂ q_{N+2} + α₁ q_{N+1} + α₀ q_N = 0")
print("Compute: if we fix α₃=1, solve for α₀,α₁,α₂ as RATIONAL FUNCTIONS of N")
print()

# For each N (with N+3 < N_max+1), solve:
# q_{N+3} + α₂(N) q_{N+2} + α₁(N) q_{N+1} + α₀(N) q_N = 0
# => α₀ = -(q_{N+3} + α₂ q_{N+2} + α₁ q_{N+1}) / q_N
# But we have 3 unknowns and 1 equation per N. Need the POLYNOMIAL structure.

# Alternative: use the 3×3 matrix to directly compute the recurrence.
# The Cayley-Hamilton-like approach:
# Let c_N = column vector that evolves as c_{N+1} = M(N)^T c_N
# Then det([c_N, c_{N+1}, c_{N+2}]) evolves predictably.

# Actually, let me just compute the 3x3 Casorati determinant for q_N.
print("Order-3 Casorati for q_N:")
for k in range(3, 20):
    D = matrix([
        [q_vals[k], q_vals[k+1], q_vals[k+2]],
        [q_vals[k+1], q_vals[k+2], q_vals[k+3]],
        [q_vals[k+2], q_vals[k+3], q_vals[k+4]]
    ])
    det_D = mp.det(D)
    if k > 3:
        D_prev = matrix([
            [q_vals[k-1], q_vals[k], q_vals[k+1]],
            [q_vals[k], q_vals[k+1], q_vals[k+2]],
            [q_vals[k+1], q_vals[k+2], q_vals[k+3]]
        ])
        det_prev = mp.det(D_prev)
        if det_prev != 0:
            ratio = det_D / det_prev
            print(f"  k={k:2d}: C_q3(k)/C_q3(k-1) = {nstr(ratio, 20)}")

# For order-3 recurrence of q_N, the Casorati ratio should equal
# (-1)^3 * α₀(k-1) / α₃(k-1) = -α₀(k-1)/α₃(k-1)
# which is a RATIONAL FUNCTION of k.

# Let's also check the order-4 Casorati to confirm order is exactly 3:
print("\nOrder-4 Casorati for q_N (should be ~0 if order is 3):")
for k in range(3, 15):
    D4 = matrix([
        [q_vals[k], q_vals[k+1], q_vals[k+2], q_vals[k+3]],
        [q_vals[k+1], q_vals[k+2], q_vals[k+3], q_vals[k+4]],
        [q_vals[k+2], q_vals[k+3], q_vals[k+4], q_vals[k+5]],
        [q_vals[k+3], q_vals[k+4], q_vals[k+5], q_vals[k+6]]
    ])
    det_D4 = mp.det(D4)
    # Compare to C3 scale
    D3 = matrix([
        [q_vals[k], q_vals[k+1], q_vals[k+2]],
        [q_vals[k+1], q_vals[k+2], q_vals[k+3]],
        [q_vals[k+2], q_vals[k+3], q_vals[k+4]]
    ])
    det_D3 = mp.det(D3)
    if det_D3 != 0:
        rel = det_D4 / det_D3
        print(f"  k={k:2d}: C4/C3 = {float(abs(rel)):.3e}")
