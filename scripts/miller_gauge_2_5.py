#!/usr/bin/env python3
"""Miller's algorithm: extract the c=-16 gauge solution by backward recursion.

The c=-16 Poincaré root is the SMALLEST absolute root (16 vs 0.47 and 543.5).
Forward recursion amplifies the dominant c=543.5 mode.
Backward recursion amplifies the c=-16 mode (making it the "dominant" one backward).

After extracting h(n), compute r(n) = h(n+1)/h(n) and identify it as a rational function.
"""
from mpmath import mp, mpf, nstr, matrix
mp.dps = 300  # High precision to handle cancellations

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

# First compute the scalar recurrence coefficients numerically.
# From the MATRIX, the scalar sequence q_N = (A·T_N)_{1,0} satisfies order-3.
# But I need the RECURRENCE COEFFICIENTS c_i(N) to run backward.

# From derive_recurrence_2_5.py, the recurrence is:
# q_{N+3} = c_2(N) q_{N+2} + c_1(N) q_{N+1} + c_0(N) q_N
# (derived by eliminating p,r from the matrix relation)

# Let me recompute c_i(N) for each N directly from the matrix.
def get_recurrence_coeffs(N):
    """Compute c_0(N), c_1(N), c_2(N) such that 
    q_{N+3} = c_2(N) q_{N+2} + c_1(N) q_{N+1} + c_0(N) q_N.
    Uses Cramer's rule on the 3-step matrix relation."""
    M_N = M_mat(N)
    M_N1 = M_mat(N+1)
    M_N2 = M_mat(N+2)
    
    m00 = M_N[0,0]; m01 = M_N[0,1]; m02 = M_N[0,2]
    m10 = M_N[1,0]; m11 = M_N[1,1]; m12 = M_N[1,2]
    m20 = M_N[2,0]; m21 = M_N[2,1]; m22 = M_N[2,2]
    m00p = M_N1[0,0]; m10p = M_N1[1,0]; m20p = M_N1[2,0]
    m01p = M_N1[0,1]; m11p = M_N1[1,1]; m21p = M_N1[2,1]
    m02p = M_N1[0,2]; m12p = M_N1[1,2]; m22p = M_N1[2,2]
    
    A11 = m10; A12 = m20
    A21 = m11*m10p + m12*m20p
    A22 = m21*m10p + m22*m20p
    det_A = A11*A22 - A12*A21
    
    pN_q0 = (-m00*A22 - (-m01*m10p-m02*m20p)*A12) / det_A
    pN_q1 = (A22 - (-m00p)*A12) / det_A
    pN_q2 = (-A12) / det_A
    rN_q0 = (A11*(-m01*m10p-m02*m20p) - A21*(-m00)) / det_A
    rN_q1 = (A11*(-m00p) - A21) / det_A
    rN_q2 = A11 / det_A
    
    pN1_q0 = m01 + pN_q0*m11 + rN_q0*m21
    pN1_q1 = pN_q1*m11 + rN_q1*m21
    pN1_q2 = pN_q2*m11 + rN_q2*m21
    rN1_q0 = m02 + pN_q0*m12 + rN_q0*m22
    rN1_q1 = pN_q1*m12 + rN_q1*m22
    rN1_q2 = pN_q2*m12 + rN_q2*m22
    
    pN2_q0 = pN1_q0*m11p + rN1_q0*m21p
    pN2_q1 = m01p + pN1_q1*m11p + rN1_q1*m21p
    pN2_q2 = pN1_q2*m11p + rN1_q2*m21p
    rN2_q0 = pN1_q0*m12p + rN1_q0*m22p
    rN2_q1 = m02p + pN1_q1*m12p + rN1_q1*m22p
    rN2_q2 = pN1_q2*m12p + rN1_q2*m22p
    
    m10pp = M_N2[1,0]; m20pp = M_N2[2,0]; m00pp = M_N2[0,0]
    
    c0 = m10pp*pN2_q0 + m20pp*rN2_q0
    c1 = m10pp*pN2_q1 + m20pp*rN2_q1
    c2 = m00pp + m10pp*pN2_q2 + m20pp*rN2_q2
    
    return c0, c1, c2

# Backward recursion from N_max:
# q_N = (q_{N+3} - c_2(N) q_{N+2} - c_1(N) q_{N+1}) / c_0(N)
N_max = 80
N_start = 3  # Start collecting from here

# Compute coefficients
coeffs = {}
for N in range(N_start, N_max + 3):
    coeffs[N] = get_recurrence_coeffs(N)

# Backward recursion with arbitrary final conditions
h = [mpf(0)] * (N_max + 4)
h[N_max + 2] = mpf(1)
h[N_max + 1] = mpf(0)
h[N_max] = mpf(0)

for N in range(N_max - 1, N_start - 1, -1):
    c0, c1, c2 = coeffs[N]
    # q_{N+3} = c2 q_{N+2} + c1 q_{N+1} + c0 q_{N}
    # → q_{N} = (q_{N+3} - c2 q_{N+2} - c1 q_{N+1}) / c0
    h[N] = (h[N+3] - c2 * h[N+2] - c1 * h[N+1]) / c0

# Compute ratios r(n) = h(n+1)/h(n)
print("=== Ratios r(n) = h(n+1)/h(n) from backward recursion ===")
ratios = []
for N in range(N_start, min(N_start + 25, N_max)):
    if h[N] != 0:
        r = h[N+1] / h[N]
        ratios.append((N, r))
        # Also compute r(n) / (-16 * n^7) — should approach 1
        normalized = r / (-16 * mpf(N)**7)
        print(f"  r({N:2d}) = {nstr(r, 15)},  r/(-16n^7) = {nstr(normalized, 12)}")

# Try to identify r(n) as a rational function
# r(n) = -16 * P(n)/Q(n) where P, Q are degree-7 polynomials with integer coefficients
# Factor: r(n) should be a product of linear factors with half-integer or integer roots.
# Try to identify the zeros and poles of r(n).
print("\n=== Testing if r(n) has specific simple zeros ===")
# Test r(n) at n = -k/2 for various k (half-integer)
for N in range(N_start + 5, N_start + 15):
    if h[N] != 0:
        r = h[N+1] / h[N]
        r_over_n7 = r / (mpf(N)**7)
        # r(n) = -16 * prod(n+ai) / prod(n+bj)
        # If n+ai = 0 then r(n) = 0, i.e., h(n+1) = 0.
        # So zeros of h are at n = -ai.
        # Check: where does h pass through zero?

print("\n=== Looking for zeros of h(n) ===")
for N in range(N_start, N_start + 20):
    if abs(h[N]) < 1e-10:
        print(f"  h({N}) ≈ 0")
    elif N > N_start and h[N] * h[N-1] < 0:
        print(f"  Sign change between h({N-1}) and h({N})")

# More useful: compute r(n)/(prod of specific factors) and check for constancy
# Try: r(n) = -16 * (2n+3)(2n+5)^2(2n+7)^3(2n+9)^2(n+1)(n+2)^2 / [(2n+5)(2n+7)^2(2n+9)^3(2n+11)^2(n+3)(n+4)^2]
# This is a guess based on the half-integer roots in c_0 and c_3.

# Actually let me try to identify r(n) by computing r(n) * n^{-7} for several n and using LLL
# to find the polynomial p(n) = r(n) / (-16).
print("\n=== r(n)/(-16) values for LLL identification ===")
for N in range(5, 20):
    if h[N] != 0:
        r = h[N+1] / h[N]
        p = r / mpf(-16)
        # p(N) should be a ratio P(N)/Q(N) of degree-7 polynomials
        # P(N) = prod (N + ai), Q(N) = prod (N + bj)
        print(f"  n={N}: r(n)/(-16) = {nstr(p, 20)}")
