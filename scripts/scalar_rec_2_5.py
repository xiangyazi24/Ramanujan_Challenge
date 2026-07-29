#!/usr/bin/env python3
"""Extract the scalar recurrence for Problem 2.5 CMF denominators.

The 3×3 CMF M(n) with initial matrix A gives denominators Q_{N,j}.
These satisfy a third-order scalar recurrence with polynomial coefficients.
We extract the recurrence by computing many terms and using the
ansatz-fitting method.
"""
from fractions import Fraction

def M(n):
    m11 = (-2*n-5)*(n+3)**2*(136*n**4+1424*n**3+5548*n**2+9551*n+6141)
    m12 = 384*n**6+6384*n**5+44168*n**4+162698*n**3+336377*n**2+369933*n+169011
    m13 = -480*n**4-4980*n**3-19210*n**2-32690*n-20730
    m21 = (n+2)**2*(n+3)**2*(4*n+10)*(48*n**3+386*n**2+1017*n+879)
    m22 = (n+2)**2*(-272*n**5-3848*n**4-21732*n**3-61184*n**2-85761*n-47808)
    m23 = (n+2)**2*(320*n**3+2540*n**2+6610*n+5640)
    m31 = (-4*n-10)*(n+2)**2*(n+3)**2*(32*n**4+302*n**3+1037*n**2+1530*n+813)
    m32 = (n+2)**2*(192*n**6+2984*n**5+19116*n**4+64452*n**3+120256*n**2+117279*n+46476)
    m33 = (n+2)**2*(-16*n**5-408*n**4-2912*n**3-8884*n**2-12254*n-6240)
    return [[m11,m12,m13],[m21,m22,m23],[m31,m32,m33]]

def matmul(A, B):
    rows, cols = len(A), len(B[0])
    inner = len(B)
    return [[sum(A[i][k]*B[k][j] for k in range(inner)) for j in range(cols)] for i in range(rows)]

A_init = [[30921, -32972, 8240],
          [33750, -36000, 9000]]

# Compute Q_{N,1} = (A · M_0 · ... · M_{N-1})[1][0]
N_max = 40
prod = [[1,0,0],[0,1,0],[0,0,1]]
Q = []
for n in range(N_max):
    AM = matmul(A_init, prod)
    Q.append(AM[1][0])
    prod = matmul(prod, M(n))
AM = matmul(A_init, prod)
Q.append(AM[1][0])

print(f"Computed {len(Q)} terms of Q_N")
print(f"Q_0 = {Q[0]}")
print(f"Q_1 = {Q[1]}")
print(f"Q_2 = {Q[2]}")

# The scalar recurrence should be: a_3(n)*Q_{n+3} + a_2(n)*Q_{n+2} + a_1(n)*Q_{n+1} + a_0(n)*Q_n = 0
# where a_i(n) are polynomials.
# The Poincaré polynomial tells us the degrees are all the same (degree d).
# From the 3×3 matrix with degree-7 entries, the scalar recurrence has
# degree pattern (28, 21, 14, 7) -- but let's just find the ratios first.

# Compute the ratios a_0/a_3, a_1/a_3, a_2/a_3 for each n
print("\n--- Extracting recurrence ratios ---")
# For a 3-term recurrence: Q_{n+3} = r_2(n)*Q_{n+2} + r_1(n)*Q_{n+1} + r_0(n)*Q_n
# This holds for n >= some n_0.
# From the recurrence: r_2 = -a_2/a_3, r_1 = -a_1/a_3, r_0 = -a_0/a_3

for n in range(min(15, len(Q)-3)):
    if Q[n] != 0 and Q[n+1] != 0 and Q[n+2] != 0 and Q[n+3] != 0:
        # Check if a simple rational relation holds
        # Q_{n+3} = r2*Q_{n+2} + r1*Q_{n+1} + r0*Q_n
        # This is underdetermined with 3 unknowns.
        # Use n, n+1, n+2 to get 3 equations for r2(n), r1(n), r0(n)
        pass

# Better approach: use the known det M(n) to get the product of Poincaré roots
# det M(n) = -8*(n+1)*(n+2)^6*(n+3)^5*(2n+3)^2*(2n+5)^3*(2n+7)^4
# The product of eigenvalues at step n is det M(n).

# The three eigenvalue sequences λ_1(n), λ_2(n), λ_3(n) of M(n) satisfy:
# tr M(n) = λ_1 + λ_2 + λ_3
# etc.

# Actually, the cleanest approach: compute the recurrence coefficients
# using the Casorati determinant method.
# The denominator Q_N and its shifts satisfy a 3rd order recurrence.
# I need 4 coefficients a_0(n), ..., a_3(n).
# For polynomial coefficients of degree d, I need enough data points.

# Let's try: assume a_i are polynomials of degree D.
# Then we need 4*(D+1) unknowns, and each value of n gives 1 equation.
# With D=28 (degree pattern (28,21,14,7)), we need 4*29 = 116 equations.
# But we only have 40 terms, which gives 37 equations.
# So this approach needs more terms.

# Alternative: use the companion matrix directly.
# The eigenvalues of M(n) are the "step-n" characteristic roots.
# The tr, sum-of-2x2-minors, and det give the symmetric functions.

print("\n--- Trace, cofactor sum, det of M(n) ---")
for n in range(5):
    Mn = M(n)
    tr = Mn[0][0] + Mn[1][1] + Mn[2][2]
    # 2x2 minors sum = M00*M11-M01*M10 + M00*M22-M02*M20 + M11*M22-M12*M21
    s2 = (Mn[0][0]*Mn[1][1] - Mn[0][1]*Mn[1][0] +
          Mn[0][0]*Mn[2][2] - Mn[0][2]*Mn[2][0] +
          Mn[1][1]*Mn[2][2] - Mn[1][2]*Mn[2][1])
    det_val = (Mn[0][0]*(Mn[1][1]*Mn[2][2]-Mn[1][2]*Mn[2][1])
              -Mn[0][1]*(Mn[1][0]*Mn[2][2]-Mn[1][2]*Mn[2][0])
              +Mn[0][2]*(Mn[1][0]*Mn[2][1]-Mn[1][1]*Mn[2][0]))
    det_formula = -8*(n+1)*(n+2)**6*(n+3)**5*(2*n+3)**2*(2*n+5)**3*(2*n+7)**4
    print(f"n={n}: tr={tr}, s2={s2}, det={det_val}")
    print(f"      det formula = {det_formula}, match = {det_val == det_formula}")

# Growth rate analysis
print("\n--- Growth rates Q_{n+1}/Q_n ---")
for n in range(min(20, len(Q)-1)):
    if Q[n] != 0:
        ratio = Q[n+1] / Q[n]
        print(f"  Q_{n+1}/Q_{n} = {ratio:.6e}")
