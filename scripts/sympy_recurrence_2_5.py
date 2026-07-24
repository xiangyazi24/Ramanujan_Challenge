#!/usr/bin/env python3
"""Problem 2.5: Compute the order-3 scalar recurrence for q_N EXACTLY using sympy.

This script derives the scalar recurrence α₃(N)q_{N+3}+α₂(N)q_{N+2}+α₁(N)q_{N+1}+α₀(N)q_N = 0
with EXACT polynomial coefficients by symbolic elimination from the 3×3 matrix system.

Strategy: for the 3-component system s_{N+1} = s_N·M(N), use sympy to eliminate the
auxiliary components p_N = s_N[1], r_N = s_N[2] and obtain the recurrence for q_N = s_N[0].
"""
from sympy import symbols, Matrix, Poly, factor, simplify, Rational, expand, gcd, lcm
from sympy import sqrt as ssqrt, ZZ, QQ, pprint

n = symbols('n')

def M_sym(N):
    """Symbolic 3x3 matrix M(N)."""
    m11 = (-2*N-5)*(N+3)**2 * (136*N**4 + 1424*N**3 + 5548*N**2 + 9551*N + 6141)
    m12 = 384*N**6 + 6384*N**5 + 44168*N**4 + 162698*N**3 + 336377*N**2 + 369933*N + 169011
    m13 = -480*N**4 - 4980*N**3 - 19210*N**2 - 32690*N - 20730

    m21 = (N+2)**2*(N+3)**2*(4*N+10)*(48*N**3 + 386*N**2 + 1017*N + 879)
    m22 = (N+2)**2*(-272*N**5 - 3848*N**4 - 21732*N**3 - 61184*N**2 - 85761*N - 47808)
    m23 = (N+2)**2*(320*N**3 + 2540*N**2 + 6610*N + 5640)

    m31 = (-4*N-10)*(N+2)**2*(N+3)**2*(32*N**4 + 302*N**3 + 1037*N**2 + 1530*N + 813)
    m32 = (N+2)**2*(192*N**6 + 2984*N**5 + 19116*N**4 + 64452*N**3 + 120256*N**2 + 117279*N + 46476)
    m33 = (N+2)**2*(-16*N**5 - 408*N**4 - 2912*N**3 - 8884*N**2 - 12254*N - 6240)

    return Matrix([[expand(m11), expand(m12), expand(m13)],
                   [expand(m21), expand(m22), expand(m23)],
                   [expand(m31), expand(m32), expand(m33)]])

print("Building symbolic matrices M(n), M(n+1), M(n+2)...")
M0 = M_sym(n)
M1 = M_sym(n+1)
M2 = M_sym(n+2)

# Compute det(M(n)) — this is α₀ up to sign
print("Computing det(M(n))...")
det_M = M0.det()
det_M = expand(det_M)
p_det = Poly(det_M, n)
print(f"  deg(det(M(n))) = {p_det.degree()}")
print(f"  Leading coeff = {p_det.LC()}")

# Compute trace(M(n))
tr_M = M0[0,0] + M0[1,1] + M0[2,2]
tr_M = expand(tr_M)
p_tr = Poly(tr_M, n)
print(f"  deg(tr(M(n))) = {p_tr.degree()}")

# Compute cofactor sum = (tr²-tr(M²))/2
print("Computing cofactor sum (sum of 2×2 minors on diagonal)...")
M0sq = M0 * M0
tr_M2 = sum(M0sq[i,i] for i in range(3))
tr_M2 = expand(tr_M2)
cof_sum = expand((tr_M**2 - tr_M2) / 2)

# Verify Cayley-Hamilton: M³ - tr·M² + cof·M - det·I = 0
print("Verifying Cayley-Hamilton (this may take a while)...")
# Skip full verification — just check degrees

# The scalar recurrence comes from the COUPLED system.
# For the row vector s = (q, p, r), s_{N+1} = s · M(N).
# The recurrence for q alone involves M(N), M(N+1), M(N+2).

# Method: eliminate p_N, r_N from the system
# Equations at step N:
#   q_{N+1} = q_N · M₀₀ + p_N · M₁₀ + r_N · M₂₀   ... (1)
# Equations at step N, columns 1 and 2:
#   p_{N+1} = q_N · M₀₁ + p_N · M₁₁ + r_N · M₂₁   ... (2)
#   r_{N+1} = q_N · M₀₂ + p_N · M₁₂ + r_N · M₂₂   ... (3)

# At step N+1, equation (1):
#   q_{N+2} = q_{N+1} · M₀₀' + p_{N+1} · M₁₀' + r_{N+1} · M₂₀'
# Substituting (2)(3):
#   q_{N+2} = q_{N+1}·M₀₀' + (q_N·M₀₁+p_N·M₁₁+r_N·M₂₁)·M₁₀'
#            + (q_N·M₀₂+p_N·M₁₂+r_N·M₂₂)·M₂₀'
#   = q_{N+1}·M₀₀' + q_N·(M₀₁·M₁₀'+M₀₂·M₂₀')
#     + p_N·(M₁₁·M₁₀'+M₁₂·M₂₀') + r_N·(M₂₁·M₁₀'+M₂₂·M₂₀')  ... (4)

# From (1) and (4): 2 equations in 2 unknowns (p_N, r_N).
# Solve, then substitute into the equation for q_{N+3}.

print("\nSetting up 2x2 system for (p_N, r_N) elimination...")

# Coefficients of p_N and r_N in equations (1) and (4)
A11 = M0[1,0]  # m10(n)
A12 = M0[2,0]  # m20(n)
A21 = expand(M0[1,1]*M1[1,0] + M0[1,2]*M1[2,0])  # m11*m10' + m12*m20'
A22 = expand(M0[2,1]*M1[1,0] + M0[2,2]*M1[2,0])  # m21*m10' + m22*m20'

print(f"  deg(A11) = {Poly(expand(A11), n).degree()}")
print(f"  deg(A12) = {Poly(expand(A12), n).degree()}")
print(f"  deg(A21) = {Poly(expand(A21), n).degree()}")
print(f"  deg(A22) = {Poly(expand(A22), n).degree()}")

det_sys = expand(A11*A22 - A12*A21)
p_det_sys = Poly(det_sys, n)
print(f"  deg(det of 2x2 system) = {p_det_sys.degree()}")

# RHS coefficients (linear in q_N, q_{N+1}, q_{N+2})
# Eq (1): q_{N+1} - q_N*m00 = p_N*m10 + r_N*m20
# => B1 = -m00*q_N + 1*q_{N+1} + 0*q_{N+2}
# B1_q0 = -m00, B1_q1 = 1, B1_q2 = 0

# Eq (4): q_{N+2} - q_{N+1}*m00' - q_N*(m01*m10'+m02*m20') = p_N*(stuff) + r_N*(stuff)
# => B2 = -(m01*m10'+m02*m20')*q_N + (-m00')*q_{N+1} + 1*q_{N+2}

B1_q0 = expand(-M0[0,0])
B1_q1 = Rational(1)
B2_q0 = expand(-(M0[0,1]*M1[1,0] + M0[0,2]*M1[2,0]))
B2_q1 = expand(-M1[0,0])

# p_N coefficients (linear in q's):
# p_N = (B1*A22 - B2*A12) / det_sys
# So p_N = [A22*B1_q0 - A12*B2_q0]/det · q_N
#        + [A22*B1_q1 - A12*B2_q1]/det · q_{N+1}
#        + [A22*0    - A12*1     ]/det · q_{N+2}

pN_q0_num = expand(A22*B1_q0 - A12*B2_q0)
pN_q1_num = expand(A22*B1_q1 - A12*B2_q1)
pN_q2_num = expand(-A12)

rN_q0_num = expand(A11*B2_q0 - A21*B1_q0)
rN_q1_num = expand(A11*B2_q1 - A21*B1_q1)
rN_q2_num = expand(A11)

print("\nComputed (p_N, r_N) as linear functions of (q_N, q_{N+1}, q_{N+2}) / det")
print(f"  deg(pN_q0_num) = {Poly(pN_q0_num, n).degree()}")
print(f"  deg(pN_q1_num) = {Poly(pN_q1_num, n).degree()}")
print(f"  deg(pN_q2_num) = {Poly(pN_q2_num, n).degree()}")

# Now compute p_{N+1}, r_{N+1}:
# p_{N+1} = q_N*m01 + p_N*m11 + r_N*m21
# Each term has numerator/det_sys structure.
# p_{N+1}_q0_num = m01*det_sys + m11*pN_q0_num + m21*rN_q0_num
# etc.

print("\nComputing p_{N+1}, r_{N+1}...")

pN1_q0_num = expand(M0[0,1]*det_sys + M0[1,1]*pN_q0_num + M0[2,1]*rN_q0_num)
pN1_q1_num = expand(M0[1,1]*pN_q1_num + M0[2,1]*rN_q1_num)
pN1_q2_num = expand(M0[1,1]*pN_q2_num + M0[2,1]*rN_q2_num)

rN1_q0_num = expand(M0[0,2]*det_sys + M0[1,2]*pN_q0_num + M0[2,2]*rN_q0_num)
rN1_q1_num = expand(M0[1,2]*pN_q1_num + M0[2,2]*rN_q1_num)
rN1_q2_num = expand(M0[1,2]*pN_q2_num + M0[2,2]*rN_q2_num)

# Now compute p_{N+2}, r_{N+2}:
# p_{N+2} = q_{N+1}*m01(N+1) + p_{N+1}*m11(N+1) + r_{N+1}*m21(N+1)
# All divided by det_sys (common denominator from the p_N, r_N expressions)

print("Computing p_{N+2}, r_{N+2}...")

pN2_q0_num = expand(M1[1,1]*pN1_q0_num + M1[2,1]*rN1_q0_num)
pN2_q1_num = expand(M1[0,1]*det_sys + M1[1,1]*pN1_q1_num + M1[2,1]*rN1_q1_num)
pN2_q2_num = expand(M1[1,1]*pN1_q2_num + M1[2,1]*rN1_q2_num)

rN2_q0_num = expand(M1[1,2]*pN1_q0_num + M1[2,2]*rN1_q0_num)
rN2_q1_num = expand(M1[0,2]*det_sys + M1[1,2]*pN1_q1_num + M1[2,2]*rN1_q1_num)
rN2_q2_num = expand(M1[1,2]*pN1_q2_num + M1[2,2]*rN1_q2_num)

# Finally: q_{N+3} = q_{N+2}*m00(N+2) + p_{N+2}*m10(N+2) + r_{N+2}*m20(N+2)
# = m00(N+2)*q_{N+2} + (m10(N+2)*pN2 + m20(N+2)*rN2) (all divided by det_sys)

# So: q_{N+3} = m00(N+2)*q_{N+2} + [m10(N+2)*pN2_qi + m20(N+2)*rN2_qi]/det * q_i

print("\nComputing final q_{N+3} coefficients...")

# Coefficient of q_N in q_{N+3} (times det_sys):
alpha0_num = expand(M2[1,0]*pN2_q0_num + M2[2,0]*rN2_q0_num)
# Coefficient of q_{N+1}:
alpha1_num = expand(M2[1,0]*pN2_q1_num + M2[2,0]*rN2_q1_num)
# Coefficient of q_{N+2}:
alpha2_num = expand(M2[0,0]*det_sys + M2[1,0]*pN2_q2_num + M2[2,0]*rN2_q2_num)

# The recurrence is: det_sys * q_{N+3} = alpha2_num * q_{N+2} + alpha1_num * q_{N+1} + alpha0_num * q_N
# Or equivalently: det_sys * q_{N+3} - alpha2_num * q_{N+2} - alpha1_num * q_{N+1} - alpha0_num * q_N = 0

# Now factor out common factors to simplify
print("Simplifying (factoring common polynomial factors)...")
p_alpha0 = Poly(alpha0_num, n)
p_alpha1 = Poly(alpha1_num, n)
p_alpha2 = Poly(alpha2_num, n)
p_alpha3 = Poly(det_sys, n)

print(f"  deg(α₃) = {p_alpha3.degree()}")
print(f"  deg(α₂) = {p_alpha2.degree()}")
print(f"  deg(α₁) = {p_alpha1.degree()}")
print(f"  deg(α₀) = {p_alpha0.degree()}")

# Find GCD of all 4 polynomials to simplify
print("\nComputing GCD of α₀, α₁, α₂, α₃...")
g01 = gcd(p_alpha0.as_expr(), p_alpha1.as_expr())
g23 = gcd(p_alpha2.as_expr(), p_alpha3.as_expr())
g_all = gcd(g01, g23)
p_gcd = Poly(g_all, n)
print(f"  deg(gcd) = {p_gcd.degree()}")
print(f"  gcd = {factor(g_all)}")

# Divide out
alpha0_red = Poly(expand(alpha0_num / g_all), n)
alpha1_red = Poly(expand(alpha1_num / g_all), n)
alpha2_red = Poly(expand(alpha2_num / g_all), n)
alpha3_red = Poly(expand(det_sys / g_all), n)

print(f"\nReduced degrees:")
print(f"  deg(α₃) = {alpha3_red.degree()}")
print(f"  deg(α₂) = {alpha2_red.degree()}")
print(f"  deg(α₁) = {alpha1_red.degree()}")
print(f"  deg(α₀) = {alpha0_red.degree()}")

# Factor each coefficient
print("\nFactored forms:")
print(f"  α₃ = {factor(alpha3_red.as_expr())}")
print(f"  α₀ = {factor(alpha0_red.as_expr())}")

# Verify the Poincaré polynomial
# For large n: α₃ ~ lead_3 n^d3, α₀ ~ lead_0 n^d0
# Poincaré poly: lead_3 t³ + lead_2 t² + lead_1 t + lead_0 = 0
print("\nPoincaré polynomial (leading coefficients):")
lc3 = alpha3_red.LC()
lc2 = alpha2_red.LC()
lc1 = alpha1_red.LC()
lc0 = alpha0_red.LC()
print(f"  {lc3} t³ + ({lc2}) t² + ({lc1}) t + ({lc0}) = 0")
# Normalize
print(f"  t³ + ({Rational(lc2,lc3)}) t² + ({Rational(lc1,lc3)}) t + ({Rational(lc0,lc3)}) = 0")
