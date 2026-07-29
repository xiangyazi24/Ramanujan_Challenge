#!/usr/bin/env python3
"""Find the rational eigenvalue of M(n) for Problem 2.5.

The Poincaré polynomial (c+16)(c²+544c+256) has rational root c=-16.
This means M(n) has an eigenvalue lambda_3(n) that is a rational function
of n, with lambda_3(n) ~ -16 n^7 for large n.

If the char poly of M(n) factors over Q(n) as
  (lambda - lambda_3(n))(quadratic) = 0,
then lambda_3(n) is the gauge: h(n) = product_{k=0}^{n-1} lambda_3(k).
"""
from sympy import symbols, Matrix, expand, factor, Poly, cancel, div, gcd, Rational
from sympy import roots as symroots, solve

n, lam = symbols('n lambda')

def M_sym(N):
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

print("Computing characteristic polynomial of M(n)...")
M = M_sym(n)
cp = M.charpoly(lam)
cp_expr = cp.as_expr()

# The char poly is lambda^3 + p(n)*lambda^2 + q(n)*lambda + r(n)
# where p = -tr, r = -det.
cp_poly = Poly(cp_expr, lam)
coeffs = cp_poly.all_coeffs()  # [1, p(n), q(n), r(n)]
print(f"Char poly coefficients as polys in n:")
for i, c in enumerate(coeffs):
    p = Poly(expand(c), n)
    print(f"  lambda^{3-i}: degree {p.degree()}, LC = {p.LC()}")

print(f"\nLC pattern: {[Poly(expand(c), n).LC() for c in coeffs]}")

# If lambda_3(n) is a rational root, then the char poly evaluated at lambda=lambda_3(n)
# is zero. Try lambda_3(n) = -16n^7 + an^6 + bn^5 + ... (polynomial ansatz)

# First try: is -16n^7 an eigenvalue? (probably not exactly, but let's check correction terms)
# Substitute lambda = -16n^7 into the char poly and see what terms survive.

print("\n=== Checking if lambda_3(n) is a polynomial in n ===")
# Try lambda_3 = -16*n^7 + a6*n^6 + ... (degree 7 polynomial)
# The char poly: lam^3 + p(n)*lam^2 + q(n)*lam + r(n) = 0

# Coefficients:
p_n = expand(coeffs[1])  # degree 7, LC = -560 (= -tr)
q_n = expand(coeffs[2])  # degree 14, LC = 8960 (sum of 2x2 minors)
r_n = expand(coeffs[3])  # degree 21, LC = -4096 (= -det)

# Check: substitute lam = -16*n^7 + a*n^6 into cp and set degree-21 terms to 0
a = symbols('a')
lam_trial = -16*n**7 + a*n**6

residual = expand(lam_trial**3 + p_n * lam_trial**2 + q_n * lam_trial + r_n)

# Get the degree-21 coefficient of the residual
res_poly = Poly(residual, n)
print(f"Residual degree: {res_poly.degree()}")
# The degree-21 coefficient should be zero (leading balance)
# For lam = -16n^7 + an^6:
# lam^3 = (-16)^3 n^21 + 3(-16)^2(a)n^20 + ...
# p(n)*lam^2 = (-560n^7+...)(256n^14+...) → top: -560·256 n^21 + ...
# q(n)*lam = (8960n^14+...)(-16n^7+...) → top: 8960·(-16) n^21 + ...
# r(n) → 4096 n^21 + ...

# Degree 21 terms:
d21_lam3 = (-16)**3  # = -4096
d21_plam2 = (-560) * (-16)**2  # = -560 * 256 = -143360
d21_qlam = 8960 * (-16)  # = -143360
d21_r = -4096  # -det LC actually... wait

# Actually, r(n) = -det(M(n)). det has LC = -4096 (from earlier).
# So r(n) has LC = -(-4096) = 4096? No — the char poly has r(n) as the constant term.
# For a 3x3 matrix: char poly = lam^3 - tr(M) lam^2 + (cofactor sum) lam - det(M)
# So r(n) = -det(M(n)). det LC = -4096, so r LC = -(-4096) = 4096.

# Hmm wait, the charpoly function in sympy returns: det(lambda I - M) = lambda^3 - tr lambda^2 + ... - det
# So: r(n) = -det(M(n)).

# Degree 21 balance:
# (-16)^3 + (-tr_LC)(-16)^2 + (cofactor_LC)(-16) + (-det_LC) = ?
# = -4096 + 560·256 + 8960·(-16) + 4096
# = -4096 + 143360 - 143360 + 4096
# = 0! ✓

print(f"\n  Degree-21 balance: (-16)^3 + 560*256 + 8960*(-16) + 4096 = {(-16)**3 + 560*256 + 8960*(-16) + 4096}")

# Good! The degree-21 terms cancel, so lam = -16n^7 + lower is consistent.
# Now find the coefficient of n^20.

# For lam = -16n^7 + a*n^6:
# lam^3 at n^20: 3*(-16)^2*a = 768a
# p(n)*lam^2 at n^20: p has degree 7, lam^2 has degree 14.
#   p at n^7: LC = -560. p at n^6: next coeff.
#   lam^2 at n^14: (-16)^2 = 256. lam^2 at n^13: 2*(-16)*a = -32a.
#   (p*lam^2) at n^20: (-560)(n^7)(256 n^14) → n^21 (already counted)
#   At n^20: either p_{n^7}*lam^2_{n^13} or p_{n^6}*lam^2_{n^14}
#   = (-560)(-32a) + (next_coeff_of_p)(256)

# This is getting complex. Let me just compute it symbolically.

# Extract the coefficient of n^20 from the residual
if res_poly.degree() >= 20:
    coeff_20 = res_poly.nth(20)
    print(f"  Coeff of n^20 in residual: {expand(coeff_20)}")
    # Solve for a
    sol_a = solve(coeff_20, a)
    print(f"  Solving for a: a = {sol_a}")
    if sol_a:
        a_val = sol_a[0]
        print(f"  So lambda_3(n) = -16n^7 + {a_val}n^6 + ...")

        # Continue: find next coefficient
        b = symbols('b')
        lam_trial2 = -16*n**7 + a_val*n**6 + b*n**5
        residual2 = expand(lam_trial2**3 + p_n * lam_trial2**2 + q_n * lam_trial2 + r_n)
        res_poly2 = Poly(residual2, n)
        if res_poly2.degree() >= 19:
            coeff_19 = res_poly2.nth(19)
            sol_b = solve(coeff_19, b)
            print(f"  Coeff of n^19: solving gives b = {sol_b}")
            if sol_b:
                b_val = sol_b[0]
                print(f"  lambda_3(n) = -16n^7 + {a_val}n^6 + {b_val}n^5 + ...")
