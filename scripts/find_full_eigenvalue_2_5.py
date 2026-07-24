#!/usr/bin/env python3
"""Find the FULL rational eigenvalue lambda_3(n) of M(n) for Problem 2.5.

lambda_3(n) is a degree-7 polynomial in n satisfying:
  lambda_3^3 + T(n)*lambda_3^2 + Q(n)*lambda_3 + D(n) = 0
where T, Q, D are the char poly coefficients.
"""
from sympy import symbols, Matrix, expand, Poly, solve, Rational

n = symbols('n')

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

M = M_sym(n)
lam = symbols('lambda')
cp = M.charpoly(lam).as_expr()
cp_poly = Poly(cp, lam)
coeffs = cp_poly.all_coeffs()  # [1, T(n), Q(n), D(n)]
T_n = expand(coeffs[1])
Q_n = expand(coeffs[2])
D_n = expand(coeffs[3])

# Ansatz: lambda_3(n) = sum_{k=0}^{7} a_k n^k
a0, a1, a2, a3, a4, a5, a6, a7 = symbols('a0 a1 a2 a3 a4 a5 a6 a7')
lam3 = a7*n**7 + a6*n**6 + a5*n**5 + a4*n**4 + a3*n**3 + a2*n**2 + a1*n + a0

# Equation: lam3^3 + T*lam3^2 + Q*lam3 + D = 0
residual = expand(lam3**3 + T_n*lam3**2 + Q_n*lam3 + D_n)
res_poly = Poly(residual, n)

# Extract coefficients of n^21 down to n^14 (8 equations for 8 unknowns)
# Actually, the residual has degree 21. We need all coefficients to be 0.
# With a7=-16 (known), we have 7 unknowns: a6,...,a0.
# Set a7 = -16 and solve iteratively from highest degree.

subs_dict = {a7: -16}
residual = residual.subs(a7, -16)
res_poly = Poly(expand(residual), n)
print(f"After a7=-16: residual degree = {res_poly.degree()}")

# Solve coefficient by coefficient, from highest to lowest
for deg in range(res_poly.degree(), 13, -1):
    coeff = res_poly.nth(deg)
    coeff_simplified = expand(coeff.subs(subs_dict))
    # Find which variable to solve for
    unknowns = [a6, a5, a4, a3, a2, a1, a0]
    for u in unknowns:
        if u not in subs_dict and coeff_simplified.has(u):
            sol = solve(coeff_simplified, u)
            if sol:
                subs_dict[u] = sol[0]
                print(f"  n^{deg}: {u} = {sol[0]}")
                # Apply to residual
                residual = residual.subs(u, sol[0])
                res_poly = Poly(expand(residual), n)
                break
    else:
        # All unknowns solved, check that this coeff is zero
        val = expand(coeff_simplified)
        print(f"  n^{deg}: consistency check = {val}")

# Print the full eigenvalue
print(f"\n=== Rational eigenvalue lambda_3(n) ===")
lam3_full = -16*n**7
for sym, name in [(a6,'a6'),(a5,'a5'),(a4,'a4'),(a3,'a3'),(a2,'a2'),(a1,'a1'),(a0,'a0')]:
    if sym in subs_dict:
        lam3_full += subs_dict[sym] * n**(int(name[1]))
lam3_full = expand(lam3_full)
print(f"  lambda_3(n) = {lam3_full}")

from sympy import factor as sym_factor
print(f"\n  Factored: {sym_factor(lam3_full)}")

# Verify: check that char poly evaluates to 0 at lambda = lambda_3(n)
print(f"\n=== Verification ===")
check = expand(lam3_full**3 + T_n*lam3_full**2 + Q_n*lam3_full + D_n)
check_poly = Poly(check, n) if check != 0 else None
if check == 0:
    print("  lambda_3^3 + T*lambda_3^2 + Q*lambda_3 + D = 0  ✓")
else:
    print(f"  Residual degree: {check_poly.degree() if check_poly else 0}")
    print(f"  Residual: {check}")

# Now compute the quotient: char_poly / (lambda - lambda_3)
print(f"\n=== Factoring char poly ===")
from sympy import quo, rem
# char_poly in lambda = lambda^3 + T*lambda^2 + Q*lambda + D
# Divide by (lambda - lambda_3)
# Result: lambda^2 + B(n)*lambda + C(n)
# where B = T + lambda_3, C = Q + lambda_3 * (T + lambda_3)

B_n = expand(T_n + lam3_full)
C_n = expand(Q_n + lam3_full * B_n)
print(f"  Quotient: lambda^2 + B(n)*lambda + C(n)")
print(f"  B(n) = {B_n}")
print(f"  C(n) = {C_n}")
pB = Poly(B_n, n)
pC = Poly(C_n, n)
print(f"  deg(B) = {pB.degree()}, LC(B) = {pB.LC()}")
print(f"  deg(C) = {pC.degree()}, LC(C) = {pC.LC()}")

# The Poincaré poly of the quadratic: t^2 + 544t + 256 = 0
# LC(B) should be 544, LC(C) should be 256
print(f"\n  Expected: LC(B) = 544, LC(C) = 256")
print(f"  Match B: {pB.LC() == 544}")
print(f"  Match C: {pC.LC() == 256}")

# Factored forms
print(f"\n  B(n) factored: {sym_factor(B_n)}")
print(f"  C(n) factored: {sym_factor(C_n)}")

# The gauge h(n) = product_{k=0}^{n-1} lambda_3(k)
print(f"\n=== Gauge h(n) = product lambda_3(k) ===")
print(f"  lambda_3(n) = {sym_factor(lam3_full)}")
print(f"  h(n+1)/h(n) = lambda_3(n)")
print(f"  First values:")
for k in range(6):
    val = lam3_full.subs(n, k)
    print(f"    lambda_3({k}) = {val}")
