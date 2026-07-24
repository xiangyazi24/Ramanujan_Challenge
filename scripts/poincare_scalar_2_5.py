#!/usr/bin/env python3
"""Compute the Poincaré polynomial for the order-3 scalar recurrence.

With degree pattern (7, 14, 21, 28) and substitution q_n ~ prod mu(k),
mu(k) ~ c*k^7, the characteristic equation is:
  LC(alpha3) * c^3 + LC(alpha2) * c^2 + LC(alpha1) * c + LC(alpha0) = 0

The roots c1,c2,c3 should be the eigenvalues of C_inf:
  -272-192*sqrt(2), -272+192*sqrt(2), -16
"""
from sympy import symbols, Matrix, Poly, factor, expand, gcd, cancel, sqrt, Rational, solve

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

# Rebuild the order-3 recurrence to get exact LCs
print("Building symbolic matrices...")
M0 = M_sym(n)
M1 = M_sym(n+1)
M2 = M_sym(n+2)

print("Deriving order-3 recurrence coefficients...")
A11 = M0[1,0]; A12 = M0[2,0]
A21 = expand(M0[1,1]*M1[1,0] + M0[1,2]*M1[2,0])
A22 = expand(M0[2,1]*M1[1,0] + M0[2,2]*M1[2,0])
det_sys = expand(A11*A22 - A12*A21)

B1_q0 = expand(-M0[0,0]); B1_q1 = Rational(1)
B2_q0 = expand(-(M0[0,1]*M1[1,0] + M0[0,2]*M1[2,0]))
B2_q1 = expand(-M1[0,0])

pN_q0_num = expand(A22*B1_q0 - A12*B2_q0)
pN_q1_num = expand(A22*B1_q1 - A12*B2_q1)
pN_q2_num = expand(-A12)
rN_q0_num = expand(A11*B2_q0 - A21*B1_q0)
rN_q1_num = expand(A11*B2_q1 - A21*B1_q1)
rN_q2_num = expand(A11)

pN1_q0_num = expand(M0[0,1]*det_sys + M0[1,1]*pN_q0_num + M0[2,1]*rN_q0_num)
pN1_q1_num = expand(M0[1,1]*pN_q1_num + M0[2,1]*rN_q1_num)
pN1_q2_num = expand(M0[1,1]*pN_q2_num + M0[2,1]*rN_q2_num)
rN1_q0_num = expand(M0[0,2]*det_sys + M0[1,2]*pN_q0_num + M0[2,2]*rN_q0_num)
rN1_q1_num = expand(M0[1,2]*pN_q1_num + M0[2,2]*rN_q1_num)
rN1_q2_num = expand(M0[1,2]*pN_q2_num + M0[2,2]*rN_q2_num)

m01p = M1[0,1]; m11p = M1[1,1]; m21p = M1[2,1]
m02p = M1[0,2]; m12p = M1[1,2]; m22p = M1[2,2]
pN2_q0_num = expand(m11p*pN1_q0_num + m21p*rN1_q0_num)
pN2_q1_num = expand(m01p*det_sys + m11p*pN1_q1_num + m21p*rN1_q1_num)
pN2_q2_num = expand(m11p*pN1_q2_num + m21p*rN1_q2_num)
rN2_q0_num = expand(m12p*pN1_q0_num + m22p*rN1_q0_num)
rN2_q1_num = expand(m02p*det_sys + m12p*pN1_q1_num + m22p*rN1_q1_num)
rN2_q2_num = expand(m12p*pN1_q2_num + m22p*rN1_q2_num)

alpha0_raw = expand(M2[1,0]*pN2_q0_num + M2[2,0]*rN2_q0_num)
alpha1_raw = expand(M2[1,0]*pN2_q1_num + M2[2,0]*rN2_q1_num)
alpha2_raw = expand(M2[0,0]*det_sys + M2[1,0]*pN2_q2_num + M2[2,0]*rN2_q2_num)
alpha3_raw = det_sys

print("Computing GCD...")
g01 = gcd(alpha0_raw, alpha1_raw)
g23 = gcd(alpha2_raw, alpha3_raw)
g_all = gcd(g01, g23)

a0 = cancel(alpha0_raw / g_all)
a1 = cancel(alpha1_raw / g_all)
a2 = cancel(alpha2_raw / g_all)
a3 = cancel(alpha3_raw / g_all)

p0 = Poly(a0, n); p1 = Poly(a1, n); p2 = Poly(a2, n); p3 = Poly(a3, n)
print(f"Degrees: a3={p3.degree()}, a2={p2.degree()}, a1={p1.degree()}, a0={p0.degree()}")

lc3 = p3.LC(); lc2 = p2.LC(); lc1 = p1.LC(); lc0 = p0.LC()
print(f"\nLeading coefficients:")
print(f"  LC(a3) = {lc3}")
print(f"  LC(a2) = {lc2}")
print(f"  LC(a1) = {lc1}")
print(f"  LC(a0) = {lc0}")

# Characteristic equation: lc3 * c^3 + lc2 * c^2 + lc1 * c + lc0 = 0
c = symbols('c')
char_poly = lc3*c**3 + lc2*c**2 + lc1*c + lc0
print(f"\nCharacteristic polynomial: {char_poly} = 0")
print(f"Factored: {factor(char_poly)}")

roots = solve(char_poly, c)
print(f"\nRoots:")
for r in roots:
    print(f"  c = {r} = {float(r.evalf()):.6f}")

# Expected roots: -272-192*sqrt(2), -272+192*sqrt(2), -16
expected = [-272-192*sqrt(2), -272+192*sqrt(2), -16]
print(f"\nExpected roots from matrix eigenvalues:")
for e in expected:
    print(f"  {e} = {float(e.evalf()):.6f}")

# Check: does the char poly vanish at each expected root?
print(f"\nChar poly at expected roots:")
for e in expected:
    val = char_poly.subs(c, e)
    val = expand(val)
    print(f"  P({e}) = {val}")

# === Now the key: Petkovsek test for the eigenvalue-1 mode (c = -16) ===
# The hypergeometric solution has h(n+1)/h(n) ~ -16 n^7 + lower order.
# Check if S - (-16 n^7) (approximately) is a right factor.
# More precisely: define r(n) as the exact eigenvalue of M(n) in the "third mode".
# The Petkovsek algorithm searches for r(n) as a rational function.

# For now, let's check the second-order leading term.
# The sub-leading correction: q_n ~ prod_{k} (c * k^7 + d * k^6 + ...)
# After removing the leading c*k^7, what's the sub-leading?

# Actually, let me just try: gauge q_n by h(n) where
# h(n+1)/h(n) = -16 * ((n+2)*(n+3))^2 * (2*n+5)*(2*n+7)^2 / something
# (reading off from det(M) / (product of other eigenvalues))

# det(M(n)) = product of eigenvalues = -8*(n+1)*(n+2)^6*(n+3)^5*(2n+3)^2*(2n+5)^3*(2n+7)^4
# The eigenvalue product ~ (-16)^3 * n^21 for the "three identical" scaling.
# But eigenvalue 3 = det(M) / (eigenvalue 1 * eigenvalue 2)

# For the matrix eigenvalue problem: if v is a right eigenvector of M with eigenvalue lambda,
# then M*v = lambda*v. The scalar recurrence q_n = e_1^T * prod M(k) * initial has
# three modes corresponding to the three eigenvalues.

# The third eigenvalue lambda_3(n) satisfies:
# lambda_1(n) * lambda_2(n) * lambda_3(n) = det M(n)
# lambda_1(n) + lambda_2(n) + lambda_3(n) = trace M(n)
# etc.

# For large n:
# lambda_1 ~ (-272-192*sqrt(2)) * n^7
# lambda_2 ~ (-272+192*sqrt(2)) * n^7
# lambda_3 ~ -16 * n^7
# lambda_1 * lambda_2 = (272^2 - 192^2*2) n^14 = (73984-73728) n^14 = 256 n^14
# lambda_1 * lambda_2 * lambda_3 = 256 * (-16) * n^21 = -4096 n^21
# Check: det LC = -4096. ✓

# So lambda_3(n) = det(M(n)) / (lambda_1(n) * lambda_2(n))
# For the gauge: h(n+1)/h(n) = lambda_3(n)
# But lambda_3(n) is NOT rational — it involves sqrt(2).

# Wait! The individual eigenvalues involve sqrt(2), but the char poly is over Q.
# So lambda_3(n) is an algebraic function of n, not rational.
# That means there's NO hypergeometric right factor with rational r(n)!

# Unless... lambda_3(n) happens to be rational despite the char poly having irrational roots for large n.

# Let's check: is the characteristic polynomial of M(n) factored over Q(sqrt(2))?
# The char poly is degree 3 in lambda. Its roots are the eigenvalues.
# For generic n, this factors as (lambda - lambda_3(n)) * (quadratic with irrational roots).
# The question is: is the linear factor (lambda - lambda_3(n)) over Q[n]?

# This would mean lambda_3(n) is a rational function of n (a polynomial divided by leading coeff).
# If so, the char poly factors as (lambda - lambda_3(n)) * (quadratic) over Q(n).

print("\n\n=== Testing if char poly of M(n) factors over Q(n) ===")
lam = symbols('lambda')
M = M_sym(n)
charpoly_expr = M.charpoly(lam)
# Get the polynomial
cp = Poly(charpoly_expr.as_expr(), lam, domain='ZZ[n]')
print(f"Char poly degree in lambda: {cp.degree()}")
print(f"Char poly coefficients (as functions of n):")
coeffs = cp.all_coeffs()
for i, c_val in enumerate(coeffs):
    p = Poly(expand(c_val), n)
    print(f"  coeff of lambda^{3-i}: degree {p.degree()}, LC = {p.LC()}")

# Try to factor the char poly over Q(n)
print("\nAttempting to factor char poly over Q(n)...")
factored = factor(charpoly_expr.as_expr(), lam)
print(f"Factored form: {factored}")
