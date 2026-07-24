#!/usr/bin/env python3
"""Compute the full Poincaré polynomial of the order-3 operator for Problem 2.5.

The reduced operator: α₃(n) q_{n+3} + α₂(n) q_{n+2} + α₁(n) q_{n+1} + α₀(n) q_n = 0
with deg(α₃)=7, deg(α₂)=14, deg(α₁)=21, deg(α₀)=28.

Poincaré polynomial: LC(α₃) t³ + LC(α₂) t² + LC(α₁) t + LC(α₀) = 0
after dividing by n^28.

But with step-7 degrees, we need to normalize differently:
The Poincaré polynomial comes from setting q_n ~ r^n * n^ρ and looking at n→∞.
With degree pattern (7, 14, 21, 28), each coeff contributes n^{deg} * r^{shift}.
So the char equation involves LC(α₃)*r³ + LC(α₂)*r² + LC(α₁)*r + LC(α₀) = 0
but weighted by the degree gaps.

Actually: the recurrence is α₃ q_{n+3} + α₂ q_{n+2} + α₁ q_{n+1} + α₀ q_n = 0.
For large n, q_n ~ C * r^n * n^ρ, and α_j(n) ~ a_j * n^{d_j}.
The dominant balance is: a₃ n^7 r³ + a₂ n^{14} r² + a₁ n^{21} r + a₀ n^{28} = 0.
This is NOT the standard Poincaré polynomial because the degrees are different.

The standard Poincaré polynomial requires that all α_j have the same degree.
If they don't, we need the BIRKHOFF-ADAMS theory instead.

But wait — the degree-7 stepping means we can do a gauge transform.
Set q_n = n^{-7} * q̃_n, then the new recurrence for q̃_n has coefficients
with degrees shifted. More specifically:

α₃(n) * (n+3)^{-7} * q̃_{n+3} + α₂(n) * (n+2)^{-7} * q̃_{n+2} + ...

Actually the right thing: divide the recurrence by n^{28} (the max degree):
(α₃/n^{28}) * r³ + (α₂/n^{28}) * r² + (α₁/n^{28}) * r + (α₀/n^{28}) = 0
→ (a₃/n^{21}) r³ + (a₂/n^{14}) r² + (a₁/n^7) r + a₀ = 0

So for finite r, only the α₀ term survives unless r ~ n^7. This means
the Poincaré roots are of order n^7, i.e., q_n grows superexponentially.

This matches the numerical check: v_k/q_k ~ k^{14} growth suggests
the ratio of consecutive terms grows polynomially.

The CORRECT Poincaré analysis: set q_n = w_n * h(n) where h(n) absorbs
the polynomial growth. With h(n) = n^{7n} or similar (superexponential gauge).

Actually, let me just compute the leading coefficients and think about it.
"""
from sympy import symbols, Matrix, Poly, factor, expand, gcd, cancel, Rational

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

print("=== Characteristic polynomial of M(n) for large n ===")
M = M_sym(n)
charpoly = M.charpoly()
print(f"Char poly of M(n) (in lambda): {charpoly}")

# The eigenvalues of M(n) for large n give the Poincaré roots
# Let's compute det(M(n))
det_M = M.det()
det_M = expand(det_M)
det_poly = Poly(det_M, n)
print(f"\ndet(M(n)): degree = {det_poly.degree()}, LC = {det_poly.LC()}")
print(f"det(M(n)) = {factor(det_M)}")

# trace(M(n))
tr_M = M.trace()
tr_poly = Poly(expand(tr_M), n)
print(f"\ntr(M(n)): degree = {tr_poly.degree()}, LC = {tr_poly.LC()}")

# Now let's look at the eigenvalue asymptotics.
# For large n, M(n) ~ n^d * A where A is the leading coefficient matrix.
# Each entry M[i,j] is a polynomial in n. Let's find the degree and LC of each.
print("\n=== Degree and leading coefficient of each M(n) entry ===")
for i in range(3):
    for j in range(3):
        p = Poly(expand(M[i,j]), n)
        print(f"  M[{i},{j}]: deg={p.degree()}, LC={p.LC()}")

# For the Poincaré analysis, we need the leading matrix.
# But the entries have DIFFERENT degrees! Let's check.
