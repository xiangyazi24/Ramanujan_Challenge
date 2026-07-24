"""
Riccati equation search for first-order right factor of P2.7 recurrence over Q(μ₀)(n).

KEY DISCOVERY: The quadratic factors of the leading (p₃) and trailing (p₀) coefficients
are related by an integer shift:
  Q_p₃(n) = n² + 149/22·n + 5430/473   (from A(n))
  Q_p₀(n) = n² + 105/22·n + 5399/946   (from D(n))
  Q_p₃(n-1) = Q_p₀(n)  ← VERIFIED

This is the ONLY integer-difference pair among the roots of p₃ and p₀.

ANSATZ for hypergeometric ratio r(n) = u_{n+1}/u_n:
  r(n) = λ₀ · Q_p₃(n) · P(n) / (Q_p₃(n-1) · P(n+1))
where P(n) is a degree-2 polynomial to be determined.

This gives σ = 2 - 2 = 0, consistent with the proved σ=0.
"""
from sage.all import *
from ore_algebra import OreAlgebra

print("=== Riccati search for first-order factor over Q(μ₀)(n) ===\n")

# Define number field K = Q(μ₀)
Rx.<x> = PolynomialRing(QQ)
K.<mu0> = NumberField(4*x^3 - 220*x^2 + 8*x - 1)
print(f"K = Q(μ₀), degree {K.degree()}")

# Poincaré root λ₀ = μ₀/64
lam0 = mu0 / K(64)
print(f"λ₀ = μ₀/64")

# Verify: 1/μ₀ = 4μ₀²-220μ₀+8
inv_mu0 = 4*mu0^2 - 220*mu0 + 8
assert inv_mu0 * mu0 == 1

# The quadratic factors
Rn.<n> = PolynomialRing(K)

Q_p3 = n^2 + K(149)/K(22)*n + K(5430)/K(473)  # from p₃ = A(n)
Q_p0 = n^2 + K(105)/K(22)*n + K(5399)/K(946)  # from p₀ = -D(n)

# Verify Q_p3(n-1) = Q_p0(n)
Q_p3_shifted = Q_p3(n-1)
print(f"\nQ_p₃(n) = {Q_p3}")
print(f"Q_p₀(n) = {Q_p0}")
print(f"Q_p₃(n-1) = {Q_p3_shifted}")
assert Q_p3_shifted == Q_p0, "SHIFT VERIFICATION FAILED"
print("✓ Q_p₃(n-1) = Q_p₀(n) verified!")

# Define the P2.7 recurrence coefficients
def A_K(z):
    return K(1024)*(2*z+5)^4*(2*z+7)^3*(2*z+9)^3*(K(946)*z^2+K(6407)*z+K(10860))

def B_K(z):
    return K(128)*(2*z+7)^3*(2*z+9)^3*(K(104060)*z^6+K(1745370)*z^5+K(12145238)*z^4+K(44886481)*z^3+K(92943995)*z^2+K(102256019)*z+K(46709052))

def C_K(z):
    return K(16)*(z+3)^4*(2*z+9)^3*(K(3784)*z^5+K(57792)*z^4+K(351019)*z^3+K(1059230)*z^2+K(1587211)*z+K(944620))

def D_K(z):
    return (z+3)^4*(z+4)^6*(K(946)*z^2+K(4515)*z+K(5399))

# The Riccati equation for r(n):
# p₃(n)·r(n)·r(n+1)·r(n+2) + p₂(n)·r(n)·r(n+1) + p₁(n)·r(n) + p₀(n) = 0
# where p₃=A, p₂=-B(n+2), p₁=C(n+1), p₀=-D(n)

p3 = Rn(A_K(n))
p2 = Rn(-B_K(n+2))
p1 = Rn(C_K(n+1))
p0 = Rn(-D_K(n))

print(f"\nDegrees: p₃={p3.degree()}, p₂={p2.degree()}, p₁={p1.degree()}, p₀={p0.degree()}")

# ANSATZ: r(n) = λ₀ · Q_p₃(n) · P(n) / (Q_p₃(n-1) · P(n+1))
# where P(n) = n² + p·n + q, p,q ∈ K
#
# Products:
# r(n)·r(n+1)·r(n+2) = λ₀³ · Q_p₃(n+2)·P(n) / (Q_p₃(n-1)·P(n+3))
# r(n)·r(n+1) = λ₀² · Q_p₃(n+1)·P(n) / (Q_p₃(n-1)·P(n+2))
# r(n) = λ₀ · Q_p₃(n)·P(n) / (Q_p₃(n-1)·P(n+1))

# Substitute into Riccati and clear denominator Q_p₃(n-1)·P(n+1)·P(n+2)·P(n+3):
# p₃·λ₀³·Q_p₃(n+2)·P(n)·P(n+1)·P(n+2)
# + p₂·λ₀²·Q_p₃(n+1)·P(n)·P(n+1)·P(n+3)
# + p₁·λ₀·Q_p₃(n)·P(n)·P(n+2)·P(n+3)
# + p₀·Q_p₃(n-1)·P(n+1)·P(n+2)·P(n+3) = 0

# Use symbolic p, q
K_pq.<pp, qq> = PolynomialRing(K)
Rn_pq.<nn> = PolynomialRing(K_pq)

# P(n) = nn^2 + pp*nn + qq
def P_pq(z):
    return z^2 + pp*z + qq

def Q_p3_pq(z):
    return z^2 + K(149)/K(22)*z + K(5430)/K(473)

# Build the Riccati polynomial
print("\nBuilding Riccati polynomial (this may take a moment)...")

term1 = Rn_pq(p3(nn)) * lam0^3 * Q_p3_pq(nn+2) * P_pq(nn) * P_pq(nn+1) * P_pq(nn+2)
term2 = Rn_pq(p2(nn)) * lam0^2 * Q_p3_pq(nn+1) * P_pq(nn) * P_pq(nn+1) * P_pq(nn+3)
term3 = Rn_pq(p1(nn)) * lam0  * Q_p3_pq(nn)   * P_pq(nn) * P_pq(nn+2) * P_pq(nn+3)
term4 = Rn_pq(p0(nn))         * Q_p3_pq(nn-1)  * P_pq(nn+1) * P_pq(nn+2) * P_pq(nn+3)

print("  Computing terms...")
F_riccati = term1 + term2 + term3 + term4
print(f"  Riccati polynomial degree in n: {F_riccati.degree()}")

# Extract coefficients of nn^k — each must be 0
coeffs_n = F_riccati.list()
print(f"  Number of coefficient equations: {len(coeffs_n)}")

# Each coefficient is a polynomial in pp, qq over K
# We need to solve the system coeffs_n[k](pp, qq) = 0 for all k

# Check the highest-degree coefficients first (they should involve the characteristic equation)
print("\n=== Checking coefficient equations ===")
for i, c in enumerate(reversed(coeffs_n)):
    if c == 0:
        print(f"  [n^{len(coeffs_n)-1-i}]: = 0 (automatically satisfied)")
    else:
        # Show degree and a few terms
        deg_p = max([e[0] for e in c.exponents()] + [0]) if c != 0 else -1
        deg_q = max([e[1] for e in c.exponents()] + [0]) if c != 0 else -1
        print(f"  [n^{len(coeffs_n)-1-i}]: nonzero, max pp-degree={deg_p}, max qq-degree={deg_q}")
    if i >= 10:
        remaining = sum(1 for c in coeffs_n if c != 0) - sum(1 for c in reversed(coeffs_n[:len(coeffs_n)-11]) if c != 0)
        print(f"  ... (showing first 11 coefficient equations)")
        break

# Try to solve: collect all nonzero coefficients as a system in pp, qq
print("\n=== Solving the system ===")
nonzero_coeffs = [c for c in coeffs_n if c != 0]
print(f"Number of nonzero equations: {len(nonzero_coeffs)}")

if len(nonzero_coeffs) == 0:
    print("All coefficients are zero — r(n) = λ₀ for ANY P(n)!")
else:
    # Try to find the ideal and solve
    I = K_pq.ideal(nonzero_coeffs)
    print(f"Ideal has {len(I.gens())} generators")

    # Compute Groebner basis
    print("Computing Groebner basis...")
    try:
        gb = I.groebner_basis()
        print(f"Groebner basis has {len(gb)} elements:")
        for g in gb:
            print(f"  {g}")

        if 1 in gb or K_pq.one() in gb:
            print("\n★ RESULT: The ideal is the whole ring — NO SOLUTION EXISTS.")
            print("  There is NO first-order right factor with this ansatz.")
        else:
            print("\n★ RESULT: Solution(s) exist!")
            # Try to find the variety
            try:
                V = I.variety()
                print(f"Variety has {len(V)} points:")
                for pt in V:
                    print(f"  pp = {pt[pp]}, qq = {pt[qq]}")
            except Exception as e:
                print(f"  (Could not compute variety: {e})")
                print(f"  Groebner basis: {gb}")
    except Exception as e:
        print(f"Groebner basis failed: {e}")

        # Fallback: try resultant
        print("\nFallback: computing resultant...")
        try:
            # Take two nonzero coefficients and compute resultant wrt pp
            c1 = nonzero_coeffs[0]
            c2 = nonzero_coeffs[1]
            res_pp = c1.resultant(c2, pp)
            print(f"Resultant wrt pp: degree in qq = {res_pp.degree(qq) if res_pp != 0 else 'zero'}")
            if res_pp == 0:
                print("Resultant is zero — equations share a common factor in pp")
            else:
                qq_roots = Rn(res_pp).roots(K)
                print(f"Roots in K: {qq_roots}")
        except Exception as e:
            print(f"Resultant failed: {e}")

# Also try the SIMPLER ansatz: r(n) = constant (sanity check)
print("\n=== Sanity check: r(n) = λ₀ ===")
F_const = p3 * lam0^3 + p2 * lam0^2 + p1 * lam0 + p0
print(f"p₃λ₀³ + p₂λ₀² + p₁λ₀ + p₀ = polynomial of degree {F_const.degree()}")
if F_const == 0:
    print("  ZERO — constant ratio works!")
else:
    print(f"  Nonzero — constant ratio fails (expected)")
    # Check leading coefficient
    lc_val = F_const.leading_coefficient()
    print(f"  Leading coefficient: {lc_val}")

print("\nDone.")
