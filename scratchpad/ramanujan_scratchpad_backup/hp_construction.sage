#!/usr/bin/env sage
"""
Hermite-Padé Type II polynomials for ζ(2)+ζ(3).

Two measures on (0,1):
  dμ₂ = -log(t) dt,   moments: ∫ t^k dμ₂ = 1/(k+1)²
  dμ₃ = ½ log²(t) dt, moments: ∫ t^k dμ₃ = 1/(k+1)³

Step-line multi-index: at step n, orthogonality constraints are
  ∫ t^i Q_n(t) dμ₂ = 0 for i=0,...,n₂-1
  ∫ t^i Q_n(t) dμ₃ = 0 for i=0,...,n₃-1
where n₂+n₃ = n (or n+1, depending on convention).

We try several conventions for the multi-index.
"""

# Moments
def mu2(i, j):
    """∫₀¹ t^{i+j} (-log t) dt = 1/(i+j+1)²"""
    return QQ(1) / QQ(i+j+1)^2

def mu3(i, j):
    """∫₀¹ t^{i+j} (½ log² t) dt = 1/(i+j+1)³"""
    return QQ(1) / QQ(i+j+1)^3

def compute_HP(n, n2, n3, normalize_at=None):
    """
    Find Q_n(x) of degree n such that:
      ∫ t^i Q_n(t) dμ₂ = 0 for i=0,...,n2-1
      ∫ t^i Q_n(t) dμ₃ = 0 for i=0,...,n3-1
    Total conditions: n2+n3. Must equal n for a unique (up to scalar) solution.
    """
    assert n2 + n3 == n, f"Need n2+n3=n, got {n2}+{n3}={n2+n3} vs n={n}"

    # Build the n × (n+1) constraint matrix
    M = matrix(QQ, n, n+1)
    for i in range(n2):
        for j in range(n+1):
            M[i, j] = mu2(i, j)
    for i in range(n3):
        for j in range(n+1):
            M[n2+i, j] = mu3(i, j)

    ker = M.right_kernel()
    if ker.dimension() != 1:
        return None, ker.dimension()

    v = ker.basis()[0]
    R_poly = PolynomialRing(QQ, 'x')
    x = R_poly.gen()
    Q_poly = sum(v[j] * x^j for j in range(n+1))

    if normalize_at is not None:
        val = Q_poly(normalize_at)
        if val != 0:
            Q_poly = Q_poly / val
    return Q_poly, 1

# Try different step-line conventions
print("="*60)
print("HP Type II for dμ₂, dμ₃")
print("="*60)

N_max = 10
results = {}

for convention in ['balanced', 'mu2_first', 'mu3_first']:
    print(f"\n--- Convention: {convention} ---")
    q_vals = []

    for n in range(1, N_max+1):
        if convention == 'balanced':
            n2 = (n+1) // 2
            n3 = n - n2
        elif convention == 'mu2_first':
            n2 = (n + 1) // 2
            n3 = n - n2
        elif convention == 'mu3_first':
            n3 = (n + 1) // 2
            n2 = n - n3

        Q_poly, dim = compute_HP(n, n2, n3)
        if Q_poly is not None:
            q0 = Q_poly(0)
            q_vals.append(q0)
            print(f"  n={n}: n2={n2}, n3={n3}, Q_n(0) = {q0}")
        else:
            q_vals.append(None)
            print(f"  n={n}: n2={n2}, n3={n3}, kernel dim = {dim}")

    results[convention] = q_vals

# Now check: does q_n(0) satisfy any recognizable recurrence?
from ore_algebra import *
from ore_algebra import guess

R.<n> = PolynomialRing(QQ)
A.<Sn> = OreAlgebra(R)

for convention, vals in results.items():
    # Remove Nones
    clean = [v for v in vals if v is not None]
    if len(clean) < 6:
        continue
    print(f"\n{'='*60}")
    print(f"Guessing recurrence for {convention} (values: {clean[:6]})")
    print(f"{'='*60}")
    try:
        L = guess(clean, A)
        print(f"  Order: {L.order()}, Degree: {max(c.degree() for c in L.list())}")
        # Print the operator
        if L.order() <= 3 and max(c.degree() for c in L.list()) <= 20:
            print(f"  L = {L}")
    except Exception as e:
        print(f"  Could not guess: {e}")

# Also try: compute remainders and check if ζ(2)+ζ(3) approximation works
print(f"\n{'='*60}")
print("Remainder analysis (balanced convention)")
print(f"{'='*60}")

for n in range(1, min(N_max+1, 8)):
    n2 = (n+1) // 2
    n3 = n - n2
    Q_poly, dim = compute_HP(n, n2, n3)
    if Q_poly is None:
        continue

    # Compute ∫ Q_n(t) dμ₂ and ∫ Q_n(t) dμ₃
    coeffs = Q_poly.list()
    int2 = sum(coeffs[j] * QQ(1)/(j+1)^2 for j in range(len(coeffs)))
    int3 = sum(coeffs[j] * QQ(1)/(j+1)^3 for j in range(len(coeffs)))

    q0 = Q_poly(0)
    if q0 != 0:
        # p_n ≈ (int2 + int3) * Q_n(0) if the approximation is good
        ratio = (int2 + int3) / q0
        print(f"  n={n}: (int2+int3)/Q(0) = {float(ratio):.15f}, "
              f"int2 = {float(int2):.10e}, int3 = {float(int3):.10e}")

print(f"\n  ζ(2)+ζ(3) ≈ {float(pi^2/6 + zeta(3)):.15f}")

# Also: HP with normalization Q_n(0) = 1
print(f"\n{'='*60}")
print("Alternative: HP with Q_n(0) normalized to 1")
print(f"{'='*60}")
for n in range(1, min(N_max+1, 8)):
    n2 = (n+1) // 2
    n3 = n - n2
    Q_poly, dim = compute_HP(n, n2, n3, normalize_at=0)
    if Q_poly is None:
        continue
    # Compute approximation
    coeffs = Q_poly.list()
    int2 = sum(coeffs[j] * QQ(1)/(j+1)^2 for j in range(len(coeffs)))
    int3 = sum(coeffs[j] * QQ(1)/(j+1)^3 for j in range(len(coeffs)))
    print(f"  n={n}: int2+int3 = {float(int2+int3):.15f}, "
          f"error = {float(int2+int3-pi^2/6-zeta(3)):.6e}")
