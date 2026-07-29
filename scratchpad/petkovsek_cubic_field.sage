"""
Petkovšek analysis for P2.7 recurrence over Q(μ₀)(n).

σ=0 removes the obstruction to first-order factoring over the cubic field.
If a hypergeometric solution exists over Q(μ₀)(n), we get explicit w⁰.

Strategy:
1. Build the order-3 recurrence operator in the shift algebra over Q(μ₀)(n)
2. Try right_factor(order=1) via ore_algebra
3. If that fails, implement Petkovšek's algorithm manually:
   - Find the universal denominator from roots of leading/trailing coefficients
   - Enumerate candidate rational functions r(n)
   - Check which ones satisfy the Riccati equation
"""
from sage.all import *
from ore_algebra import OreAlgebra

print("=== Petkovšek analysis over Q(μ₀)(n) ===\n")

# Step 1: Define the number field K = Q(μ₀)
Rx_qq.<x_qq> = PolynomialRing(QQ)
# Characteristic polynomial: 4μ³-220μ²+8μ-1 = 0
# Use the monic form: μ³ - 55μ² + 2μ - 1/4 = 0
# Or equivalently: work with ν = μ₀ root of 4x³-220x²+8x-1
K.<mu0> = NumberField(4*x_qq^3 - 220*x_qq^2 + 8*x_qq - 1)
print(f"Number field K = Q(μ₀) defined, degree {K.degree()}")
print(f"Minimal polynomial: {K.defining_polynomial()}")

# Verify μ₀ is the correct root (dominant, ≈ 54.96)
embs = K.embeddings(RealField(100))
real_embs = [e for e in embs if e(mu0) > 1]  # dominant root > 1
if real_embs:
    print(f"μ₀ ≈ {float(real_embs[0](mu0)):.10f}")
else:
    print("WARNING: no real embedding > 1 found")
    for e in K.embeddings(ComplexField(100)):
        print(f"  embedding: {e(mu0)}")

# Step 2: Build the recurrence operator over K(n)
Rn_K.<n> = PolynomialRing(K)
Kn = Rn_K.fraction_field()

# The shift algebra K(n)[S]
OS.<Sn> = OreAlgebra(Kn, 'Sn')

# Define the P2.7 coefficients over K
def A_K(z):
    return K(1024)*(2*z+5)^4*(2*z+7)^3*(2*z+9)^3*(K(946)*z^2+K(6407)*z+K(10860))

def B_K(z):
    return K(128)*(2*z+7)^3*(2*z+9)^3*(K(104060)*z^6+K(1745370)*z^5+K(12145238)*z^4+K(44886481)*z^3+K(92943995)*z^2+K(102256019)*z+K(46709052))

def C_K(z):
    return K(16)*(z+3)^4*(2*z+9)^3*(K(3784)*z^5+K(57792)*z^4+K(351019)*z^3+K(1059230)*z^2+K(1587211)*z+K(944620))

def D_K(z):
    return (z+3)^4*(z+4)^6*(K(946)*z^2+K(4515)*z+K(5399))

# Build the operator: A(n)S³ - B(n+2)S² + C(n+1)S - D(n)
# L u_n = A(n) u_{n+3} - B(n+2) u_{n+2} + C(n+1) u_{n+1} - D(n) u_n = 0
# In shift operator form:
L_rec = Kn(A_K(n))*Sn^3 - Kn(B_K(n+2))*Sn^2 + Kn(C_K(n+1))*Sn - Kn(D_K(n))

print(f"\nRecurrence operator L: order = {L_rec.order()}")

# Step 3: Try ore_algebra's built-in methods
print("\n=== Attempting factorization via ore_algebra ===")

# Method 1: right_factor of order 1
print("\nTrying right_factor(1)...")
try:
    rf = L_rec.right_factor(1)
    print(f"  SUCCESS! Right factor of order 1: {rf}")
    # If found, r(n) = -rf[0]/rf[1] gives the ratio u_{n+1}/u_n
    r_n = -rf[0]/rf[1]
    print(f"  r(n) = u_{{n+1}}/u_n = {r_n}")
except Exception as e:
    print(f"  Failed: {e}")

# Method 2: Try to factor the full operator
print("\nTrying factor()...")
try:
    facts = L_rec.factor()
    print(f"  SUCCESS! Factorization: {facts}")
except Exception as e:
    print(f"  Failed: {e}")

# Method 3: right_factors (enumerate all)
print("\nTrying right_factors(1)...")
try:
    rfs = L_rec.right_factors(1)
    print(f"  Found {len(rfs)} right factors of order 1:")
    for i, rf in enumerate(rfs):
        print(f"    [{i}]: {rf}")
except Exception as e:
    print(f"  Failed: {e}")

# Step 4: Manual Petkovšek analysis
print("\n=== Manual Petkovšek analysis ===")

# Get the polynomial (non-rational) form of the recurrence
# L = p₃(n)S³ + p₂(n)S² + p₁(n)S + p₀(n)
# where p₃(n) = A(n), p₀(n) = -D(n)

# Leading coefficient (at S³): A(n)
# Trailing coefficient (at S⁰): -D(n)

# Factor these over K
p3_poly = Rn_K(A_K(n))
p0_poly = Rn_K(-D_K(n))

print(f"\np₃(n) = A(n), degree = {p3_poly.degree()}")
print(f"p₀(n) = -D(n), degree = {p0_poly.degree()}")

# Factor p₃ and p₀ over K[n]
print("\nFactoring p₃(n) = A(n) over K[n]:")
p3_factors = p3_poly.factor()
print(f"  {p3_factors}")

print("\nFactoring p₀(n) = -D(n) over K[n]:")
p0_factors = p0_poly.factor()
print(f"  {p0_factors}")

# Roots of p₃ and p₀
print("\nRoots of p₃(n) over K:")
for f, m in p3_factors:
    if f.degree() == 1:
        root = -f[0]/f[1]
        print(f"  n = {root} (multiplicity {m})")
    else:
        print(f"  irreducible factor of degree {f.degree()} (mult {m}): {f}")

print("\nRoots of p₀(n) over K:")
for f, m in p0_factors:
    if f.degree() == 1:
        root = -f[0]/f[1]
        print(f"  n = {root} (multiplicity {m})")
    else:
        print(f"  irreducible factor of degree {f.degree()} (mult {m}): {f}")

# Petkovšek's universal denominator:
# For r(n) = u_{n+1}/u_n, the possible poles of r(n) come from:
# - roots α of p₀(n): r(n) can have a pole at n = α
# - roots β of p₃(n): r(n) can have a zero at n = β (but shifted)
# The key constraint is on INTEGER DIFFERENCES between roots of p₃ and p₀.

print("\n=== Integer difference analysis ===")
# Collect all roots
roots_p3 = []
for f, m in p3_factors:
    if f.degree() == 1:
        roots_p3.append((-f[0]/f[1], m))

roots_p0 = []
for f, m in p0_factors:
    if f.degree() == 1:
        roots_p0.append((-f[0]/f[1], m))

print(f"\nRoots of p₃: {[r for r,m in roots_p3]}")
print(f"Roots of p₀: {[r for r,m in roots_p0]}")

# Check for integer differences α - β where α ∈ roots(p₀), β ∈ roots(p₃)
print("\nInteger differences (root_of_p₀ - root_of_p₃):")
for r0, m0 in roots_p0:
    for r3, m3 in roots_p3:
        diff = r0 - r3
        if diff in QQ:
            d = QQ(diff)
            if d in ZZ:
                print(f"  {r0} - {r3} = {d} (INTEGER)")
            else:
                pass  # non-integer rational, skip

# Step 5: Characteristic polynomial at infinity
# The Poincaré roots are solutions of 4μ³-220μ²+8μ-1 = 0
# Over K, this factors as (μ - μ₀) · q(μ) for some quadratic q
print("\n=== Characteristic polynomial over K ===")
chi_poly = K['mu'](4*K['mu'].gen()^3 - 220*K['mu'].gen()^2 + 8*K['mu'].gen() - 1)
print(f"χ(μ) = {chi_poly}")
chi_factors = chi_poly.factor()
print(f"Factored: {chi_factors}")

# The Poincaré root λ₀ = μ₀/64 (since the recurrence has degree 18 in n)
lambda0 = mu0 / K(64)
print(f"\nλ₀ = μ₀/64 = {lambda0}")

# The adjoint has multipliers 1/λ₀, so the slow adjoint base is 64/μ₀
# In K: 64/μ₀
# Since 4μ₀³-220μ₀²+8μ₀-1 = 0, we have μ₀(4μ₀²-220μ₀+8) = 1
# So 1/μ₀ = 4μ₀²-220μ₀+8
inv_mu0 = 4*mu0^2 - 220*mu0 + 8
print(f"1/μ₀ = 4μ₀²-220μ₀+8 = {inv_mu0}")
assert inv_mu0 * mu0 == 1, "Inverse check failed!"
print("  Verified: (4μ₀²-220μ₀+8)·μ₀ = 1 ✓")

print(f"\n64/μ₀ = {K(64)*inv_mu0}")

# Step 6: Build the ADJOINT recurrence and try to factor it
# Adjoint: L̃ = -D(n)S³ + C(n)S² - B(n)S + A(n-3)
print("\n=== Adjoint recurrence ===")
L_adj = -Kn(D_K(n))*Sn^3 + Kn(C_K(n))*Sn^2 - Kn(B_K(n))*Sn + Kn(A_K(n-3))
print(f"Adjoint operator L̃: order = {L_adj.order()}")

print("\nTrying right_factor(1) on the ADJOINT...")
try:
    rf_adj = L_adj.right_factor(1)
    print(f"  SUCCESS! Right factor: {rf_adj}")
except Exception as e:
    print(f"  Failed: {e}")

print("\nTrying factor() on the ADJOINT...")
try:
    facts_adj = L_adj.factor()
    print(f"  SUCCESS! Factorization: {facts_adj}")
except Exception as e:
    print(f"  Failed: {e}")

# Step 7: Try the "desingularized" form
# Divide through by the leading coefficient to get monic form,
# then try Petkovšek's method via the Riccati substitution
print("\n=== Riccati substitution check ===")
# For the forward recurrence in monic form:
# u_{n+3} + (p₂/p₃)(n) u_{n+2} + (p₁/p₃)(n) u_{n+1} + (p₀/p₃)(n) u_n = 0
# Substitute u_n = λ₀ⁿ · v_n where v_{n+1}/v_n → 1
# Then v satisfies: v_{n+3} + (p₂/p₃)(n)/λ₀ · v_{n+2} + (p₁/p₃)(n)/λ₀² · v_{n+1} + (p₀/p₃)(n)/λ₀³ · v_n = 0

# For σ=0, the Birkhoff theory says v_n → const, so w_n = v_n satisfies a recurrence
# with Poincaré root 1 (dominant) and σ=0.

# Equivalently, conjugate L by λ₀:
# (S - λ₀) divides the transformed operator iff r(n) = λ₀ for all n,
# which would mean u_n = c·λ₀ⁿ exactly. This won't happen (the solution isn't geometric).
# But a generalized Pochhammer ratio might work.

print("\nDone with initial analysis.")
print("\nSummary of Petkovšek prospects:")
print("- The recurrence is order 3 over K(n)")
print("- p₃(n) = A(n) has roots: n ∈ {-5/2, -7/2, -9/2} (half-integers) + roots of 946n²+6407n+10860")
print("- p₀(n) = -D(n) has roots: n ∈ {-3, -4} + roots of 946n²+4515n+5399")
print("- σ=0 means no n^σ correction (Petkovšek over K is viable)")
print("- If ore_algebra methods fail, need manual Riccati equation construction")
