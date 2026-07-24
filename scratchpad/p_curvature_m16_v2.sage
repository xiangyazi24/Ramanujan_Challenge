"""
P-curvature computation for M₁₆ (v2 — fixed OreAlgebra variable naming).

Computes D^p mod M₁₆ in F_p(t)[D] for good primes p.
If Ψ_p = 0 for all tested primes → global nilpotence → geometric origin.
"""
from sage.all import *
from ore_algebra import OreAlgebra

# Build M₁₆ over Q(t)
Rn.<n> = PolynomialRing(QQ)
Kn = Rn.fraction_field()
Rt.<t> = PolynomialRing(QQ)
Kt = Rt.fraction_field()
OD.<Dt> = OreAlgebra(Kt, 'Dt')
theta = Kt(t)*Dt

def eval_poly_at_op(f, T):
    out = OD.zero()
    for a in reversed(Rn(f).list()):
        out = out*T + Kt(a)
    return out

def primitive_coefficients(c):
    c = [Kn(a) for a in c]
    den = Rn.one()
    for a in c:
        den = lcm(den, Rn(a.denominator()))
    nums = [Rn(den*a) for a in c]
    g = nums[0]
    for a in nums[1:]:
        g = gcd(g, a)
    nums = [a.quo_rem(g)[0] for a in nums]
    scl = QQ(1)/QQ(nums[-1].leading_coefficient())
    return [Rn(scl*a) for a in nums]

# Original P2.7 rational recurrence
def A(z): return 1024*(2*z+5)^4*(2*z+7)^3*(2*z+9)^3*(946*z^2+6407*z+10860)
def B(z): return 128*(2*z+7)^3*(2*z+9)^3*(104060*z^6+1745370*z^5+12145238*z^4+44886481*z^3+92943995*z^2+102256019*z+46709052)
def C(z): return 16*(z+3)^4*(2*z+9)^3*(3784*z^5+57792*z^4+351019*z^3+1059230*z^2+1587211*z+944620)
def D_coeff(z): return (z+3)^4*(z+4)^6*(946*z^2+4515*z+5399)

crat = [Kn(-D_coeff(n)/A(n)), Kn(C(n+1)/A(n+1)), Kn(-B(n+2)/A(n+2)), Kn(1)]
p_coeff = primitive_coefficients(crat)

# Build Euler quotient M₁₆
Rx.<x> = PolynomialRing(QQ)
shifted = [Rx(p_coeff[j](x-j)) for j in range(4)]
gshift = shifted[0]
for f in shifted[1:]:
    gshift = gcd(gshift, f)
gshift = gshift.monic()
assert gshift == x^2 + QQ(105)/22*x + QQ(5399)/946

Rpoly = []
for f in shifted:
    q, rem = f.quo_rem(gshift)
    assert rem == 0
    Rpoly.append(q)

M16 = OD.zero()
for j in range(4):
    M16 += Kt(t^j)*eval_poly_at_op(Rpoly[j], -theta-1)
assert M16.order() == 16
print(f"M16 order = {M16.order()}")

# Extract polynomial coefficients and clear denominators
coeffs_M16 = [M16[k] for k in range(17)]
denom = Rt.one()
for c in coeffs_M16:
    c = Kt(c)
    denom = lcm(denom, Rt(c.denominator()))
int_coeffs = [Rt(denom * Kt(c)) for c in coeffs_M16]

# Scale to make everything integer-coefficient
# Find the overall lcm of all denominators of all polynomial coefficients
from fractions import Fraction
overall_denom = ZZ(1)
for poly in int_coeffs:
    for coeff in poly.list():
        if coeff != 0:
            overall_denom = lcm(overall_denom, ZZ(QQ(coeff).denominator()))
int_coeffs_Z = [Rt([ZZ(QQ(c)*overall_denom) for c in poly.list()]) for poly in int_coeffs]
print(f"Cleared denominators (scale = {overall_denom})")
print(f"Coefficient degrees: {[c.degree() for c in int_coeffs_Z]}")

# Verify leading coefficient
lc = int_coeffs_Z[16]
print(f"Leading coeff degree = {lc.degree()}")

# P-curvature computation
# Work directly over F_p[t] with variable name 't' (not 'tp')
def p_curvature_mod(int_coeffs_Z, p_prime, verbose=True):
    """
    Compute D^p mod L in F_p(t)[D].
    Returns 0 if p-curvature vanishes, else the order of D^p mod L.
    """
    Fp = GF(p_prime)

    # Use a fresh polynomial ring over Fp with variable 's' to avoid name conflicts
    Rp = PolynomialRing(Fp, 's')
    s = Rp.gen()
    Kp = Rp.fraction_field()
    ODp = OreAlgebra(Kp, ('Ds', {}, {s: 1}))
    Ds = ODp.gen()

    d = len(int_coeffs_Z) - 1  # order

    # Reduce coefficients mod p, substituting t -> s
    L_mod = ODp.zero()
    for k in range(d+1):
        coeffs_list = Rt(int_coeffs_Z[k]).list()
        poly_mod = sum(Fp(ZZ(c)) * s^i for i, c in enumerate(coeffs_list))
        L_mod += Kp(poly_mod) * Ds^k

    if L_mod == 0:
        if verbose:
            print(f"  p={p_prime}: operator ZERO mod p (bad prime)")
        return -1

    if L_mod[d] == 0:
        if verbose:
            print(f"  p={p_prime}: leading coeff vanishes (bad prime)")
        return -1

    # Compute D^p mod L using repeated squaring
    # D^p in the Ore algebra: use the algebra's power method
    if verbose:
        print(f"  p={p_prime}: computing Ds^{p_prime} mod L (order {d})...")

    # Method: iterative squaring of Ds, reducing mod L at each step
    # Represent power as binary
    result = ODp.one()  # Ds^0
    base = Ds           # Ds^1
    exp = p_prime

    while exp > 0:
        if exp % 2 == 1:
            result = result * base
            if result.order() >= d:
                _, result = result.quo_rem(L_mod)
        exp //= 2
        if exp > 0:
            base = base * base
            if base.order() >= d:
                _, base = base.quo_rem(L_mod)

    if result == 0:
        if verbose:
            print(f"  p={p_prime}: Ψ_p = 0 ★ GLOBALLY NILPOTENT at this prime!")
        return 0
    else:
        pc_order = result.order()
        if verbose:
            print(f"  p={p_prime}: Ψ_p ≠ 0, order = {pc_order}")
        return pc_order

print("\n=== P-curvature of M16 ===")
for pp in [5, 7, 11, 13, 17, 19, 23, 29, 31]:
    print(f"\nPrime p = {pp}:")
    try:
        rank = p_curvature_mod(int_coeffs_Z, pp)
    except Exception as e:
        print(f"  FAILED: {e}")
        import traceback
        traceback.print_exc()

print("\nDone.")
