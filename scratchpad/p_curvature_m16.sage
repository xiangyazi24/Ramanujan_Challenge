"""
P-curvature computation for M₁₆.

If the p-curvature is 0 for all good primes, M₁₆ is globally nilpotent
(comes from geometry), which would strongly support the period interpretation.

Algorithm: Compute D^p mod M₁₆ in F_p(t)[D].
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
def D(z): return (z+3)^4*(z+4)^6*(946*z^2+4515*z+5399)

crat = [Kn(-D(n)/A(n)), Kn(C(n+1)/A(n+1)), Kn(-B(n+2)/A(n+2)), Kn(1)]
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
print(f"M16 leading term degree = {Rt(M16[16].numerator()).degree()}")

# Extract polynomial coefficients of M16 = Σ c_k(t) Dt^k
# Clear denominators to work over Z[t]
coeffs_M16 = [M16[k] for k in range(17)]  # c_0(t) through c_16(t)

# Find common denominator
denom = Rt.one()
for c in coeffs_M16:
    c = Kt(c)
    denom = lcm(denom, Rt(c.denominator()))

int_coeffs = [Rt(denom * Kt(c)) for c in coeffs_M16]
print(f"\nCommon denominator degree = {denom.degree()}")
print(f"Integer coefficient degrees = {[c.degree() for c in int_coeffs]}")

# P-curvature: compute D^p mod L in F_p(t)[D]
# For L = sum c_k(t) D^k of order d, the p-curvature is the remainder of D^p mod L.

def p_curvature_mod(L_int_coeffs, p_prime, verbose=True):
    """
    Compute the p-curvature of L (given by integer polynomial coefficients)
    modulo a prime p.

    Returns the rank of the p-curvature matrix (order of D^p mod L).
    If result is 0, the p-curvature vanishes at this prime.
    """
    Fp = GF(p_prime)
    Rp.<tp> = PolynomialRing(Fp)
    Kp = Rp.fraction_field()
    ODp.<Dp> = OreAlgebra(Kp, 'Dp')

    d = len(L_int_coeffs) - 1  # order of L

    # Reduce coefficients mod p
    L_mod = ODp.zero()
    for k in range(d+1):
        poly_mod = Rp([Fp(c) for c in Rt(L_int_coeffs[k]).list()])
        L_mod += Kp(poly_mod) * Dp^k

    if L_mod == 0:
        if verbose:
            print(f"  p={p_prime}: operator is ZERO mod p (bad prime)")
        return -1

    # Check if leading coefficient vanishes
    if L_mod[d] == 0:
        if verbose:
            print(f"  p={p_prime}: leading coefficient vanishes (bad prime)")
        return -1

    # Compute D^p mod L
    # Start with D^1, repeatedly multiply by D and reduce
    # In ore_algebra: D * f(t) = f(t)*D + f'(t) (Leibniz rule)

    # Actually, use the right_quo_rem from ore_algebra
    Dp_power = Dp  # D^1
    for i in range(2, p_prime + 1):
        Dp_power = Dp_power * Dp
        # Reduce mod L
        if Dp_power.order() >= d:
            _, Dp_power = Dp_power.quo_rem(L_mod)
        if i % max(1, p_prime // 10) == 0 and verbose:
            print(f"    step {i}/{p_prime}, current order = {Dp_power.order()}")

    # The p-curvature is D^p mod L
    pc_order = Dp_power.order() if Dp_power != 0 else -1

    if verbose:
        if Dp_power == 0:
            print(f"  p={p_prime}: p-curvature = 0 (GLOBALLY NILPOTENT at this prime!)")
        else:
            print(f"  p={p_prime}: p-curvature has order {pc_order}")

    return pc_order if Dp_power != 0 else 0

# Check which primes are "good"
# The leading coefficient is t^16 * chi(t) where chi = t^3 - 55/64*t^2 + ...
# Clearing denominators of chi: 2^20 * chi(t) = 2^20*t^3 - 55*2^14*t^2 + 2^10*t - 1
print(f"\nLeading coefficient involves 2^20 in denominators")
print(f"Bad primes for chi: 2 (denominator factor)")

# Compute p-curvature for small primes
# Note: p must be > order for meaningful results, and computation of D^p is O(p^2)
# For order 16, smallest useful primes are 17, 19, 23, ...
# But even p=5 can be informative (just need p > 0)

print("\n=== P-curvature of M16 ===")
# Start with small primes for speed
for pp in [5, 7, 11, 13, 17, 19, 23]:
    print(f"\nPrime p = {pp}:")
    try:
        rank = p_curvature_mod(int_coeffs, pp)
    except Exception as e:
        print(f"  FAILED: {e}")

print("\nDone.")
