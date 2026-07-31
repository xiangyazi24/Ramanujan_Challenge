#!/usr/bin/env sage
"""
Mellin-adjoint ODE for P2.7.

For recurrence sum_j p_j(n) u_{n+j} = 0 with moments u_n = int t^n phi(t) dt,
the density phi satisfies:
  M_dual = sum_j t^j * p_j(-theta-j-1) = 0
where theta = t*d/dt.

This gives an ODE of order <= max(deg p_j) = 18.
"""
from ore_algebra import *

Rn.<n> = PolynomialRing(QQ)
Rt.<t> = PolynomialRing(QQ)
Qt = Rt.fraction_field()
Dt_alg.<Dt> = OreAlgebra(Qt)

# P2.7 coefficient polynomials
A_fn = Rn(1024*(2*n+5)**4*(2*n+7)**3*(2*n+9)**3*(946*n**2+6407*n+10860))
B_fn = Rn(128*(2*n+7)**3*(2*n+9)**3*(104060*n**6+1745370*n**5+12145238*n**4+44886481*n**3+92943995*n**2+102256019*n+46709052))
C_fn = Rn(16*(n+3)**4*(2*n+9)**3*(3784*n**5+57792*n**4+351019*n**3+1059230*n**2+1587211*n+944620))
D_fn = Rn((n+3)**4*(n+4)**6*(946*n**2+4515*n+5399))

# Build the monic operator's primitive polynomial coefficients
Qn = Rn.fraction_field()
A_rec.<Sn> = OreAlgebra(Qn, 'Sn')

L_monic = (Sn**3
           - B_fn(n=n+2)/A_fn(n=n+2) * Sn**2
           + C_fn(n=n+1)/A_fn(n=n+1) * Sn
           - D_fn(n=n)/A_fn(n=n))

coeffs_monic = [Qn(L_monic[i]) for i in range(4)]
den = Rn.one()
for c in coeffs_monic:
    den = lcm(den, Rn(c.denominator()))
nums = [Rn(den * c) for c in coeffs_monic]
g = nums[0]
for c in nums[1:]:
    g = gcd(g, c)
pnums = [c // g for c in nums]

print("Primitive polynomial coefficients (degree 18 each):")
for j in range(4):
    print(f"  p{j}: degree {pnums[j].degree()}")

# Build theta = t*Dt as an operator in the Dt algebra
# Key: theta^k = (t*Dt)^k
# In ore_algebra: t*Dt is NOT directly theta, but we can use the relation
# Dt * t = t*Dt + 1, so t*Dt = Dt*t - 1, i.e. theta = Dt*t - 1
# More useful: t*Dt = theta, where theta*f = t*f'
# In the ODE, we need to evaluate p_j(-theta-j-1) which is a polynomial in theta

# Approach: directly compute the operator in the Dt algebra
# p_j(x) = sum_k c_{j,k} x^k
# p_j(-theta-j-1) = sum_k c_{j,k} (-theta-j-1)^k
# Each (-theta-j-1)^k expands to a polynomial in theta
# theta^k in the Dt algebra: theta = t*Dt, so theta^2 = t*Dt*t*Dt = t*(t*Dt+1)*Dt = t^2*Dt^2 + t*Dt

# Let me build theta as an operator and compute powers
# In ore_algebra with standard differential algebra R[t]<Dt>,
# the Euler operator is theta = t*Dt
# We compute (t*Dt)^k iteratively

def theta_power(k, alg, t_var, Dt_var):
    """Compute (t*Dt)^k in the given OreAlgebra."""
    if k == 0:
        return alg.one()
    result = t_var * Dt_var
    for _ in range(k - 1):
        # Multiply by theta = t*Dt on the LEFT
        result = t_var * Dt_var * result
    return result

print("\nBuilding Mellin adjoint operator...")

# For each j, compute p_j(-theta-j-1) * t^j
# Then sum over j

M_dual = Dt_alg.zero()

for j in range(4):
    # p_j(x) = sum_k pnums[j][k] * x^k
    # p_j(-theta-j-1) = sum_k pnums[j][k] * (-theta-j-1)^k
    # (-theta-j-1)^k = sum binomial terms of theta^m

    # First build (-theta-j-1)^k for each k
    # Let alpha = -(j+1)
    # (-theta-j-1) = -(theta + j + 1) = -(theta - alpha) where alpha = -(j+1)

    shift = -(j + 1)

    # Build the polynomial p_j evaluated at (shift - theta)
    # where shift = -(j+1) and the argument is -theta - j - 1 = -(theta + j + 1)

    # p_j(-(theta+j+1)) = sum_k pnums[j][k] * (-(theta+j+1))^k

    # Compute (-(theta+j+1))^k iteratively
    # Start with the operator -(theta+j+1) = -(t*Dt + j + 1)
    neg_theta_shift = -(Qt(t) * Dt + QQ(j + 1))

    # Build p_j(neg_theta_shift) = sum_k c_k * (neg_theta_shift)^k
    pj_eval = Dt_alg.zero()
    power = Dt_alg.one()  # (neg_theta_shift)^0 = 1

    deg_j = pnums[j].degree()
    for k in range(deg_j + 1):
        c_k = QQ(pnums[j][k])
        if c_k != 0:
            pj_eval += c_k * power
        if k < deg_j:
            power = power * neg_theta_shift

    # Multiply by t^j
    contribution = Qt(t)**j * pj_eval
    M_dual += contribution

    print(f"  j={j}: p_{j}(-theta-{j+1}) * t^{j} computed (order {pj_eval.order()})")

print(f"\nMellin adjoint operator M_dual:")
print(f"  Order: {M_dual.order()}")
deg = max(c.degree() for c in M_dual.list() if c != 0)
print(f"  Max coefficient degree: {deg}")

# Leading coefficient
lc = M_dual.leading_coefficient()
print(f"  Leading coefficient: degree {lc.degree()}")
print(f"  Leading coefficient factored: {lc.factor()}")

# Trailing coefficient (Dt^0)
tc = M_dual.list()[0]
print(f"  Trailing coefficient: degree {tc.degree()}")
if tc != 0:
    print(f"  Trailing coefficient factored: {tc.factor()}")

# Singular points
print(f"\nSingular points of M_dual:")
sing = lc.roots(QQbar)
for rt, mult in sing:
    print(f"  t = {CC(rt):.6f} (mult {mult})")

# Also check t = infinity
print(f"\n  t = infinity: order {M_dual.order()}, need indicial analysis")

# Indicial polynomial at t = 0
print(f"\nIndicial polynomial at t = 0:")
try:
    Ra.<alpha> = PolynomialRing(QQ)
    ind_0 = M_dual.indicial_polynomial(t, alpha)
    print(f"  Degree: {ind_0.degree()}")
    rts = ind_0.roots(QQbar)
    for rt, mult in sorted(rts, key=lambda x: CC(x[0]).real()):
        print(f"  alpha = {CC(rt):.6f} (mult {mult})")
except Exception as e:
    print(f"  Error: {e}")

# Verify: does the Mellin dual relate to the to_D operator?
# The generating function ODE and the Mellin dual ODE should be related
# by the substitution z ↔ 1/t (or similar)

print("\n--- Checking relationship to GF ODE ---")
Rz.<z> = PolynomialRing(QQ)
D_z_alg.<Dz> = OreAlgebra(Rz)
L_prim = sum((Qn(pnums[i]) * Sn**i for i in range(4)), A_rec.zero())
try:
    L_toD = L_prim.to_D(D_z_alg)
    print(f"GF ODE order: {L_toD.order()}")
    print(f"Mellin dual order: {M_dual.order()}")
    print(f"(Should be same or related by substitution z = 1/t)")
except Exception as e:
    print(f"to_D error: {e}")

print("\nDone.")
