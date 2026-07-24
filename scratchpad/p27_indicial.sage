#!/usr/bin/env sage
"""
Compute indicial polynomials at finite singularities of the P2.7 ODE,
and try to find right factors to reduce the order.
"""
from ore_algebra import *
from ore_algebra import guess

Rn.<n> = PolynomialRing(QQ)

A_fn = Rn(1024*(2*n+5)**4*(2*n+7)**3*(2*n+9)**3*(946*n**2+6407*n+10860))
B_fn = Rn(128*(2*n+7)**3*(2*n+9)**3*(104060*n**6+1745370*n**5+12145238*n**4+44886481*n**3+92943995*n**2+102256019*n+46709052))
C_fn = Rn(16*(n+3)**4*(2*n+9)**3*(3784*n**5+57792*n**4+351019*n**3+1059230*n**2+1587211*n+944620))
D_fn = Rn((n+3)**4*(n+4)**6*(946*n**2+4515*n+5399))

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

L_prim = sum((Qn(pnums[i]) * Sn**i for i in range(4)), A_rec.zero())

R2.<z> = PolynomialRing(QQ)
D_alg.<Dz> = OreAlgebra(R2)

print("Converting to ODE...")
L_ode = L_prim.to_D(D_alg)
print(f"Order: {L_ode.order()}, max degree: {max(c.degree() for c in L_ode.list())}")

lc = L_ode.leading_coefficient()
print(f"\nLeading coefficient: {lc.factor()}")

# The cubic factor
cubic = R2(z**3 - 512*z**2 + 901120*z - 1048576)
print(f"\nCubic singularity polynomial: {cubic}")
print(f"Roots (numerical):")
for rt in cubic.roots(CC):
    print(f"  z = {rt[0]} (mult {rt[1]})")

# Indicial polynomial at each finite singular point
print("\n" + "="*60)
print("Indicial polynomials at finite singular points")
print("="*60)

Ra.<alpha> = PolynomialRing(QQ)

# At z = 0
print("\n--- z = 0 ---")
try:
    ind_0 = L_ode.indicial_polynomial(z, alpha)
    print(f"Indicial polynomial degree: {ind_0.degree()}")
    print(f"Roots:")
    for rt, mult in ind_0.roots(QQbar):
        print(f"  alpha = {rt} (mult {mult})")
except Exception as e:
    print(f"Error: {e}")

# At z = infinity (substitute z -> 1/w)
print("\n--- z = infinity ---")
try:
    ind_inf = L_ode.indicial_polynomial(~z, alpha)
    print(f"Indicial polynomial degree: {ind_inf.degree()}")
    print(f"Roots:")
    for rt, mult in ind_inf.roots(QQbar):
        print(f"  alpha = {rt} (mult {mult})")
except Exception as e:
    print(f"Error: {e}")

# At the dominant singularity z = z0 (root of cubic)
# We work over the splitting field
print("\n--- z = z0 (dominant real root of cubic) ---")
Rz0.<z0> = QQ.extension(cubic)
print(f"z0 is a root of {cubic}")
try:
    R2ext.<zz> = PolynomialRing(Rz0)
    D_ext.<Dzz> = OreAlgebra(R2ext)
    L_ext = D_ext.zero()
    for i, c in enumerate(L_ode.list()):
        L_ext += R2ext(c.map_coefficients(Rz0)) * Dzz**i
    ind_z0 = L_ext.indicial_polynomial(zz - z0, alpha)
    print(f"Indicial polynomial degree: {ind_z0.degree()}")
    rts = ind_z0.roots(QQbar)
    print(f"Roots:")
    for rt, mult in sorted(rts, key=lambda x: (CC(x[0]).real(), CC(x[0]).imag())):
        print(f"  alpha = {CC(rt):.6f} (mult {mult})")
except Exception as e:
    print(f"Error computing indicial at z0: {e}")
    import traceback
    traceback.print_exc()

# Try numerical approach for indicial at dominant singularity
print("\n--- Numerical indicial at dominant singularity ---")
z0_num = cubic.roots(CC, multiplicities=False)
z0_num.sort(key=lambda x: abs(x.imag()))
z0_real = z0_num[0]
print(f"z0 ≈ {z0_real}")
print(f"1/z0 ≈ {1/z0_real} (Poincaré root μ₀)")

# Check: is 1/z0 close to the known Poincaré root?
poincare = R2(1048576*z**3 - 901120*z**2 + 512*z - 1)
mu_roots = poincare.roots(CC, multiplicities=False)
mu_roots.sort(key=lambda x: -abs(x))
print(f"Poincaré roots: {[str(r) for r in mu_roots]}")
print(f"Dominant μ₀ ≈ {mu_roots[0]}")
print(f"64/μ₀ = 1/λ₀ ≈ {64/mu_roots[0]} (should = z0)")

# Try to find a right factor of order 3 or 4
print("\n" + "="*60)
print("Searching for right factors...")
print("="*60)

for ord in [3, 4, 6]:
    print(f"\nSearching for right factor of order {ord}...")
    try:
        factors = L_ode.right_factors(order=ord)
        if factors:
            print(f"  Found {len(factors)} right factor(s)!")
            for i, f in enumerate(factors):
                print(f"  Factor {i}: order={f.order()}, deg={max(c.degree() for c in f.list())}")
        else:
            print(f"  No right factor of order {ord}")
    except Exception as e:
        print(f"  Error: {e}")

print("\nDone.")
