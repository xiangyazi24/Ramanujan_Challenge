#!/usr/bin/env sage
"""
Verify Q5117's key finding: M† = M₁₆ · Q₂₀₉(θ)/946
Build Mellin adjoint and check the Q₂₀₉(θ) right factor.
"""
from ore_algebra import *

Rn.<n> = PolynomialRing(QQ)
Rt.<t> = PolynomialRing(QQ)
Qt = Rt.fraction_field()
Dt_alg.<Dt> = OreAlgebra(Qt)

# P2.7 primitive polynomial coefficients
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

# Build Mellin adjoint: M† = sum_j t^j * p_j(-theta-j-1)
# where theta = t*Dt

# Helper: compute (-theta-j-1)^k as an operator in Dt_alg
# theta = t*Dt
# (-theta-j-1)^k = (-(t*Dt) - j - 1)^k

print("Building Mellin adjoint M† = Σ_j t^j p_j(-θ-j-1)...")

M_dual = Dt_alg.zero()
for j in range(4):
    base_op = -(Qt(t) * Dt + QQ(j + 1))  # -(theta + j + 1) = -theta - j - 1
    # Build p_j(base_op) = sum_k pnums[j][k] * base_op^k
    pj_eval = Dt_alg.zero()
    power = Dt_alg.one()
    for k in range(pnums[j].degree() + 1):
        c_k = QQ(pnums[j][k])
        if c_k != 0:
            pj_eval += c_k * power
        if k < pnums[j].degree():
            power = power * base_op
    M_dual += Qt(t)**j * pj_eval
    print(f"  j={j} done")

print(f"\nM† order: {M_dual.order()}")

# Check leading coefficient
lc = Rt(M_dual.leading_coefficient().numerator())
print(f"Leading coeff degree: {lc.degree()}")
print(f"Leading coeff factor: {lc.factor()}")

# Now check for Q₂₀₉(θ) right factor
# Q₂₀₉(θ) = 946*θ² - 2623*θ + 1830
# In terms of Dt: θ = t*Dt, θ² = t*Dt*t*Dt = t*(t*Dt+1)*Dt = t²*Dt² + t*Dt
# So Q₂₀₉(θ) = 946*(t²*Dt² + t*Dt) - 2623*(t*Dt) + 1830
#             = 946*t²*Dt² + (946-2623)*t*Dt + 1830
#             = 946*t²*Dt² - 1677*t*Dt + 1830

Q209_theta = 946*Qt(t)**2*Dt**2 - 1677*Qt(t)*Dt + 1830
print(f"\nQ₂₀₉(θ) = {Q209_theta}")
print(f"Q₂₀₉(θ) order: {Q209_theta.order()}")

# Check if Q₂₀₉(θ) right-divides M†
print("\nChecking if Q₂₀₉(θ)/946 right-divides M†...")
Q209_normalized = Qt(t)**2*Dt**2 - QQ(1677)/QQ(946)*Qt(t)*Dt + QQ(1830)/QQ(946)
try:
    Q_rem, R_rem = M_dual.quo_rem(Q209_normalized)
    if R_rem == 0:
        print("*** YES — Q₂₀₉(θ)/946 is an exact right factor! ***")
        print(f"M₁₆ = M† / (Q₂₀₉(θ)/946)")
        print(f"M₁₆ order: {Q_rem.order()}")
    else:
        print(f"No — remainder order: {R_rem.order()}")
        # Try with the unnormalized version
        Q_rem2, R_rem2 = M_dual.quo_rem(Q209_theta)
        if R_rem2 == 0:
            print("But Q₂₀₉(θ) (unnormalized) IS a right factor!")
            print(f"Quotient order: {Q_rem2.order()}")
        else:
            print(f"Unnormalized remainder order: {R_rem2.order()}")
except Exception as e:
    print(f"Division error: {e}")
    import traceback
    traceback.print_exc()

# Alternative: try right_factors
print("\n--- Searching for Euler right factors ---")
try:
    rf2 = M_dual.right_factors(order=2)
    if rf2:
        print(f"Found {len(rf2)} right factor(s) of order 2:")
        for f in rf2:
            print(f"  {f}")
    else:
        print("No Euler right factor of order 2 found by right_factors()")
except Exception as e:
    print(f"right_factors error: {e}")

# Try: verify Q209(theta) kills t^rho for rho = roots of Q209
print("\n--- Verifying Q₂₀₉(θ) has Euler solutions ---")
RX.<X> = PolynomialRing(QQ)
Q209_poly = 946*X**2 - 2623*X + 1830
disc = 2623**2 - 4*946*1830
print(f"Q₂₀₉ discriminant: {disc}")
print(f"Factor: {ZZ(abs(disc)).factor()}")
rho_roots = Q209_poly.roots(QQbar)
for rt, mult in rho_roots:
    print(f"  ρ = {CC(rt):.10f} (mult {mult})")

print("\nDone.")
