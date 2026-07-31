#!/usr/bin/env sage
"""
Factor the order-16 quotient M₁₆ from the Mellin adjoint.
M† = M₁₆ · Q₂₀₉(θ)/946 (verified)
Now: what is the structure of M₁₆?
"""
from ore_algebra import *

Rn.<n> = PolynomialRing(QQ)
Rt.<t> = PolynomialRing(QQ)
Qt = Rt.fraction_field()
Dt_alg.<Dt> = OreAlgebra(Qt)

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

# Build Mellin adjoint
M_dual = Dt_alg.zero()
for j in range(4):
    base_op = -(Qt(t) * Dt + QQ(j + 1))
    pj_eval = Dt_alg.zero()
    power = Dt_alg.one()
    for k in range(pnums[j].degree() + 1):
        c_k = QQ(pnums[j][k])
        if c_k != 0:
            pj_eval += c_k * power
        if k < pnums[j].degree():
            power = power * base_op
    M_dual += Qt(t)**j * pj_eval

print(f"M† order: {M_dual.order()}")

# Divide by Q209(theta)/946
Q209_norm = Qt(t)**2*Dt**2 - QQ(1677)/QQ(946)*Qt(t)*Dt + QQ(1830)/QQ(946)
M16, rem = M_dual.quo_rem(Q209_norm)
assert rem == 0, "Q209 factor failed!"
print(f"M₁₆ order: {M16.order()}")

# Analyze M₁₆ structure
lc16 = M16.leading_coefficient()
print(f"\nM₁₆ leading coefficient:")
# Convert to polynomial if possible
try:
    lc16_num = Rt(lc16.numerator())
    lc16_den = Rt(lc16.denominator())
    print(f"  numerator degree: {lc16_num.degree()}")
    print(f"  denominator degree: {lc16_den.degree()}")
    print(f"  numerator factored: {lc16_num.factor()}")
except:
    print(f"  {lc16}")

tc16 = M16.list()[0]
print(f"\nM₁₆ trailing coefficient (Dt^0):")
try:
    tc16_num = Rt(tc16.numerator())
    tc16_den = Rt(tc16.denominator())
    print(f"  numerator degree: {tc16_num.degree()}")
    print(f"  denominator degree: {tc16_den.degree()}")
    print(f"  numerator factored: {tc16_num.factor()}")
except:
    print(f"  {tc16}")

# Try to find additional Euler right factors (polynomials in theta = t*Dt alone)
print("\n--- Checking for additional Euler factors ---")

# Check theta-a for various rational a
for a_num in range(-5, 20):
    for a_den in [1, 2]:
        a = QQ(a_num) / QQ(a_den)
        euler1 = Qt(t)*Dt - a
        Q_test, R_test = M16.quo_rem(euler1)
        if R_test == 0:
            print(f"  θ - {a} is a right factor!")

# Check quadratic Euler factors with integer coefficients
print("Checking quadratic Euler factors...")
for c1 in range(-20, 20):
    for c0 in range(-20, 20):
        euler2 = Qt(t)**2*Dt**2 + QQ(c1)*Qt(t)*Dt + QQ(c0)
        try:
            Q_test, R_test = M16.quo_rem(euler2)
            if R_test == 0:
                print(f"  θ² + {c1}θ + {c0} is a right factor!")
        except:
            pass

# Check for right factors of small order
print("\n--- Checking for right factors of order 1-4 ---")
for r in [1, 2, 3, 4]:
    print(f"Order {r}:")
    try:
        rfs = M16.right_factors(order=r)
        if rfs:
            for f in rfs:
                print(f"  Found: order={f.order()}")
                flc = f.leading_coefficient()
                print(f"    leading: {flc}")
        else:
            print(f"  None found")
    except Exception as e:
        print(f"  Error: {e}")

# Singular points of M₁₆
print("\n--- Singular points of M₁₆ ---")
# Get the leading coefficient
try:
    lc_num = Rt(lc16.numerator())
    roots = lc_num.roots(QQbar)
    print(f"Finite singular points ({len(roots)}):")
    for rt, mult in roots:
        print(f"  t = {CC(rt):.10f} (mult {mult})")
except Exception as e:
    print(f"Error: {e}")

# Indicial polynomial at t=0 for M₁₆
print("\n--- Indicial polynomial of M₁₆ at t=0 ---")
try:
    Ra.<alpha> = PolynomialRing(QQ)
    ind_0 = M16.indicial_polynomial(t, alpha)
    print(f"Degree: {ind_0.degree()}")
    rts = ind_0.roots(QQbar)
    print("Roots:")
    for rt, mult in sorted(rts, key=lambda x: CC(x[0]).real()):
        print(f"  α = {CC(rt):.6f} (mult {mult})")
except Exception as e:
    print(f"Error: {e}")

# Check if M₁₆ has a symmetric structure
# The coboundary h_n involves half-integer Pochhammer symbols
# In the Mellin variable, these correspond to t^{1/2} factors
# Check if M₁₆ has the substitution t -> t² as a symmetry

print("\n--- Checking t -> t² pullback structure ---")
# If M₁₆ = pullback of a lower-order operator through t -> t²,
# then M₁₆ would factor as P(t²,t²Dt²) where Dt² = 1/(2t) Dt
# This would halve the effective order

# For now, just check the parity of the leading coefficient
print("Leading coefficient parity check:")
try:
    lc_poly = Rt(lc16.numerator())
    even_part = Rt(sum(lc_poly[2*k]*t**(2*k) for k in range((lc_poly.degree()//2)+1)))
    odd_part = Rt(sum(lc_poly[2*k+1]*t**(2*k+1) for k in range(lc_poly.degree()//2+1)))
    print(f"  Even part nonzero: {even_part != 0}")
    print(f"  Odd part nonzero: {odd_part != 0}")
except Exception as e:
    print(f"  Error: {e}")

# Also try: GCRD with the AESZ #209 Mellin dual
print("\n--- Building AESZ #209 Mellin dual ---")
# AESZ inner recurrence: a_3(n)u_{n+3} + a_2(n)u_{n+2} + a_1(n)u_{n+1} + a_0(n)u_n = 0
# From proof.tex: c_3(n) = (n+3)^4(946n^2+3053n+2475), c_0(n) = -(n+1)^2(2n+1)(2n+3)(946n^2+4945n+6474)
# Need the full recurrence. Let's guess it from the sequence.

inner = []
for nn in range(60):
    inner.append(ZZ(sum(binomial(nn,k)**2 * binomial(nn+k,nn) * binomial(nn+2*k,nn)
                        for k in range(nn+1))))
inner_full = [ZZ(binomial(2*nn, nn)) * inner[nn] for nn in range(60)]

Sn_alg = OreAlgebra(Qn, 'Sn')
Sn_gen = Sn_alg.gen()

from ore_algebra import guess
try:
    L_aesz_rec = guess(inner_full[:50], Sn_alg, order=3)
    print(f"AESZ recurrence: order {L_aesz_rec.order()}")
    aesz_lc = L_aesz_rec.leading_coefficient()
    print(f"Leading: {Rn(aesz_lc.numerator()).factor()}")
    aesz_tc = L_aesz_rec.list()[0]
    print(f"Trailing: {Rn(aesz_tc.numerator()).factor()}")

    # Build AESZ Mellin dual
    aesz_pnums = []
    aesz_coeffs = [Qn(L_aesz_rec[i]) for i in range(L_aesz_rec.order() + 1)]
    aesz_den = Rn.one()
    for c in aesz_coeffs:
        aesz_den = lcm(aesz_den, Rn(c.denominator()))
    aesz_nums = [Rn(aesz_den * c) for c in aesz_coeffs]
    aesz_g = aesz_nums[0]
    for c in aesz_nums[1:]:
        aesz_g = gcd(aesz_g, c)
    aesz_pnums = [c // aesz_g for c in aesz_nums]

    aesz_deg = max(p.degree() for p in aesz_pnums)
    print(f"AESZ primitive degree: {aesz_deg}")

    # Build AESZ Mellin adjoint
    M_aesz = Dt_alg.zero()
    for j in range(len(aesz_pnums)):
        base_op = -(Qt(t) * Dt + QQ(j + 1))
        pj_eval = Dt_alg.zero()
        power = Dt_alg.one()
        for k in range(aesz_pnums[j].degree() + 1):
            c_k = QQ(aesz_pnums[j][k])
            if c_k != 0:
                pj_eval += c_k * power
            if k < aesz_pnums[j].degree():
                power = power * base_op
        M_aesz += Qt(t)**j * pj_eval

    print(f"AESZ Mellin dual order: {M_aesz.order()}")

    # Check if Q209(theta) right-divides AESZ Mellin dual too!
    Q_test_aesz, R_test_aesz = M_aesz.quo_rem(Q209_norm)
    if R_test_aesz == 0:
        print("*** Q₂₀₉(θ) ALSO right-divides AESZ Mellin dual! ***")
        print(f"AESZ quotient order: {Q_test_aesz.order()}")
    else:
        print("Q₂₀₉(θ) does NOT right-divide AESZ Mellin dual")
        # Check for any order-2 Euler right factor
        for a_num in range(-5, 10):
            for a_den in [1, 2, 3]:
                a = QQ(a_num) / QQ(a_den)
                for b_num in range(-5, 10):
                    for b_den in [1, 2, 3]:
                        b = QQ(b_num) / QQ(b_den)
                        euler2 = Qt(t)**2*Dt**2 + a*Qt(t)*Dt + b
                        try:
                            _, R_test = M_aesz.quo_rem(euler2)
                            if R_test == 0:
                                print(f"  AESZ has Euler right factor: θ² + {a}θ + {b}")
                        except:
                            pass

    # GCRD of M₁₆ and AESZ Mellin dual
    print("\nGCRD(M₁₆, M_aesz):")
    try:
        G = M16.gcrd(M_aesz)
        print(f"  Order: {G.order()}")
        if G.order() > 0:
            print("  *** NONTRIVIAL GCRD! ***")
    except Exception as e:
        print(f"  Error: {e}")

except Exception as e:
    print(f"AESZ computation error: {e}")
    import traceback
    traceback.print_exc()

print("\nDone.")
