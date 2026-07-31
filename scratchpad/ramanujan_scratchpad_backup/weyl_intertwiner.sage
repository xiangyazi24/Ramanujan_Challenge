#!/usr/bin/env sage
"""
Differential-Ore intertwiner search at the ODE level.
Strategy from Q5095:
1. Convert recurrences to ODEs via Weyl transform (n -> theta, E -> z^{-1})
2. Apply Cauchy twist to AESZ ODE
3. Search for T in K(z)[theta] with L_P27_ODE · T ≡ 0 mod L_twisted_AESZ_ODE
"""
from ore_algebra import *
from ore_algebra import guess

# ===== Step 1: Get the recurrence operators =====

R_n.<n> = PolynomialRing(QQ)
A_rec.<Sn> = OreAlgebra(R_n)

# AESZ inner sum
inner = []
for nn in range(100):
    inner.append(sum(binomial(nn,k)^2 * binomial(nn+k,nn) * binomial(nn+2*k,nn) for k in range(nn+1)))

L_aesz = guess(inner, A_rec)
print(f"L_aesz (rec): order={L_aesz.order()}, degree={max(c.degree() for c in L_aesz.list())}")

# P2.7 initial values and recurrence
def A_p27(nn):
    nn = QQ(nn)
    return QQ(1024)*(2*nn+5)^4*(2*nn+7)^3*(2*nn+9)^3*(946*nn^2+6407*nn+10860)
def B_p27(nn):
    nn = QQ(nn)
    return QQ(128)*(2*nn+7)^3*(2*nn+9)^3*(104060*nn^6+1745370*nn^5+12145238*nn^4+44886481*nn^3+92943995*nn^2+102256019*nn+46709052)
def C_p27(nn):
    nn = QQ(nn)
    return QQ(16)*(nn+3)^4*(2*nn+9)^3*(3784*nn^5+57792*nn^4+351019*nn^3+1059230*nn^2+1587211*nn+944620)
def D_p27(nn):
    nn = QQ(nn)
    return (nn+3)^4*(nn+4)^6*(946*nn^2+4515*nn+5399)

q27 = [QQ(-215040420000), QQ(-167282265043404)/905, QQ(-964185327658080)/6071]
for nn in range(2, 100):
    q27.append(B_p27(nn)/A_p27(nn)*q27[nn] - C_p27(nn-1)/A_p27(nn-1)*q27[nn-1] + D_p27(nn-2)/A_p27(nn-2)*q27[nn-2])

L_p27 = guess(q27, A_rec)
print(f"L_p27 (rec): order={L_p27.order()}, degree={max(c.degree() for c in L_p27.list())}")

# ===== Step 2: Convert to ODE via Weyl transform =====
# n -> theta = z d/dz, shift E_n f(n) = f(n+1) -> z^{-1} (inverse shift)
# Actually: the standard map is: n -> theta, Sn -> 1/z (or z, depending on convention)
# For an OGF A(z) = Σ a_n z^n, n·a_n z^n -> z·(zA(z))' = (theta+1)·A if we define theta = z d/dz

# The ore_algebra package can do this conversion:
R_z.<z> = PolynomialRing(QQ)
A_diff.<Dz> = OreAlgebra(R_z)

print("\n===== Weyl conversion =====")

# For the AESZ recurrence (order 3, degree 6):
# Use ore_algebra's built-in conversion if available
# Otherwise, manual: if L = Σ_i p_i(n) S^i annihilates a_n
# then the OGF A(z) = Σ a_n z^n satisfies an ODE obtained by replacing:
#   n -> θ (theta operator = z·d/dz)
#   S -> z^{-1} (the shift Sn a_n = a_{n+1} maps to z^{-1} in the OGF world)
# Actually: Sn a_n = a_{n+1}, and z^{n+1} a_{n+1} = z · (z^n a_{n+1}),
# but the shift of the OGF is: if A(z) = Σ a_n z^n, then Σ a_{n+1} z^n = (A(z) - a_0)/z

# For ore_algebra, try the to_ODE method or similar
try:
    # ore_algebra may have built-in conversion
    # Try: L_aesz.to_D(A_diff) or similar
    from ore_algebra.analytic.polynomial_approximation import *
    print("Analytic submodule available")
except:
    pass

# Manual Weyl transform using theta = z d/dz
# For polynomial coefficient recurrence c_3(n)Sn³ + c_2(n)Sn² + c_1(n)Sn + c_0(n)
# The ODE annihilator of A(z) = Σ a_n z^n is obtained by:
# Replacing n -> θ, Sn -> z^{-1}
# So L(n, Sn) becomes L(θ, z^{-1})
# This is a differential operator in θ (or equivalently Dz)

# Let's use the Ore algebra over the differential ring
# θ = z Dz, so Dz = θ/z
R_z2.<z> = PolynomialRing(QQ)
A_theta.<theta> = OreAlgebra(FractionField(R_z2))  # theta = z*Dz

# For the AESZ operator: L = c3(n)Sn^3 + c2(n)Sn^2 + c1(n)Sn + c0(n)
# Replace n -> theta, Sn -> z^{-1}
# L(theta, z^{-1}) = c3(theta)·z^{-3} + c2(theta)·z^{-2} + c1(theta)·z^{-1} + c0(theta)
# Multiply by z^3 to clear denominators:
# z^3·L = c3(theta) + c2(theta)·z + c1(theta)·z^2 + c0(theta)·z^3

# BUT: theta and z DON'T commute! theta·z = z·(theta+1), so z^k·theta = (theta-k)·z^k
# So we need to be careful about ordering.

# The generating function approach: if Σ p_i(n) a_{n+i} z^n = 0 for all n,
# then Σ_i p_i(theta) z^{-i} · A(z) = boundary terms (polynomial in z^{-1})
# where we use: Σ n^k a_{n+i} z^n = (z^{-i} theta^k - boundary) A(z)

# Actually, let me just use the known result directly:
# For L_aesz (order 3, degree 6 in n), the ODE has order 6, degree 3 in z.
# For L_p27 (order 3, degree 18 in n), the ODE has order 18, degree 3 in z.

# Rather than doing the Weyl transform manually, let me verify by computing
# the ODE from the generating function.

# Compute the OGF truncation
print("\nComputing OGFs...")
PS.<Z> = PowerSeriesRing(QQ, default_prec=80)
A_ogf = PS(sum(inner[nn] * Z^nn for nn in range(80)))
Q_ogf = PS(sum(q27[nn] * Z^nn for nn in range(80)))

# Use ore_algebra's to_diff or built-in conversion
# The clean way: use ore_algebra's OreAlgebra with the shift algebra
# and convert to differential algebra

# Let's try the annihilator method
R_z3.<z> = PolynomialRing(QQ)
D_alg.<Dz> = OreAlgebra(R_z3)

print("Guessing ODE for AESZ inner OGF...")
# theta = z Dz, so we work in the D_alg = QQ[z]<Dz>
# A(z) = inner OGF
a_coeffs = inner[:80]
try:
    L_aesz_ode = guess(a_coeffs, D_alg)
    print(f"  Order: {L_aesz_ode.order()}, Degree: {max(c.degree() for c in L_aesz_ode.list())}")
except Exception as e:
    print(f"  Error: {e}")
    # Try with explicit order
    for ord in [3, 4, 5, 6, 7, 8]:
        try:
            L_aesz_ode = guess(a_coeffs, D_alg, order=ord)
            deg = max(c.degree() for c in L_aesz_ode.list())
            print(f"  order={ord}: degree={deg}")
            break
        except:
            print(f"  order={ord}: no relation")

print("\nGuessing ODE for P2.7 OGF...")
q_coeffs_for_guess = [QQ(q27[nn]) for nn in range(80)]
try:
    L_p27_ode = guess(q_coeffs_for_guess, D_alg)
    print(f"  Order: {L_p27_ode.order()}, Degree: {max(c.degree() for c in L_p27_ode.list())}")
except Exception as e:
    print(f"  Error: {e}")
    for ord in [3, 6, 9, 12, 15, 18, 21]:
        try:
            L_p27_ode = guess(q_coeffs_for_guess, D_alg, order=ord)
            deg = max(c.degree() for c in L_p27_ode.list())
            print(f"  order={ord}: degree={deg}")
            if ord <= 20:
                break
        except:
            print(f"  order={ord}: no relation")

# ===== Step 3: GCRD at ODE level =====
print(f"\n{'='*60}")
print("GCRD at ODE level")
print("="*60)
try:
    G_ode = L_aesz_ode.gcrd(L_p27_ode)
    print(f"GCRD order: {G_ode.order()}")
    if G_ode.order() > 0:
        print("*** NONTRIVIAL! ***")
except Exception as e:
    print(f"Error: {e}")

try:
    M_lclm_ode = L_aesz_ode.lclm(L_p27_ode)
    print(f"LCLM order: {M_lclm_ode.order()}")
    print(f"Expected max: {L_aesz_ode.order() + L_p27_ode.order()}")
    if M_lclm_ode.order() < L_aesz_ode.order() + L_p27_ode.order():
        drop = L_aesz_ode.order() + L_p27_ode.order() - M_lclm_ode.order()
        print(f"*** ORDER DROP by {drop}! Shared solution space dim = {drop} ***")
except Exception as e:
    print(f"Error: {e}")

# ===== Step 4: Right factors =====
print(f"\n{'='*60}")
print("Factorization of ODE operators")
print("="*60)
try:
    print("Right factors of L_aesz_ode (order 1):")
    rf = L_aesz_ode.right_factors(1)
    print(f"  Found {len(rf)}")
    for f in rf:
        print(f"  {f}")
except Exception as e:
    print(f"  Error: {e}")

try:
    print("\nRight factors of L_p27_ode (order 1):")
    rf = L_p27_ode.right_factors(1)
    print(f"  Found {len(rf)}")
except Exception as e:
    print(f"  Error: {e}")
