#!/usr/bin/env sage
"""
Compute the P2.7 ODE via explicit Weyl transform of the recurrence.
Order-3 degree-18 recurrence -> Order-18 degree-3 ODE.
Then compute GCRD with the AESZ inner ODE.
"""
from ore_algebra import *
from ore_algebra import guess

# ===== Get the AESZ inner ODE =====
R.<z> = PolynomialRing(QQ)
D.<Dz> = OreAlgebra(R)

inner = []
for nn in range(100):
    inner.append(sum(binomial(nn,k)^2 * binomial(nn+k,nn) * binomial(nn+2*k,nn) for k in range(nn+1)))

print("Guessing AESZ inner ODE...")
L_a = guess(inner[:100], D, order=4)
print(f"L_a: order={L_a.order()}, degree={max(c.degree() for c in L_a.list())}")

# Check: is this the full annihilator or can we find order 3?
try:
    L_a3 = guess(inner[:100], D, order=3)
    print(f"  Also found order 3: degree={max(c.degree() for c in L_a3.list())}")
except:
    print("  No order-3 ODE exists.")

# ===== Weyl transform of P2.7 recurrence =====
# P2.7: A(n) u_{n+1} - B(n) u_n + C(n-1) u_{n-1} - D(n-2) u_{n-2} = 0
# Wait, the cleared form from ore_analysis2 is:
# c3 Sn^3 + c2 Sn^2 + c1 Sn + c0 = 0
# where c3 = A0 A1 A2, c2 = -A0 A1 B2, c1 = A0 A2 C1, c0 = -A1 A2 D0

# From the recurrence: for OGF F(z) = Σ u_n z^n
# The Weyl transform maps: n -> θ = z d/dz, Sn -> z^{-1}
# But more precisely: if L = Σ_i p_i(n) Sn^i annihilates u_n
# then F(z) is annihilated by: Σ_i z^{-i} p_i(θ+i) (or similar)
# Let me use the ore_algebra conversion

R2.<n> = PolynomialRing(QQ)
A_rec.<Sn> = OreAlgebra(R2)

# Build the P2.7 recurrence operator
A_fn = R2(1024*(2*n+5)^4*(2*n+7)^3*(2*n+9)^3*(946*n^2+6407*n+10860))
B_fn = R2(128*(2*n+7)^3*(2*n+9)^3*(104060*n^6+1745370*n^5+12145238*n^4+44886481*n^3+92943995*n^2+102256019*n+46709052))
C_fn = R2(16*(n+3)^4*(2*n+9)^3*(3784*n^5+57792*n^4+351019*n^3+1059230*n^2+1587211*n+944620))
D_fn = R2((n+3)^4*(n+4)^6*(946*n^2+4515*n+5399))

# The standard clearing: at index n, the recurrence is
# A(n) u_{n+1} = B(n) u_n - C(n-1) u_{n-1} + D(n-2) u_{n-2}
# Shift so that it reads: p3(n) u_{n+3} + p2(n) u_{n+2} + p1(n) u_{n+1} + p0(n) u_n = 0
# Set m = n-2: A(m+2) u_{m+3} = B(m+2) u_{m+2} - C(m+1) u_{m+1} + D(m) u_m
# So: A(n+2) u_{n+3} - B(n+2) u_{n+2} + C(n+1) u_{n+1} - D(n) u_n = 0
# i.e., p3 = A(n+2), p2 = -B(n+2), p1 = C(n+1), p0 = -D(n)

p3 = A_fn(n=n+2)
p2 = -B_fn(n=n+2)
p1 = C_fn(n=n+1)
p0 = -D_fn(n=n)

L_p27_rec = p3*Sn^3 + p2*Sn^2 + p1*Sn + p0
print(f"\nL_p27 (rec): order={L_p27_rec.order()}, degree={max(c.degree() for c in L_p27_rec.list())}")

# Verify it annihilates q_n
q27 = [QQ(-215040420000), QQ(-167282265043404)/905, QQ(-964185327658080)/6071]
def Af(nn):
    nn = QQ(nn)
    return QQ(1024)*(2*nn+5)^4*(2*nn+7)^3*(2*nn+9)^3*(946*nn^2+6407*nn+10860)
def Bf(nn):
    nn = QQ(nn)
    return QQ(128)*(2*nn+7)^3*(2*nn+9)^3*(104060*nn^6+1745370*nn^5+12145238*nn^4+44886481*nn^3+92943995*nn^2+102256019*nn+46709052)
def Cf(nn):
    nn = QQ(nn)
    return QQ(16)*(nn+3)^4*(2*nn+9)^3*(3784*nn^5+57792*nn^4+351019*nn^3+1059230*nn^2+1587211*nn+944620)
def Df(nn):
    nn = QQ(nn)
    return (nn+3)^4*(nn+4)^6*(946*nn^2+4515*nn+5399)

for nn in range(2, 50):
    q27.append(Bf(nn)/Af(nn)*q27[nn] - Cf(nn-1)/Af(nn-1)*q27[nn-1] + Df(nn-2)/Af(nn-2)*q27[nn-2])

# Verify
ok = True
for nn in range(47):
    val = sum(QQ(L_p27_rec.list()[i](n=nn)) * q27[nn+i] for i in range(4))
    if val != 0:
        ok = False
        if nn < 3:
            print(f"  Residual at n={nn}: {float(val)}")
print(f"L_p27 annihilates q: {ok}")

# ===== Use ore_algebra's to_D if available =====
# Try converting the recurrence to ODE
print("\nAttempting Weyl conversion via ore_algebra...")
try:
    # ore_algebra might have a method to_D or annihilator_of_associate
    # or we can use the dfinite module
    from ore_algebra.dfinite_function import *
    print("  DFiniteFunction module available")
except:
    print("  DFiniteFunction not available")

# Try the manual Weyl transform
# For the OGF F(z) = Σ u_n z^n, if Σ p_i(n) u_{n+i} = 0 for all n >= 0,
# then: [Σ_i p_i(θ+i) z^{-i}] F(z) = polynomial (boundary terms)
# Multiply by z^ord to clear:
# [Σ_i p_i(θ+i) z^{ord-i}] F(z) = z^ord × polynomial
# = p_0(θ+0) z^3 + p_1(θ+1) z^2 + p_2(θ+2) z + p_3(θ+3)
# acts on F(z) to give boundary terms

# The annihilator is: p_0(θ) z^3 + p_1(θ+1) z^2 + p_2(θ+2) z + p_3(θ+3)
# Wait, I need to be more careful.
# If the recurrence is p_3 u_{n+3} + p_2 u_{n+2} + p_1 u_{n+1} + p_0 u_n = 0
# and F(z) = Σ u_n z^n, then:
# Σ_n [p_3(n) u_{n+3} z^n + p_2(n) u_{n+2} z^n + p_1(n) u_{n+1} z^n + p_0(n) u_n z^n] = 0
#
# Now: Σ p_i(n) u_{n+i} z^n = Σ p_i(n) u_{n+i} z^{n+i} z^{-i}
# And Σ f(n) u_n z^n = f(θ) F(z) where θ = z d/dz
# So: Σ p_i(n) u_{n+i} z^n = z^{-i} p_i(θ-i) Σ u_m z^m = z^{-i} p_i(θ-i) F(z) + boundary

# Actually no. Let me be more precise.
# Σ_n p(n) a_{n+k} z^n = ???
# Let m = n+k, so n = m-k:
# = Σ_m p(m-k) a_m z^{m-k} = z^{-k} Σ_m p(m-k) a_m z^m
# = z^{-k} p(θ-k) F(z)  where θ acts as z d/dz

# Therefore: Σ_n [p_3(n) a_{n+3} + p_2(n) a_{n+2} + p_1(n) a_{n+1} + p_0(n) a_n] z^n
# = z^{-3} p_3(θ-3) F + z^{-2} p_2(θ-2) F + z^{-1} p_1(θ-1) F + p_0(θ) F = 0
# (modulo boundary terms)

# Multiply by z^3:
# p_3(θ-3) F + z p_2(θ-2) F + z^2 p_1(θ-1) F + z^3 p_0(θ) F = boundary

# But z commutes with scalar multiplication: z · g(θ) = g(θ-1) · z
# Actually, z · θ = (θ+1) · z, so z · g(θ) = g(θ+1) · z
# Therefore: z^k · g(θ) = g(θ+k) · z^k

# Wait no: (z θ) f = z · (z f') = z² f' = θ(z f) - z f? No...
# θ = z d/dz. θ(z f) = z d/dz(z f) = z(f + z f') = z f + z² f' = z f + θ(f) · z
# Hmm, this is getting confused. Let me use the standard identity:
# θ z^k = z^k (θ + k)
# So: z^k θ = (θ + k) z^k - k z^k = θ z^k  ... that's circular
# The correct identity: [θ, z] = z, so θ z = z θ + z = z(θ+1)
# More generally: θ z^k = z^k (θ + k)
# And: z^k θ^m ... this requires normal ordering.
# For a polynomial g(θ): z^k g(θ) = g(θ-k) z^k ??? Let me check:
# z · g(θ) f(z) = z · g(z d/dz) f(z)
# vs g(θ-1) z f(z) = g(z d/dz - 1) (z f(z))
# For g(θ) = θ: z θ f = z · z f' = z² f'
# (θ-1) z f = z d/dz (z f) - z f = z(f + z f') - z f = z² f' ✓

# So: z g(θ) = g(θ-1) z and z^k g(θ) = g(θ-k) z^k

# Therefore: p_3(θ-3) + z p_2(θ-2) + z^2 p_1(θ-1) + z^3 p_0(θ)
# = p_3(θ-3) + p_2(θ-2-1) z + ... wait, the z in front doesn't shift θ further.
#
# NO: the operator is already in normal form. Let me write it as:
# L_ODE = p_3(θ-3) + z p_2(θ-2) + z² p_1(θ-1) + z³ p_0(θ)
# where each term is polynomial_in_theta * power_of_z
# This is in z-θ normal form: z^j come AFTER theta polynomials

# Wait, actually the ordering matters. In ore_algebra, the convention is:
# Dz · z = z · Dz + 1 (in the Dz algebra)
# θ = z Dz, so θ · z = z Dz · z = z(z Dz + 1) = z θ + z

# For the ODE annihilator of F(z):
# From the recurrence, we derived:
# [p_3(θ-3) + z · p_2(θ-2) + z² · p_1(θ-1) + z³ · p_0(θ)] · F(z) = boundary
# The operator in brackets is the ODE annihilator (ignoring boundary).

# To express this in the Dz algebra:
# θ = z Dz
# So p_i(θ-k) is a polynomial in z Dz, which is a polynomial in Dz with z-dependent coefficients.

# For the AESZ inner (degree 6 in n, order 3):
# L_ODE has terms: p_3(θ-3), z·p_2(θ-2), z²·p_1(θ-1), z³·p_0(θ)
# Each p_i has degree 6 in θ, so order 6 in Dz.
# But the highest z power is z³.
# So L_ODE has z-degree 3 and Dz-order 6.

# For P2.7 (degree 18 in n, order 3):
# L_ODE has Dz-order 18 and z-degree 3.

# This is too large to compute symbolically. Let me verify with the AESZ case first.

# ===== Verify Weyl transform for AESZ inner =====
print("\n===== Verifying Weyl transform for AESZ =====")

# Get the AESZ inner recurrence coefficients
L_inner = guess(inner[:100], A_rec)
print(f"L_inner (rec): order={L_inner.order()}, degree={max(c.degree() for c in L_inner.list())}")

# The recurrence: L = c3(n) Sn³ + c2(n) Sn² + c1(n) Sn + c0(n)
# meaning: c3(n) a_{n+3} + c2(n) a_{n+2} + c1(n) a_{n+1} + c0(n) a_n = 0
# Weyl transform ODE: c3(θ-3) + z c2(θ-2) + z² c1(θ-1) + z³ c0(θ)

# But wait, I need to get the recurrence in the form p_i(n) a_{n+i}
# The ore_algebra convention: L = c0(n) + c1(n) Sn + c2(n) Sn² + c3(n) Sn³
# meaning: c0(n) a_n + c1(n) a_{n+1} + c2(n) a_{n+2} + c3(n) a_{n+3} = 0
# So p_0 = c0, p_1 = c1, p_2 = c2, p_3 = c3

coeffs_inner = L_inner.list()
# Weyl transform: L_ODE = p_3(θ-3) + z p_2(θ-2) + z² p_1(θ-1) + z³ p_0(θ)
# But we need to express θ = z Dz in terms of Dz

# Let's just construct the ODE numerically and verify
# For a truncated power series check:
# Apply the ODE operator to A(z) and check it's zero (mod boundary terms)

# Actually, let me try ore_algebra's built-in conversion
try:
    # Some versions of ore_algebra have to_D()
    L_ode_from_rec = L_inner.to_D(D)
    print(f"Converted to ODE: order={L_ode_from_rec.order()}, degree={max(c.degree() for c in L_ode_from_rec.list())}")
except AttributeError as e:
    print(f"to_D not available: {e}")
    # Try alternative methods
    try:
        from ore_algebra.analytic.differential_operator import DifferentialOperator
        print("DifferentialOperator available")
    except:
        pass

# ===== Alternative: use annihilator method =====
# If L_inner annihilates a_n, then the same sequence as an element of the
# D-finite algebra can give us the ODE. Let's use ore_algebra's DFiniteFunction
# or similar.

# Actually, for the AESZ inner, we already guessed L_a (order 4, degree 8).
# Let me check if the order-6 Weyl transform factors through L_a.

# The Weyl ODE should have order 6 (= max degree of recurrence coefficients).
# But guess found order 4. This means either:
# (a) The Weyl ODE factors, and A(z) satisfies a factor of order 4
# (b) The minimal ODE has order 4 (lower than the Weyl transform)

# For P2.7: Weyl gives order 18. But the minimal ODE might be lower.
# Let me try guessing with high terms and more data.

# First, normalize q_n to clear denominators
# q_n has denominators: 1, 905, 6071, ...
# Let d_n = lcm of denominators up to n
# Use a normalized sequence q_n * d_n

# Actually, a smarter approach: compute q_n * N! for large N to clear ALL denominators
# Or use the common denominator structure.

# The P2.7 recurrence generates q_n from q_0, q_1, q_2 which are:
# q_0 = -215040420000 (integer)
# q_1 = -167282265043404/905
# q_2 = -964185327658080/6071

# The denominators 905 and 6071 factor as:
# 905 = 5 × 181
# 6071 = 11 × 521? No: 6071/11 = 551.9...
# 6071 = 7 × 867 + 2... Let me compute
print(f"\nDenominator analysis:")
print(f"  q_1 denominator: {QQ(-167282265043404)/QQ(905)}")
print(f"  q_2 denominator: {QQ(-964185327658080)/QQ(6071)}")
print(f"  905 = {factor(905)}")
print(f"  6071 = {factor(6071)}")

# Try normalizing: multiply by a common factor to make integer sequence
# Find lcm of all denominators for first 50 terms
from sage.all import lcm as lcm_fn
denoms = [q27[nn].denominator() for nn in range(50)]
common = lcm_fn(denoms)
print(f"  Common denominator (first 50): {common}")
print(f"  # digits: {len(str(common))}")

# If the common denominator is manageable, try guessing with the integer sequence
q_int = [QQ(q27[nn] * common) for nn in range(50)]
print(f"  q_int[0] = {q_int[0]}")
print(f"  q_int[1] = {q_int[1]}")

# The ODE for q_n·D (where D is the common denominator) is the SAME ODE
# (just with a different initial condition).
# So let me try guessing the ODE from this integer sequence.

print("\nGuessing P2.7 ODE from normalized integer sequence...")
for ord in range(3, 25):
    try:
        L_q_ode = guess([ZZ(x) for x in q_int[:50]], D, order=ord)
        deg = max(c.degree() for c in L_q_ode.list())
        print(f"  order={ord}: degree={deg}")
        break
    except:
        if ord % 3 == 0 or ord <= 6:
            print(f"  order={ord}: no")

if 'L_q_ode' in dir():
    print(f"\nL_q ODE: order={L_q_ode.order()}, degree={max(c.degree() for c in L_q_ode.list())}")

    # GCRD with AESZ inner ODE
    print(f"\nGCRD(L_a, L_q) at ODE level:")
    try:
        G = L_a.gcrd(L_q_ode)
        print(f"  Order: {G.order()}")
        if G.order() > 0:
            print("  *** NONTRIVIAL GCRD! ***")
            print(f"  Degree: {max(c.degree() for c in G.list())}")
    except Exception as e:
        print(f"  Error: {e}")

    # LCLM
    try:
        M = L_a.lclm(L_q_ode)
        print(f"\nLCLM order: {M.order()}")
        expected = L_a.order() + L_q_ode.order()
        if M.order() < expected:
            drop = expected - M.order()
            print(f"  *** ORDER DROP by {drop}! ***")
    except Exception as e:
        print(f"  LCLM error: {e}")

    # Singular loci
    print(f"\nSingular loci:")
    lc = L_q_ode.list()[-1]
    tc = L_q_ode.list()[0]
    print(f"  Leading: {lc.factor()}")
    print(f"  Trailing: {tc.factor()}")
