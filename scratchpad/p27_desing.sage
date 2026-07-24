#!/usr/bin/env sage
"""
Try to desingularize the order-21 P2.7 ODE and find a smaller annihilator.
Also try LCLM factorization approaches.
"""
from ore_algebra import *

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
print(f"Original ODE: order {L_ode.order()}")

# Try desingularization
print("\n--- Desingularization ---")
try:
    L_desing = L_ode.desingularize()
    print(f"Desingularized: order {L_desing.order()}")
    print(f"Leading: {L_desing.leading_coefficient().factor()}")
except Exception as e:
    print(f"desingularize() failed: {e}")

# Try to find generalized series solutions at z=0
print("\n--- Generalized series at z=0 ---")
try:
    sols = L_ode.generalized_series_solutions(10)
    print(f"Found {len(sols)} solutions")
    for i, s in enumerate(sols[:5]):
        print(f"  sol {i}: leading term ≈ {str(s)[:80]}...")
except Exception as e:
    print(f"generalized_series_solutions failed: {e}")

# Try to reduce the recurrence first, then convert
print("\n--- Trying alternative: reduce recurrence THEN convert ---")

# Factor the polynomial gcd structure
p3 = pnums[3]
p0 = pnums[0]
print(f"p3 = {p3.factor()}")
print(f"p0 = {p0.factor()}")

# Check if the recurrence has common factors we can cancel
print(f"\ngcd(p3, p0) = {gcd(p3, p0).factor()}")
print(f"gcd(p3, p2) = {gcd(pnums[3], pnums[2]).factor()}")
print(f"gcd(all) = {gcd(gcd(pnums[0], pnums[1]), gcd(pnums[2], pnums[3])).factor()}")

# Alternative approach: work with the Euler operator theta = z*Dz
# In terms of theta, n corresponds to theta, and Sn to z^{-1}
# So the recurrence sum p_j(n) u_{n+j} = 0 becomes
# sum p_j(theta) z^{-j} f(z) = 0, i.e.
# z^{-3}[p3(theta) + p2(theta)*z + p1(theta)*z^2 + p0(theta)*z^3] f(z) = 0

print("\n--- Euler operator approach ---")
# theta = z*Dz
# Build the operator in terms of theta and z directly
# p_j(theta) * z^j summed, then multiply by z^3

R3.<z3> = PolynomialRing(QQ)
theta_alg.<theta> = OreAlgebra(R3)

# Actually let's just try to see the structure
# The operator is L = z^3 * sum_{j=0}^3 p_j(theta) * z^{-j}
# = p3(theta) + p2(theta)*z + p1(theta)*z^2 + p0(theta)*z^3

# In terms of Dz = theta/z, we have theta = z*Dz
# Let's compute this directly in the Weyl algebra

print("\n--- Computing via Weyl algebra ---")
# In W = Q<z, Dz>, theta = z*Dz
# p_j(theta) is a polynomial in z*Dz
# The operator is sum_{j=0}^3 z^j * p_j(z*Dz)

# Actually, the correct correspondence for sequences to GF is:
# If sum p_j(n) a_{n+j} = 0, then the GF f(z) = sum a_n z^n satisfies
# sum_{j=0}^r p_j(theta) z^{-j} f = 0
# where theta = z*d/dz
# Multiplying by z^r: sum_{j=0}^r p_j(theta) z^{r-j} f = 0

# For our recurrence (r=3):
# p_0(theta)*z^3 + p_1(theta)*z^2 + p_2(theta)*z + p_3(theta) = 0 (applied to f)

# Let's compute this step by step
# theta^k = (z*Dz)^k in the Weyl algebra

# First, let's just confirm the order by counting
# p_j has degree 18, so p_j(theta) = sum c_{j,k} theta^k, k=0..18
# theta^k = sum Stirling(k,m) z^m Dz^m (in the Weyl algebra)
# So p_j(theta) has ODE order up to 18
# The full operator has ODE order up to 18 + 3 (from the z^{3-j} shifts? no, z is a multiplication, not a differentiation)

# Actually z^j multiplied by theta^k = z^j (z Dz)^k, which is z^{j+k} Dz^k + lower order terms
# So p_j(theta) * z^{3-j} has max order = deg(p_j) = 18
# The full operator is a sum of 4 such terms, each of order ≤ 18
# The leading term (Dz^18 coefficient) of p_j(theta)*z^{3-j} comes from the leading coefficient of p_j times z^{3-j+18}
# So the leading order term in Dz^18 has z-coefficient = sum_j lc(p_j) * z^{21-j}
# This is a polynomial of degree 21 in z times Dz^18

# Wait, that gives order 18, not 21. Let me reconsider.
# theta^k = (z Dz)^k. This is z^k Dz^k + lower order in Dz.
# p_j(theta) = sum_{k=0}^{18} c_{j,k} z^k Dz^k + (lower order Dz terms)
# p_j(theta) * z^{3-j} = sum_k c_{j,k} z^{k+3-j} Dz^k + ...

# Hmm, so the highest Dz order from each p_j(theta) is 18.
# And the leading coefficient of Dz^18 in p_j(theta)*z^{3-j} is lc(p_j) * z^{18+3-j} = lc(p_j) * z^{21-j}

# So the full operator has Dz^18 with coefficient sum_j lc(p_j) * z^{21-j}
# But wait, theta^k ≠ z^k Dz^k in general! Theta^2 = z Dz z Dz = z(z Dz^2 + Dz) = z^2 Dz^2 + z Dz
# So theta^k has leading term z^k Dz^k, but also has lower order terms.

# So the ORDER is indeed at most 18. The to_D() giving 21 means it's NOT computing the minimal annihilator.
# Let me try to compute the Euler representation directly.

# Key insight: the Euler representation should give order AT MOST deg(p_j) = 18.
# In fact, it should be exactly 18 (since deg is the same for all p_j).
# The order-21 from to_D() is an artifact of the implementation.

print("The Euler representation gives ODE of order ≤ max(deg p_j) = 18")
print("The to_D() result of order 21 is NOT minimal.")
print()

# Let's try to compute the minimal operator via the Euler approach directly
# Build the operator sum_{j=0}^3 z^{3-j} * p_j(z*Dz)
# in the Weyl algebra Q[z]<Dz>

# Helper: compute (z*Dz)^k in terms of z^i * Dz^j
# Using the normal ordering: (zDz)^k = sum_{j=0}^k S(k,j) z^j Dz^j
# where S(k,j) are Stirling numbers of the second kind

from sage.combinat.combinat import stirling_number2

MAX_K = 19  # degree 18, so we need k=0..18

# Build z^{3-j} * p_j(zDz) for j=0,1,2,3
# p_j(theta) = sum_{k=0}^{18} p_j[k] * theta^k
# theta^k = sum_{m=0}^k S(k,m) z^m Dz^m
# z^{3-j} * theta^k = sum_m S(k,m) z^{3-j+m} Dz^m

# Full operator: sum_j z^{3-j} p_j(theta) = sum_j sum_k p_j[k] sum_m S(k,m) z^{3-j+m} Dz^m
# = sum_m [sum_j sum_k p_j[k] S(k,m) z^{3-j+m}] Dz^m
# = sum_m z^m [sum_j sum_k p_j[k] S(k,m) z^{3-j}] Dz^m

# So coefficient of Dz^m is: sum_{j=0}^3 z^{3-j+m} * sum_{k=m}^{18} p_j[k] * S(k,m)
# = z^m * sum_{j=0}^3 z^{3-j} * sum_{k=m}^{18} p_j[k] * S(k,m)

# Let's compute this

# First extract polynomial coefficients
pcoeffs = []
for j in range(4):
    pcoeffs.append([QQ(pnums[j][k]) for k in range(pnums[j].degree() + 1)]
                   + [QQ(0)] * (18 - pnums[j].degree()))

# Precompute Stirling numbers
stirl = {}
for k in range(MAX_K):
    for m in range(k + 1):
        stirl[(k, m)] = ZZ(stirling_number2(k, m))

# Build ODE coefficients
# The operator is sum_{m=0}^{18} C_m(z) * Dz^m
# where C_m(z) = sum_{j=0}^3 z^{3-j+m} * (sum_{k=m}^{18} pcoeffs[j][k] * S(k,m))

ode_coeffs_z = {}
for m in range(MAX_K):
    C_m = R2(0)
    for j in range(4):
        inner = QQ(0)
        for k in range(m, min(19, len(pcoeffs[j]))):
            if k < len(pcoeffs[j]):
                inner += pcoeffs[j][k] * stirl.get((k, m), 0)
        if inner != 0:
            C_m += inner * z**(3 - j + m)
    if C_m != 0:
        ode_coeffs_z[m] = C_m

max_order = max(ode_coeffs_z.keys())
print(f"Euler operator has order {max_order}")

# Build the ODE in ore_algebra
L_euler = D_alg.zero()
for m in range(max_order + 1):
    if m in ode_coeffs_z:
        L_euler += ode_coeffs_z[m] * Dz**m

print(f"ODE order: {L_euler.order()}")
deg = max(c.degree() for c in L_euler.list())
print(f"Max coefficient degree: {deg}")
print(f"Leading coefficient: {L_euler.leading_coefficient().factor()}")

# Verify this ODE is correct by checking it annihilates q_n
# If L f(z) = 0 where f(z) = sum q_n z^n, then
# [z^n] L f = 0 for all n
# We can check this directly from the recurrence

# Actually let's verify by checking a few coefficients
def Af(nn):
    nn = QQ(nn)
    return QQ(1024)*(2*nn+5)**4*(2*nn+7)**3*(2*nn+9)**3*(946*nn**2+6407*nn+10860)
def Bf(nn):
    nn = QQ(nn)
    return QQ(128)*(2*nn+7)**3*(2*nn+9)**3*(104060*nn**6+1745370*nn**5+12145238*nn**4+44886481*nn**3+92943995*nn**2+102256019*nn+46709052)
def Cf(nn):
    nn = QQ(nn)
    return QQ(16)*(nn+3)**4*(2*nn+9)**3*(3784*nn**5+57792*nn**4+351019*nn**3+1059230*nn**2+1587211*nn+944620)
def Df(nn):
    nn = QQ(nn)
    return (nn+3)**4*(nn+4)**6*(946*nn**2+4515*nn+5399)

q = [QQ(0)] * 100
q[0] = QQ(-215040420000)
q[1] = QQ(-167282265043404) / QQ(905)
q[2] = QQ(-964185327658080) / QQ(6071)
for nn in range(2, 99):
    q[nn+1] = Bf(nn)/Af(nn)*q[nn] - Cf(nn-1)/Af(nn-1)*q[nn-1] + Df(nn-2)/Af(nn-2)*q[nn-2]

# Check: the action of L_euler on f(z) = sum q_n z^n
# Dz^m f(z) = sum_{n>=m} n!/(n-m)! q_n z^{n-m}
# z^j Dz^m f(z) = sum_{n>=m} n!/(n-m)! q_n z^{n-m+j}
# [z^N] z^j Dz^m f(z) = (N-j+m)!/(N-j)! q_{N-j+m}  [when N-j+m >= m, i.e. N >= j]
# = falling(N-j+m, m) q_{N-j+m}

from sage.arith.misc import falling_factorial

print("\nVerifying Euler ODE annihilates q_n:")
ok = True
for N in range(50):
    val = QQ(0)
    for m in range(max_order + 1):
        if m not in ode_coeffs_z:
            continue
        # C_m(z) = sum_d c_{m,d} z^d
        Cm = ode_coeffs_z[m]
        for d in range(Cm.degree() + 1):
            c_md = QQ(Cm[d])
            if c_md == 0:
                continue
            # [z^N] c_md z^d Dz^m f = c_md * falling(N-d+m, m) * q[N-d+m]
            idx = N - d + m
            if 0 <= idx < len(q) and idx >= m:
                ff = QQ(falling_factorial(idx, m))
                val += c_md * ff * q[idx]
    if val != 0:
        ok = False
        if N < 5:
            print(f"  N={N}: NONZERO residual")
if ok:
    print("  *** ALL ZERO — Euler ODE is correct! ***")
else:
    print("  SOME NONZERO (Euler ODE may need correction)")

# Factor analysis
print(f"\nLeading coefficient factorization:")
print(f"  {L_euler.leading_coefficient().factor()}")
print(f"\nTrailing coefficient (Dz^0):")
if 0 in ode_coeffs_z:
    print(f"  {ode_coeffs_z[0].factor()}")

# Try to get a right factor
print("\n--- Searching for right factors ---")
for r in [3, 4]:
    try:
        rf = L_euler.right_factors(order=r)
        if rf:
            print(f"Right factor of order {r}: FOUND!")
            for f in rf:
                print(f"  order={f.order()}, deg={max(c.degree() for c in f.list())}")
                print(f"  leading: {f.leading_coefficient().factor()}")
        else:
            print(f"No right factor of order {r}")
    except Exception as e:
        print(f"Error searching for order-{r} right factors: {e}")

print("\nDone.")
