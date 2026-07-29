#!/usr/bin/env sage
"""
Build the P2.7 ODE via direct Euler operator construction.

For recurrence sum_{j=0}^3 p_j(n) u_{n+j} = 0,
the GF f(z) = sum u_n z^n satisfies:
  sum_{j=0}^3 z^{3-j} p_j(z*Dz) f(z) = 0

Build this directly in Q[z]<Dz>.
"""
from ore_algebra import *
from sage.combinat.combinat import stirling_number2

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

# All have degree 18
DEG = 18
for j in range(4):
    assert pnums[j].degree() == DEG, f"p{j} degree {pnums[j].degree()}"

# Build the ODE: sum_{j=0}^3 z^{3-j} * p_j(theta) where theta = z*Dz
# theta^k in normal ordering: (z*Dz)^k = sum_{m=0}^k S(k,m) z^m Dz^m
# So z^{3-j} * p_j(theta) = z^{3-j} * sum_{k=0}^{18} p_j[k] * sum_{m=0}^k S(k,m) z^m Dz^m
#   = sum_{k=0}^{18} p_j[k] * sum_{m=0}^k S(k,m) z^{3-j+m} Dz^m

# Collect by Dz^m:
# coeff of Dz^m = sum_{j=0}^3 sum_{k=m}^{18} p_j[k] * S(k,m) * z^{3-j+m}

R2.<z> = PolynomialRing(QQ)
D_alg.<Dz> = OreAlgebra(R2)

# Precompute Stirling numbers of the second kind
MAX_K = DEG + 1
stirl = {}
for k in range(MAX_K):
    for m in range(k + 1):
        stirl[(k, m)] = ZZ(stirling_number2(k, m))

# Extract polynomial coefficients
pcoeffs = []
for j in range(4):
    pc = []
    for k in range(DEG + 1):
        pc.append(QQ(pnums[j][k]))
    pcoeffs.append(pc)

# Build ODE coefficients for each Dz^m
print("Building Euler ODE operator...")
ode_terms = {}  # m -> polynomial in z (coefficient of Dz^m)

for m in range(DEG + 1):
    C_m = R2(0)
    for j in range(4):
        for k in range(m, DEG + 1):
            s = stirl.get((k, m), ZZ(0))
            if s == 0 or pcoeffs[j][k] == 0:
                continue
            C_m += pcoeffs[j][k] * s * z**(3 - j + m)
    if C_m != 0:
        ode_terms[m] = C_m

max_order = max(ode_terms.keys())
print(f"Max ODE order: {max_order}")

# Build the operator
L_euler = D_alg.zero()
for m in sorted(ode_terms.keys()):
    L_euler += ode_terms[m] * Dz**m

print(f"ODE operator: order {L_euler.order()}")
deg_ode = max(c.degree() for c in L_euler.list())
print(f"Max coefficient degree: {deg_ode}")

# Check leading coefficient
print(f"\nLeading coefficient (of Dz^{max_order}):")
lc = L_euler.leading_coefficient()
print(f"  degree: {lc.degree()}")
print(f"  factored: {lc.factor()}")

# If leading coefficient has z factors, we can reduce!
# Extract z power from leading coefficient
v = 0
while lc[v] == 0:
    v += 1
print(f"  z-valuation of leading coeff: {v}")

if v > 0:
    # The operator has a z^v factor in the leading term
    # This means z=0 is an irregular singular point or there are apparent solutions
    print(f"  Leading coefficient has z^{v} factor")

# Trailing coefficient
tc = L_euler.list()[0]
print(f"\nTrailing coefficient (of Dz^0):")
print(f"  degree: {tc.degree()}")
print(f"  factored: {tc.factor()}")

# Verify by checking annihilation on q_n
print("\n--- Verifying on q_n ---")

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

q = [QQ(0)] * 80
q[0] = QQ(-215040420000)
q[1] = QQ(-167282265043404) / QQ(905)
q[2] = QQ(-964185327658080) / QQ(6071)
for nn in range(2, 79):
    q[nn+1] = Bf(nn)/Af(nn)*q[nn] - Cf(nn-1)/Af(nn-1)*q[nn-1] + Df(nn-2)/Af(nn-2)*q[nn-2]

from sage.arith.misc import falling_factorial

ok = True
for N in range(40):
    val = QQ(0)
    for m in sorted(ode_terms.keys()):
        Cm = ode_terms[m]
        for d in range(Cm.degree() + 1):
            c_md = QQ(Cm[d])
            if c_md == 0:
                continue
            idx = N - d + m
            if 0 <= idx < len(q) and idx >= m:
                ff = QQ(falling_factorial(idx, m))
                val += c_md * ff * q[idx]
    if val != 0:
        ok = False
        if N < 3:
            print(f"  N={N}: NONZERO")

if ok:
    print("  ALL ZERO — Euler ODE verified!")
else:
    # Maybe we need to account for boundary terms at n=0,1,2
    print("  Some nonzero residuals (checking if only at boundaries)")
    # The Euler relation may not hold for n < 3 due to initial conditions
    ok2 = True
    for N in range(3, 40):
        val = QQ(0)
        for m in sorted(ode_terms.keys()):
            Cm = ode_terms[m]
            for d in range(Cm.degree() + 1):
                c_md = QQ(Cm[d])
                if c_md == 0:
                    continue
                idx = N - d + m
                if 0 <= idx < len(q) and idx >= m:
                    ff = QQ(falling_factorial(idx, m))
                    val += c_md * ff * q[idx]
        if val != 0:
            ok2 = False
            if N < 10:
                print(f"    N={N}: NONZERO")
    if ok2:
        print("  Zero for N >= 3 (boundary effect at N < 3)")
    else:
        print("  Still nonzero — Euler approach needs checking")

# Compare with to_D() result
print("\n--- Comparing with to_D() result ---")
L_prim = sum((Qn(pnums[i]) * Sn**i for i in range(4)), A_rec.zero())
L_toD = L_prim.to_D(D_alg)
print(f"to_D order: {L_toD.order()}")
print(f"Euler order: {L_euler.order()}")

# Check if one divides the other
try:
    G = L_toD.gcrd(L_euler)
    print(f"GCRD order: {G.order()}")
except Exception as e:
    print(f"GCRD computation: {e}")

# Check if L_toD = L_euler * (something)
try:
    Q, R = L_toD.quo_rem(L_euler)
    if R == 0:
        print(f"to_D = Euler * (order {Q.order()} quotient) — Euler is a RIGHT FACTOR!")
    else:
        print(f"Euler does NOT right-divide to_D (nonzero remainder)")
except Exception as e:
    print(f"Division: {e}")

# Try left division too
try:
    Q2, R2_rem = L_toD.quo_rem(L_euler)
    print(f"Left division remainder order: {R2_rem.order() if R2_rem != 0 else 'zero'}")
except:
    pass

print("\nDone.")
