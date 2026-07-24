#!/usr/bin/env sage
"""
Find the MINIMAL-ORDER ODE annihilating the P2.7 generating function.

The to_D() conversion gave order 21 which is inflated. The true minimal
order should be much lower. Strategy:
1. Compute 150+ exact q_n terms
2. Try guessing a D-finite ODE of increasing order
3. Verify annihilation
"""
from ore_algebra import *
from ore_algebra import guess

# Recurrence coefficients as functions of QQ
Rn.<n> = PolynomialRing(QQ)

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

# Compute q_n with exact rationals
N_TERMS = 160
q = [QQ(0)] * N_TERMS
q[0] = QQ(-215040420000)
q[1] = QQ(-167282265043404) / QQ(905)
q[2] = QQ(-964185327658080) / QQ(6071)

for nn in range(2, N_TERMS - 1):
    q[nn+1] = Bf(nn)/Af(nn)*q[nn] - Cf(nn-1)/Af(nn-1)*q[nn-1] + Df(nn-2)/Af(nn-2)*q[nn-2]

print(f"Computed {N_TERMS} exact q_n values")
print(f"q[0] = {q[0]}")
print(f"q[1] denom digits: {len(str(q[1].denominator()))}")
print(f"q[50] numer digits: {len(str(q[50].numerator()))}")
print(f"q[50] denom digits: {len(str(q[50].denominator()))}")
print(f"q[100] numer digits: {len(str(q[100].numerator()))}")
print(f"q[100] denom digits: {len(str(q[100].denominator()))}")

# Try to guess the minimal ODE
R2.<z> = PolynomialRing(QQ)
D_alg.<Dz> = OreAlgebra(R2)

# Also set up the Sn algebra for guessing recurrence
Qn2 = Rn.fraction_field()
Sn_alg.<Sn> = OreAlgebra(Qn2, 'Sn')

print("\n" + "="*60)
print("Guessing minimal ODE from sequence")
print("="*60)

# Try increasing orders
for n_terms in [80, 120, 150]:
    seq = q[:n_terms]
    for ord in [3, 4, 5, 6, 9, 12]:
        for max_deg in [6, 8, 10, 12, 15, 18]:
            try:
                L = guess(seq, D_alg, order=ord, degree=max_deg)
                deg = max(c.degree() for c in L.list())
                # Verify on remaining terms
                ver_ok = True
                for k in range(n_terms, min(n_terms + 10, N_TERMS)):
                    val = sum(QQ(L.list()[i](z=0)) * q[k] for i in range(ord + 1))
                    if val != 0:
                        ver_ok = False
                        break
                status = "VERIFIED" if ver_ok else "failed verification"
                print(f"  order={ord}, deg={max_deg}, terms={n_terms}: "
                      f"FOUND (actual deg={deg}) [{status}]")
                if ver_ok:
                    print(f"    Leading coeff: {L.leading_coefficient().factor()}")
                    print(f"    Full operator:")
                    for i in range(ord + 1):
                        print(f"      Dz^{i}: {L.list()[i].factor()}")
                    # Done!
                    raise StopIteration
            except StopIteration:
                raise
            except:
                pass
    print(f"  No ODE found with {n_terms} terms for orders 3-12, degrees 6-18")

print("\n" + "="*60)
print("Trying to guess with Sn algebra (recurrence) as cross-check")
print("="*60)

for ord in [3]:
    for max_deg in [12, 15, 18, 20]:
        try:
            L = guess(q[:80], Sn_alg, order=ord, degree=max_deg)
            deg = max(c.degree() for c in L.list())
            print(f"  Recurrence found: order={ord}, degree={deg}")
            print(f"  Leading: {L.leading_coefficient().factor()}")
        except:
            pass

print("\nDone.")
