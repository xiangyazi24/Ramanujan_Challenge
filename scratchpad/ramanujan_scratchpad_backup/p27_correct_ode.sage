#!/usr/bin/env sage
"""
P2.7 ODE via CORRECT operator (Q5107 monic form).

The STANDARD form A(n+2)u(n+3)-B(n+2)u(n+2)+C(n+1)u(n+1)-D(n)u(n) = 0
does NOT annihilate the P2.7 sequences (residuals ~10^25).

The CORRECT monic form (proof.tex eq 2.1, Q5107):
  u(n+1) = B(n)/A(n)*u(n) - C(n-1)/A(n-1)*u(n-1) + D(n-2)/A(n-2)*u(n-2)

As shift operator at index n:
  L = Sn^3 - B(n+2)/A(n+2)*Sn^2 + C(n+1)/A(n+1)*Sn - D(n)/A(n)

Clearing denominators (multiply by A(n)*A(n+1)*A(n+2)):
  p3 = A(n)*A(n+1)*A(n+2)
  p2 = -B(n+2)*A(n)*A(n+1)
  p1 = C(n+1)*A(n)*A(n+2)    # NOTE: A(n+2), not A(n+1)!
  p0 = -D(n)*A(n+1)*A(n+2)
"""
from ore_algebra import *
from ore_algebra import guess

# ============================================================
# Coefficient polynomials
# ============================================================
Rn.<n> = PolynomialRing(QQ)

A_fn = Rn(1024*(2*n+5)^4*(2*n+7)^3*(2*n+9)^3*(946*n^2+6407*n+10860))
B_fn = Rn(128*(2*n+7)^3*(2*n+9)^3*(104060*n^6+1745370*n^5+12145238*n^4+44886481*n^3+92943995*n^2+102256019*n+46709052))
C_fn = Rn(16*(n+3)^4*(2*n+9)^3*(3784*n^5+57792*n^4+351019*n^3+1059230*n^2+1587211*n+944620))
D_fn = Rn((n+3)^4*(n+4)^6*(946*n^2+4515*n+5399))

# ============================================================
# CORRECT operator (Q5107 form)
# ============================================================
Qn = Rn.fraction_field()
A_rec.<Sn> = OreAlgebra(Qn, 'Sn')

L_monic = (Sn^3
           - B_fn(n=n+2)/A_fn(n=n+2) * Sn^2
           + C_fn(n=n+1)/A_fn(n=n+1) * Sn
           - D_fn(n=n)/A_fn(n=n))

print("Monic operator constructed (Q5107 form)")
print(f"Order: {L_monic.order()}")

# Clear denominators and take primitive part
coeffs_monic = [Qn(L_monic[i]) for i in range(4)]
den = Rn.one()
for c in coeffs_monic:
    den = lcm(den, Rn(c.denominator()))
nums = [Rn(den * c) for c in coeffs_monic]
g = nums[0]
for c in nums[1:]:
    g = gcd(g, c)
pnums = [c // g for c in nums]

print(f"\nPrimitive polynomial coefficients:")
for i in range(4):
    print(f"  p{i} degree: {pnums[i].degree()}")

L_prim = sum((Qn(pnums[i]) * Sn^i for i in range(4)), A_rec.zero())

# ============================================================
# Verify that this operator annihilates q_n
# ============================================================
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

q = [QQ(-215040420000), QQ(-167282265043404)/905, QQ(-964185327658080)/6071]
for nn in range(2, 55):
    q.append(Bf(nn)/Af(nn)*q[nn] - Cf(nn-1)/Af(nn-1)*q[nn-1] + Df(nn-2)/Af(nn-2)*q[nn-2])

print("\nVerifying operator annihilates q_n:")
ok = True
for nn in range(52):
    val = sum(QQ(pnums[i](n=nn)) * q[nn+i] for i in range(4))
    if val != 0:
        ok = False
        if nn < 5:
            print(f"  n={nn}: residual = {float(val):.3e}")
if ok:
    print("  *** ALL ZERO — operator is correct! ***")
else:
    print("  SOME NONZERO RESIDUALS")

# ============================================================
# Factor the primitive leading and trailing coefficients
# ============================================================
print(f"\nLeading coefficient (p3) factorization:")
print(f"  {pnums[3].factor()}")
print(f"\nTrailing coefficient (p0) factorization:")
print(f"  {pnums[0].factor()}")

# ============================================================
# Try to convert to ODE using ore_algebra
# ============================================================
print("\n" + "="*60)
print("Converting to differential equation...")
print("="*60)

try:
    # ore_algebra's to_D() converts a shift operator to differential
    R2.<z> = PolynomialRing(QQ)
    D_alg.<Dz> = OreAlgebra(R2)

    L_ode = L_prim.to_D(D_alg)
    print(f"\nDifferential operator obtained!")
    print(f"Order: {L_ode.order()}")
    deg = max(c.degree() for c in L_ode.list())
    print(f"Max coefficient degree: {deg}")

    # Singular points
    lc = L_ode.leading_coefficient()
    tc = L_ode.list()[0]
    print(f"\nLeading coefficient factorization:")
    print(f"  {lc.factor()}")
    print(f"\nTrailing coefficient factorization:")
    print(f"  {tc.factor()}")

    # Local exponents at singular points
    print("\n" + "="*60)
    print("Local exponents at singular points")
    print("="*60)

    sing_roots = lc.roots(QQbar)
    print(f"Number of singular points (finite): {len(sing_roots)}")
    for rt, mult in sing_roots:
        print(f"  z = {rt} (multiplicity {mult})")
        if rt.abs() < 10:
            try:
                indicial = L_ode.indicial_polynomial(z - QQ(rt) if rt in QQ else None)
                print(f"    Indicial polynomial: {indicial}")
            except:
                pass

except Exception as e:
    print(f"to_D conversion failed: {e}")
    import traceback
    traceback.print_exc()

    # Fallback: try guessing ODE from numerical data
    print("\nFallback: guessing ODE from q_n sequence...")
    R2.<z> = PolynomialRing(QQ)
    D_alg.<Dz> = OreAlgebra(R2)

    # Use ore_algebra's guess with exact rational q_n
    for ord in [3, 4, 6, 9, 12, 15, 18]:
        try:
            L_guess = guess(q[:50], D_alg, order=ord)
            deg = max(c.degree() for c in L_guess.list())
            print(f"  ODE found: order={ord}, degree={deg}")
            print(f"  Leading: {L_guess.leading_coefficient().factor()}")
            break
        except:
            print(f"  order={ord}: no ODE found")

# ============================================================
# Also compute AESZ #209 ODE for comparison
# ============================================================
print("\n" + "="*60)
print("AESZ #209 for comparison")
print("="*60)

inner = []
for nn in range(60):
    inner.append(sum(binomial(nn,k)^2 * binomial(nn+k,nn) * binomial(nn+2*k,nn)
                     for k in range(nn+1)))

R2.<z> = PolynomialRing(QQ)
D_alg.<Dz> = OreAlgebra(R2)

try:
    L_aesz = guess(inner[:60], D_alg, order=4)
    print(f"AESZ ODE: order={L_aesz.order()}, degree={max(c.degree() for c in L_aesz.list())}")
    print(f"Leading: {L_aesz.leading_coefficient().factor()}")

    if 'L_ode' in dir():
        print(f"\nGCRD(L_p27_ode, L_aesz):")
        try:
            G = L_ode.gcrd(L_aesz)
            print(f"  Order: {G.order()}")
            if G.order() > 0:
                print("  *** NONTRIVIAL GCRD! ***")
            else:
                print("  GCRD = 1 (no common factor)")
        except Exception as e:
            print(f"  Error: {e}")
except Exception as e:
    print(f"AESZ guess failed: {e}")

print("\nDone.")
