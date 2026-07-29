"""
Ore factorization of the P2.7 recurrence over K = Q(mu_0).

If a first-order right factor S - r(n) exists with r in K(n),
then the cubic-trace hypothesis is confirmed and the P2.7 module
decomposes into rank-one pieces over K.

Based on Q5102 §6, §10.
"""
from sage.all import *
from ore_algebra import OreAlgebra

# Number field K = Q(mu) where mu is a root of the Poincare cubic
Qx = PolynomialRing(QQ, 'x')
x = Qx.gen()
P = 4*x^3 - 220*x^2 + 8*x - 1
assert P.is_irreducible()
print("Poincare cubic:", P)
print("Discriminant:", P.discriminant())

K = NumberField(P, 'mu')
mu = K.gen()
print("Field K = Q(mu), disc =", K.discriminant())
print("Signature:", K.signature())

# Verify mu satisfies the polynomial
assert 4*mu^3 - 220*mu^2 + 8*mu - 1 == 0
print("mu satisfies P: OK")

# Set up Ore algebra K(n)<S>
Rn = PolynomialRing(K, 'n')
n = Rn.gen()
Fn = Rn.fraction_field()
Ore = OreAlgebra(Fn, 'S')
S = Ore.gen()

# P2.7 coefficient polynomials (in n)
def A(z):
    return Fn(1024*(2*z+5)^4*(2*z+7)^3*(2*z+9)^3*(946*z^2+6407*z+10860))

def B(z):
    return Fn(128*(2*z+7)^3*(2*z+9)^3*(104060*z^6+1745370*z^5+12145238*z^4+44886481*z^3+92943995*z^2+102256019*z+46709052))

def C(z):
    return Fn(16*(z+3)^4*(2*z+9)^3*(3784*z^5+57792*z^4+351019*z^3+1059230*z^2+1587211*z+944620))

def D(z):
    return Fn((z+3)^4*(z+4)^6*(946*z^2+4515*z+5399))

# Forward operator: A(n+2) u_{n+3} - B(n+2) u_{n+2} + C(n+1) u_{n+1} - D(n) u_n = 0
L = A(n+2)*S^3 - B(n+2)*S^2 + C(n+1)*S - D(n)

print("\nOperator L constructed.")
print("Order:", L.order())

# Check leading coefficient
lc = L.leading_coefficient()
print("Leading coefficient degree:", Rn(lc.numerator()).degree())

# Attempt factorization
print("\n" + "="*60)
print("Attempting factorization of L over K(n)<S>...")
print("="*60)
import sys
sys.stdout.flush()

try:
    factors = L.factor()
    print("\nFactorization result:")
    print(factors)

    # Check if there's a first-order factor
    for f, mult in factors:
        print(f"\n  Factor of order {f.order()}, multiplicity {mult}")
        if f.order() == 1:
            # Extract r(n) from S - r(n)
            print("  FIRST-ORDER FACTOR FOUND!")
            print("  Factor:", f)
except Exception as e:
    print(f"\nFactorization raised an exception: {type(e).__name__}: {e}")
    print("\nTrying right_factor search instead...")
    try:
        rf = L.right_factor(1)
        print("Right factor of order 1:", rf)
    except Exception as e2:
        print(f"Right factor search also failed: {type(e2).__name__}: {e2}")

print("\n" + "="*60)
print("Also checking: does L factor over Q(n)?")
print("="*60)

# Same operator over Q(n) for comparison
Rn_Q = PolynomialRing(QQ, 'n')
n_Q = Rn_Q.gen()
Fn_Q = Rn_Q.fraction_field()
Ore_Q = OreAlgebra(Fn_Q, 'S')
S_Q = Ore_Q.gen()

def AQ(z): return Fn_Q(1024*(2*z+5)^4*(2*z+7)^3*(2*z+9)^3*(946*z^2+6407*z+10860))
def BQ(z): return Fn_Q(128*(2*z+7)^3*(2*z+9)^3*(104060*z^6+1745370*z^5+12145238*z^4+44886481*z^3+92943995*z^2+102256019*z+46709052))
def CQ(z): return Fn_Q(16*(z+3)^4*(2*z+9)^3*(3784*z^5+57792*z^4+351019*z^3+1059230*z^2+1587211*z+944620))
def DQ(z): return Fn_Q((z+3)^4*(z+4)^6*(946*z^2+4515*z+5399))

L_Q = AQ(n_Q+2)*S_Q^3 - BQ(n_Q+2)*S_Q^2 + CQ(n_Q+1)*S_Q - DQ(n_Q)

try:
    factors_Q = L_Q.factor()
    print("Factorization over Q(n):", factors_Q)
except Exception as e:
    print(f"Factorization over Q(n) failed: {e}")

print("\nDone.")
