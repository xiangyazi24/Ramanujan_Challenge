"""
Ore factorization of P2.7 over K = Q(mu_0).
Fixed OreAlgebra construction.
"""
from sage.all import *
from ore_algebra import OreAlgebra

# Number field
Qx = PolynomialRing(QQ, 'x')
x = Qx.gen()
P = 4*x^3 - 220*x^2 + 8*x - 1
K = NumberField(P, 'mu')
mu = K.gen()
print("K = Q(mu), disc =", K.discriminant(), ", sig =", K.signature())

# Polynomial ring and Ore algebra over K
R = PolynomialRing(K, 'n')
n = R.gen()
A_ore = OreAlgebra(R, 'Sn')
Sn = A_ore.gen()
print("Ore algebra constructed:", A_ore)

# P2.7 coefficients as elements of K[n]
def Ac(z):
    return R(1024*(2*z+5)^4*(2*z+7)^3*(2*z+9)^3*(946*z^2+6407*z+10860))

def Bc(z):
    return R(128*(2*z+7)^3*(2*z+9)^3*(104060*z^6+1745370*z^5+12145238*z^4+44886481*z^3+92943995*z^2+102256019*z+46709052))

def Cc(z):
    return R(16*(z+3)^4*(2*z+9)^3*(3784*z^5+57792*z^4+351019*z^3+1059230*z^2+1587211*z+944620))

def Dc(z):
    return R((z+3)^4*(z+4)^6*(946*z^2+4515*z+5399))

# Forward operator: A(n+2) u_{n+3} - B(n+2) u_{n+2} + C(n+1) u_{n+1} - D(n) u_n = 0
L = Ac(n+2)*Sn^3 - Bc(n+2)*Sn^2 + Cc(n+1)*Sn - Dc(n)
print("Operator L: order", L.order())

# First try factoring over Q(n)
print("\n--- Factoring over Q(n) first ---")
R_Q = PolynomialRing(QQ, 'n')
n_Q = R_Q.gen()
A_Q = OreAlgebra(R_Q, 'Sn')
Sn_Q = A_Q.gen()

def AQ(z): return R_Q(1024*(2*z+5)^4*(2*z+7)^3*(2*z+9)^3*(946*z^2+6407*z+10860))
def BQ(z): return R_Q(128*(2*z+7)^3*(2*z+9)^3*(104060*z^6+1745370*z^5+12145238*z^4+44886481*z^3+92943995*z^2+102256019*z+46709052))
def CQ(z): return R_Q(16*(z+3)^4*(2*z+9)^3*(3784*z^5+57792*z^4+351019*z^3+1059230*z^2+1587211*z+944620))
def DQ(z): return R_Q((z+3)^4*(z+4)^6*(946*z^2+4515*z+5399))

L_Q = AQ(n_Q+2)*Sn_Q^3 - BQ(n_Q+2)*Sn_Q^2 + CQ(n_Q+1)*Sn_Q - DQ(n_Q)

import sys
sys.stdout.flush()

try:
    print("Factoring L over Q(n)...")
    sys.stdout.flush()
    fQ = L_Q.factor()
    print("Result:", fQ)
except Exception as e:
    print(f"  Exception: {type(e).__name__}: {e}")

sys.stdout.flush()

# Now factor over K(n)
print("\n--- Factoring over K(n) ---")
sys.stdout.flush()
try:
    print("Factoring L over K(n)...")
    sys.stdout.flush()
    fK = L.factor()
    print("Result:", fK)
    for f, m in fK:
        print(f"  order {f.order()}, mult {m}: {f}")
except Exception as e:
    print(f"  Exception: {type(e).__name__}: {e}")

    # Try searching for right factors directly
    print("\nSearching for right factor of order 1...")
    sys.stdout.flush()
    try:
        rf = L.right_factor(1)
        print("Right factor found:", rf)
    except Exception as e2:
        print(f"  No order-1 right factor: {type(e2).__name__}: {e2}")

    print("\nSearching for right factor of order 2...")
    sys.stdout.flush()
    try:
        rf2 = L.right_factor(2)
        print("Right factor found:", rf2)
    except Exception as e3:
        print(f"  No order-2 right factor: {type(e3).__name__}: {e3}")

print("\nDone.")
