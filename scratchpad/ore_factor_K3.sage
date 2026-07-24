"""
Ore factorization of P2.7 using right_factors().
"""
from sage.all import *
from ore_algebra import OreAlgebra
import sys

# Number field K = Q(mu)
Qx = PolynomialRing(QQ, 'x')
x = Qx.gen()
P = 4*x^3 - 220*x^2 + 8*x - 1
K = NumberField(P, 'mu')
mu = K.gen()
print("K = Q(mu), disc =", K.discriminant())

# First test: factor over Q(n)
print("\n=== Testing over Q(n) ===")
R_Q = PolynomialRing(QQ, 'n')
n_Q = R_Q.gen()
A_Q = OreAlgebra(R_Q, 'Sn')
Sn_Q = A_Q.gen()

def AQ(z): return R_Q(1024*(2*z+5)^4*(2*z+7)^3*(2*z+9)^3*(946*z^2+6407*z+10860))
def BQ(z): return R_Q(128*(2*z+7)^3*(2*z+9)^3*(104060*z^6+1745370*z^5+12145238*z^4+44886481*z^3+92943995*z^2+102256019*z+46709052))
def CQ(z): return R_Q(16*(z+3)^4*(2*z+9)^3*(3784*z^5+57792*z^4+351019*z^3+1059230*z^2+1587211*z+944620))
def DQ(z): return R_Q((z+3)^4*(z+4)^6*(946*z^2+4515*z+5399))

L_Q = AQ(n_Q+2)*Sn_Q^3 - BQ(n_Q+2)*Sn_Q^2 + CQ(n_Q+1)*Sn_Q - DQ(n_Q)
print("L over Q: order", L_Q.order())

print("Searching for right factors of order 1 over Q(n)...")
sys.stdout.flush()
try:
    rf1_Q = L_Q.right_factors(1)
    print("  Right factors of order 1:", rf1_Q)
except Exception as e:
    print(f"  Exception: {e}")

print("Searching for right factors of order 2 over Q(n)...")
sys.stdout.flush()
try:
    rf2_Q = L_Q.right_factors(2)
    print("  Right factors of order 2:", rf2_Q)
except Exception as e:
    print(f"  Exception: {e}")

# Now over K(n)
print("\n=== Testing over K(n) ===")
R_K = PolynomialRing(K, 'n')
n_K = R_K.gen()
A_K = OreAlgebra(R_K, 'Sn')
Sn_K = A_K.gen()

def AK(z): return R_K(1024*(2*z+5)^4*(2*z+7)^3*(2*z+9)^3*(946*z^2+6407*z+10860))
def BK(z): return R_K(128*(2*z+7)^3*(2*z+9)^3*(104060*z^6+1745370*z^5+12145238*z^4+44886481*z^3+92943995*z^2+102256019*z+46709052))
def CK(z): return R_K(16*(z+3)^4*(2*z+9)^3*(3784*z^5+57792*z^4+351019*z^3+1059230*z^2+1587211*z+944620))
def DK(z): return R_K((z+3)^4*(z+4)^6*(946*z^2+4515*z+5399))

L_K = AK(n_K+2)*Sn_K^3 - BK(n_K+2)*Sn_K^2 + CK(n_K+1)*Sn_K - DK(n_K)
print("L over K: order", L_K.order())

print("Searching for right factors of order 1 over K(n)...")
sys.stdout.flush()
try:
    rf1_K = L_K.right_factors(1)
    print("  Right factors of order 1:", rf1_K)
    if rf1_K:
        print("  *** FIRST-ORDER RIGHT FACTOR FOUND! ***")
        for rf in rf1_K:
            print("    Factor:", rf)
            # Extract r(n) = -coefficient of Sn^0 / coefficient of Sn
            coeffs = rf.to_list()
            if len(coeffs) == 2:
                r_n = -coeffs[0] / coeffs[1]
                print("    r(n) = -p0/p1 =", r_n)
    else:
        print("  No order-1 right factors exist over K(n).")
except Exception as e:
    print(f"  Exception: {type(e).__name__}: {e}")

print("\nSearching for right factors of order 2 over K(n)...")
sys.stdout.flush()
try:
    rf2_K = L_K.right_factors(2)
    print("  Right factors of order 2:", rf2_K)
    if rf2_K:
        print("  *** ORDER-2 RIGHT FACTOR FOUND! ***")
except Exception as e:
    print(f"  Exception: {type(e).__name__}: {e}")

print("\nDone.")
