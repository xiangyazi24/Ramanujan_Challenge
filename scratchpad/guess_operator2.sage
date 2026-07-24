"""
Guess the minimal polynomial-coefficient recurrence operator
for the ACTUAL P2.7 sequences (proof.tex convention).

Uses 200 terms to handle degree-36 coefficients.
Need at least order*(degree+1)+1 = 3*37+1 = 112 terms for order 3, degree 36.
"""
from sage.all import *
from ore_algebra import OreAlgebra, guess
import sys

# P2.7 coefficient polynomials
def AA(n):
    return (1024*(2*n+5)**4*(2*n+7)**3*(2*n+9)**3
            *(946*n**2 + 6407*n + 10860))

def BB(n):
    return (128*(2*n+7)**3*(2*n+9)**3
            *(104060*n**6 + 1745370*n**5 + 12145238*n**4
              + 44886481*n**3 + 92943995*n**2
              + 102256019*n + 46709052))

def CC(n):
    return (16*(n+3)**4*(2*n+9)**3
            *(3784*n**5 + 57792*n**4 + 351019*n**3
              + 1059230*n**2 + 1587211*n + 944620))

def DD(n):
    return ((n+3)**4*(n+4)**6
            *(946*n**2 + 4515*n + 5399))

# Correct monic recurrence (proof.tex)
def alpha(n): return QQ(BB(n)) / QQ(AA(n))
def beta_func(n): return -QQ(CC(n-1)) / QQ(AA(n-1))
def gamma_func(n): return QQ(DD(n-2)) / QQ(AA(n-2))

# Compute exact rational q_n
N = 200
print(f"Computing {N} exact rational q_n...")
sys.stdout.flush()
q = [QQ(0)] * (N + 3)
q[0] = QQ(-215040420000)
q[1] = QQ(-167282265043404) / QQ(905)
q[2] = QQ(-964185327658080) / QQ(6071)

for n in range(2, N + 2):
    q[n+1] = alpha(n)*q[n] + beta_func(n)*q[n-1] + gamma_func(n)*q[n-2]

print(f"q[0] = {q[0]}")
print(f"q[10] = {float(q[10]):.6e}")
print(f"q[100] = {float(q[100]):.6e}")
print(f"q[199] = {float(q[199]):.6e}")

# Set up Ore algebra
R = PolynomialRing(QQ, 'n')
n_var = R.gen()
A_ore = OreAlgebra(R, 'Sn')
Sn = A_ore.gen()

seq_data = q[:N]

# Try with increasing degree bounds
for deg in [12, 16, 20, 24, 28, 32, 36, 40]:
    needed = 3*(deg+1)+1
    if needed > N:
        print(f"\nSkipping degree {deg}: need {needed} terms, have {N}")
        continue
    print(f"\nGuessing order=3, degree={deg} (need {needed} terms)...")
    sys.stdout.flush()
    try:
        ops = guess(seq_data, A_ore, order=3, degree=deg)
        if ops:
            print(f"  FOUND!")
            for L in ops:
                print(f"  Order: {L.order()}")
                coeffs = L.to_list()
                for i, c in enumerate(coeffs):
                    print(f"  coeff[Sn^{i}]: degree {c.degree()}, lc = {c.leading_coefficient()}")

                # Check GCD of coefficients
                from sage.arith.misc import gcd as sage_gcd
                g = coeffs[0]
                for c in coeffs[1:]:
                    g = g.gcd(c)
                print(f"  GCD of coefficients: degree {g.degree()}")
                if g.degree() > 0:
                    print(f"  GCD = {g.factor()}")
                    # Divide out GCD
                    reduced = [c // g for c in coeffs]
                    for i, c in enumerate(reduced):
                        print(f"  reduced coeff[Sn^{i}]: degree {c.degree()}")

                # Verify
                print("  Verifying...")
                for k in range(min(N-4, 100)):
                    val = sum(coeffs[j](n_var=k) * q[k+j] for j in range(len(coeffs)))
                    if val != 0:
                        print(f"    FAILS at k={k}: residual = {float(val):.6e}")
                        break
                else:
                    print(f"    Verified for k=0..{min(N-4, 100)-1}")
            break
        else:
            print(f"  No relations found")
    except Exception as e:
        print(f"  Exception: {type(e).__name__}: {e}")

# Also try higher order
for ord_val in [4, 5, 6]:
    for deg in [12, 16]:
        needed = ord_val*(deg+1)+1
        if needed > N:
            continue
        print(f"\nGuessing order={ord_val}, degree={deg} (need {needed} terms)...")
        sys.stdout.flush()
        try:
            ops = guess(seq_data, A_ore, order=ord_val, degree=deg)
            if ops:
                print(f"  FOUND!")
                for L in ops:
                    print(f"  Order: {L.order()}")
                    coeffs = L.to_list()
                    for i, c in enumerate(coeffs):
                        print(f"  coeff[Sn^{i}]: degree {c.degree()}")
                break
            else:
                print(f"  No relations found")
        except Exception as e:
            print(f"  Exception: {type(e).__name__}: {e}")

print("\nDone.")
