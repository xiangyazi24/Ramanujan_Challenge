"""
Guess the minimal polynomial-coefficient recurrence operator
for the ACTUAL P2.7 sequences (proof.tex convention).

Strategy: compute q_n as exact rationals, then use ore_algebra
to guess the annihilating operator.
"""
from sage.all import *
from ore_algebra import OreAlgebra, guess

# P2.7 coefficient polynomials
Rx = PolynomialRing(QQ, 'x')
x = Rx.gen()

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
def gamma(n): return QQ(DD(n-2)) / QQ(AA(n-2))

# Compute exact rational q_n
print("Computing exact rational q_n...")
N = 60
q = [QQ(0)] * (N + 3)
q[0] = QQ(-215040420000)
q[1] = QQ(-167282265043404) / QQ(905)
q[2] = QQ(-964185327658080) / QQ(6071)

for n in range(2, N + 2):
    q[n+1] = alpha(n)*q[n] + beta_func(n)*q[n-1] + gamma(n)*q[n-2]

print(f"q[0] = {q[0]}")
print(f"q[1] = {q[1]}")
print(f"q[3] = {q[3]}")
print(f"q[5] = {q[5]}")

# Set up Ore algebra for guessing
R = PolynomialRing(QQ, 'n')
n_var = R.gen()
A_ore = OreAlgebra(R, 'Sn')
Sn = A_ore.gen()

# Guess the annihilating operator
print("\nGuessing operator from q_n data...")
import sys
sys.stdout.flush()

# Try with enough terms for a degree-18 order-3 operator
# Need at least order*(degree+1) + order + 1 terms
# For order 3, degree 18: 3*19 + 3 + 1 = 61 terms
seq_data = q[:N]

try:
    ops = guess(seq_data, A_ore, order=3, degree=20)
    if ops:
        print(f"Found operator(s)!")
        for L in ops:
            print(f"  Order: {L.order()}")
            coeffs = L.to_list()
            for i, c in enumerate(coeffs):
                print(f"  coeff[Sn^{i}]: degree {c.degree()}")

            # Verify it annihilates the sequence
            print("  Verifying...")
            for k in range(min(N-4, 50)):
                val = sum(coeffs[j] * q[k+j] for j in range(len(coeffs)))
                if val(n_var=k) != 0:
                    print(f"    FAILS at k={k}!")
                    break
            else:
                print("    Verified for k=0..{0}".format(min(N-4, 50)-1))
    else:
        print("No operator found at order=3, degree=20")
except Exception as e:
    print(f"Exception: {type(e).__name__}: {e}")

# Also try order 3 with higher degree
print("\nTrying order=3, degree=36...")
sys.stdout.flush()
try:
    ops36 = guess(seq_data, A_ore, order=3, degree=40)
    if ops36:
        for L in ops36:
            print(f"  Order: {L.order()}")
            coeffs = L.to_list()
            for i, c in enumerate(coeffs):
                print(f"  coeff[Sn^{i}]: degree {c.degree()}")
    else:
        print("  No operator found")
except Exception as e:
    print(f"  Exception: {type(e).__name__}: {e}")

# Compare with eq 2.1 operator
print("\n--- Checking if eq 2.1 annihilates the correct q_n ---")
def p3(k): return AA(k+2)
def p2(k): return -BB(k+2)
def p1(k): return CC(k+1)
def p0(k): return -DD(k)

# Check if the eq 2.1 operator kills our sequence
for k in range(5):
    val = p3(k)*q[k+3] + p2(k)*q[k+2] + p1(k)*q[k+1] + p0(k)*q[k]
    if val == 0:
        print(f"  k={k}: eq 2.1 residual = 0 (EXACT)")
    else:
        print(f"  k={k}: eq 2.1 residual = {float(val):.6e} (NONZERO)")

print("\nDone.")
