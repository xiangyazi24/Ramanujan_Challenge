"""
Try to match P2.7 q_n to a hypergeometric sum.

Key structural clue from proof.tex: the coboundary h_n has Pochhammer
structure (3)_n^4 (4)_n^6 / [(5/2)_n^4 (7/2)_n^3 (9/2)_n^3],
which belongs to the Whipple 4F3 <-> 5F4 family (Dauguet-Zudilin).

Also, AESZ #209 period is:
  A_n = C(2n,n) * sum_{k=0}^{n} C(n,k)^2 * C(n+k,n) * C(n+2k,n)

Strategy: compute q_n for several values, try to express as a sum
involving factorials/Pochhammer symbols.
"""
from mpmath import mp, mpf, fac, binomial, zeta, pi as mppi
from fractions import Fraction

mp.dps = 50

# P2.7 recurrence coefficients (exact rational)
def A(n):
    return 1024*(2*n+5)**4*(2*n+7)**3*(2*n+9)**3*(946*n**2 + 6407*n + 10860)

def B(n):
    return 128*(2*n+7)**3*(2*n+9)**3*(104060*n**6 + 1745370*n**5 + 12145238*n**4
                + 44886481*n**3 + 92943995*n**2 + 102256019*n + 46709052)

def C(n):
    return 16*(n+3)**4*(2*n+9)**3*(3784*n**5 + 57792*n**4 + 351019*n**3
                + 1059230*n**2 + 1587211*n + 944620)

def D(n):
    return (n+3)**4*(n+4)**6*(946*n**2 + 4515*n + 5399)

def alpha(n): return Fraction(B(n), A(n))
def beta(n): return Fraction(-C(n-1), A(n-1))
def gamma(n): return Fraction(D(n-2), A(n-2))

# Compute q_n as exact rationals
from fractions import Fraction
q = [Fraction(0)] * 15
p = [Fraction(0)] * 15

q[0] = Fraction(-215040420000)
q[1] = Fraction(-167282265043404, 905)
q[2] = Fraction(-964185327658080, 6071)

p[0] = Fraction(-612218384750)
p[1] = Fraction(-9525021973931919, 18100)
p[2] = Fraction(-29561828382772029, 65380)

for n in range(2, 12):
    a = alpha(n)
    b = beta(n)
    g = gamma(n)
    q[n+1] = a * q[n] + b * q[n-1] + g * q[n-2]
    p[n+1] = a * p[n] + b * p[n-1] + g * p[n-2]

print("=== P2.7 q_n values (exact rational) ===")
for n in range(8):
    print(f"q[{n}] = {q[n]}")
    print(f"  numerator = {q[n].numerator}")
    print(f"  denominator = {q[n].denominator}")
    print()

# Factor q[0]
import sympy
q0_int = -215040420000
print(f"q[0] = {q0_int}")
print(f"factorization: {sympy.factorint(abs(q0_int))}")

# Factor p[0]
p0_int = -612218384750
print(f"\np[0] = {p0_int}")
print(f"factorization: {sympy.factorint(abs(p0_int))}")

# Check AESZ #209 period for comparison
print("\n=== AESZ #209 period A_n ===")
def aesz209(n):
    """A_n = C(2n,n) * sum_{k=0}^{n} C(n,k)^2 * C(n+k,n) * C(n+2k,n)"""
    n = int(n)
    result = 0
    c2n_n = sympy.binomial(2*n, n)
    for k in range(n+1):
        term = sympy.binomial(n,k)**2 * sympy.binomial(n+k,n) * sympy.binomial(n+2*k,n)
        result += term
    return c2n_n * result

for n in range(8):
    print(f"A[{n}] = {aesz209(n)}")

# Check ratios between q_n and AESZ period
print("\n=== Ratio q_n / A_n (check for proportionality) ===")
for n in range(8):
    a_n = aesz209(n)
    if a_n != 0:
        ratio = Fraction(q[n].numerator, q[n].denominator * a_n)
        print(f"q[{n}] / A[{n}] = {float(ratio):.6e}")

# The coboundary hypergeometric term
print("\n=== Coboundary h_n / h_0 ===")
# h_n = 2^{-20n} * (3)_n^4 * (4)_n^6 / [(5/2)_n^4 * (7/2)_n^3 * (9/2)_n^3]
def pochhammer(a, n):
    """Rising factorial (a)_n = a(a+1)...(a+n-1)"""
    result = Fraction(1)
    for j in range(n):
        result *= Fraction(a).limit_denominator(1000) + j
    return result

def pochhammer_exact(a_num, a_den, n):
    """Exact Pochhammer for rational a = a_num/a_den."""
    result = Fraction(1)
    for j in range(n):
        result *= Fraction(a_num + j * a_den, a_den)
    return result

def h_ratio(n):
    """h_n / h_0"""
    top = pochhammer_exact(3, 1, n)**4 * pochhammer_exact(4, 1, n)**6
    bot = pochhammer_exact(5, 2, n)**4 * pochhammer_exact(7, 2, n)**3 * pochhammer_exact(9, 2, n)**3
    return Fraction(1, 2**(20*n)) * top / bot

for n in range(6):
    h = h_ratio(n)
    print(f"h[{n}]/h[0] = {float(h):.10e}")

# Product of h_0..h_{n-1}
print("\n=== Product of h ratios (det C(0)*...*C(n-1)) ===")
prod = Fraction(1)
for n in range(1, 6):
    prod *= h_ratio(n)
    print(f"prod h[0..{n}] = {float(prod):.10e}")

# Check: is det C(n) = D(n)/A(n) equal to [R(n)/R(n+1)] * delta(n)?
print("\n=== Verify det C(n) = D(n)/A(n) ===")
for n in range(5):
    det_c = Fraction(D(n), A(n))
    R_n = 946*n**2 + 4515*n + 5399
    R_n1 = 946*(n+1)**2 + 4515*(n+1) + 5399
    delta_n = h_ratio(n+1) / h_ratio(n) if n > 0 else h_ratio(1)
    Rr = Fraction(R_n, R_n1)
    # should have det_c = R(n)/R(n+1) * delta(n)
    # delta(n) = h(n+1)/h(n) = h_ratio(n+1)/h_ratio(n)
    if n == 0:
        delta_n = h_ratio(1) / Fraction(1)
    else:
        delta_n = h_ratio(n+1) / h_ratio(n)
    expected = Rr * delta_n
    print(f"n={n}: det C = {float(det_c):.10e}, R/R' * delta = {float(expected):.10e}, match = {abs(float(det_c) - float(expected)) < 1e-20}")

print("\nDone.")
