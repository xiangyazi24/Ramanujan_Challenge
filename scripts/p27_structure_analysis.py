#!/usr/bin/env python3
"""P2.7: Analyze the structure of q_n and p_n sequences.
Factor the initial values, look at denominator patterns, find normalization."""
from fractions import Fraction
from math import gcd
import sympy

def A(n): return 1024*(2*n+5)**4*(2*n+7)**3*(2*n+9)**3*(946*n**2+6407*n+10860)
def B(n): return 128*(2*n+7)**3*(2*n+9)**3*(104060*n**6+1745370*n**5+12145238*n**4+44886481*n**3+92943995*n**2+102256019*n+46709052)
def C(n): return 16*(n+3)**4*(2*n+9)**3*(3784*n**5+57792*n**4+351019*n**3+1059230*n**2+1587211*n+944620)
def D(n): return (n+3)**4*(n+4)**6*(946*n**2+4515*n+5399)

q = [None]*250
p = [None]*250
q[0] = Fraction(-215040420000)
q[1] = Fraction(-167282265043404, 905)
q[2] = Fraction(-964185327658080, 6071)
p[0] = Fraction(-612218384750)
p[1] = Fraction(-9525021973931919, 18100)
p[2] = Fraction(-29561828382772029, 65380)

print("=== Initial value factorizations ===")
for name, val in [("q0", -215040420000), ("p0", -612218384750)]:
    fac = sympy.factorint(abs(val))
    print(f"  |{name}| = {abs(val)} = {' * '.join(f'{p}^{e}' if e > 1 else str(p) for p, e in sorted(fac.items()))}")

print(f"\n  q1 denom = 905 = {sympy.factorint(905)}")
print(f"  q2 denom = 6071 = {sympy.factorint(6071)}")
print(f"  p1 denom = 18100 = {sympy.factorint(18100)}")
print(f"  p2 denom = 65380 = {sympy.factorint(65380)}")

# Compute more terms
print("\n=== Computing q_n, p_n for n=0..100 ===", flush=True)
for n in range(2, 100):
    An = A(n)
    Bn = B(n)
    Cn_1 = C(n-1)
    An_1 = A(n-1)
    Dn_2 = D(n-2)
    An_2 = A(n-2)
    q[n+1] = Fraction(Bn, An) * q[n] - Fraction(Cn_1, An_1) * q[n-1] + Fraction(Dn_2, An_2) * q[n-2]
    p[n+1] = Fraction(Bn, An) * p[n] - Fraction(Cn_1, An_1) * p[n-1] + Fraction(Dn_2, An_2) * p[n-2]
    if n % 20 == 0:
        print(f"  n={n+1} done, q denom size = {len(str(q[n+1].denominator))} digits", flush=True)

# Analyze denominators
print("\n=== Denominator analysis ===")
for n in range(11):
    d = q[n].denominator
    if d > 1:
        fac = sympy.factorint(d) if d < 10**15 else f"({len(str(d))} digits)"
        print(f"  q_{n} denom = {d} = {fac}")
    else:
        print(f"  q_{n} is integer: {q[n]}")

# Check: is there a common denominator pattern?
# Look at lcm of denominators
print("\n=== LCM pattern ===")
from functools import reduce
for N in [5, 10, 15, 20]:
    denoms = [q[n].denominator for n in range(N+1)]
    lcm_val = reduce(lambda a, b: a * b // gcd(a, b), denoms)
    print(f"  lcm(denom(q_0),...,denom(q_{N})) has {len(str(lcm_val))} digits")

# Check normalization: multiply all q_n by some common factor
# Factor A(0), A(1), A(2) to understand the denominator accumulation
print("\n=== Coefficient factorization ===")
for n in range(5):
    a = A(n)
    print(f"  A({n}) = {a}")
    print(f"       = {sympy.factorint(a)}")

print()
for n in range(5):
    d = D(n)
    print(f"  D({n}) = {d}")
    print(f"       = {sympy.factorint(d)}")

# The determinant ratio D(n)/A(n) governs the CMF
print("\n=== D(n)/A(n) ratios ===")
for n in range(5):
    r = Fraction(D(n), A(n))
    print(f"  D({n})/A({n}) = {r}")

# Compute q_n * lcm_denom and check if we get integers
print("\n=== Trying to find integer normalization ===")
# Strategy: compute product of A(j) / gcd(A(j), B(j)*C(j)*D(j)) for j=0..n
# to understand what denominators accumulate

# Let's look at the "normalized" recurrence: define Q_n = q_n * prod(A(j), j=0..n-1)
# or similar

# Actually, let's look at: q_n * denom(q_n) = numerator(q_n)
print("  Numerators of q_n (first terms):")
for n in range(8):
    print(f"    q_{n} num = {q[n].numerator}")

# Check: the Apery numbers for ζ(3) are a_n = Σ C(n,k)²C(n+k,k)², all integers.
# For P2.7, q_n are NOT integers in general. What if we normalize?

# Try: define hat_q_n = q_n * Π_{j=0}^{n-1} A(j) / (something)
# The recurrence is u_{n+1} = (B/A) u - (C'/A') u' + (D''/A'') u''
# Writing v_n = u_n * Π A(j): doesn't simplify nicely.

# Alternative: check if there's a "gauge transform" s_n such that s_n * q_n is integer
# s_n = lcm of all denominators appearing in q_0,...,q_n

# Check the denominator growth
print("\n=== Denominator digit growth ===")
for n in range(0, 101, 5):
    if q[n] is not None:
        d = q[n].denominator
        nd = q[n].numerator
        print(f"  n={n:3d}: num {len(str(abs(nd))):4d} digits, den {len(str(d)):4d} digits")

# Let's also look at the "error" e_n = p_n - L*q_n to understand the structure
from mpmath import mp, mpf, zeta
mp.dps = 200
L = mpf(zeta(2)) + mpf(zeta(3))
print(f"\n=== Error sequence analysis ===")
print(f"  L = ζ(2)+ζ(3) = {mp.nstr(L, 50)}")

for n in [0, 1, 2, 5, 10, 20, 30]:
    if q[n] is not None:
        pn = mpf(p[n].numerator) / mpf(p[n].denominator)
        qn = mpf(q[n].numerator) / mpf(q[n].denominator)
        en = pn - L * qn
        print(f"  n={n:3d}: |e_n| ≈ {mp.nstr(abs(en), 10)}")

print("\nDone.", flush=True)
