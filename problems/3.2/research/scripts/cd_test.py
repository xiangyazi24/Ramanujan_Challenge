#!/usr/bin/env python3
"""c_d = -5*B_{d-1}(-d) + B_{d-2}(-d) over Q, d=2..30. Nonzero? closed form?"""
from fractions import Fraction as F
import math

def coefs(x):
    # a(x) = P(x)/(x+1)^3, beta(x) = -x^3/(x+1)^3, exact rational; x integer (negative ok, x != -1)
    den = F((x+1)**3)
    return F((2*x+1)*(17*x*x+17*x+5), 1)/den, F(-(x**3),1)/den

def B_at(k, base):
    """B_k evaluated at r=base (exact rational), via recursion in k."""
    Aprev, Bprev = F(1), F(0)          # k=0
    if k == 0: return F(0)
    a1, b1 = coefs(base)
    Acur, Bcur = a1, b1                 # k=1
    for j in range(1, k):
        aj, bj = coefs(base + j)
        Acur, Aprev = aj*Acur + bj*Aprev, Acur
        Bcur, Bprev = aj*Bcur + bj*Bprev, Bcur
    return Bcur

print(f"{'d':>3} {'c_d':>40} {'num factored?':>20}")
for d in range(2, 31):
    # poles of B_{d-1}, B_{d-2} at r in {-1..-(d-1)},{-1..-(d-2)}: r=-d is regular, safe
    cd = F(-5)*B_at(d-1, -d) + B_at(d-2, -d)
    num, den = cd.numerator, cd.denominator
    print(f"{d:>3} {str(cd)[:40]:>40} num={num}")

# verify closed form c_d = d^3 * b_{d-1}
b = [1,5]
for n in range(1, 32):
    num = (2*n+1)*(17*n*n+17*n+5)*b[n] - n**3*b[n-1]
    q, rem = divmod(num, (n+1)**3); assert rem == 0
    b.append(q)
allok = True
for d in range(2, 31):
    cd = F(-5)*B_at(d-1, -d) + B_at(d-2, -d)
    if cd != F(d**3 * b[d-1]):
        allok = False; print("MISMATCH at d =", d)
print("c_d = d^3 * b_{d-1} for d=2..30:", "VERIFIED" if allok else "FAILED")
