#!/usr/bin/env python3
"""Wide sweep: prime factors of G(d1,d2)-type gcds for d1 < 20, d2 <= 50.
Theorem 4's exceptional set is fed only by primes dividing these; verify all prime factors are O(d2)."""
from fractions import Fraction as F
from math import gcd
import sys

def coefs(x):
    den = F((x+1)**3)
    return F((2*x+1)*(17*x*x+17*x+5),1)/den, F(-(x**3),1)/den

def ABat(k, base):
    if k == 0: return F(1), F(0)
    Ap, Bp = F(1), F(0)
    a1, b1 = coefs(base); Ac, Bc = a1, b1
    for j in range(1, k):
        aj, bj = coefs(base + j)
        Ac, Ap = aj*Ac + bj*Ap, Ac
        Bc, Bp = aj*Bc + bj*Bp, Bc
    return Ac, Bc

b = [1,5]
for n in range(1, 60):
    num = (2*n+1)*(17*n*n+17*n+5)*b[n] - n**3*b[n-1]
    q, rm = divmod(num, (n+1)**3); assert rm == 0
    b.append(q)

def largest_pf(n, cap=10**6):
    if n <= 1: return 1
    lp = 1; d = 2
    while d*d <= n and d < cap:
        while n % d == 0: lp = d; n //= d
        d += 1
    if n > 1: lp = max(lp, n if n < cap*cap else -n)  # negative flags huge residual factor
    return lp

worst = 0; worst_pair = None; flagged = []
for d1 in range(1, 20):
    for d2 in range(d1+1, 51):
        A1v, B1v = ABat(d1, -d2)
        K = F(d2**3)*((1 - A1v)*F(b[d2-1]) - F(b[d2])*B1v)
        Kn = abs(K.numerator)
        delta = d2 - d1
        for other in (abs(b[delta]-1), abs(b[d1-1]) if d1 >= 1 else 1):
            g = gcd(Kn, other)
            if g > 1:
                lp = largest_pf(g)
                if lp < 0: flagged.append((d1,d2,g)); continue
                if lp > worst: worst, worst_pair = lp, (d1, d2, g)
print("largest prime factor over sweep:", worst, "at (d1,d2,gcd) =", worst_pair)
print("huge-residual flags:", flagged[:10], "count:", len(flagged))
