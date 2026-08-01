#!/usr/bin/env python3
"""Refined conjecture: gcd(b_delta-1, numerator(K/d2^3)) bounded by absolute constant.
Sweep d1 < 16, d2 <= 22."""
from fractions import Fraction as F
from math import gcd

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
for n in range(1, 50):
    num = (2*n+1)*(17*n*n+17*n+5)*b[n] - n**3*b[n-1]
    q, rm = divmod(num, (n+1)**3); assert rm == 0
    b.append(q)

mx = 1; arg = None; vals = []
for d1 in range(1, 16):
    for d2 in range(d1+1, 23):
        A1v, B1v = ABat(d1, -d2)
        bracket = (1 - A1v)*F(b[d2-1]) - F(b[d2])*B1v   # K = d2^3 * bracket
        Bn = abs(bracket.numerator)
        if Bn == 0: print("bracket=0 at", d1, d2); continue
        g = gcd(Bn, abs(b[d2-d1]-1))
        vals.append(g)
        if g > mx: mx, arg = g, (d1, d2)
print("max gcd(b_delta-1, num(K/d2^3)) =", mx, "at", arg)
import collections
print("gcd distribution:", dict(collections.Counter(vals)))
def fac(n):
    f={}; d=2
    while d*d<=n:
        while n%d==0: f[d]=f.get(d,0)+1; n//=d
        d+=1
    if n>1: f[n]=f.get(n,0)+1
    return f
print("factorization of max:", fac(mx))
