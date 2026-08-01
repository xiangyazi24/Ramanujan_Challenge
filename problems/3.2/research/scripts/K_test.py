#!/usr/bin/env python3
"""K(d1,d2) = d2^3*[(1-A_{d1}(-d2))*b_{d2-1} - b_{d2}*B_{d1}(-d2)] over Q.
Check: K != 0? gcd/divisibility by (b_{d2-d1}-1)? Verify K is the pole-3 coeff of Psi at -d2."""
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
for n in range(1, 40):
    num = (2*n+1)*(17*n*n+17*n+5)*b[n] - n**3*b[n-1]
    q, rm = divmod(num, (n+1)**3); assert rm == 0
    b.append(q)

# numeric verification of K as pole coeff: K ?= lim (x+d2)^3 Psi(x) at x=-d2, test via x close (exact rational Laurent by evaluating at x=-d2+t for symbolic t? do numeric check with two sample x offsets using exact rationals and Richardson-free direct: multiply and cancel)
import sympy as sp
X = sp.symbols('x')
def AB_sym(k):
    a = sp.expand((2*X+1)*(17*X**2+17*X+5))/(X+1)**3
    beta = -X**3/(X+1)**3
    Ap, Bp = sp.Integer(1), sp.Integer(0)
    Ac, Bc = a, beta
    for j in range(1, k):
        aj = a.subs(X, X+j); bj = beta.subs(X, X+j)
        Ac, Ap = sp.cancel(aj*Ac + bj*Ap), Ac
        Bc, Bp = sp.cancel(aj*Bc + bj*Bp), Bc
    return sp.cancel(Ac), sp.cancel(Bc)

print(f"{'d1':>3}{'d2':>4} {'K':>28} {'b_delta-1':>16} {'K % (b_delta-1)==0?':>20} {'gcd':>12}")
for (d1,d2) in [(1,2),(1,3),(2,3),(1,4),(2,4),(3,4),(2,5),(3,5),(2,6),(3,7)]:
    A1v, B1v = ABat(d1, -d2)
    K = F(d2**3) * ((1 - A1v)*F(b[d2-1]) - F(b[d2])*B1v)
    delta = d2 - d1
    bd1 = b[delta] - 1
    # K rational; clear denominator
    Knum, Kden = K.numerator, K.denominator
    divis = (Knum % bd1 == 0) if bd1 != 0 else None
    g = gcd(abs(Knum), abs(bd1))
    print(f"{d1:>3}{d2:>4} {str(K)[:28]:>28} {bd1:>16} {str(divis):>20} {g:>12}")
    # cross-check with symbolic pole coefficient for small cases
    if d2 <= 4:
        Ad1s, Bd1s = AB_sym(d1); Ad2s, Bd2s = AB_sym(d2)
        Psi = sp.cancel((1-Ad1s)*Bd2s - (1-Ad2s)*Bd1s)
        pole = sp.limit(sp.cancel(Psi*(X+d2)**3), X, -d2)
        match = sp.nsimplify(pole) == sp.nsimplify(sp.Rational(Knum, Kden))
        print(f"      symbolic pole coeff = {pole}  matches K: {match}")

print("\nLarger sweep: gcd(b_delta-1, Knum) and gcd(b_{d1-1}, Knum)")
mx1 = mx2 = 1
for d1 in range(1, 13):
    for d2 in range(d1+1, 16):
        A1v, B1v = ABat(d1, -d2)
        K = F(d2**3) * ((1 - A1v)*F(b[d2-1]) - F(b[d2])*B1v)
        Kn = abs(K.numerator)
        if Kn == 0: print("K=0 at", d1, d2); continue
        g1 = gcd(Kn, abs(b[d2-d1]-1)); g2 = gcd(Kn, abs(b[d1-1])) if d1>=1 else 0
        mx1 = max(mx1, g1); mx2 = max(mx2, g2)
        if g1 > 1000 or g2 > 1000: print(f"  big gcd at ({d1},{d2}): g1={g1} g2={g2}")
print("max gcd(K, b_delta-1) =", mx1, "; max gcd(K, b_{d1-1}) =", mx2, "over d1<13,d2<16")
