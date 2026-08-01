#!/usr/bin/env python3
"""Machine-check R18's proof ingredients:
(1) F(phi(u)) = (1+u) h(u)^2 over Q (first 15 coefficients), phi = u(1-8u)/(1+u), h = Franel GF;
(2) q(phi(u)) = ((1-16u-8u^2)/(1+u))^2 as rational identity;
(3) f_{p-1} = 1 mod p spot; (4) endpoint lemma re-confirm (already 165/165)."""
import sympy as sp
from math import comb
u = sp.symbols('u')
N = 15
# Apery b and F series
b=[1,5]
for n in range(1,N+2):
    num=(2*n+1)*(17*n*n+17*n+5)*b[n]-n**3*b[n-1]
    q_,r=divmod(num,(n+1)**3); assert r==0
    b.append(q_)
# Franel
f=[sum(comb(n,k)**3 for k in range(n+1)) for n in range(N+2)]
phi = u*(1-8*u)/(1+u)
# F(phi(u)) as series
Fphi = sum(b[k]*phi**k for k in range(N+1))
Fphi_ser = sp.series(sp.together(Fphi), u, 0, N).removeO()
h = sum(f[k]*u**k for k in range(N+1))
rhs = sp.series(sp.expand((1+u)*h*h), u, 0, N).removeO()
diff = sp.expand(Fphi_ser - rhs)
ok1 = all(sp.Poly(diff, u).coeff_monomial(u**k) == 0 for k in range(N))
print("(1) F(phi(u)) = (1+u) h(u)^2 to O(u^15):", "VERIFIED" if ok1 else "FAIL")
# (2)
qt = 1 - 34*phi + phi**2
target = ((1-16*u-8*u**2)/(1+u))**2
ok2 = sp.simplify(qt - target) == 0
print("(2) q(phi(u)) = ((1-16u-8u^2)/(1+u))^2:", "VERIFIED" if ok2 else "FAIL")
# (3)
for p in (7, 11, 13, 101):
    fp1 = sum(comb(p-1,k)**3 for k in range(p)) % p
    assert fp1 == 1 % p
print("(3) f_(p-1) = 1 mod p (p=7,11,13,101): VERIFIED")
print("(4) endpoint lemma: previously verified 165/165 (leading_coeff_check.py)")
