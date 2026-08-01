#!/usr/bin/env python3
"""Premise check for the alpha=1/2 route: are gcd(numB_{d1}, numB_{d2}) over Q of bounded degree?
B_d(x) = -x^3/(x+1)^3 A_{d-1}(x+1). Common roots mod p of two B-numerators (beyond classified
linear factors) require p | Res — the Frobenius-ledger structure. Also factor A_d over Q: one big
irreducible like the Psi cores?"""
import sympy as sp
from functools import lru_cache
x = sp.symbols('x')

@lru_cache(None)
def AB(d):
    a = sp.expand((2*x+1)*(17*x**2+17*x+5))/(x+1)**3
    beta = -x**3/(x+1)**3
    Ap, Bp = sp.Integer(1), sp.Integer(0)
    Ac, Bc = a, beta
    for j in range(1, d):
        aj = a.subs(x, x+j); bj = beta.subs(x, x+j)
        Ac, Ap = sp.cancel(aj*Ac + bj*Ap), Ac
        Bc, Bp = sp.cancel(aj*Bc + bj*Bp), Bc
    return Ac, Bc

nums = {}
for d in range(1, 13):
    _, Bd = AB(d)
    n = sp.fraction(sp.together(Bd))[0]
    nums[d] = sp.Poly(sp.expand(n), x).primitive()[1]

print("A) pairwise gcd degrees over Q (d1<d2<=12):")
mx = 0
for d1 in range(1, 13):
    for d2 in range(d1+1, 13):
        g = sp.gcd(nums[d1].as_expr(), nums[d2].as_expr())
        dg = sp.Poly(g, x).degree() if g != 1 else 0
        mx = max(mx, dg)
print("   max pairwise gcd degree:", mx)

print("B) factorization shape of B-numerators (linear-classified vs big core):")
for d in (4, 7, 10, 12):
    fl = sp.factor_list(nums[d].as_expr())
    shape = sorted(sp.Poly(b, x).degree() for b, e in fl[1])
    print(f"   d={d}: deg={nums[d].degree()}  factor degrees={shape}")
