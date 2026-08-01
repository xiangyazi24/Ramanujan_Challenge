#!/usr/bin/env python3
"""Refined ledger: strip ALL linear factors (x+c), c integer (singular/forced/classified),
then D_core = sum of remaining radical degrees. Decides the fork on the CORE."""
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
        Bc, Bp = sp.cancel(bj*Bp + aj*Bc), Bc
    return Ac, Bc

HMAX = 10
Dcore = 0; Jcore = 0; lin_deg = 0; rows = []
for k in range(2, HMAX+1):
    for h in range(1, k):
        Ah, Bh = AB(h); Ak, Bk = AB(k)
        psi = sp.cancel((1-Ah)*Bk - (1-Ak)*Bh)
        num = sp.fraction(sp.together(psi))[0]
        fl = sp.factor_list(sp.expand(num))
        core_rad = 0; core_fac = 0; lin = 0
        for b, e in fl[1]:
            pb = sp.Poly(b, x)
            if pb.degree() == 1:
                lin += 1  # classified: mirror/type-II/singular linear family
            else:
                core_rad += pb.degree(); core_fac += 1
        Dcore += core_rad; Jcore += core_fac; lin_deg += lin
        rows.append((h, k, core_rad, core_fac))
print(f"H={HMAX}: D_core={Dcore}  D_core/H^3={Dcore/HMAX**3:.3f}  D_core/H^2={Dcore/HMAX**2:.3f}  "
      f"J_core={Jcore}  J_core/pairs={Jcore/len(rows):.2f}  linear-factor count={lin_deg}")
print("sample (h,k,core_rad,#core_fac):", rows[:10])
print("largest core:", sorted(rows, key=lambda r: -r[2])[:5])
