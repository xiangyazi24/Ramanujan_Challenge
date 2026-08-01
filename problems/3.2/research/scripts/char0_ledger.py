#!/usr/bin/env python3
"""Characteristic-zero squarefree-degree ledger (R12b's decisive fork):
D_H = sum deg rad(Phi_{h,k}) over Q after removing mirror factor (2x+h+k+1) and known
singular factors. D ~ H^3 => all-prime (RES) false, need scale-coupled Frobenius sieve;
D << H^2 => (RES) trivially true for ALL primes."""
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

HMAX = 10
D = 0; J = 0; tot_deg = 0; details = []
for k in range(2, HMAX+1):
    for h in range(1, k):
        Ah, Bh = AB(h); Ak, Bk = AB(k)
        psi = sp.cancel((1-Ah)*Bk - (1-Ak)*Bh)
        num = sp.fraction(sp.together(psi))[0]
        f = sp.Poly(sp.expand(num), x)
        # remove content
        f = f.primitive()[1]
        if (h-k) % 2 == 0:
            mirror = sp.Poly(2*x + h + k + 1, x)
            q, r = sp.div(f, mirror)
            assert r.is_zero, (h,k)
            f = q
        # remove even-lag type-II centers: for even d in {h,k}, factor (2x+d+1)? type-II center r=-(d+1)/2:
        for d in (h, k):
            if d % 2 == 0:
                lin = sp.Poly(2*x + d + 1, x)
                q, r = sp.div(f, lin)
                if r.is_zero: f = q
        tot_deg += f.degree()
        sf = sp.Poly(sp.factor_list(f.as_expr())[1] and sp.prod([b for b,e in sp.factor_list(f.as_expr())[1]]) or 1, x) if f.degree()>0 else f
        fl = sp.factor_list(f.as_expr())
        raddeg = sum(sp.Poly(b, x).degree() for b, e in fl[1])
        nfac = len(fl[1])
        D += raddeg; J += nfac
        details.append((h,k,f.degree(),raddeg,nfac))
print(f"H={HMAX}: total residual degree={tot_deg}, D (rad degree)={D}, D/H^3={D/HMAX**3:.3f}, "
      f"J (#factors)={J}, J/H^2={J/HMAX**2:.3f}")
print("sample (h,k,deg,raddeg,#fac):", details[:8], "...", details[-4:])
sq_drop = sum(1 for _,_,dg,rd,_ in details if rd < dg)
print(f"pairs with nonsquarefree residual (rad < deg): {sq_drop} / {len(details)}")
