#!/usr/bin/env python3
"""CRON_Mh_galois.py — Galois maximality certificates for the symmetrized gap polynomials M_h.

Pipeline (appendix AH.5):
  N_h = gap continuant (N_1=1, N_2=A_1, N_{d+1}=A_d N_d - B_d N_{d-1};
        A_d = 34(X+d)^3+51(X+d)^2+27(X+d)+5, B_d = (X+d)^6).
  N-hat_h = primitive part (remove forced factor (2X+h+1) for even h).
  Reflection N_h(-h-1-X) = ±N_h(X)  ==>  N-hat_h(X) = c * M_h(T^2)|_{T=2X+h+1}
  (checked: substitution X=(T-h-1)/2 yields an even polynomial in T).
  M_h in Z[Y] primitive = the true Galois object (Q6463/AH.1 Capelli framing).

Certificates that Gal(M_h) = S_n (n = deg M_h):
  - h=2..5 (n<=6): sympy galois_group exact  -> S_1, S_3, S_4, S_6.
  - h=6..11 (n=7..15): Dedekind cycle types mod primes:
      (a) irreducible over Q  ==> transitive;
      (b) squarefree factorization type [q,1,...] with q prime > n/2  ==> primitive
          (classical lemma: transitive + prime cycle of length > n/2);
      (c) type with EXACTLY ONE even part equal to 2 (rest odd)  ==> an odd power
          is a transposition;
      (d) primitive + transposition  ==> S_n (Jordan).
Run: python3 CRON_Mh_galois.py   (expected: S_n verdict for every h=2..11)
"""
import sys
from math import gcd
from sympy import Poly, symbols, ZZ, div, factor_list, expand, Rational, isprime

X, T, Y = symbols('X T Y')

def P_at(e):
    t = X + e
    return 34*t**3 + 51*t**2 + 27*t + 5

def build_N(hmax):
    N = {1: Poly(1, X, domain=ZZ), 2: Poly(P_at(1), X, domain=ZZ)}
    for d in range(2, hmax):
        N[d+1] = Poly(P_at(d), X, domain=ZZ)*N[d] - Poly((X+d)**6, X, domain=ZZ)*N[d-1]
    return N

def M_of(N, h):
    f = N[h]
    if h % 2 == 0:
        q, r = div(f, Poly(2*X + h + 1, X, domain=ZZ), domain='QQ')
        assert r.is_zero, f"forced factor missing at h={h}"
        f = Poly(q, X)
    g = expand(f.as_expr().subs(X, (T - h - 1)/Rational(2)))
    gp = Poly(g, T)
    assert all(c == 0 for (m,), c in gp.terms() if m % 2 == 1), f"h={h}: not even in T"
    coeffs = {m//2: c for (m,), c in gp.terms() if m % 2 == 0}
    cs = [coeffs.get(k, 0) for k in range(max(coeffs), -1, -1)]
    L = 1
    for c in cs:
        if hasattr(c, 'q'): L = L*c.q//gcd(L, c.q)
    csz = [int(c*L) for c in cs]
    cont = 0
    for c in csz: cont = gcd(cont, abs(c))
    return Poly([c//cont for c in csz], Y, domain=ZZ)

def certify(Mz, pmax=6000):
    n = Mz.degree()
    prim = transp = None
    p = 3
    while p < pmax and not (prim and transp):
        if isprime(p) and Mz.all_coeffs()[0] % p != 0:
            try:
                _, facs = factor_list(Mz.as_expr(), modulus=p)
                if all(m == 1 for _, m in facs):
                    degs = sorted(Poly(g, Y).degree() for g, _ in facs)
                    if [d for d in degs if d % 2 == 0] == [2]:
                        transp = transp or (p, tuple(degs))
                    for d in degs:
                        if isprime(d) and d > n/2 and degs.count(d) == 1 \
                           and all(dd == 1 or dd == d for dd in degs):
                            prim = prim or (p, d)
            except Exception:
                pass
        p += 2
    return prim, transp

def main():
    N = build_N(14)
    fails = 0
    for h in range(2, 12):
        Mz = M_of(N, h); n = Mz.degree()
        _, facs = factor_list(Mz.as_expr())
        irr = (len(facs) == 1 and facs[0][1] == 1)
        if n <= 6:
            if n == 1:
                verdict = "S_1 (linear)"
            else:
                from sympy.polys.numberfields.galoisgroups import galois_group
                G, _ = galois_group(Mz, by_name=True)
                verdict = str(G).split('.')[-1]
            ok = irr
        else:
            prim, transp = certify(Mz)
            ok = irr and prim and transp
            verdict = f"S_{n} [prim-cycle@{prim}, transp-type@{transp}]" if ok else "UNDECIDED"
        print(f"h={h:2d} deg M_h={n:2d} irred={irr} => {verdict}")
        if not ok: fails += 1
    print("ALL S_n" if fails == 0 else f"{fails} UNDECIDED")
    return fails

if __name__ == "__main__":
    sys.exit(main())
