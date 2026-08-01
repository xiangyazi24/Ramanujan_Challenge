#!/usr/bin/env python3
"""Exact ramification of f_h and exact analysis of X_{h,k} for small (h,k):
   - critical points/values of f_h = N_h(r)/D_h(r), D_h = prod (r+j)^3
   - reducibility of F_{h,k}(r,r') = N_h(r) D_k(r') - N_k(r') D_h(r) over Q
   - genus estimate via singularities if irreducible (sympy genus not available;
     use RH with TRUE ramification data of both maps)"""
import sympy as sp
r, s = sp.symbols('r s')

def Ppoly(u): return 34*u**3 + 51*u**2 + 27*u + 5
def Npoly(h):
    if h == 1: return sp.Integer(1)
    N1, N2 = sp.Integer(1), sp.expand(Ppoly(r+1))
    if h == 2: return N2
    a, b = N1, N2
    for m in range(2, h):
        a, b = b, sp.expand(Ppoly(r+m)*b - (r+m)**6*a)
    return b
def Dpoly(h): return sp.prod([(r+j)**3 for j in range(1, h+1)])

for h in [2,3,4,5]:
    N = Npoly(h); D = Dpoly(h)
    f = N/D
    # derivative numerator
    num = sp.expand(sp.diff(N, r)*D - N*sp.diff(D, r))
    num = sp.cancel(num / sp.gcd(num, D))  # remove pole factors if any
    num = sp.expand(sp.simplify(num))
    fac = sp.factor_list(num)
    # critical points = roots of num that are not poles; get their f-values structure
    polys = [(pl, m) for pl, m in fac[1]]
    print(f"h={h}: deg N={sp.degree(N,r)}, crit-numerator deg={sp.degree(num,r)}")
    print(f"   factors of f' numerator: {[(sp.degree(pl,r), m) for pl,m in polys]}")
    # discriminant of N (simple zeros?)
    print(f"   N squarefree: {sp.gcd(N, sp.diff(N,r)) == 1}")

# reducibility of F_{h,k} over Q
for (h,k) in [(2,3),(2,4),(3,4),(2,5),(3,5)]:
    Nh, Dh = Npoly(h), Dpoly(h)
    Nk = Npoly(k).subs(r, s); Dk = Dpoly(k).subs(r, s)
    F = sp.expand(Nh*Dk - Nk*Dh)
    fl = sp.factor_list(F)
    degs = [(sp.degree(pl, r), sp.degree(pl, s)) for pl, _ in fl[1]]
    print(f"X_{h},{k}: factors bidegrees = {degs}")
