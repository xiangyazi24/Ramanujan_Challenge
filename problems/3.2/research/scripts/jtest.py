#!/usr/bin/env python3
"""J-test: does M(r) preserve a symmetric bilinear form conformally?
M~ = [[P(r), -r^3],[(r+1)^3, 0]] (denominator-cleared transfer matrix).
Seek symmetric J(r) with polynomial entries deg <= DMAX and sign eps in {+1,-1}:
   M~^T J(r+1) M~ = eps * r^3 (r+1)^3 J(r).
(Gauge + determinant argument forces lambda = eps*det, det J(r+1) = eps^2 det J(r).)
Solve the linear system for J's coefficients over QQ via sympy nullspace.
"""
import sympy as sp

r = sp.symbols('r')
P = sp.expand((2*r+1)*(17*r**2+17*r+5))
Mt = sp.Matrix([[P, -r**3],[(r+1)**3, 0]])
DMAX = 12

def run(eps):
    # J(r) = [[A(r), B(r)],[B(r), C(r)]]
    coeffs = []
    A = sum(sp.symbols(f'a{i}')*r**i for i in range(DMAX+1))
    B = sum(sp.symbols(f'b{i}')*r**i for i in range(DMAX+1))
    C = sum(sp.symbols(f'c{i}')*r**i for i in range(DMAX+1))
    unks = [sp.symbols(f'{ch}{i}') for ch in 'abc' for i in range(DMAX+1)]
    J  = sp.Matrix([[A, B],[B, C]])
    J1 = sp.Matrix([[A.subs(r, r+1), B.subs(r, r+1)],[B.subs(r, r+1), C.subs(r, r+1)]])
    E = sp.expand(Mt.T * J1 * Mt - eps * r**3 * (r+1)**3 * J)
    eqs = []
    for expr in (E[0,0], E[0,1], E[1,1]):
        poly = sp.Poly(sp.expand(expr), r)
        eqs.extend(poly.all_coeffs())
    Amat, _ = sp.linear_eq_to_matrix(eqs, unks)
    ns = Amat.nullspace()
    print(f"eps={eps:+d}: nullspace dim = {len(ns)}")
    for v in ns:
        sol = {u: val for u, val in zip(unks, v)}
        JA = sp.simplify(A.subs(sol)); JB = sp.simplify(B.subs(sol)); JC = sp.simplify(C.subs(sol))
        print("  J =", [sp.factor(JA), sp.factor(JB), sp.factor(JC)])
    return ns

for eps in (1, -1):
    run(eps)
