#!/usr/bin/env python3
"""C_h(T) = Res_x(N_h(x) - T*q_h(x)^3, A_h(x)), A_h = N_h' q_h - 3 N_h q_h'.
Check irreducibility over Q for h=2..10; report leading/constant factorizations
and reduction patterns mod small primes (Eisenstein hunting)."""
import sympy as sp
x, T = sp.symbols('x T')
def Ppoly(u): return 34*u**3 + 51*u**2 + 27*u + 5
def Npoly(h):
    if h == 1: return sp.Integer(1)
    a, b = sp.Integer(1), sp.expand(Ppoly(x+1))
    for m in range(2, h): a, b = b, sp.expand(Ppoly(x+m)*b - (x+m)**6*a)
    return b if h >= 2 else a
for h in range(2, 11):
    N = Npoly(h)
    q = sp.prod([(x+j) for j in range(1, h+1)])
    A = sp.expand(sp.diff(N,x)*q - 3*N*sp.diff(q,x))
    C = sp.resultant(sp.expand(N - T*q**3), A, x)
    C = sp.Poly(sp.expand(C), T)
    cont = sp.gcd(list(C.all_coeffs()))
    Cp = sp.Poly([c//cont for c in C.all_coeffs()], T)
    fl = sp.factor_list(Cp.as_expr())
    degs = sorted(sp.degree(f, T) for f,_ in fl[1])
    lead = sp.factorint(abs(Cp.all_coeffs()[0])); const = sp.factorint(abs(Cp.all_coeffs()[-1]))
    print(f"h={h}: degC={C.degree()} factor_degs={degs} lead={dict(lead)} const={dict(const)}")
