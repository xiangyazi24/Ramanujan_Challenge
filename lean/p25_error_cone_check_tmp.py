#!/usr/bin/env python3
"""Exact SymPy verification of the moving Catalan-error cone."""

import sympy as s

n = s.symbols("n")
Q = s.Rational


def positive_matrix():
    return s.Matrix([
        [(2*n+5)*(n+3)**2*(136*n**4+1424*n**3+5548*n**2+9551*n+6141),
         384*n**6+6384*n**5+44168*n**4+162698*n**3+336377*n**2+369933*n+169011,
         480*n**4+4980*n**3+19210*n**2+32690*n+20730],
        [(n+2)**2*(n+3)**2*(4*n+10)*(48*n**3+386*n**2+1017*n+879),
         (n+2)**2*(272*n**5+3848*n**4+21732*n**3+61184*n**2+85761*n+47808),
         (n+2)**2*(320*n**3+2540*n**2+6610*n+5640)],
        [(4*n+10)*(n+2)**2*(n+3)**2*(32*n**4+302*n**3+1037*n**2+1530*n+813),
         (n+2)**2*(192*n**6+2984*n**5+19116*n**4+64452*n**3+120256*n**2+117279*n+46476),
         (n+2)**2*(16*n**5+408*n**4+2912*n**3+8884*n**2+12254*n+6240)],
    ])


A = positive_matrix()
h = (n + 1)**2 * (4*n + 10) / (4*n + 3)
hn = (n + 2)**2 * (4*n + 14) / (4*n + 7)
lx = Q(5, 4) - Q(1, 4)/(n + 1)
ux = Q(5, 4) - Q(1, 16)/(n + 1)
ly = 2 - 1/(n + 1)**2
uy = 2 - Q(1, 4)/(n + 1)**2
lxn, uxn, lyn, uyn = [v.subs(n, n + 1) for v in (lx, ux, ly, uy)]


def shifted_nonnegative(expr, start=1):
    num, den = s.cancel(expr).as_numer_denom()
    poly = s.Poly(s.expand(num.subs(n, n + start)), n)
    return all(c >= 0 for c in poly.all_coeffs()), s.factor(num), s.factor(den), poly


for corner, (x, y) in enumerate((
        (lx, ly), (lx, uy), (ux, ly), (ux, uy))):
    t = [s.cancel(h*A[0, j] - x*A[1, j] - y*A[2, j]) for j in range(3)]
    tests = (
        ("den", t[0]),
        ("xlo", hn*(-t[1]) - lxn*t[0]),
        ("xhi", uxn*t[0] - hn*(-t[1])),
        ("ylo", hn*(-t[2]) - lyn*t[0]),
        ("yhi", uyn*t[0] - hn*(-t[2])),
    )
    print("corner", corner)
    for label, expr in tests:
        ok, num, den, poly = shifted_nonnegative(expr, start=1)
        print(label, ok, "deg", poly.degree(), "min", min(poly.all_coeffs()),
              "den", den)
        if not ok:
            bad = [(poly.degree() - i, c) for i, c in enumerate(poly.all_coeffs()) if c < 0]
            print(" bad", bad[:10], "factor", num)
