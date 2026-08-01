#!/usr/bin/env python3
"""Temporary exact search for a scaled projective error box."""

import sympy as s

n = s.symbols("n", nonnegative=True, integer=True)
P = s.Matrix([
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


def eventually_nonnegative(poly, start=0):
    poly = s.Poly(s.cancel(poly), n)
    shifted = s.Poly(poly.as_expr().subs(n, n + start), n)
    return all(c >= 0 for c in shifted.all_coeffs())


def test_box(lx, ux, ly, uy, start=0, verbose=False):
    bounds = []
    for X in (lx, ux):
        for Y in (ly, uy):
            x = X/(n+2)**2
            y = Y/(n+2)**2
            d = P[0, 0]-x*P[1, 0]-y*P[2, 0]
            e1 = -P[0, 1]+x*P[1, 1]+y*P[2, 1]
            e2 = -P[0, 2]+x*P[1, 2]+y*P[2, 2]
            tests = [d, (n+3)**2*e1-lx*d, ux*d-(n+3)**2*e1,
                     (n+3)**2*e2-ly*d, uy*d-(n+3)**2*e2]
            flags = [eventually_nonnegative(t, start) for t in tests]
            bounds.append(flags)
            if verbose:
                print(X, Y, flags)
                for t, ok in zip(tests, flags):
                    if not ok:
                        num = s.cancel(t).as_numer_denom()[0]
                        print(" bad", s.factor(num.subs(n, n+start)))
    return all(all(row) for row in bounds)


candidates_x = [s.Rational(5,4), s.Rational(6,5), s.Rational(19,16),
                s.Rational(25,16), s.Rational(8,5), s.Rational(5,3)]
candidates_y = [s.Rational(2), s.Rational(19,10), s.Rational(15,7),
                s.Rational(13,6), s.Rational(11,5), s.Rational(9,4)]

for start in (0, 1, 2, 4, 8, 16, 32):
    for lx in candidates_x[:3]:
        for ux in candidates_x[3:]:
            if lx >= ux:
                continue
            for ly in candidates_y[:2]:
                for uy in candidates_y[2:]:
                    if ly >= uy:
                        continue
                    if test_box(lx, ux, ly, uy, start):
                        print("FOUND", start, lx, ux, ly, uy, flush=True)
                        raise SystemExit

print("candidate diagnostic", flush=True)
test_box(s.Rational(5,4), s.Rational(25,16), s.Rational(2),
         s.Rational(15,7), 0, True)
