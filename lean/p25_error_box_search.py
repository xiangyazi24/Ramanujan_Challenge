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
PP = [[s.Poly(P[i, j], n, domain=s.QQ) for j in range(3)] for i in range(3)]
SCALE = s.Poly((n + 2)**2, n, domain=s.QQ)
NEXT_SCALE = s.Poly((n + 3)**2, n, domain=s.QQ)


def eventually_nonnegative(poly, start=0):
    shifted = poly.shift(start)
    return all(c >= 0 for c in shifted.all_coeffs())


def test_box(lx, ux, ly, uy, start=0, verbose=False):
    bounds = []
    for X in (lx, ux):
        for Y in (ly, uy):
            d = SCALE*PP[0][0]-X*PP[1][0]-Y*PP[2][0]
            e1 = -SCALE*PP[0][1]+X*PP[1][1]+Y*PP[2][1]
            e2 = -SCALE*PP[0][2]+X*PP[1][2]+Y*PP[2][2]
            tests = [d, NEXT_SCALE*e1-lx*d, ux*d-NEXT_SCALE*e1,
                     NEXT_SCALE*e2-ly*d, uy*d-NEXT_SCALE*e2]
            flags = [eventually_nonnegative(t, start) for t in tests]
            bounds.append(flags)
            if verbose:
                print(X, Y, flags)
                for t, ok in zip(tests, flags):
                    if not ok:
                        print(" bad", s.factor(t.shift(start).as_expr()))
    return all(all(row) for row in bounds)


def test_az_box(la, ua, lz, uz, start=0, verbose=False):
    """Use X=5/4+A/h, Y=2+(8/3)A/h+Z/h^2, h=n+2."""
    h = s.Poly(n + 2, n, domain=s.QQ)
    hp = s.Poly(n + 3, n, domain=s.QQ)
    rows = []
    for A in (la, ua):
        for Z in (lz, uz):
            xnum = s.Rational(5, 4)*h + A
            ynum = 2*h*h + s.Rational(8, 3)*A*h + Z
            d = h**4*PP[0][0] - h*xnum*PP[1][0] - ynum*PP[2][0]
            e1 = -h**4*PP[0][1] + h*xnum*PP[1][1] + ynum*PP[2][1]
            e2 = -h**4*PP[0][2] + h*xnum*PP[1][2] + ynum*PP[2][2]
            na = hp**3*e1 - s.Rational(5, 4)*hp*d
            nz = hp**2*(hp**2*(e2-s.Rational(8, 3)*e1)
                        + s.Rational(4, 3)*d)
            tests = [d, na-la*d, ua*d-na, nz-lz*d, uz*d-nz]
            flags = [eventually_nonnegative(t, start) for t in tests]
            rows.append(flags)
            if verbose:
                print("AZ", A, Z, flags)
                for t, ok in zip(tests, flags):
                    if not ok:
                        print(" bad", s.factor(t.shift(start).as_expr()))
    return all(all(row) for row in rows)


candidates_x = [s.Rational(5,4), s.Rational(6,5), s.Rational(19,16),
                s.Rational(25,16), s.Rational(8,5), s.Rational(5,3)]
candidates_y = [s.Rational(2), s.Rational(19,10), s.Rational(15,7),
                s.Rational(13,6), s.Rational(11,5), s.Rational(9,4)]

for start in (2, 3, 4, 8, 16):
    for la in (s.Rational(1, 6), s.Rational(3, 16), s.Rational(9, 50)):
        for ua in (s.Rational(7, 32), s.Rational(1, 4), s.Rational(1, 3)):
            for lz in (-s.Rational(3, 5), -s.Rational(1, 2), -s.Rational(9, 20)):
                for uz in (-s.Rational(1, 4), -s.Rational(1, 5), 0):
                    if test_az_box(la, ua, lz, uz, start):
                        print("FOUND AZ", start, la, ua, lz, uz, flush=True)
                        raise SystemExit

print("AZ candidate diagnostic", flush=True)
test_az_box(s.Rational(3,16), s.Rational(7,32), -s.Rational(1,2),
            -s.Rational(1,4), 3, True)
