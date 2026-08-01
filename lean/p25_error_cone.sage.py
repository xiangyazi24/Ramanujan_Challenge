#!/usr/bin/env sage-python
"""Temporary exact corner checks for the Catalan-error projective cone."""

from sage.all import QQ, PolynomialRing, matrix

R = PolynomialRing(QQ, "n")
n = R.gen()


def positive_matrix(n):
    return matrix(R, 3, 3, [
        (2*n+5)*(n+3)**2*(136*n**4+1424*n**3+5548*n**2+9551*n+6141),
        384*n**6+6384*n**5+44168*n**4+162698*n**3+336377*n**2+369933*n+169011,
        480*n**4+4980*n**3+19210*n**2+32690*n+20730,
        (n+2)**2*(n+3)**2*(4*n+10)*(48*n**3+386*n**2+1017*n+879),
        (n+2)**2*(272*n**5+3848*n**4+21732*n**3+61184*n**2+85761*n+47808),
        (n+2)**2*(320*n**3+2540*n**2+6610*n+5640),
        (4*n+10)*(n+2)**2*(n+3)**2*(32*n**4+302*n**3+1037*n**2+1530*n+813),
        (n+2)**2*(192*n**6+2984*n**5+19116*n**4+64452*n**3+120256*n**2+117279*n+46476),
        (n+2)**2*(16*n**5+408*n**4+2912*n**3+8884*n**2+12254*n+6240),
    ])


A = positive_matrix(n)
# A shrinking rectangle suggested by the first hundred exact projective
# iterates.  The constants are deliberately loose.
lx = QQ(5)/4 - QQ(1)/4/(n+1)
ux = QQ(5)/4 - QQ(1)/16/(n+1)
ly = 2 - 1/(n+1)**2
uy = 2 - QQ(1)/4/(n+1)**2
lx_next, ux_next, ly_next, uy_next = [
    value(n=n+1) for value in [lx, ux, ly, uy]
]
vertices = [(x, y) for x in [lx, ux] for y in [ly, uy]]


def t(j, x, y):
    h = (n + 1)**2 * (4*n + 10) / (4*n + 3)
    return h * A[0, j] - x * A[1, j] - y * A[2, j]


h_next = (n + 2)**2 * (4*n + 14) / (4*n + 7)


tests = {
    "denominator": lambda x, y: t(0, x, y),
    "x_lower": lambda x, y: h_next*(-t(1,x,y))-lx_next*t(0,x,y),
    "x_upper": lambda x, y: ux_next*t(0,x,y)-h_next*(-t(1,x,y)),
    "y_lower": lambda x, y: h_next*(-t(2,x,y))-ly_next*t(0,x,y),
    "y_upper": lambda x, y: uy_next*t(0,x,y)-h_next*(-t(2,x,y)),
}

for name, expression in tests.items():
    print(name)
    for x, y in vertices:
        rational = expression(x, y)
        polynomial = R(rational.numerator())
        shifted = polynomial(n=n+2)
        coefficients = polynomial.list()
        print(x, y, "degree", polynomial.degree(),
              "shifted_coefficient_min", min(shifted.list()),
              "all_shifted_nonnegative", all(c >= 0 for c in shifted.list()),
              "factor", polynomial.factor())
