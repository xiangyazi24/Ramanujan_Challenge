#!/usr/bin/env sage-python
"""Eliminate the six-dimensional Lima partial-sum comparison systems."""

from sage.all import QQ, PolynomialRing, matrix, vector


R = PolynomialRing(QQ, "n")
n = R.gen()
K = R.fraction_field()


def challenge():
    return matrix(K, 3, 3, [
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


def shift(value, amount):
    return value.apply_map(lambda entry: entry(n=n+amount))


def lima_ratio(m):
    return -(m+1)**3*(3*m+5)/((2*m+3)**3*(3*m+2))


def transition(kind):
    m = 2*n+4
    a = lima_ratio(m)
    s = a*lima_ratio(m+1)
    forcing = 1+a if kind == "lower" else -(a+s)
    C = challenge()
    zero = matrix(K, 3, 3, 0)
    return C.augment(zero).stack((forcing*C).augment(s*C))


def recurrence(kind, coordinate):
    seed = vector(K, [1 if i == coordinate else 0 for i in range(6)])
    columns = [seed]
    product = matrix.identity(K, 6)
    for step in range(6):
        product *= shift(transition(kind), step)
        columns.append(product*seed)
    basis = matrix(K, 6, 6, lambda row, col: columns[col][row])
    coefficients = basis.solve_right(columns[6])
    print("kind", kind, "coordinate", coordinate, flush=True)
    print("det", basis.det().factor(), flush=True)
    for index, coefficient in enumerate(coefficients):
        print(index, coefficient.factor(), flush=True)


recurrence("lower", 2)
recurrence("upper", 0)
