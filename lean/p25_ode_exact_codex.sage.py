#!/usr/bin/env sage-python
"""Exact recurrence-to-differential conversion for the normalized P2.5 module."""

import sys

from sage.all import QQ, PolynomialRing, matrix, vector

sys.path.insert(0, "/Users/huangx/Library/SageMath-10-9/lib/python3.14/site-packages")
from ore_algebra import OreAlgebra


Rn = PolynomialRing(QQ, "n")
n = Rn.gen()
Kn = Rn.fraction_field()
RS = OreAlgebra(Kn, names=("Sn",))


def sigma(value, amount=1):
    return value(n=n + amount)


def normalized_challenge():
    entries = [
        (-2*n-5)*(n+3)**2*(136*n**4+1424*n**3+5548*n**2+9551*n+6141),
        384*n**6+6384*n**5+44168*n**4+162698*n**3+336377*n**2+369933*n+169011,
        -(480*n**4+4980*n**3+19210*n**2+32690*n+20730),
        (n+2)**2*(n+3)**2*(4*n+10)*(48*n**3+386*n**2+1017*n+879),
        (n+2)**2*(-272*n**5-3848*n**4-21732*n**3-61184*n**2-85761*n-47808),
        (n+2)**2*(320*n**3+2540*n**2+6610*n+5640),
        (-4*n-10)*(n+2)**2*(n+3)**2*(32*n**4+302*n**3+1037*n**2+1530*n+813),
        (n+2)**2*(192*n**6+2984*n**5+19116*n**4+64452*n**3+120256*n**2+117279*n+46476),
        (n+2)**2*(-16*n**5-408*n**4-2912*n**3-8884*n**2-12254*n-6240),
    ]
    delta = -2*(n+2)**2*(n+3)**2*(2*n+5)*(2*n+7)**2
    return matrix(Kn, 3, 3, entries) / delta


def scalar_operator(column=0):
    transition = normalized_challenge()
    basis = vector(Kn, [1 if i == column else 0 for i in range(3)])
    columns = [basis]
    product = matrix.identity(Kn, 3)
    for shift in range(3):
        product *= transition.apply_map(lambda value: sigma(value, shift))
        columns.append(product*basis)
    cyclic = matrix(Kn, 3, 3, lambda row, col: columns[col][row])
    relation = cyclic.solve_right(-columns[3])
    return RS(list(relation)+[1]).normalize()


rec = scalar_operator()
print("REC", flush=True)
for coefficient in rec:
    print(coefficient.factor(), flush=True)

Rz = PolynomialRing(QQ, "z")
z = Rz.gen()
Kz = Rz.fraction_field()
DO = OreAlgebra(Kz, names=("Dz",))
Dz = DO.gen()
theta = DO(z) * Dz


def eval_at_operator(rational, operator):
    denominator = Rn(rational.denominator())
    if denominator.degree() != 0:
        raise ValueError("nonconstant recurrence denominator")
    polynomial = Rn(rational.numerator())
    value = DO.zero()
    for coefficient in reversed(polynomial.list()):
        value = value * operator + DO(coefficient)
    return DO(Kz(QQ.one() / QQ(denominator))) * value


L = DO.zero()
for shift, coefficient in enumerate(rec):
    L += DO(z**(3-shift)) * eval_at_operator(coefficient, theta - shift)
print("DIFF order", L.order(), "degree", L.degree(), flush=True)
print("leading", L[L.order()].factor(), flush=True)
for derivative_order, coefficient in enumerate(L):
    print("D", derivative_order, coefficient.factor(), flush=True)

for side, method in [("right", L.right_factors)]:
    for order in range(1, 4):
        print(side, "factors", order, flush=True)
        try:
            factors = method(order=order, infolevel=1)
            print(factors, flush=True)
        except Exception as error:
            print("factor error", type(error).__name__, error, flush=True)
