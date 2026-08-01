#!/usr/bin/env sage-python
"""Temporary exact scalar-operator factorization for P2.5."""

import sys

from sage.all import QQ, PolynomialRing, matrix, vector

sys.path.insert(0, "/Users/huangx/Library/SageMath-10-9/lib/python3.14/site-packages")
from ore_algebra import OreAlgebra


R = PolynomialRing(QQ, "n")
n = R.gen()
K = R.fraction_field()
OA = OreAlgebra(K, names=("Sn",))
Sn = OA.gen()


def sigma(value, amount=1):
    return value(n=n+amount)


def challenge():
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
    return matrix(K, 3, 3, entries)/delta


def scalar_operator(column):
    transition = challenge()
    basis = vector(K, [1 if i == column else 0 for i in range(3)])
    columns = [basis]
    product = matrix.identity(K, 3)
    for shift in range(3):
        product *= transition.apply_map(lambda value: sigma(value, shift))
        columns.append(product*basis)
    cyclic = matrix(K, 3, 3, lambda row, col: columns[col][row])
    relation = cyclic.solve_right(-columns[3])
    operator = OA(list(relation)+[1]).normalize()
    return operator


for column in range(1):
    print("column", column, flush=True)
    operator = scalar_operator(column)
    print("operator order", operator.order(), "degree", operator.degree(), flush=True)
    print("coefficients", flush=True)
    for coefficient in operator:
        print(coefficient.factor(), flush=True)
    for side, method in [("left", operator.left_factors),
                         ("right", operator.right_factors)]:
        for order in (1, 2):
            print(side, "factors", order, flush=True)
            try:
                print(method(order=order, infolevel=1), flush=True)
            except Exception as error:
                print("factor error", type(error).__name__, error, flush=True)
    print("hypergeometric solutions", flush=True)
    try:
        print(operator.generalized_series_solutions(2), flush=True)
    except Exception as error:
        print("series error", type(error).__name__, error, flush=True)
