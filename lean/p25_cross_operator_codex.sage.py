#!/usr/bin/env sage-python
"""Construct exact scalar operators for challenge/Wilson tensor coordinates."""

import sys

from sage.all import QQ, PolynomialRing, matrix, vector

sys.path.insert(0, "/Users/huangx/Library/SageMath-10-9/lib/python3.14/site-packages")
from ore_algebra import OreAlgebra


R = PolynomialRing(QQ, "n")
n = R.gen()
K = R.fraction_field()
A = OreAlgebra(K, "Sn")


def sigma(value, amount=1):
    return value(n=n+amount)


def challenge_matrix():
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
    return matrix(K, 3, 3, entries) / delta


def wilson_matrix():
    m = n + 2
    a4 = 12265 + 29296*m + 26176*m**2 + 10368*m**3 + 1536*m**4
    c4 = 313 + 1904*m + 4288*m**2 + 4224*m**3 + 1536*m**4
    b10 = (
        111992515 + 1144683736*m + 5147619352*m**2
        + 13412393984*m**3 + 22433518592*m**4
        + 25185342464*m**5 + 19235018752*m**6
        + 9876373504*m**7 + 3265527808*m**8
        + 628359168*m**9 + 53477376*m**10
    )
    a = 4*(m+1)**2*(4*m+1)**2*(4*m+3)**2*a4
    c = 4*(m+2)**2*(4*m+5)**2*(4*m+7)**2*c4
    return matrix(K, [[0, -a/c], [1, b10/c]])


transition = challenge_matrix().tensor_product(wilson_matrix())


def coordinate_operator(coordinate):
    column = vector(K, [1 if index == coordinate else 0 for index in range(6)])
    columns = [column]
    product = matrix.identity(K, 6)
    for shift in range(6):
        product *= transition.apply_map(lambda value: sigma(value, shift))
        columns.append(product * column)
    basis = matrix(K, 6, 6, lambda row, col: columns[col][row])
    relation = basis.solve_right(-columns[6])
    operator = A(list(relation) + [1]).normalize()
    return operator


for coordinate in (0, 2, 4):
    print("coordinate", coordinate, flush=True)
    operator = coordinate_operator(coordinate)
    print("order", operator.order(), flush=True)
    for shift in range(operator.order()+1):
        coefficient = operator[shift]
        numerator = coefficient.numerator()
        denominator = coefficient.denominator()
        print(
            "c", shift,
            "numdeg", numerator.degree(),
            "dendeg", denominator.degree(),
            "numfactor", numerator.factor(),
            "denfactor", denominator.factor(),
            flush=True,
        )
    try:
        print("factor", operator.factor(), flush=True)
    except Exception as error:
        print("factor error", repr(error), flush=True)
