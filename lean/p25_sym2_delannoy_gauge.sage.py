#!/usr/bin/env sage-python
"""Recover a rational row-module gauge from Sym^2 Delannoy to P2.5."""

import sys

from sage.all import QQ, PolynomialRing, matrix, vector

sys.path.insert(0, "/Users/huangx/Library/SageMath-10-9/lib/python3.14/site-packages")
from ore_algebra import OreAlgebra


R = PolynomialRing(QQ, "n")
n = R.gen()
K = R.fraction_field()
OA = OreAlgebra(K, names=("Sn",))


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


def symmetric_square_row(value):
    p, q, r, s = value[0, 0], value[0, 1], value[1, 0], value[1, 1]
    return matrix(K, [
        [p*p, p*q, q*q],
        [2*p*r, p*s+q*r, 2*q*s],
        [r*r, r*s, s*s],
    ])


def delannoy():
    return matrix(K, [[0, -(n+1)/(n+2)], [1, 3*(2*n+3)/(n+2)]])


def flatten(value):
    return vector(K, [value[i, j] for i in range(3) for j in range(3)])


def recover(left, right):
    basis = []
    for index in range(9):
        value = matrix(K, 3, 3, 0)
        value[index // 3, index % 3] = 1
        basis.append(value)
    transition = matrix(
        K, 9, 9,
        lambda i, j: flatten(left.inverse()*basis[j]*right)[i],
    )
    for seed_index in range(9):
        print("seed", seed_index, flush=True)
        seed = [K.zero()]*9
        seed[seed_index] = K.one()
        rows = [vector(K, seed)]
        for _ in range(9):
            rows.append(vector(K, [sigma(x) for x in rows[-1]])*transition)
        krylov = matrix(K, rows[:9])
        if krylov.det() == 0:
            print("noncyclic", flush=True)
            continue
        relation = krylov.transpose().solve_right(-rows[9])
        operator = OA(list(relation)+[1]).normalize()
        print("operator", operator.order(), operator.degree(), flush=True)
        solutions = operator.rational_solutions()
        print("solutions", len(solutions), flush=True)
        for solution in solutions:
            scalar = solution[0]
            rhs = vector(K, [sigma(scalar, shift) for shift in range(9)])
            gauge = matrix(K, 3, 3, list(krylov.solve_right(rhs)))
            assert gauge.apply_map(sigma) == left.inverse()*gauge*right
            print("rank", gauge.rank(), flush=True)
            for row in gauge.rows():
                print([entry.factor() for entry in row], flush=True)
            if gauge.rank() == 3:
                return gauge
    return None


source = symmetric_square_row(delannoy())
target = challenge()
print("det", source.det().factor(), target.det().factor(), flush=True)
answer = recover(source, target)
print("FOUND", answer is not None, flush=True)
