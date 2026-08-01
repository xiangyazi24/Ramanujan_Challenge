#!/usr/bin/env sage-python
"""Test the rational difference-module equivalence with Sym^2 Wilson--Pade."""

import sys

from sage.all import QQ, PolynomialRing, matrix, vector

sys.path.insert(0, "/Users/huangx/Library/SageMath-10-9/lib/python3.14/site-packages")
from ore_algebra import OreAlgebra


R = PolynomialRing(QQ, "n")
n = R.gen()
K = R.fraction_field()
OA = OreAlgebra(K, names=("Sn",))


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
    # Transition from the two Wilson bounds at index m=n+2 to index m+1.
    d0 = 8*(n+3)**2*(4*n+11)**2
    d1 = d0*(4*n+13)**2
    return matrix(K, [
        [
            (4*n+9)**2*(40*n**2+228*n+325)/d0,
            (1536*n**4+16512*n**3+66496*n**2+118896*n+79641)/d0,
        ],
        [
            (4*n+9)**2*(1536*n**4+18432*n**3+82880*n**2+165504*n+123841)/d1,
            (59392*n**6+1012736*n**5+7184384*n**4+27140352*n**3+
             57583336*n**2+65059404*n+30580677)/d1,
        ],
    ])


def symmetric_square(value):
    a, b, c, d = value[0, 0], value[0, 1], value[1, 0], value[1, 1]
    return matrix(K, [
        [a*a, 2*a*b, b*b],
        [a*c, a*d+b*c, b*d],
        [c*c, 2*c*d, d*d],
    ])


def sigma(value, amount=1):
    return value(n=n+amount)


def flatten(value):
    return vector(K, [value[i, j] for i in range(value.nrows())
                      for j in range(value.ncols())])


def recover(left, right):
    dimension = left.nrows()*right.nrows()
    basis = []
    for index in range(dimension):
        value = matrix(K, left.nrows(), right.nrows(), 0)
        value[index // right.nrows(), index % right.nrows()] = 1
        basis.append(value)
    transition = matrix(
        K, dimension, dimension,
        lambda i, j: flatten(left.inverse()*basis[j]*right)[i],
    )
    for seed_index in range(dimension):
        print("seed", seed_index, flush=True)
        seed = [K.zero()]*dimension
        seed[seed_index] = K.one()
        rows = [vector(K, seed)]
        for _ in range(dimension):
            rows.append(vector(K, [sigma(x) for x in rows[-1]])*transition)
        krylov = matrix(K, rows[:dimension])
        if krylov.det() == 0:
            print("noncyclic", flush=True)
            continue
        relation = krylov.transpose().solve_right(-rows[dimension])
        operator = OA(list(relation)+[1]).normalize()
        print("operator", operator.order(), operator.degree(), flush=True)
        solutions = operator.rational_solutions()
        print("solutions", len(solutions), flush=True)
        for solution in solutions:
            scalar = solution[0]
            rhs = vector(K, [sigma(scalar, shift) for shift in range(dimension)])
            gauge = matrix(K, 3, 3, list(krylov.solve_right(rhs)))
            assert gauge.apply_map(sigma) == left.inverse()*gauge*right
            print("rank", gauge.rank(), flush=True)
            for row in gauge.rows():
                print([entry.factor() for entry in row], flush=True)
            if gauge.rank() == 3:
                return gauge
    return None


target = challenge_matrix().transpose()
source = symmetric_square(wilson_matrix())
print("determinants", target.det().factor(), source.det().factor(), flush=True)
answer = recover(target, source)
print("FOUND", answer is not None, flush=True)
