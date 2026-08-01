#!/usr/bin/env sage-python
"""Exact rational quotient search from the P2.5 module to the Wilson module."""

from sage.all import QQ, PolynomialRing, matrix, vector

import sys
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


def sigma(value, amount=1):
    return value(n=n+amount)


def flatten(value):
    return vector(K, [value[i, j] for i in range(value.nrows())
                      for j in range(value.ncols())])


def recover(left, right):
    rows, columns = left.nrows(), right.nrows()
    dimension = rows*columns
    basis = []
    for index in range(dimension):
        value = matrix(K, rows, columns, 0)
        value[index // columns, index % columns] = 1
        basis.append(value)
    transition = matrix(
        K, dimension, dimension,
        lambda i, j: flatten(left.inverse()*basis[j]*right)[i],
    )
    for seed_index in range(dimension):
        print("seed", seed_index, flush=True)
        seed = [K.zero()]*dimension
        seed[seed_index] = K.one()
        krylov_rows = [vector(K, seed)]
        for _ in range(dimension):
            krylov_rows.append(
                vector(K, [sigma(x) for x in krylov_rows[-1]])*transition
            )
        krylov = matrix(K, krylov_rows[:dimension])
        if krylov.det() == 0:
            print("noncyclic", flush=True)
            continue
        relation = krylov.transpose().solve_right(-krylov_rows[dimension])
        operator = OA(list(relation)+[1]).normalize()
        print("operator order", operator.order(),
              "degrees", [c.numerator().degree() for c in operator], flush=True)
        solutions = operator.rational_solutions()
        print("solution count", len(solutions), flush=True)
        for solution in solutions:
            scalar = solution[0]
            rhs = vector(K, [sigma(scalar, shift) for shift in range(dimension)])
            gauge_vector = krylov.solve_right(rhs)
            gauge = matrix(K, rows, columns, list(gauge_vector))
            assert gauge.apply_map(sigma) == left.inverse()*gauge*right
            print("rank", gauge.rank(), flush=True)
            for row in gauge.rows():
                print([entry.factor() for entry in row], flush=True)
            if gauge.rank() == columns:
                return gauge
    return None


target = challenge_matrix()
source = wilson_matrix()
print("quotient 3x2", flush=True)
quotient = recover(target, source)
print("quotient found", quotient is not None, flush=True)
print("embedding 2x3", flush=True)
embedding = recover(source, target)
print("embedding found", embedding is not None, flush=True)
