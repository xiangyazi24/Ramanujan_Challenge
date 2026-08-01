#!/usr/bin/env sage-python
"""Search a rational quotient after folding two challenge steps."""

import sys

from sage.all import QQ, PolynomialRing, matrix, vector

sys.path.insert(0, "/Users/huangx/Library/SageMath-10-9/lib/python3.14/site-packages")
from ore_algebra import OreAlgebra


R = PolynomialRing(QQ, "r")
r = R.gen()
K = R.fraction_field()
OA = OreAlgebra(K, names=("Sr",))


def challenge_at(x):
    entries = [
        (-2*x-5)*(x+3)**2*(136*x**4+1424*x**3+5548*x**2+9551*x+6141),
        384*x**6+6384*x**5+44168*x**4+162698*x**3+336377*x**2+369933*x+169011,
        -(480*x**4+4980*x**3+19210*x**2+32690*x+20730),
        (x+2)**2*(x+3)**2*(4*x+10)*(48*x**3+386*x**2+1017*x+879),
        (x+2)**2*(-272*x**5-3848*x**4-21732*x**3-61184*x**2-85761*x-47808),
        (x+2)**2*(320*x**3+2540*x**2+6610*x+5640),
        (-4*x-10)*(x+2)**2*(x+3)**2*(32*x**4+302*x**3+1037*x**2+1530*x+813),
        (x+2)**2*(192*x**6+2984*x**5+19116*x**4+64452*x**3+120256*x**2+117279*x+46476),
        (x+2)**2*(-16*x**5-408*x**4-2912*x**3-8884*x**2-12254*x-6240),
    ]
    delta = -2*(x+2)**2*(x+3)**2*(2*x+5)*(2*x+7)**2
    return matrix(K, 3, 3, entries) / delta


def wilson_at(x):
    d0 = 8*(x+3)**2*(4*x+11)**2
    d1 = d0*(4*x+13)**2
    return matrix(K, [
        [
            (4*x+9)**2*(40*x**2+228*x+325)/d0,
            (1536*x**4+16512*x**3+66496*x**2+118896*x+79641)/d0,
        ],
        [
            (4*x+9)**2*(1536*x**4+18432*x**3+82880*x**2+165504*x+123841)/d1,
            (59392*x**6+1012736*x**5+7184384*x**4+27140352*x**3+
             57583336*x**2+65059404*x+30580677)/d1,
        ],
    ])


def sigma(value, amount=1):
    return value(r=r+amount)


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
        print("operator", operator.order(),
              [c.numerator().degree() for c in operator], flush=True)
        solutions = operator.rational_solutions()
        print("solutions", len(solutions), flush=True)
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


for phase in (0, 1):
    left = challenge_at(2*r+phase)*challenge_at(2*r+phase+1)
    # Wilson interval index is r for even and r+1 for odd.
    right = wilson_at(r+phase).transpose()
    print("PHASE", phase, "det ratio", (left.det()/right.det()).factor(), flush=True)
    answer = recover(left, right)
    print("FOUND", answer is not None, flush=True)
