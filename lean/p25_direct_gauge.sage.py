#!/usr/bin/env sage-python
"""Temporary exact recovery of a direct-Catalan 3F2 gauge for P2.5."""

import sys

from sage.all import QQ, PolynomialRing, matrix, vector

sys.path.insert(0, "/Users/huangx/Library/SageMath-10-9/lib/python3.14/site-packages")
from ore_algebra import OreAlgebra

R = PolynomialRing(QQ, "n")
n = R.gen()
K = R.fraction_field()
OA = OreAlgebra(K, names=("Sn",))


def source_matrix():
    # Exact RamanujanTools trajectory matrix for the determinant-compatible
    # spectral hit (-2,+2,0;0,0), based at
    # 3F2(1/2,1/2,1;3/2,3/2;-1) = Catalan's constant.
    position = [K(1)/2 - 2*n, K(1)/2 + 2*n, 1, K(3)/2, K(3)/2]

    def theta_matrix(pos):
        x0, x1, x2, y0, y1 = pos
        # Companion matrix of the monic differential polynomial
        # theta(theta+y0-1)(theta+y1-1) + prod(theta+xi).
        c2 = ((y0 - 1) + (y1 - 1) + x0 + x1 + x2) / 2
        c1 = ((y0 - 1)*(y1 - 1) + x0*x1 + x0*x2 + x1*x2) / 2
        c0 = x0*x1*x2 / 2
        return matrix(K, [[0, 0, -c0], [1, 0, -c1], [0, 1, -c2]])

    def x_positive(pos, axis):
        return matrix.identity(K, 3) + theta_matrix(pos) / pos[axis]

    def x_negative(pos, axis):
        shifted = list(pos)
        shifted[axis] -= 1
        return x_positive(shifted, axis).inverse()

    def diagonal_step(pos):
        pos = list(pos)
        result = matrix.identity(K, 3)
        # RamanujanTools' first valid path is reversed(x0,x1): x1 then x0.
        result *= x_positive(pos, 1)
        pos[1] += 1
        result *= x_negative(pos, 0)
        pos[0] -= 1
        return result

    first = diagonal_step(position)
    second_position = list(position)
    second_position[0] -= 1
    second_position[1] += 1
    return first * diagonal_step(second_position)


def target_matrix():
    g = 2*(n+2)**2*(n+3)**2*(2*n+5)*(2*n+7)**2
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
    ]) / g


def sigma(value, amount=1):
    return value(n=n + amount)


def flatten(value):
    return vector(K, [value[row, column]
                      for row in range(value.nrows())
                      for column in range(value.ncols())])


left = source_matrix()
right = target_matrix()
dimension = 9
basis = []
for index in range(dimension):
    value = matrix(K, 3, 3, 0)
    value[index // 3, index % 3] = 1
    basis.append(value)
transition = matrix(
    K, dimension, dimension,
    lambda row, column: flatten(left.inverse() * basis[column] * right)[row],
)

for seed_index in range(dimension):
    print("seed", seed_index, flush=True)
    seed = [K.zero()] * dimension
    seed[seed_index] = K.one()
    rows = [vector(K, seed)]
    for _ in range(dimension):
        rows.append(vector(K, [sigma(entry) for entry in rows[-1]]) * transition)
    krylov = matrix(K, rows[:dimension])
    if krylov.det() == 0:
        print("noncyclic", flush=True)
        continue
    relation = krylov.transpose().solve_right(-rows[dimension])
    operator = OA(list(relation) + [1]).normalize()
    print("operator order", operator.order(), "degree", operator.degree(), flush=True)
    solutions = operator.rational_solutions()
    print("solution count", len(solutions), flush=True)
    for solution in solutions:
        scalar = solution[0]
        rhs = vector(K, [sigma(scalar, shift) for shift in range(dimension)])
        gauge_vector = krylov.solve_right(rhs)
        gauge = matrix(K, 3, 3, list(gauge_vector))
        assert gauge.apply_map(sigma) == left.inverse() * gauge * right
        print("rank", gauge.rank(), flush=True)
        for row in gauge.rows():
            print([entry.factor() for entry in row], flush=True)
        if gauge.rank() == 3:
            raise SystemExit(0)

print("no full-rank gauge", flush=True)
