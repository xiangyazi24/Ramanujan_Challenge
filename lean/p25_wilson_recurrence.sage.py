#!/usr/bin/env sage-python
"""Temporary exact recurrence/gauge calculations for the r=2 Wilson-Pade family."""

import sys

from sage.all import QQ, PolynomialRing, binomial, matrix, vector

sys.path.insert(0, "/Users/huangx/Library/SageMath-10-9/lib/python3.14/site-packages")
from ore_algebra import OreAlgebra


def generalized_binomial(value, count):
    answer = QQ.one()
    for index in range(count):
        answer *= value - index
        answer /= index + 1
    return answer


def wilson_u(index):
    a = QQ(4 * index - 1) / 2
    return sum(binomial(index, j) * generalized_binomial(a, j)
               * generalized_binomial(a + j, j)
               for j in range(index + 1))


R = PolynomialRing(QQ, "n")
n = R.gen()
values = [wilson_u(index) for index in range(48)]
degree = 10
rows = []
for index in range(40):
    rows.append([values[index + shift] * QQ(index) ** power
                 for shift in range(3) for power in range(degree + 1)])
kernel = matrix(QQ, rows).right_kernel_matrix()
print("kernel dimension", kernel.nrows(), flush=True)
coefficients = list(kernel[0])
scale = coefficients[-1] ** -1
coefficients = [value * scale for value in coefficients]
polynomials = [sum(coefficients[shift * (degree + 1) + power] * n ** power
                   for power in range(degree + 1)) for shift in range(3)]
for index, polynomial in enumerate(polynomials):
    print("a" + str(index), polynomial.factor(), flush=True)
for index in range(40, 46):
    assert sum(polynomials[shift](index) * values[index + shift]
               for shift in range(3)) == 0
print("verified", flush=True)


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


source = matrix(K, [[0, -polynomials[0] / polynomials[2]],
                    [1, -polynomials[1] / polynomials[2]]])
target = challenge_matrix()


def sigma(value, amount=1):
    return value(n=n + amount)


def flatten(value):
    return vector(K, [value[row, column]
                      for row in range(value.nrows())
                      for column in range(value.ncols())])


def recover_rectangular(left, right):
    row_count = left.nrows()
    column_count = right.nrows()
    dimension = row_count * column_count
    basis = []
    for index in range(dimension):
        value = matrix(K, row_count, column_count, 0)
        value[index // column_count, index % column_count] = 1
        basis.append(value)
    left_inverse = left.inverse()
    transition = matrix(
        K, dimension, dimension,
        lambda row, column: flatten(left_inverse * basis[column] * right)[row],
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
        print("operator", operator, flush=True)
        solutions = operator.rational_solutions()
        print("solutions", solutions, flush=True)
        for solution in solutions:
            scalar = solution[0]
            rhs = vector(K, [sigma(scalar, shift) for shift in range(dimension)])
            gauge_vector = krylov.solve_right(rhs)
            gauge = matrix(K, row_count, column_count, list(gauge_vector))
            assert gauge.apply_map(sigma) == left_inverse * gauge * right
            print("rank", gauge.rank(), flush=True)
            for row in gauge.rows():
                print([entry.factor() for entry in row], flush=True)
            if gauge.rank() == row_count:
                return gauge
    return None


print("source companion", source, flush=True)
recovered = recover_rectangular(source, target)
print("recovered", recovered is not None, flush=True)
