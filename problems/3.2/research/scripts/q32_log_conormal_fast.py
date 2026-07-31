#!/usr/bin/env python3
"""Exact first-conormal audit for the Apéry logarithmic reduction.

This is a lightweight Python/SymPy implementation of the two-jet
calculation in ``/tmp/p32_local_conormal_connection.sage``.  It avoids
polynomial Groebner bases: the completed local ring is represented
directly in the coordinates

    u = y-z,  v = x-2z-1,  w = z^2-1/2

through total (u,v,w)-degree two.  The script constructs logarithmic
fields of weights u,v,w, computes their induced first-conormal
connection, and tests the six-dimensional tangent gauge space on a
deterministic basis.  The raw determinant is not expected to be
gauge-invariant; the audit records that failure without dumping the
very large determinant polynomials.
"""

from collections import defaultdict
from fractions import Fraction
from hashlib import sha256
from itertools import permutations
from math import comb

import sympy as sp


def clean(value):
    return {key: coefficient for key, coefficient in value.items() if coefficient}


def add(left, right):
    out = defaultdict(lambda: Fraction(0))
    for key, coefficient in left.items():
        out[key] += coefficient
    for key, coefficient in right.items():
        out[key] += coefficient
    return clean(out)


def scale(scalar, value):
    scalar = Fraction(scalar)
    return clean({key: scalar * coefficient for key, coefficient in value.items()})


def monomial(u=0, v=0, w=0, z=0, coefficient=1):
    """Return a reduced jet monomial, using z^2=1/2+w."""

    out = defaultdict(lambda: Fraction(0))
    quotient, parity = divmod(z, 2)
    for extra_w in range(quotient + 1):
        if u + v + w + extra_w > 2:
            continue
        out[(u, v, w + extra_w, parity)] += (
            Fraction(coefficient)
            * comb(quotient, extra_w)
            * Fraction(1, 2) ** (quotient - extra_w)
        )
    return clean(out)


def multiply(left, right):
    out = {}
    for left_key, left_coefficient in left.items():
        for right_key, right_coefficient in right.items():
            out = add(
                out,
                monomial(
                    left_key[0] + right_key[0],
                    left_key[1] + right_key[1],
                    left_key[2] + right_key[2],
                    left_key[3] + right_key[3],
                    left_coefficient * right_coefficient,
                ),
            )
    return out


ONE = monomial()
Z = monomial(z=1)
U = monomial(u=1)
V = monomial(v=1)
W = monomial(w=1)
X = add(add(ONE, scale(2, Z)), V)
Y = add(Z, U)


def power(value, exponent):
    out = ONE
    for _ in range(exponent):
        out = multiply(out, value)
    return out


def derivative(value, coordinate):
    """Apply theta_x, theta_y, or theta_z in local coordinates."""

    actions = {
        0: ({}, {}, X, {}),
        1: ({}, Y, {}, {}),
        2: (Z, scale(-1, Z), scale(-2, Z), add(ONE, scale(2, W))),
    }[coordinate]
    out = {}
    for key, coefficient in value.items():
        degrees = (key[3], key[0], key[1], key[2])  # z,u,v,w
        generators = (Z, U, V, W)
        for generator_index, degree in enumerate(degrees):
            if not degree or not actions[generator_index]:
                continue
            term = monomial(coefficient=coefficient * degree)
            for index, generator in enumerate(generators):
                exponent = degrees[index] - (index == generator_index)
                term = multiply(term, power(generator, exponent))
            out = add(out, multiply(term, actions[generator_index]))
    return out


def inverse_unit(value):
    base = {
        key: coefficient
        for key, coefficient in value.items()
        if sum(key[:3]) == 0
    }
    scalar = base.get((0, 0, 0, 0), Fraction(0))
    z_coefficient = base.get((0, 0, 0, 1), Fraction(0))
    norm = scalar**2 - z_coefficient**2 / 2
    base_inverse = add(
        monomial(coefficient=scalar / norm),
        monomial(z=1, coefficient=-z_coefficient / norm),
    )
    delta = add(value, scale(-1, base))
    epsilon = multiply(base_inverse, delta)
    return multiply(
        base_inverse,
        add(add(ONE, scale(-1, epsilon)), multiply(epsilon, epsilon)),
    )


F = multiply(
    add(ONE, X),
    multiply(
        add(ONE, Y),
        multiply(
            add(ONE, Z),
            add(multiply(add(ONE, Y), add(ONE, Z)), multiply(X, multiply(Y, Z))),
        ),
    ),
)
F_INV = inverse_unit(F)
CRITICAL = [
    multiply(add(derivative(F, coordinate), scale(-1, F)), F_INV)
    for coordinate in range(3)
]

BASIS_H = (
    ONE,
    Z,
    U,
    monomial(u=1, z=1),
    V,
    monomial(v=1, z=1),
    W,
    monomial(w=1, z=1),
)
POSITIVE_KEYS = tuple(
    (iu, iv, total - iu - iv, ez)
    for total in (1, 2)
    for iu in range(total + 1)
    for iv in range(total - iu + 1)
    for ez in (0, 1)
)


def sympy_fraction(value):
    if isinstance(value, Fraction):
        return sp.Rational(value.numerator, value.denominator)
    return sp.sympify(value)


columns = []
for coordinate in range(3):
    for basis_value in BASIS_H:
        contribution = multiply(basis_value, CRITICAL[coordinate])
        columns.append(
            [sympy_fraction(contribution.get(key, Fraction(0))) for key in POSITIVE_KEYS]
        )
WEIGHT_MATRIX = sp.Matrix(18, 24, lambda row, col: columns[col][row])
assert WEIGHT_MATRIX.rank() == 18


def solve_weight(target):
    rhs = sp.Matrix(
        [sympy_fraction(target.get(key, Fraction(0))) for key in POSITIVE_KEYS]
    )
    solution, parameters = WEIGHT_MATRIX.gauss_jordan_solve(rhs)
    solution = solution.subs({parameter: 0 for parameter in parameters})
    fields = []
    for coordinate in range(3):
        value = {}
        for index, basis_value in enumerate(BASIS_H):
            coefficient = solution[8 * coordinate + index]
            value = add(
                value,
                scale(
                    Fraction(int(coefficient.p), int(coefficient.q)),
                    basis_value,
                ),
            )
        fields.append(value)
    residual = {}
    for coordinate in range(3):
        residual = add(residual, multiply(fields[coordinate], CRITICAL[coordinate]))
    assert add(residual, scale(-1, target)) == {}
    return fields


FIELDS = [solve_weight(target) for target in (U, V, W)]
KERNEL = WEIGHT_MATRIX.nullspace()
assert len(KERNEL) == 6

a, b, c, M = sp.symbols("a b c M")
EXPONENTS = (a, b, c)


def first_jet_coefficients(value):
    remainder = (
        sympy_fraction(value.get((0, 0, 0, 0), Fraction(0))),
        sympy_fraction(value.get((0, 0, 0, 1), Fraction(0))),
    )
    coefficients = []
    for key0 in ((1, 0, 0), (0, 1, 0), (0, 0, 1)):
        coefficients.append(
            (
                sympy_fraction(value.get(key0 + (0,), Fraction(0))),
                sympy_fraction(value.get(key0 + (1,), Fraction(0))),
            )
        )
    return remainder, tuple(coefficients)


def pair_add(left, right):
    return (sp.expand(left[0] + right[0]), sp.expand(left[1] + right[1]))


def pair_mul(left, right):
    return (
        sp.expand(left[0] * right[0] + left[1] * right[1] / 2),
        sp.expand(left[0] * right[1] + left[1] * right[0]),
    )


def pair_scale(scalar, value):
    return (sp.expand(scalar * value[0]), sp.expand(scalar * value[1]))


def determinant_pair(matrix):
    out = (sp.Integer(0), sp.Integer(0))
    for permutation in permutations(range(3)):
        inversions = sum(
            permutation[i] > permutation[j]
            for i in range(3)
            for j in range(i + 1, 3)
        )
        term = (sp.Integer(1), sp.Integer(0))
        for row, column in enumerate(permutation):
            term = pair_mul(term, matrix[row][column])
        out = pair_add(out, pair_scale((-1) ** inversions, term))
    return out


def connection_for_fields(selected_fields):
    connection = [[None] * 3 for _ in range(3)]
    base_remainders = []
    for column, field in enumerate(selected_fields):
        beta = {}
        value = {}
        for coordinate in range(3):
            beta = add(beta, derivative(field[coordinate], coordinate))
            for key, coefficient in field[coordinate].items():
                value[key] = value.get(key, 0) + coefficient * EXPONENTS[coordinate]
        value = add(value, beta)
        remainder, jets = first_jet_coefficients(value)
        base_remainders.append(remainder)
        for row in range(3):
            connection[row][column] = jets[row]
    block = [[connection[row][column] for column in range(3)] for row in range(3)]
    for index in range(3):
        block[index][index] = pair_add(block[index][index], (M, 0))
    det = determinant_pair(block)
    norm = sp.factor(det[0] ** 2 - det[1] ** 2 / 2)
    return connection, base_remainders, det, norm


def field_from_vector(vector):
    fields = []
    for coordinate in range(3):
        value = {}
        for index, basis_value in enumerate(BASIS_H):
            coefficient = vector[8 * coordinate + index]
            value = add(
                value,
                scale(Fraction(int(coefficient.p), int(coefficient.q)), basis_value),
            )
        fields.append(value)
    return fields


CONNECTION, BASE_REMAINDERS, DET_PAIR, NORM = connection_for_fields(FIELDS)
print("WEIGHT_MATRIX_RANK", WEIGHT_MATRIX.rank())
print("WEIGHT_GAUGE_DIMENSION_Q", len(KERNEL))


def polynomial_summary(value):
    value = sp.Poly(sp.expand(value), M, a, b, c)
    payload = str(value.as_expr()).encode()
    return {
        "total_degree": value.total_degree(),
        "terms": len(value.terms()),
        "sha256": sha256(payload).hexdigest()[:16],
    }


print("BASE_REMAINDER_COUNT", len(BASE_REMAINDERS))
print("RAW_NORM_SUMMARY", polynomial_summary(NORM))

gauge_norms = set()
for field_index in range(3):
    for kernel_vector in KERNEL:
        modified = [[dict(component) for component in field] for field in FIELDS]
        gauge = field_from_vector(kernel_vector)
        modified[field_index] = [
            add(modified[field_index][coordinate], gauge[coordinate])
            for coordinate in range(3)
        ]
        gauge_norms.add(str(connection_for_fields(modified)[3]))
print("GAUGE_NORM_COUNT", len(gauge_norms))
print(
    "GAUGE_NORM_DIGESTS",
    sorted(sha256(norm.encode()).hexdigest()[:16] for norm in gauge_norms),
)
assert len(gauge_norms) > 1
print("RAW_CONORMAL_DETERMINANT_GAUGE_INVARIANT", False)
print("Q32_LOG_CONORMAL_FAST=PASS")
