#!/usr/bin/env sage
"""Audit the exact first-order factor of the doubled-period operator.

The certified order-three operator P for J_n has a rational first-order
right factor.  This script proves the factorization, identifies its
rational kernel solution, checks that P has no common scalar solution
with the Apéry recurrence, and gives a finite-field determinant
certificate excluding the simplest degree-at-most-30 polynomial gauge
for the distinguished sequence J.

The bounded-degree gauge exclusion is deliberately not promoted to a
general no-gauge theorem.
"""

import ast
from pathlib import Path

from ore_algebra import *


R = QQ["n"]
n = R.gen()
K = R.fraction_field()
OA = OreAlgebra(R, names=("Sn",))
Sn = OA.gen()


def stored_candidate():
    candidate_path = Path(__file__).with_name(
        "q32_doubled_period_recurrence_guess.py"
    )
    tree = ast.parse(candidate_path.read_text(encoding="utf-8"))
    for statement in tree.body:
        if not isinstance(statement, ast.Assign):
            continue
        if any(
            isinstance(target, ast.Name) and target.id == "RECURRENCE"
            for target in statement.targets
        ):
            return ast.literal_eval(statement.value)

    raise AssertionError("RECURRENCE was not found")


RECURRENCE = stored_candidate()
P = sum(
    sum(
        ZZ(coefficient) * n^degree
        for degree, coefficient in enumerate(row)
    )
    * Sn^shift
    for shift, row in enumerate(RECURRENCE)
)
assert P.order() == 3

A4 = 4 * n^5 + 12 * n^4 + 17 * n^3 + 78 * n^2 + 63 * n - 54
B4 = 4 * n^5 + 20 * n^4 + 57 * n^3 + 118 * n^2 + 107 * n + 30
right_factor = A4 / 4 * Sn - B4 / 4
left_factor, remainder = P.quo_rem(right_factor)
assert remainder == 0
assert left_factor.order() == 2

algorithmic_factors = P.right_factors(
    order=1,
    early_termination=False,
    infolevel=0,
)
assert len(algorithmic_factors) == 1
assert len(algorithmic_factors[0]) == 1
assert algorithmic_factors[0][0] == right_factor

assert A4 == 4 * (
    (n - 1 / 2)
    * (n + 3 / 2)
    * (n + 3)
    * (n^2 - n + 6)
)
assert B4 == 4 * (
    (n + 1 / 2)
    * (n + 1)
    * (n + 5 / 2)
    * (n^2 + n + 6)
)

y = K(
    (1 - 2 * n)
    * (2 * n + 3)
    * (n^2 - n + 6)
    / (9 * (n + 1) * (n + 2))
)
assert K(A4 / 4 * y(n=n + 1) - B4 / 4 * y) == 0
assert y(n=0) == 1

# The shifted Apéry recurrence annihilating b_n.
apery_operator = (
    (n + 1)^3
    - (2 * n + 3)
    * (17 * (n + 1)^2 + 17 * (n + 1) + 5)
    * Sn
    + (n + 2)^3 * Sn^2
)
assert apery_operator.order() == 2
assert P.gcrd(apery_operator).order() == 0
assert left_factor.gcrd(apery_operator).order() == 0


def polynomial_value(coefficients, index, modulus):
    out = 0
    for coefficient in reversed(coefficients):
        out = (out * index + coefficient) % modulus

    return out


modulus = 1_000_000_007
field = GF(modulus)
maximum_index = 63

# Generate b_n modulo the certificate prime.
apery = [1, 5]
for index in range(1, maximum_index):
    numerator = (
        (
            34 * index^3
            + 51 * index^2
            + 27 * index
            + 5
        )
        * apery[index]
        - index^3 * apery[index - 1]
    ) % modulus
    denominator = (index + 1)^3 % modulus
    apery.append(numerator * inverse_mod(denominator, modulus) % modulus)

# Generate J_n from the now-certified recurrence.
doubled_period = [45, 225, 3465]
for index in range(maximum_index - 2):
    coefficients = [
        polynomial_value(RECURRENCE[shift], index, modulus)
        for shift in range(4)
    ]
    assert coefficients[3] != 0
    next_value = -sum(
        coefficients[shift] * doubled_period[index + shift]
        for shift in range(3)
    )
    doubled_period.append(
        next_value * inverse_mod(coefficients[3], modulus) % modulus
    )

transformed = []
for index in range(maximum_index):
    a4 = ZZ(A4(n=index)) % modulus
    b4 = ZZ(B4(n=index)) % modulus
    transformed.append(
        (
            a4 * doubled_period[index + 1]
            - b4 * doubled_period[index]
        )
        % modulus
    )

# If
#
#   (A4(n)J_{n+1}-B4(n)J_n)
#       = U(n)b_n + V(n)b_{n+1}
#
# held over QQ with deg U,deg V <= 30, the following 63 integral
# columns would be dependent over QQ.  A nonzero determinant modulo one
# prime proves their integral determinant is nonzero.
degree_bound = 30
rows = []
for index in range(2 * degree_bound + 3):
    rows.append(
        [
            field(apery[index] * index^degree)
            for degree in range(degree_bound + 1)
        ]
        + [
            field(apery[index + 1] * index^degree)
            for degree in range(degree_bound + 1)
        ]
        + [field(transformed[index])]
    )

gauge_matrix = matrix(field, rows)
assert gauge_matrix.nrows() == gauge_matrix.ncols() == 63
determinant_residue = gauge_matrix.det()
assert determinant_residue != 0

print("Q32_DOUBLED_PERIOD_FACTOR_AUDIT=PASS")
print("OPERATOR_ORDER", P.order())
print("RIGHT_FACTOR_ORDER", right_factor.order())
print("LEFT_FACTOR_ORDER", left_factor.order())
print("APERY_GCRD_ORDER", P.gcrd(apery_operator).order())
print(
    "LEFT_FACTOR_APERY_GCRD_ORDER",
    left_factor.gcrd(apery_operator).order(),
)
print("RATIONAL_KERNEL", y)
print("POLYNOMIAL_GAUGE_DEGREE_BOUND", degree_bound)
print("GAUGE_DETERMINANT_MODULUS", modulus)
print("GAUGE_DETERMINANT_RESIDUE", determinant_residue)
