#!/usr/bin/env python3
"""Verify the normalized moving-gap Green-kernel identities exactly."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json

import sympy as sp


EXPECTED_SHA256 = "8b27017ed37134278861b3460f86814e796916a432138233bad49359fa5a65ca"


def p_value(value):
    return (2 * value + 1) * (17 * value**2 + 17 * value + 5)


def rational_sequence(a_value: int, length: int) -> list[Fraction]:
    previous = Fraction(0)
    current = Fraction(1)
    values = [current]
    for index in range(length - 1):
        following = Fraction(
            p_value(a_value + index) * current
            - (a_value + index) ** 3 * previous,
            (a_value + index + 1) ** 3,
        )
        previous, current = current, following
        values.append(current)
    return values


def apery_rows(length: int) -> tuple[list[Fraction], list[Fraction]]:
    apery = [Fraction(1), Fraction(5)]
    companion = [Fraction(0), Fraction(1)]
    for index in range(1, length - 1):
        for row in (apery, companion):
            row.append(
                Fraction(
                    p_value(index) * row[index]
                    - index**3 * row[index - 1],
                    (index + 1) ** 3,
                )
            )
    return apery, companion


def main() -> None:
    a = sp.symbols("a")
    theta = sp.symbols("theta")
    z = sp.symbols("z")

    symbolic_u = [sp.Integer(1)]
    previous_u = sp.Integer(0)
    current_u = sp.Integer(1)
    for index in range(8):
        following_u = sp.expand(
            p_value(a + index) * current_u
            - (a + index) ** 6 * previous_u
        )
        previous_u, current_u = current_u, following_u
        symbolic_u.append(current_u)

    symbolic_v = []
    for index, value in enumerate(symbolic_u):
        denominator = sp.prod(a + step for step in range(1, index + 1)) ** 3
        symbolic_v.append(sp.cancel(value / denominator))

    for index in range(7):
        residual = sp.cancel(
            (a + index + 1) ** 3 * symbolic_v[index + 1]
            - p_value(a + index) * symbolic_v[index]
            + (a + index) ** 3
            * (sp.Integer(0) if index == 0 else symbolic_v[index - 1])
        )
        assert residual == 0

    # Coefficients of D_a F: a^3 at z^0 and zero thereafter.
    for degree in range(8):
        coefficient = ((a + degree) ** 3) * symbolic_v[degree]
        if degree >= 1:
            coefficient -= p_value(a + degree - 1) * symbolic_v[degree - 1]
        if degree >= 2:
            coefficient += (a + degree - 1) ** 3 * symbolic_v[degree - 2]
        expected = a**3 if degree == 0 else 0
        assert sp.cancel(coefficient - expected) == 0

    # The discarded eigenfunction equation already fails in degree one.
    wrong_degree_one = sp.cancel(
        ((a + 1) ** 3 - a**3) * symbolic_v[1] - p_value(a)
    )
    assert sp.cancel(
        wrong_degree_one + a**3 * p_value(a) / (a + 1) ** 3
    ) == 0
    assert wrong_degree_one != 0

    operator = (
        (theta + a) ** 3
        - z * p_value(theta + a)
        + z**2 * (theta + a + 1) ** 3
    )
    assert sp.expand(operator).coeff(theta, 3) == 1 - 34 * z + z**2

    apery, companion = apery_rows(24)
    rows = []
    for a_value in range(1, 11):
        values = rational_sequence(a_value, 11)
        for index, value in enumerate(values):
            kernel = a_value**3 * (
                apery[a_value - 1] * companion[a_value + index]
                - companion[a_value - 1] * apery[a_value + index]
            )
            assert value == kernel
            rows.append(
                [a_value, index, value.numerator, value.denominator]
            )

    assert rational_sequence(0, 12) == apery[:12]
    assert rational_sequence(1, 12) == companion[1:13]

    encoded = json.dumps(rows, separators=(",", ":")).encode("ascii")
    digest = hashlib.sha256(encoded).hexdigest()
    if EXPECTED_SHA256:
        assert digest == EXPECTED_SHA256
    print(
        "GAP_KERNEL_GREEN_VERIFY"
        " symbolic_depth=8 integer_a=1..10 coefficient_index=0..10"
    )
    print(f"sha256={digest}")


if __name__ == "__main__":
    main()
