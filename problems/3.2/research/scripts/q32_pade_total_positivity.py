#!/usr/bin/env python3
"""Audit (and refute) total positivity in the signed Newton--Padé matrix.

Let A_n be the Apéry numbers, c_k=Delta^k A_0, and

    M[k,l] = binom(k,l) Delta^(k-l) A_l
           = binom(k,l) sum_(t=0)^l binom(l,t)c_(k-l+t).

The signed Padé equations use consecutive bottom rows of M.  There is an
exact factorization after row reversal:

    M[k,k-r] = (B diag(c) B)[k,r],

where B[k,r]=binom(k,r) is the Pascal matrix.  Hence the row-reversed
triangle is totally nonnegative.  This does *not* make M totally
nonnegative: the consecutive order-six minor on rows 1,...,6 and columns
0,...,5 is negative.

This script verifies the factorization, audits small minors, records the
first structured counterexample, and retains the limited same-sign audit
for denominator degrees at most six.
"""

from __future__ import annotations

from itertools import combinations
from math import comb, gcd
from random import Random


APERY_LIMIT = 90


def apery_values(limit: int) -> list[int]:
    values = [1, 5]
    for n in range(1, limit - 1):
        polynomial = 34 * n**3 + 51 * n**2 + 27 * n + 5
        numerator = polynomial * values[n] - n**3 * values[n - 1]
        values.append(numerator // (n + 1) ** 3)
    return values


def newton_coefficients(values: list[int]) -> list[int]:
    row = values[:]
    result = []
    while row:
        result.append(row[0])
        row = [row[index + 1] - row[index] for index in range(len(row) - 1)]
    return result


def entry(k: int, ell: int, coefficients: list[int]) -> int:
    if ell > k:
        return 0
    shifted_difference = sum(
        comb(ell, t) * coefficients[k - ell + t]
        for t in range(ell + 1)
    )
    return comb(k, ell) * shifted_difference


def bareiss_determinant(matrix: list[list[int]]) -> int:
    """Exact fraction-free determinant."""

    size = len(matrix)
    if size == 0:
        return 1
    work = [row[:] for row in matrix]
    sign = 1
    previous = 1

    for pivot_index in range(size - 1):
        if work[pivot_index][pivot_index] == 0:
            swap = next(
                (
                    row
                    for row in range(pivot_index + 1, size)
                    if work[row][pivot_index]
                ),
                None,
            )
            if swap is None:
                return 0
            work[pivot_index], work[swap] = (
                work[swap],
                work[pivot_index],
            )
            sign = -sign

        pivot = work[pivot_index][pivot_index]
        for row in range(pivot_index + 1, size):
            for column in range(pivot_index + 1, size):
                numerator = (
                    work[row][column] * pivot
                    - work[row][pivot_index]
                    * work[pivot_index][column]
                )
                work[row][column] = numerator // previous
        previous = pivot

    return sign * work[-1][-1]


def verify_pascal_factorization(coefficients: list[int], limit: int) -> None:
    for k in range(limit):
        for r in range(k + 1):
            product_entry = sum(
                comb(k, middle)
                * coefficients[middle]
                * comb(middle, r)
                for middle in range(r, k + 1)
            )
            assert entry(k, k - r, coefficients) == product_entry


def exhaustive_minor_audit(
    coefficients: list[int],
    row_limit: int,
    column_limit: int,
    max_order: int,
) -> int:
    checked = 0
    for order in range(1, max_order + 1):
        for rows in combinations(range(row_limit), order):
            for columns in combinations(range(column_limit), order):
                determinant = bareiss_determinant(
                    [
                        [entry(row, column, coefficients) for column in columns]
                        for row in rows
                    ]
                )
                assert determinant >= 0, (rows, columns, determinant)
                checked += 1
    return checked


def random_minor_audit(
    coefficients: list[int],
    trials_per_order: int = 250,
) -> int:
    random = Random(762)
    checked = 0
    for order in range(2, 10):
        for _ in range(trials_per_order):
            rows = sorted(random.sample(range(55), order))
            columns = sorted(random.sample(range(18), order))
            determinant = bareiss_determinant(
                [
                    [entry(row, column, coefficients) for column in columns]
                    for row in rows
                ]
            )
            assert determinant >= 0, (rows, columns, determinant)
            checked += 1
    return checked


def primitive_pade_kernel(
    height: int,
    denominator_degree: int,
    coefficients: list[int],
) -> list[int]:
    numerator_degree = height - denominator_degree
    rows = list(range(numerator_degree + 1, height + 1))
    columns = list(range(denominator_degree + 1))

    kernel = []
    for deleted_column in columns:
        minor_columns = [
            column for column in columns if column != deleted_column
        ]
        determinant = bareiss_determinant(
            [
                [entry(row, column, coefficients) for column in minor_columns]
                for row in rows
            ]
        )
        kernel.append((-1) ** deleted_column * determinant)

    common = abs(gcd(*kernel))
    return [value // common for value in kernel]


def verify_pade_signs(coefficients: list[int]) -> int:
    """Check only the originally tested range; this is not a theorem."""
    checked = 0
    for denominator_degree in range(1, 7):
        for height in range(max(denominator_degree + 2, 6), 41):
            numerator_degree = height - denominator_degree
            kernel = primitive_pade_kernel(
                height, denominator_degree, coefficients
            )
            numerator_coefficients = [
                sum(
                    kernel[ell] * entry(k, ell, coefficients)
                    for ell in range(min(denominator_degree, k) + 1)
                )
                for k in range(numerator_degree + 1)
            ]
            signs = {value > 0 for value in numerator_coefficients}
            assert 0 not in numerator_coefficients
            assert len(signs) == 1, (
                height,
                denominator_degree,
                numerator_coefficients,
            )
            checked += 1
    return checked


def verify_total_positivity_counterexample(
    coefficients: list[int],
) -> tuple[int, list[int], list[int]]:
    """Return the first consecutive flag counterexample found at H=6, b=5."""

    rows = list(range(1, 7))
    columns = list(range(6))
    determinant = bareiss_determinant(
        [
            [entry(row, column, coefficients) for column in columns]
            for row in rows
        ]
    )
    assert determinant == -2421161406987687811000

    kernel = primitive_pade_kernel(6, 5, coefficients)
    numerator_coefficients = [
        sum(
            kernel[ell] * entry(k, ell, coefficients)
            for ell in range(min(5, k) + 1)
        )
        for k in range(2)
    ]
    assert numerator_coefficients == [
        1267742334817618036530,
        -12105807034938439055,
    ]
    return determinant, kernel, numerator_coefficients


def main() -> None:
    values = apery_values(APERY_LIMIT)
    coefficients = newton_coefficients(values)

    verify_pascal_factorization(coefficients, 40)
    exhaustive = exhaustive_minor_audit(
        coefficients,
        row_limit=14,
        column_limit=8,
        max_order=3,
    )
    random = random_minor_audit(coefficients)
    pade = verify_pade_signs(coefficients)
    counterexample, _, numerator = verify_total_positivity_counterexample(
        coefficients
    )

    print(
        "Padé multiplication-matrix audit: "
        f"factorization_rows=40, exhaustive_minors={exhaustive}, "
        f"random_minors={random}, restricted_pade_pairs={pade}, "
        f"negative_order_6_minor={counterexample}, "
        f"numerator_coefficients={numerator}"
    )


if __name__ == "__main__":
    main()
