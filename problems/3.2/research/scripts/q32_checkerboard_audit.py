#!/usr/bin/env python3
"""Exact and asymptotic audit of the Q803 checkerboard formula.

For a+b=H, let Q be the rational kernel of the signed Newton--Padé
conditions, normalized by q_b=1, and let Phat(n) be the corresponding
numerator value.  Q803 proposes the exact identity

    |Phat(n)| = (n)_(H+1)/b!
        * |Z^-(b+1, w_H/(n-s)) / Z^-(b, w_H)|,

where w_H(s)=A_s/(s!(H-s)!).  This script checks that identity against
the original Padé kernel, without using the partition formula to build
the kernel.

It also compares exact fixed-b values with the asymptotic claimed in
Q803.  The numerical comparison is evidence only; it does not supply
the uniform cofactor estimates omitted from that answer.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations
from math import comb, exp, factorial, log, prod, sqrt

from q32_pade_total_positivity import (
    apery_values,
    bareiss_determinant,
    entry,
    newton_coefficients,
)


def vandermonde(nodes: tuple[int, ...]) -> int:
    return prod(
        nodes[j] - nodes[i]
        for i in range(len(nodes))
        for j in range(i + 1, len(nodes))
    )


def zminus(
    height: int,
    particles: int,
    weight,
) -> Fraction:
    total = Fraction(0)
    for nodes in combinations(range(height + 1), particles):
        term = Fraction((-1) ** sum(nodes) * vandermonde(nodes) ** 2)
        for node in nodes:
            term *= weight(node)
        total += term
    return total


def rational_kernel(
    height: int,
    denominator_degree: int,
    coefficients: list[int],
) -> list[Fraction]:
    numerator_degree = height - denominator_degree
    rows = list(range(numerator_degree + 1, height + 1))
    columns = list(range(denominator_degree + 1))
    cofactors = []
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
        cofactors.append((-1) ** deleted_column * determinant)
    assert cofactors[-1]
    return [Fraction(value, cofactors[-1]) for value in cofactors]


def direct_projective_value(
    height: int,
    denominator_degree: int,
    evaluation_point: int,
    coefficients: list[int],
) -> Fraction:
    numerator_degree = height - denominator_degree
    kernel = rational_kernel(height, denominator_degree, coefficients)
    numerator = [
        sum(
            kernel[ell] * entry(index, ell, coefficients)
            for ell in range(min(denominator_degree, index) + 1)
        )
        for index in range(numerator_degree + 1)
    ]
    return sum(
        numerator[index] * comb(evaluation_point, index)
        for index in range(numerator_degree + 1)
    )


def partition_projective_value(
    height: int,
    denominator_degree: int,
    evaluation_point: int,
    apery: list[int],
) -> Fraction:
    weight = lambda node: Fraction(
        apery[node],
        factorial(node) * factorial(height - node),
    )
    shifted_weight = lambda node: weight(node) / (evaluation_point - node)
    falling = prod(evaluation_point - node for node in range(height + 1))
    return (
        Fraction(falling, factorial(denominator_degree))
        * zminus(height, denominator_degree + 1, shifted_weight)
        / zminus(height, denominator_degree, weight)
    )


def verify_exact_formula(max_height: int = 10) -> int:
    apery = apery_values(max_height + 2)
    coefficients = newton_coefficients(apery)
    checked = 0
    for height in range(1, max_height + 1):
        evaluation_point = 3 * height + 1
        for denominator_degree in range(height + 1):
            direct = direct_projective_value(
                height,
                denominator_degree,
                evaluation_point,
                coefficients,
            )
            partition = partition_projective_value(
                height,
                denominator_degree,
                evaluation_point,
                apery,
            )
            assert abs(direct) == abs(partition), (
                height,
                denominator_degree,
                direct,
                partition,
            )
            checked += 1
    return checked


def reciprocal_binomial_moment(
    height: int,
    power: int,
    apery: list[int],
) -> Fraction:
    return sum(
        Fraction(
            (-1) ** node * comb(height, node) * node**power,
            apery[node],
        )
        for node in range(height + 1)
    )


def verify_reciprocal_a_one(max_height: int = 10) -> int:
    """Check the explicit fixed-numerator-degree a=1 reduction."""

    apery = apery_values(max_height + 2)
    coefficients = newton_coefficients(apery)
    checked = 0
    for height in range(2, max_height + 1):
        evaluation_point = 3 * height + 1
        moments = [
            reciprocal_binomial_moment(height, power, apery)
            for power in range(3)
        ]
        explicit = (
            height
            * (evaluation_point * moments[0] - moments[1])
            / (moments[0] * moments[2] - moments[1] ** 2)
        )
        direct = direct_projective_value(
            height,
            height - 1,
            evaluation_point,
            coefficients,
        )
        assert abs(explicit) == abs(direct), (
            height,
            explicit,
            direct,
        )
        checked += 1
    return checked


def q803_fixed_strip_log_prediction(
    height: int,
    denominator_degree: int,
    coefficients: list[int],
) -> float:
    lam = 17 + 12 * sqrt(2)
    gamma = lam - 1
    numerator_degree = height - denominator_degree
    return (
        denominator_degree * log(lam / gamma)
        + log(coefficients[numerator_degree])
        + log(comb(3 * height + 1, numerator_degree))
        - (denominator_degree + 1) * log(1 - 1 / (2 * gamma))
    )


def fixed_strip_ratios(
    heights: tuple[int, ...] = (40, 80, 120, 180),
    max_denominator_degree: int = 8,
) -> list[tuple[int, int, float, float]]:
    apery = apery_values(max(heights) + 2)
    coefficients = newton_coefficients(apery)
    records = []
    for height in heights:
        for denominator_degree in range(max_denominator_degree + 1):
            exact = direct_projective_value(
                height,
                denominator_degree,
                3 * height + 1,
                coefficients,
            )
            log_prediction = q803_fixed_strip_log_prediction(
                height,
                denominator_degree,
                coefficients,
            )
            log_exact = (
                log(abs(exact.numerator))
                - log(exact.denominator)
            )
            ratio = exp(log_exact - log_prediction)
            rate = log_exact / (3 * height + 1)
            records.append((height, denominator_degree, ratio, rate))
    return records


def main() -> None:
    checked = verify_exact_formula()
    print(f"exact checkerboard numerator/denominator ratios checked={checked}")
    reciprocal_checked = verify_reciprocal_a_one()
    print(f"exact reciprocal a=1 reductions checked={reciprocal_checked}")
    for height, degree, ratio, rate in fixed_strip_ratios():
        print(
            f"H={height:3d} b={degree} "
            f"exact_over_Q803_prediction={ratio:.12f} "
            f"log_abs_Phat_over_n={rate:.12f}"
        )


if __name__ == "__main__":
    main()
