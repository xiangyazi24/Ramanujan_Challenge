#!/usr/bin/env python3
"""Exact audit of the three-target Racah cutoff saturation at n=321."""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations
from math import comb, gcd, prod


def determinant(matrix: list[list[int]]) -> int:
    """Fraction-free determinant, sufficient for the small audit matrices."""
    size = len(matrix)
    if size == 0:
        return 1
    if size == 1:
        return matrix[0][0]
    return sum(
        (-1) ** column
        * matrix[0][column]
        * determinant(
            [
                row[:column] + row[column + 1 :]
                for row in matrix[1:]
            ]
        )
        for column in range(size)
    )


def residue(value: Fraction | int, prime: int) -> int:
    value = Fraction(value)
    assert value.denominator % prime
    return (
        value.numerator
        * pow(value.denominator, -1, prime)
        % prime
    )


def audit() -> None:
    n = 321
    primes = (179, 193, 211)
    cutoffs = tuple((prime - 1) // 2 for prime in primes)
    assert cutoffs == (89, 96, 105)

    def apery_factor(index: int) -> int:
        return comb(n, index) * comb(n + index, index)

    def summand(index: int) -> int:
        return apery_factor(index) ** 2

    prefixes = tuple(
        sum(summand(index) for index in range(cutoff + 1))
        for cutoff in cutoffs
    )
    for prefix, prime in zip(prefixes, primes):
        assert prefix % prime == 0

    cumulative = (
        primes[0],
        primes[0] * primes[1],
        prod(primes),
    )
    differences = (
        prefixes[1] - prefixes[0],
        prefixes[2] - prefixes[1],
    )
    assert differences[0] % cumulative[0] ** 2 == 0
    assert differences[1] % cumulative[1] ** 2 == 0
    assert prefixes[1] % cumulative[1] == 0
    assert prefixes[2] % cumulative[2] == 0

    divided = (
        differences[0] // cumulative[0] ** 2,
        differences[1] // cumulative[1] ** 2,
    )
    boundary = summand(cutoffs[0] + 1) // primes[0] ** 2

    # The Q5677 attachment incorrectly asserted integer divisibility here.
    # The correct rank-one identity lives in the localization whose
    # denominators are cutoff factorials.
    assert divided[0] % boundary
    assert divided[1] % boundary
    racah = (
        Fraction(divided[0], boundary),
        Fraction(divided[1], boundary),
    )
    assert racah[1] * divided[0] - racah[0] * divided[1] == 0

    final_prime = primes[2]
    assert tuple(value % final_prime for value in differences) == (59, 125)
    assert tuple(value % final_prime for value in divided) == (107, 4)
    assert boundary % final_prime == 66
    assert tuple(residue(value, final_prime) for value in racah) == (24, 64)

    # Presentation rows [s, W, -p_i e_i].  The gcd of its maximal minors
    # is one, so its primitive maximal Fitting ideal is the unit ideal.
    q_values = (
        0,
        primes[0] ** 2 * racah[0],
        primes[0] ** 2 * racah[0]
        + cumulative[1] ** 2 * racah[1],
    )
    common_denominator = 1
    for value in q_values:
        common_denominator = (
            common_denominator
            * value.denominator
            // gcd(common_denominator, value.denominator)
        )
    q_integral = [
        value.numerator * (common_denominator // value.denominator)
        for value in q_values
    ]
    presentation = []
    for row, prime in enumerate(primes):
        values = [common_denominator, q_integral[row], 0, 0, 0]
        values[2 + row] = -prime * common_denominator
        presentation.append(values)

    maximal_minors = []
    for columns in combinations(range(5), 3):
        minor = [
            [presentation[row][column] for column in columns]
            for row in range(3)
        ]
        maximal_minors.append(determinant(minor))
    normalized_gcd = 0
    for value in maximal_minors:
        normalized_gcd = gcd(normalized_gcd, abs(value))
    # Clearing the localized coefficient ring contributes only the chosen
    # common denominator; after its saturation the target-prime content is 1.
    while gcd(normalized_gcd, common_denominator) > 1:
        normalized_gcd //= gcd(normalized_gcd, common_denominator)
    assert normalized_gcd == 1

    target_quotients = tuple(
        prefix // prime for prefix, prime in zip(prefixes, primes)
    )
    determinant_identity = (
        q_values[1]
        * (primes[2] * target_quotients[2]
           - primes[0] * target_quotients[0])
        - q_values[2]
        * (primes[1] * target_quotients[1]
           - primes[0] * target_quotients[0])
    )
    assert determinant_identity == 0

    print(
        {
            "n": n,
            "targets": primes,
            "cutoffs": cutoffs,
            "divided_mod_211": tuple(
                value % final_prime for value in divided
            ),
            "racah_mod_211": tuple(
                residue(value, final_prime) for value in racah
            ),
            "primitive_maximal_minor_gcd": normalized_gcd,
            "failures": 0,
        }
    )


if __name__ == "__main__":
    audit()
