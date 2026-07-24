#!/usr/bin/env python3
"""Audit a multi-hit Siegel-lattice idea using the Apéry companion.

For hit pairs (j_i,p_i), P=prod p_i, and D=lcm(1,...,L)^3, form

    c_i=(P/p_i) D a_(j_i),   B_i=b_(j_i)/p_i.

An integer relation c.z=0 cancels the rational companion coordinate, while
the remaining Apéry coordinate is P*(B.z).  The image ideal of B on ker(c)
is computed from the 2-by-K Smith minors.  For the two known q=1 triples it
is respectively 1 and 85, so the remaining integer ideal is P and 85P.

This exact saturation does not provide a short relation: for K=3 the
minimum-supnorm solution of c.z=0, B.z=image is enormous.
"""

from __future__ import annotations

from fractions import Fraction
from math import ceil, floor, gcd, lcm, prod

from q32_newton import apery_numbers


TRIPLES = (
    (321, (179, 193, 211), "reflected"),
    (11576, (8893, 9319, 11437), "direct"),
)


def extended_gcd(first: int, second: int) -> tuple[int, int, int]:
    first_sign = -1 if first < 0 else 1
    second_sign = -1 if second < 0 else 1
    old_remainder, remainder = abs(first), abs(second)
    old_first, first_coefficient = 1, 0
    old_second, second_coefficient = 0, 1
    while remainder:
        quotient = old_remainder // remainder
        old_remainder, remainder = (
            remainder,
            old_remainder - quotient * remainder,
        )
        old_first, first_coefficient = (
            first_coefficient,
            old_first - quotient * first_coefficient,
        )
        old_second, second_coefficient = (
            second_coefficient,
            old_second - quotient * second_coefficient,
        )
    return (
        old_remainder,
        old_first * first_sign,
        old_second * second_sign,
    )


def companion_numbers(limit: int) -> list[Fraction]:
    values = [Fraction(0), Fraction(6)]
    for index in range(1, limit):
        values.append(
            (
                (
                    34 * index**3
                    + 51 * index**2
                    + 27 * index
                    + 5
                )
                * values[index]
                - index**3 * values[index - 1]
            )
            / (index + 1) ** 3
        )
    return values


def minimum_three_variable_solution(
    companion_row: list[int],
    apery_row: list[int],
) -> tuple[int, list[int]]:
    """Solve c.z=0, B.z=image and minimize supnorm along the rank-one fiber."""

    first, second, third = companion_row
    first_pair_gcd = gcd(first, second)
    row_gcd = gcd(first_pair_gcd, third)
    _, bezout_first, bezout_second = extended_gcd(
        first // first_pair_gcd,
        second // first_pair_gcd,
    )
    first_basis = (
        second // first_pair_gcd,
        -first // first_pair_gcd,
        0,
    )
    second_basis = (
        -(third // row_gcd) * bezout_first,
        -(third // row_gcd) * bezout_second,
        first_pair_gcd // row_gcd,
    )
    first_image = sum(
        apery_row[index] * first_basis[index]
        for index in range(3)
    )
    second_image = sum(
        apery_row[index] * second_basis[index]
        for index in range(3)
    )
    image, first_weight, second_weight = extended_gcd(
        first_image, second_image
    )
    solution = [
        first_weight * first_basis[index]
        + second_weight * second_basis[index]
        for index in range(3)
    ]
    if sum(
        apery_row[index] * solution[index]
        for index in range(3)
    ) == -image:
        solution = [-value for value in solution]

    kernel = [
        companion_row[1] * apery_row[2]
        - companion_row[2] * apery_row[1],
        companion_row[2] * apery_row[0]
        - companion_row[0] * apery_row[2],
        companion_row[0] * apery_row[1]
        - companion_row[1] * apery_row[0],
    ]
    kernel_gcd = 0
    for value in kernel:
        kernel_gcd = gcd(kernel_gcd, value)
    kernel = [value // kernel_gcd for value in kernel]

    candidates = {0}
    for index in range(3):
        if kernel[index]:
            crossing = Fraction(
                -solution[index], kernel[index]
            )
            candidates.update(
                range(floor(crossing) - 2, ceil(crossing) + 3)
            )
    for first_index in range(3):
        for second_index in range(first_index + 1, 3):
            for sign in (1, -1):
                denominator = (
                    kernel[first_index]
                    - sign * kernel[second_index]
                )
                if denominator:
                    crossing = Fraction(
                        sign * solution[second_index]
                        - solution[first_index],
                        denominator,
                    )
                    candidates.update(
                        range(
                            floor(crossing) - 2,
                            ceil(crossing) + 3,
                        )
                    )
    _, best_shift = min(
        (
            max(
                abs(
                    solution[index]
                    + shift * kernel[index]
                )
                for index in range(3)
            ),
            shift,
        )
        for shift in candidates
    )
    solution = [
        solution[index] + best_shift * kernel[index]
        for index in range(3)
    ]
    assert sum(
        companion_row[index] * solution[index]
        for index in range(3)
    ) == 0
    assert sum(
        apery_row[index] * solution[index]
        for index in range(3)
    ) == image
    return image, solution


def main() -> None:
    for n, raw_primes, branch in TRIPLES:
        if branch == "direct":
            pairs = sorted((n - prime, prime) for prime in raw_primes)
        else:
            pairs = sorted(
                (2 * prime - 1 - n, prime)
                for prime in raw_primes
            )
        maximum_index = max(index for index, _ in pairs)
        apery = apery_numbers(maximum_index)
        companion = companion_numbers(maximum_index)
        common_denominator = 1
        for value in range(1, maximum_index + 1):
            common_denominator = lcm(
                common_denominator, value
            )
        common_denominator **= 3
        prime_product = prod(prime for _, prime in pairs)

        companion_row = [
            prime_product
            // prime
            * int(common_denominator * companion[index])
            for index, prime in pairs
        ]
        apery_row = [
            apery[index] // prime
            for index, prime in pairs
        ]
        image, solution = minimum_three_variable_solution(
            companion_row, apery_row
        )
        print(
            f"n={n} branch={branch} "
            f"indices={[index for index, _ in pairs]} "
            f"image_after_P={image} "
            f"minimum_solution_bits="
            f"{[abs(value).bit_length() for value in solution]}"
        )


if __name__ == "__main__":
    main()
