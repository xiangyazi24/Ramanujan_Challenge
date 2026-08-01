#!/usr/bin/env python3
"""Machine checks for the honest horizontal objects in the rank-two attack.

This script verifies four logically distinct facts.

1. Additive orthogonality gives the exact two-prime zero correlation.
2. Averaging that correlation over a full CRT period factors exactly.
3. A full *character-period* correlation of linear Mellin transforms has an
   exact gcd(p-1,q-1)-term formula (a genuine linearized power saving).
4. Exponentiating a Jacobi/character sum does not move the additive character
   inside the sum.  Thus the zero detector is not a Deligne sum merely because
   the coefficient being tested is one.

The last distinction is the obstruction that survives the explicit rank-two
formula.  All computations use only Python's standard library.
"""

from __future__ import annotations

import cmath
from math import gcd, lcm, pi

from codex_hgk_coefficients import (
    apery_coefficients,
    branch_for_prime,
    branch_values_from_pullback,
)


def zero_set(prime: int) -> set[int]:
    return {
        index for index, value in enumerate(apery_coefficients(prime)) if value == 0
    }


def pair_count(prime: int, other: int, start: int, length: int) -> int:
    first_zeros = zero_set(prime)
    second_zeros = zero_set(other)
    return sum(
        index % prime in first_zeros and index % other in second_zeros
        for index in range(start, start + length)
    )


def additive_pair_count(prime: int, other: int, start: int, length: int) -> complex:
    first_values = apery_coefficients(prime)
    second_values = apery_coefficients(other)
    first_root = cmath.exp(2j * pi / prime)
    second_root = cmath.exp(2j * pi / other)
    result = 0j
    for first_mode in range(prime):
        for second_mode in range(other):
            result += sum(
                first_root ** (first_mode * first_values[index % prime])
                * second_root ** (second_mode * second_values[index % other])
                for index in range(start, start + length)
            )
    return result / (prime * other)


def verify_crt_averages() -> None:
    prime, other = 17, 19
    first_zeros = zero_set(prime)
    second_zeros = zero_set(other)
    assert first_zeros == {3, 13}
    assert second_zeros == {8, 10}

    period = prime * other
    complete = pair_count(prime, other, 0, period)
    assert complete == len(first_zeros) * len(second_zeros)

    length = 31
    sliding_total = sum(pair_count(prime, other, start, length) for start in range(period))
    assert sliding_total == length * len(first_zeros) * len(second_zeros)

    start = 23
    direct = pair_count(prime, other, start, length)
    orthogonal = additive_pair_count(prime, other, start, length)
    assert abs(orthogonal.imag) < 1e-8
    assert abs(orthogonal.real - direct) < 1e-8

    # Four-prime version: over the full CRT period every compatible zero
    # quadruple occurs exactly once.
    primes = (5, 11, 17, 19)
    four_period = 1
    for value in primes:
        four_period *= value
    expected = 1
    zero_sets = []
    for value in primes:
        zeros = zero_set(value)
        zero_sets.append(zeros)
        expected *= len(zeros)
    actual = sum(
        all(index % value in zeros for value, zeros in zip(primes, zero_sets))
        for index in range(four_period)
    )
    assert actual == expected

    print(
        "CRT averages: pair and four-prime complete periods factor exactly; "
        "sliding-interval mean and additive zero detector VERIFIED"
    )


def primitive_root(prime: int) -> int:
    order = prime - 1
    factors = []
    remaining = order
    divisor = 2
    while divisor * divisor <= remaining:
        if remaining % divisor == 0:
            factors.append(divisor)
            while remaining % divisor == 0:
                remaining //= divisor
        divisor += 1
    if remaining > 1:
        factors.append(remaining)
    return next(
        candidate
        for candidate in range(2, prime)
        if all(pow(candidate, order // factor, prime) != 1 for factor in factors)
    )


def multiplicative_ordering(values: list[int], prime: int) -> list[complex]:
    generator = primitive_root(prime)
    return [complex(values[pow(generator, exponent, prime)]) for exponent in range(prime - 1)]


def mellin_transform(values: list[complex], index: int) -> complex:
    order = len(values)
    root = cmath.exp(-2j * pi / order)
    return -sum(value * root ** (index * exponent) for exponent, value in enumerate(values))


def verify_linear_mellin_correlation() -> None:
    prime, other = 13, 29
    first_values = branch_values_from_pullback(prime, branch_for_prime(prime))
    second_values = branch_values_from_pullback(other, branch_for_prime(other))
    first = multiplicative_ordering(first_values, prime)
    second = multiplicative_ordering(second_values, other)

    first_order = prime - 1
    second_order = other - 1
    common = gcd(first_order, second_order)
    period = lcm(first_order, second_order)

    left = sum(
        mellin_transform(first, index)
        * mellin_transform(second, index).conjugate()
        for index in range(period)
    )
    right = period * sum(
        first[first_order * residue // common]
        * second[second_order * residue // common].conjugate()
        for residue in range(common)
    )
    assert abs(left - right) < 1e-7 * max(1.0, abs(right))

    print(
        "linear Mellin full-period identity: "
        f"p={prime}, q={other}, gcd(p-1,q-1)={common}, period={period}, VERIFIED"
    )


def verify_character_of_sum_obstruction() -> None:
    prime = 5
    root = cmath.exp(2j * pi / prime)
    values = (1, 2)
    character_after_sum = root ** (sum(values) % prime)
    sum_after_character = sum(root**value for value in values)
    assert abs(character_after_sum - sum_after_character) > 1e-3

    # This is the exact invalid interchange that would be needed to turn
    # e_p(a * (a Jacobi sum)) into a standard complete sum in its Jacobi
    # variables.
    print("nonlinear zero-detector interchange: explicit F_5 counterexample VERIFIED")


def verify_character_period_zero_correlation() -> None:
    prime, other = 17, 19
    first_period = prime - 1
    second_period = other - 1
    common = gcd(first_period, second_period)
    period = lcm(first_period, second_period)
    first_zeros = zero_set(prime) & set(range(first_period))
    second_zeros = zero_set(other) & set(range(second_period))

    direct = sum(
        index % first_period in first_zeros
        and index % second_period in second_zeros
        for index in range(period)
    )
    stratified = sum(
        sum(index % common == residue for index in first_zeros)
        * sum(index % common == residue for index in second_zeros)
        for residue in range(common)
    )
    assert direct == stratified
    print("zero events over lcm(p-1,q-1): gcd-stratified identity VERIFIED")


def main() -> None:
    verify_crt_averages()
    verify_linear_mellin_correlation()
    verify_character_of_sum_obstruction()
    verify_character_period_zero_correlation()


if __name__ == "__main__":
    main()
