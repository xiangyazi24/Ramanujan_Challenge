#!/usr/bin/env python3
"""Audit the Bernoulli-free quartic and quintic target laws.

For a target p | b_r, n=p+r, put s=p-1-r and x=b_n/p.  The direct
and reflected exact block decompositions, together with the rank-one
endpoint defect, give two normalized residues

    d == x (1 + 7 p^3 e_p / 5)  (mod p^4),
    z == x (1 - 8 p^3 e_p / 5)  (mod p^4).

Consequently 15x == 8d+7z (mod p^4), with the Bernoulli scalar e_p
eliminated.  The stronger endpoint law modulo p^5 replaces p^3 e_p
by Delta_p=b_{p-1}-1 and gives the same cancellation modulo p^5.
This script independently constructs the required degree-five
recurrence jets and checks both statements.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from fractions import Fraction
from math import comb, gcd, isqrt


def primes_at_most(limit: int) -> list[int]:
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[:2] = b"\x00\x00"
    for prime in range(2, isqrt(limit) + 1):
        if sieve[prime]:
            sieve[prime * prime : limit + 1 : prime] = b"\x00" * (
                (limit - prime * prime) // prime + 1
            )
    return [candidate for candidate in range(2, limit + 1) if sieve[candidate]]


ORDER = 6


def recurrence_polynomial(index: int) -> int:
    return 34 * index**3 + 51 * index**2 + 27 * index + 5


def apery_numbers(limit: int) -> list[int]:
    values = [1, 5]
    for index in range(1, limit):
        numerator = (
            recurrence_polynomial(index) * values[index]
            - index**3 * values[index - 1]
        )
        denominator = (index + 1) ** 3
        assert numerator % denominator == 0
        values.append(numerator // denominator)
    return values


def add(left: list[Fraction], right: list[Fraction]) -> list[Fraction]:
    return [left[index] + right[index] for index in range(ORDER)]


def scale(polynomial: list[Fraction], scalar: int) -> list[Fraction]:
    return [scalar * coefficient for coefficient in polynomial]


def multiply(
    left: list[Fraction], right: list[Fraction]
) -> list[Fraction]:
    result = [Fraction(0) for _ in range(ORDER)]
    for left_index, left_coefficient in enumerate(left):
        for right_index, right_coefficient in enumerate(right):
            if left_index + right_index >= ORDER:
                break
            result[left_index + right_index] += (
                left_coefficient * right_coefficient
            )
    return result


def inverse_shift_cube(shift: int) -> list[Fraction]:
    return [
        Fraction((-1) ** degree * comb(degree + 2, 2), shift ** (degree + 3))
        for degree in range(ORDER)
    ]


def shifted_fundamental_solutions(
    limit: int,
) -> tuple[list[list[Fraction]], list[list[Fraction]]]:
    zero = [Fraction(0) for _ in range(ORDER)]
    one = zero.copy()
    one[0] = Fraction(1)

    b_at_x = [
        Fraction(5),
        Fraction(27),
        Fraction(51),
        Fraction(34),
        Fraction(0),
        Fraction(0),
    ]
    direct_one = multiply(b_at_x, inverse_shift_cube(1))

    direct = [one, direct_one]
    companion = [zero, one]

    for index in range(1, limit):
        shifted_b = [
            Fraction(recurrence_polynomial(index)),
            Fraction(102 * index**2 + 102 * index + 27),
            Fraction(102 * index + 51),
            Fraction(34),
            Fraction(0),
            Fraction(0),
        ]
        shifted_cube = [
            Fraction(index**3),
            Fraction(3 * index**2),
            Fraction(3 * index),
            Fraction(1),
            Fraction(0),
        ]
        inverse_denominator = inverse_shift_cube(index + 1)

        for sequence in (direct, companion):
            numerator = add(
                multiply(shifted_b, sequence[index]),
                scale(multiply(shifted_cube, sequence[index - 1]), -1),
            )
            sequence.append(multiply(numerator, inverse_denominator))

    return direct, companion


def fraction_mod(value: Fraction, modulus: int) -> int:
    assert gcd(value.denominator, modulus) == 1
    return value.numerator * pow(value.denominator, -1, modulus) % modulus


def evaluate_mod(
    polynomial: list[Fraction], argument: int, modulus: int
) -> int:
    return sum(
        fraction_mod(coefficient, modulus) * pow(argument, degree, modulus)
        for degree, coefficient in enumerate(polynomial)
    ) % modulus


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prime-limit", type=int, default=499)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    primes = [
        prime
        for prime in primes_at_most(args.prime_limit)
        if prime >= 7
    ]
    maximum_prime = max(primes)
    values = apery_numbers(2 * maximum_prime)
    direct, companion = shifted_fundamental_solutions(maximum_prime - 1)
    checks: Counter[str] = Counter()
    target_rows: dict[int, list[int]] = defaultdict(list)
    local_combined: dict[tuple[int, int], int] = {}

    for prime in primes:
        modulus3 = prime**3
        modulus4 = prime**4
        modulus5 = prime**5
        modulus6 = prime**6
        endpoint_difference = values[prime - 1] - 1
        assert endpoint_difference % modulus3 == 0
        endpoint_scalar = endpoint_difference // modulus3 % prime
        inverse_five = pow(5, -1, modulus4)

        for index in range(prime):
            if values[index] % prime != 0:
                continue
            reflected_index = prime - 1 - index
            assert values[reflected_index] % prime == 0
            upper_index = prime + index
            assert values[upper_index] % prime == 0
            quotient = values[upper_index] // prime % modulus4
            checks["target_digits"] += 1
            target_rows[upper_index].append(prime)

            direct_a = companion[index][0]
            direct_c = companion[index][1] - 3 * direct_a
            direct_stripped = (
                5 * evaluate_mod(direct[index], prime, modulus5)
                - prime**3
                * fraction_mod(direct_a + prime * direct_c, modulus5)
            ) % modulus5
            assert direct_stripped % prime == 0
            normalized_direct = direct_stripped // prime % modulus4

            direct_full = (
                direct_stripped
                - 7
                * prime**3
                * endpoint_scalar
                * (
                    values[index]
                    + prime * fraction_mod(direct[index][1], modulus5)
                )
            ) % modulus5
            assert direct_full == values[upper_index] % modulus5
            assert normalized_direct == (
                quotient
                * (
                    1
                    + 7
                    * prime**3
                    * endpoint_scalar
                    * inverse_five
                )
            ) % modulus4
            checks["direct_quartic"] += 1

            reflected_a = companion[reflected_index][0]
            reflected_c = (
                companion[reflected_index][1] - 3 * reflected_a
            )
            reflected_stripped = (
                5
                * evaluate_mod(
                    direct[reflected_index], -2 * prime, modulus5
                )
                + 584
                * prime**3
                * fraction_mod(
                    reflected_a - 2 * prime * reflected_c,
                    modulus5,
                )
            ) % modulus5
            assert reflected_stripped % prime == 0
            normalized_reflected = reflected_stripped // prime % modulus4

            reflected_full = (
                reflected_stripped
                + 8
                * prime**3
                * endpoint_scalar
                * (
                    values[reflected_index]
                    - 2
                    * prime
                    * fraction_mod(
                        direct[reflected_index][1], modulus5
                    )
                )
            ) % modulus5
            assert reflected_full == values[upper_index] % modulus5
            assert normalized_reflected == (
                quotient
                * (
                    1
                    - 8
                    * prime**3
                    * endpoint_scalar
                    * inverse_five
                )
            ) % modulus4
            checks["reflected_quartic"] += 1

            combined = (
                8 * normalized_direct + 7 * normalized_reflected
            ) % modulus4
            assert combined == 15 * quotient % modulus4
            local_combined[(upper_index, prime)] = combined
            checks["bernoulli_elimination"] += 1

            direct_d = companion[index][2] - 3 * companion[index][1]
            direct_d += 6 * direct_a
            direct_stripped5 = (
                5 * evaluate_mod(direct[index], prime, modulus6)
                - prime**3
                * fraction_mod(
                    direct_a
                    + prime * direct_c
                    + prime**2 * direct_d,
                    modulus6,
                )
            ) % modulus6
            assert direct_stripped5 % prime == 0
            normalized_direct5 = direct_stripped5 // prime % modulus5
            direct_full5 = (
                direct_stripped5
                - 7
                * endpoint_difference
                * evaluate_mod(direct[index], prime, modulus6)
            ) % modulus6
            assert direct_full5 == values[upper_index] % modulus6
            assert normalized_direct5 == (
                (values[upper_index] // prime)
                * (
                    1
                    + 7
                    * endpoint_difference
                    * pow(5, -1, modulus5)
                )
            ) % modulus5
            checks["direct_quintic"] += 1

            reflected_d = (
                companion[reflected_index][2]
                - 3 * companion[reflected_index][1]
                + 6 * reflected_a
            )
            reflected_stripped5 = (
                5
                * evaluate_mod(
                    direct[reflected_index], -2 * prime, modulus6
                )
                + 584
                * prime**3
                * fraction_mod(
                    reflected_a
                    - 2 * prime * reflected_c
                    + 4 * prime**2 * reflected_d,
                    modulus6,
                )
            ) % modulus6
            assert reflected_stripped5 % prime == 0
            normalized_reflected5 = (
                reflected_stripped5 // prime % modulus5
            )
            reflected_full5 = (
                reflected_stripped5
                + 8
                * endpoint_difference
                * evaluate_mod(
                    direct[reflected_index], -2 * prime, modulus6
                )
            ) % modulus6
            assert reflected_full5 == values[upper_index] % modulus6
            assert normalized_reflected5 == (
                (values[upper_index] // prime)
                * (
                    1
                    - 8
                    * endpoint_difference
                    * pow(5, -1, modulus5)
                )
            ) % modulus5
            checks["reflected_quintic"] += 1

            combined5 = (
                8 * normalized_direct5 + 7 * normalized_reflected5
            ) % modulus5
            assert (
                combined5
                == 15 * (values[upper_index] // prime) % modulus5
            )
            local_combined[(upper_index, -prime)] = combined5
            checks["rank_one_elimination_mod_p5"] += 1

    for upper_index, targets in target_rows.items():
        radical = 1
        for prime in targets:
            radical *= prime
        assert values[upper_index] % radical == 0
        common_quotient = values[upper_index] // radical

        for prime in targets:
            modulus4 = prime**4
            radical_without_prime = radical // prime
            combined = local_combined[(upper_index, prime)]
            recovered = (
                combined
                * pow(15 * radical_without_prime, -1, modulus4)
            ) % modulus4
            assert recovered == common_quotient % modulus4
            checks["common_quotient_mod_p4"] += 1

            modulus5 = prime**5
            combined5 = local_combined[(upper_index, -prime)]
            recovered5 = (
                combined5
                * pow(15 * radical_without_prime, -1, modulus5)
            ) % modulus5
            assert recovered5 == common_quotient % modulus5
            checks["common_quotient_mod_p5"] += 1

    print(f"prime_limit={args.prime_limit}")
    for key in sorted(checks):
        print(f"{key}_checks={checks[key]}")
    print("first_failure=None")
    print("failures=0")


if __name__ == "__main__":
    main()
