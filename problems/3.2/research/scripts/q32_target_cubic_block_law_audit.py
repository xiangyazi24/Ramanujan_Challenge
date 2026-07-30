#!/usr/bin/env python3
"""Audit the target-selective cubic Apéry block law modulo p^4."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from fractions import Fraction
from math import gcd, isqrt


def primes_at_most(limit: int) -> list[int]:
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[:2] = b"\x00\x00"
    for prime in range(2, isqrt(limit) + 1):
        if sieve[prime]:
            sieve[prime * prime : limit + 1 : prime] = b"\x00" * (
                (limit - prime * prime) // prime + 1
            )
    return [candidate for candidate in range(2, limit + 1) if sieve[candidate]]


def recurrence_polynomial(index: int) -> int:
    return 34 * index**3 + 51 * index**2 + 27 * index + 5


def apery_first_solution(limit: int) -> list[int]:
    values = [1, 5]
    for index in range(1, limit):
        numerator = (
            recurrence_polynomial(index) * values[-1]
            - index**3 * values[-2]
        )
        denominator = (index + 1) ** 3
        assert numerator % denominator == 0
        values.append(numerator // denominator)
    return values


def apery_second_solution(limit: int) -> list[Fraction]:
    values = [Fraction(0), Fraction(6)]
    for index in range(1, limit):
        values.append(
            (
                recurrence_polynomial(index) * values[-1]
                - index**3 * values[-2]
            )
            / (index + 1) ** 3
        )
    return values


def cubic_jets(
    limit: int,
) -> tuple[
    list[Fraction],
    list[Fraction],
    list[Fraction],
    list[Fraction],
]:
    """Return b,G,H,K with T_r=b_r+XG_r+X^2H_r+X^3K_r."""

    values = [Fraction(1), Fraction(5)]
    first = [Fraction(0), Fraction(12)]
    second = [Fraction(0), Fraction(0)]
    third = [Fraction(0), Fraction(-7)]

    for index in range(1, limit):
        coefficient = recurrence_polynomial(index)
        derivative = 102 * index**2 + 102 * index + 27
        half_second_derivative = 102 * index + 51
        cubic_coefficient = 34
        denominator = (index + 1) ** 3

        next_value = (
            coefficient * values[index]
            - index**3 * values[index - 1]
        ) / denominator
        next_first = (
            coefficient * first[index]
            + derivative * values[index]
            - index**3 * first[index - 1]
            - 3 * index**2 * values[index - 1]
            - 3 * (index + 1) ** 2 * next_value
        ) / denominator
        next_second = (
            coefficient * second[index]
            + derivative * first[index]
            + half_second_derivative * values[index]
            - index**3 * second[index - 1]
            - 3 * index**2 * first[index - 1]
            - 3 * index * values[index - 1]
            - 3 * (index + 1) ** 2 * next_first
            - 3 * (index + 1) * next_value
        ) / denominator
        next_third = (
            coefficient * third[index]
            + derivative * second[index]
            + half_second_derivative * first[index]
            + cubic_coefficient * values[index]
            - index**3 * third[index - 1]
            - 3 * index**2 * second[index - 1]
            - 3 * index * first[index - 1]
            - values[index - 1]
            - 3 * (index + 1) ** 2 * next_second
            - 3 * (index + 1) * next_first
            - next_value
        ) / denominator

        values.append(next_value)
        first.append(next_first)
        second.append(next_second)
        third.append(next_third)

    return values, first, second, third


def fraction_mod(value: Fraction, modulus: int) -> int:
    assert gcd(value.denominator, modulus) == 1
    return (
        value.numerator
        * pow(value.denominator, -1, modulus)
        % modulus
    )


def jet_mod(
    jets: tuple[
        list[Fraction],
        list[Fraction],
        list[Fraction],
        list[Fraction],
    ],
    index: int,
    argument: int,
    modulus: int,
) -> int:
    return sum(
        (
            argument**order
            * fraction_mod(sequence[index], modulus)
            for order, sequence in enumerate(jets)
        ),
        0,
    ) % modulus


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prime-limit", type=int, default=499)
    parser.add_argument("--quotient-limit", type=int, default=3)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    assert args.prime_limit >= 5
    assert args.quotient_limit >= 1

    primes = [
        prime
        for prime in primes_at_most(args.prime_limit)
        if prime >= 5
    ]
    maximum_prime = max(primes)
    exact_limit = (args.quotient_limit + 1) * maximum_prime
    values = apery_first_solution(exact_limit)
    companion = apery_second_solution(maximum_prime - 1)
    jets = cubic_jets(maximum_prime - 1)

    checks: Counter[str] = Counter()
    target_rows: dict[int, list[int]] = defaultdict(list)

    for index in range(maximum_prime):
        assert jets[0][index].denominator == 1
        assert jets[0][index].numerator == values[index]
        checks["integer_recurrence"] += 1

    for prime in primes:
        modulus3 = prime**3
        modulus4 = prime**4

        for quotient in range(args.quotient_limit + 1):
            assert (
                values[quotient * prime] - values[quotient]
            ) % modulus3 == 0
            assert (
                values[(quotient + 1) * prime - 1]
                - values[quotient]
            ) % modulus3 == 0
            checks["endpoint_supercongruence"] += 2

        for index in range(prime):
            if values[index] % prime != 0:
                continue
            checks["target_digits"] += 1
            if prime + index <= args.prime_limit:
                target_rows[prime + index].append(prime)

            companion_residue = fraction_mod(
                companion[index] / 6, modulus4
            )

            for quotient in range(1, args.quotient_limit + 1):
                argument = quotient * prime
                direct_rhs = (
                    values[quotient]
                    * jet_mod(jets, index, argument, modulus4)
                    - quotient**3
                    * prime**3
                    * values[quotient - 1]
                    * companion_residue
                ) % modulus4
                assert (
                    values[quotient * prime + index] % modulus4
                    == direct_rhs
                )
                checks["direct_target_cubic"] += 1

            for quotient in range(args.quotient_limit):
                argument = -(quotient + 1) * prime
                reflected_rhs = (
                    values[quotient]
                    * jet_mod(jets, index, argument, modulus4)
                    + (quotient + 1) ** 3
                    * prime**3
                    * values[quotient + 1]
                    * companion_residue
                ) % modulus4
                assert (
                    values[
                        (quotient + 1) * prime - 1 - index
                    ]
                    % modulus4
                    == reflected_rhs
                )
                checks["reflected_target_cubic"] += 1

            top_direct = (
                5 * jet_mod(jets, index, prime, modulus4)
                - prime**3 * companion_residue
            ) % modulus4
            assert values[prime + index] % modulus4 == top_direct

            lower_reflected = (
                jet_mod(jets, index, -prime, modulus4)
                + 5 * prime**3 * companion_residue
            ) % modulus4
            assert (
                values[prime - 1 - index] % modulus4
                == lower_reflected
            )

            upper_reflected = (
                5 * jet_mod(jets, index, -2 * prime, modulus4)
                + 584 * prime**3 * companion_residue
            ) % modulus4
            assert (
                values[2 * prime - 1 - index] % modulus4
                == upper_reflected
            )
            checks["top_specializations"] += 3

    for n, targets in target_rows.items():
        radical = 1
        for prime in targets:
            radical *= prime
        assert values[n] % radical == 0
        common_quotient = values[n] // radical
        checks["common_quotient_rows"] += 1

        for prime in targets:
            index = n - prime
            modulus3 = prime**3
            local_quotient = (
                5
                * (
                    values[index] // prime
                    + fraction_mod(jets[1][index], modulus3)
                    + prime
                    * fraction_mod(jets[2][index], modulus3)
                    + prime**2
                    * fraction_mod(jets[3][index], modulus3)
                )
                - prime**2
                * fraction_mod(companion[index] / 6, modulus3)
            ) % modulus3
            radical_without_prime = radical // prime
            normalized = (
                local_quotient
                * pow(radical_without_prime, -1, modulus3)
            ) % modulus3
            assert normalized == common_quotient % modulus3
            checks["common_quotient_mod_p3"] += 1

    print(f"prime_limit={args.prime_limit}")
    print(f"quotient_limit={args.quotient_limit}")
    for key in sorted(checks):
        print(f"{key}_checks={checks[key]}")
    print("first_failure=None")
    print("failures=0")


if __name__ == "__main__":
    main()
