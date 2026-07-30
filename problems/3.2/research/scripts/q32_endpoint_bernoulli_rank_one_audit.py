#!/usr/bin/env python3
"""Audit the rank-one mod-p^4 defect at Apéry block endpoints.

For the zeta(3) Apéry numbers b_n, put

    e_p = (b_{p-1} - 1) / p^3  (mod p).

The proved mod-p^4 closed forms are

    b_{mp} - b_m
      == p^3 e_p * m^3 (b_{m-1} - 17 b_m) / 12       (mod p^4),

    b_{mp-1} - b_{m-1}
      == p^3 e_p * m^3 (17 b_{m-1} - b_m) / 12       (mod p^4).

The same scalar is Bernoulli carried:

    e_p == -S_p/3 == 2 B_{p-3}/3                     (mod p),
    S_p = sum_{1 <= k <= (p-1)/2} k^(-3)             (mod p).

The script also audits the stronger experimental lift, valid for every
tested p>=7,

    b_{mp} - b_m
      == (b_{p-1} - 1) * m^3 (b_{m-1} - 17 b_m) / 12  (mod p^5),

and its reflected companion.  That final one-order lift is not yet
used as a theorem.
"""

from __future__ import annotations

import argparse
from collections import Counter
from math import isqrt


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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prime-limit", type=int, default=997)
    parser.add_argument("--quotient-limit", type=int, default=12)
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
    values = apery_numbers(args.quotient_limit * maximum_prime)
    checks: Counter[str] = Counter()

    for quotient in range(1, args.quotient_limit + 1):
        direct_numerator = quotient**3 * (
            values[quotient - 1] - 17 * values[quotient]
        )
        reflected_numerator = quotient**3 * (
            17 * values[quotient - 1] - values[quotient]
        )
        assert direct_numerator % 12 == 0
        assert reflected_numerator % 12 == 0
        checks["integral_closed_form"] += 2

    for prime in primes:
        modulus3 = prime**3
        modulus4 = prime**4
        modulus5 = prime**5
        endpoint_difference = values[prime - 1] - 1
        assert endpoint_difference % modulus3 == 0
        endpoint_scalar = endpoint_difference // modulus3 % prime

        half_cubic_sum = sum(
            pow(index, -3, prime)
            for index in range(1, (prime - 1) // 2 + 1)
        ) % prime
        assert (3 * endpoint_scalar + half_cubic_sum) % prime == 0
        checks["bernoulli_carrier"] += 1

        for quotient in range(1, args.quotient_limit + 1):
            direct_coefficient = (
                quotient**3
                * (values[quotient - 1] - 17 * values[quotient])
                // 12
            )
            reflected_coefficient = (
                quotient**3
                * (17 * values[quotient - 1] - values[quotient])
                // 12
            )

            direct_difference = (
                values[quotient * prime] - values[quotient]
            )
            reflected_difference = (
                values[quotient * prime - 1]
                - values[quotient - 1]
            )
            assert direct_difference % modulus3 == 0
            assert reflected_difference % modulus3 == 0
            assert (
                direct_difference
                - modulus3 * endpoint_scalar * direct_coefficient
            ) % modulus4 == 0
            assert (
                reflected_difference
                - modulus3 * endpoint_scalar * reflected_coefficient
            ) % modulus4 == 0
            checks["direct_rank_one"] += 1
            checks["reflected_rank_one"] += 1

            if prime >= 7:
                assert (
                    direct_difference
                    - endpoint_difference * direct_coefficient
                ) % modulus5 == 0
                assert (
                    reflected_difference
                    - endpoint_difference * reflected_coefficient
                ) % modulus5 == 0
                checks["direct_rank_one_mod_p5"] += 1
                checks["reflected_rank_one_mod_p5"] += 1

    print(f"prime_limit={args.prime_limit}")
    print(f"quotient_limit={args.quotient_limit}")
    for key in sorted(checks):
        print(f"{key}_checks={checks[key]}")
    print("first_failure=None")
    print("failures=0")


if __name__ == "__main__":
    main()
