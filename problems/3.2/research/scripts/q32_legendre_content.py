#!/usr/bin/env python3
"""Polynomial-content audit for the Legendre--Euler truncation family.

For J=floor((n-1)/3), Q_n(t)=sum_k L(n,k)t^k, and
F_i=sum_a binom(i,a)^3, put

    T_n(c) = sum_{m=0}^J [y^m]Q_n(c+y)
             * sum_{i=0}^m binom(m,i)(-c)^(m-i) F_i.

This script expands T_n(c) exactly.  For every top-half prime p and folded
index j it verifies the stronger polynomial congruence

    T_n(c) = A_j * (1 + 2*c^p)  (mod p).

Thus p divides the content (the gcd of all monomial coefficients) exactly
when p divides A_j.  The construction is a clean global package, but the
coefficient heights remain exponential.
"""

from __future__ import annotations

from math import comb, gcd, log

from q32_newton import apery_numbers


LIMIT = 300


def primes_up_to(limit: int) -> list[int]:
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[:2] = b"\x00\x00"
    for prime in range(2, int(limit**0.5) + 1):
        if sieve[prime]:
            sieve[prime * prime : limit + 1 : prime] = b"\x00" * (
                (limit - prime * prime) // prime + 1
            )
    return [value for value, is_prime in enumerate(sieve) if is_prime]


def franel_numbers(limit: int) -> list[int]:
    return [
        sum(comb(index, part) ** 3 for part in range(index + 1))
        for index in range(limit + 1)
    ]


def truncation_coefficients(n: int, franel: list[int]) -> list[int]:
    """Return the monomial coefficients of T_n(c)."""

    cutoff = (n - 1) // 3
    coefficients = [
        sum(
            comb(n, k) * comb(n + k, k) * franel[k]
            for k in range(cutoff + 1)
        )
    ] + [0] * n

    for degree in range(1, n + 1):
        lower = max(0, cutoff + 1 - degree)
        upper = min(cutoff, n - degree)
        coefficients[degree] = sum(
            (-1) ** (cutoff - index)
            * comb(n, index + degree)
            * comb(n + index + degree, index + degree)
            * comb(index + degree, index)
            * comb(degree - 1, cutoff - index)
            * franel[index]
            for index in range(lower, upper + 1)
        )
    return coefficients


def direct_transformed_value(
    n: int, center: int, franel: list[int]
) -> int:
    cutoff = (n - 1) // 3
    result = 0
    for m in range(cutoff + 1):
        kernel = sum(
            comb(n, k)
            * comb(n + k, k)
            * comb(k, m)
            * center ** (k - m)
            for k in range(m, n + 1)
        )
        shifted_franel = sum(
            comb(m, i) * (-center) ** (m - i) * franel[i]
            for i in range(m + 1)
        )
        result += kernel * shifted_franel
    return result


def evaluate(coefficients: list[int], value: int) -> int:
    result = 0
    for coefficient in reversed(coefficients):
        result = result * value + coefficient
    return result


def main() -> None:
    apery = apery_numbers(LIMIT)
    franel = franel_numbers(LIMIT)
    primes = primes_up_to(LIMIT)
    records: list[tuple[int, int, float]] = []

    for n in range(3, LIMIT + 1):
        coefficients = truncation_coefficients(n, franel)
        if n <= 40:
            for center in (-1, 0, 1, 2, 3):
                assert evaluate(coefficients, center) == direct_transformed_value(
                    n, center, franel
                )

        content = 0
        for coefficient in coefficients:
            content = gcd(content, coefficient)

        for prime in primes:
            if prime <= n // 2:
                continue
            if prime > n:
                break
            raw_index = n - prime
            folded_index = min(raw_index, prime - 1 - raw_index)
            lower_value = apery[folded_index] % prime
            for degree, coefficient in enumerate(coefficients):
                expected = lower_value if degree == 0 else 0
                if degree == prime:
                    expected = 2 * lower_value
                assert coefficient % prime == expected % prime
            assert (content % prime == 0) == (lower_value == 0)

        rate = log(content) / n if content > 1 else 0.0
        records.append((n, content, rate))

    lower = 10
    while lower < LIMIT:
        upper = min(2 * lower, LIMIT)
        winner = max(
            (record for record in records if lower < record[0] <= upper),
            key=lambda record: record[2],
        )
        print(
            f"({lower},{upper}] max_log_content_over_n={winner[2]:.9f} "
            f"at_n={winner[0]} content={winner[1]}"
        )
        lower *= 2


if __name__ == "__main__":
    main()
