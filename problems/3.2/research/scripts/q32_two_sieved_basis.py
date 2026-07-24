#!/usr/bin/env python3
"""Pure-Python audit of Q515's two-sieved central-binomial basis.

For epsilon in {0,1}, expand

    A_(2m+epsilon) = sum_{k=0}^m r_k^(epsilon)
                     * binom(2m+epsilon,m-k).

The inverse triangular transform below computes the integral residual
sequences.  A top-half prime p, with r=n-p and d=m-r, sees the moving tail

    sum_{h=0}^r binom(r,h) r_(d+h)^(epsilon) = A_n = 5*A_r (mod p).

This verifies that the basis is an exact change of coordinates, not a fixed
small certificate.
"""

from __future__ import annotations

from fractions import Fraction
from math import comb

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


def inverse_even(k: int, m: int) -> Fraction:
    if k == 0:
        return Fraction(int(m == 0))
    return (
        (-1) ** (k - m)
        * Fraction(2 * k, k + m)
        * comb(k + m, k - m)
    )


def inverse_odd(k: int, m: int) -> Fraction:
    return (
        (-1) ** (k - m)
        * Fraction(2 * k + 1, k + m + 1)
        * comb(k + m + 1, k - m)
    )


def residuals(apery: list[int], limit: int) -> tuple[list[int], list[int]]:
    even: list[int] = []
    odd: list[int] = []
    for k in range(limit + 1):
        even_value = sum(
            inverse_even(k, m) * apery[2 * m] for m in range(k + 1)
        )
        odd_value = sum(
            inverse_odd(k, m) * apery[2 * m + 1] for m in range(k + 1)
        )
        assert even_value.denominator == odd_value.denominator == 1
        even.append(even_value.numerator)
        odd.append(odd_value.numerator)
    return even, odd


def main() -> None:
    apery = apery_numbers(LIMIT + 1)
    half_limit = LIMIT // 2
    even, odd = residuals(apery, half_limit)
    primes = primes_up_to(LIMIT)

    assert even[:7] == [
        1,
        71,
        32711,
        21263474,
        16196884679,
        13494506759471,
        11910357240848882,
    ]
    assert odd[:6] == [
        5,
        1430,
        811805,
        578594525,
        463454152550,
        398546130989165,
    ]

    for m in range(half_limit + 1):
        assert apery[2 * m] == sum(
            even[k] * comb(2 * m, m - k) for k in range(m + 1)
        )
        if 2 * m + 1 <= LIMIT:
            assert apery[2 * m + 1] == sum(
                odd[k] * comb(2 * m + 1, m - k) for k in range(m + 1)
            )

    for epsilon, coefficients in ((0, even), (1, odd)):
        for m in range(1, half_limit + 1):
            n = 2 * m + epsilon
            if n > LIMIT:
                break
            for prime in primes:
                if prime <= n // 2:
                    continue
                if prime > n:
                    break
                raw_index = n - prime
                start = m - raw_index
                tail = sum(
                    comb(raw_index, h) * coefficients[start + h]
                    for h in range(raw_index + 1)
                )
                assert (tail - apery[n]) % prime == 0
                assert (tail - 5 * apery[raw_index]) % prime == 0

    # The dominant residual root is the largest real root of this reciprocal
    # quartic.  Consecutive ratios converge to about 1151.998..., not to 1.
    for coefficients in (even, odd):
        ratio = coefficients[-1] / coefficients[-2]
        assert 1100 < ratio < 1200
        print(f"last_ratio={ratio:.12f}")


if __name__ == "__main__":
    main()
