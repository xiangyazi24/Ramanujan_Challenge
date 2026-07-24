#!/usr/bin/env python3
"""Generalized Legendre--Euler content for a fixed quotient slice.

For fixed odd q, truncate the shifted Legendre--Euler expansion at

    J=floor((n-q)/(2*q+1)).

If floor(n/p)=q and j=min(n-q*p,p-1-(n-q*p)), then coefficientwise

    T_(n,q)(c) = A_j * Q_q(c^p) (mod p),
    Q_q(t)=sum_k binom(q,k)binom(q+k,k)t^k.

Thus every bad prime in the fixed-q slice divides the coefficient content.
The script verifies this identity and measures the content for q=3.
"""

from __future__ import annotations

from math import comb, gcd, isqrt, log

from q32_legendre_content import franel_numbers, primes_up_to
from q32_newton import apery_numbers


LIMIT = 400
QUOTIENT = 3


def truncation_coefficients(
    n: int, quotient: int, franel: list[int]
) -> list[int]:
    cutoff = (n - quotient) // (2 * quotient + 1)
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


def truncation_content(
    n: int,
    quotient: int,
    franel: list[int],
    strip_primes: list[int] | None = None,
) -> int:
    """Compute the content without constructing every enormous coefficient.

    After the constant coefficient, each new coefficient is needed only
    modulo the gcd accumulated so far.
    """

    cutoff = (n - quotient) // (2 * quotient + 1)
    result = abs(
        sum(
            comb(n, k) * comb(n + k, k) * franel[k]
            for k in range(cutoff + 1)
        )
    )
    if strip_primes is not None:
        for prime in strip_primes:
            if prime > isqrt(n):
                break
            while result % prime == 0:
                result //= prime
    for degree in range(1, n + 1):
        if result == 1:
            break
        lower = max(0, cutoff + 1 - degree)
        upper = min(cutoff, n - degree)
        residue = 0
        for index in range(lower, upper + 1):
            term = (
                comb(n, index + degree)
                * comb(n + index + degree, index + degree)
                * comb(index + degree, index)
                * comb(degree - 1, cutoff - index)
                * franel[index]
            )
            residue += (-1) ** (cutoff - index) * (term % result)
        result = gcd(result, residue % result)
    return result


def main() -> None:
    assert QUOTIENT > 0 and QUOTIENT % 2 == 1
    apery = apery_numbers(LIMIT)
    franel = franel_numbers(LIMIT)
    primes = primes_up_to(LIMIT)
    records: list[tuple[int, int, float]] = []

    for n in range(QUOTIENT, LIMIT + 1):
        content = truncation_content(n, QUOTIENT, franel)
        coefficients = None
        if n <= 160:
            coefficients = truncation_coefficients(n, QUOTIENT, franel)
            exact_content = 0
            for coefficient in coefficients:
                exact_content = gcd(exact_content, coefficient)
            assert content == exact_content

        for prime in primes:
            quotient, residue = divmod(n, prime)
            if quotient != QUOTIENT:
                continue
            folded = min(residue, prime - 1 - residue)
            folded_value = apery[folded] % prime
            if coefficients is not None:
                for degree, coefficient in enumerate(coefficients):
                    expected = 0
                    if degree % prime == 0:
                        digit = degree // prime
                        if digit <= QUOTIENT:
                            expected = (
                                comb(QUOTIENT, digit)
                                * comb(QUOTIENT + digit, digit)
                                * folded_value
                            )
                    assert coefficient % prime == expected % prime
            assert (content % prime == 0) == (folded_value == 0)

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
