#!/usr/bin/env python3
"""Audit folded truncation certificates for one fixed quotient slice.

Let q be odd and n=q*p+r with 0<=r<p.  Apéry reflection folds r to

    j=min(r,p-1-r) <= floor((n-q)/(2*q+1)).

The direct Apéry prefix and the Strehl--Franel prefix through this cutoff
both reduce to A_j modulo p.  Their gcd therefore contains every bad prime
in the fixed-q slice.  This script measures that gcd and verifies the local
congruences; it does not prove a subexponential bound.
"""

from __future__ import annotations

from math import comb, factorial, gcd, log

from q32_newton import apery_numbers
from q32_strehl_gcd import franel_numbers, primes_up_to, valuation


LIMIT = 1200
QUOTIENT = 3


def certificates(
    n: int, quotient: int, franel: list[int]
) -> tuple[int, int]:
    cutoff = (n - quotient) // (2 * quotient + 1)
    direct = 0
    strehl = 0
    for k in range(cutoff + 1):
        kernel = comb(n, k) * comb(n + k, k)
        direct += kernel * kernel
        strehl += kernel * franel[k]
    return direct, strehl


def main() -> None:
    assert QUOTIENT > 0 and QUOTIENT % 2 == 1
    apery = apery_numbers(LIMIT)
    franel = franel_numbers(LIMIT)
    primes = primes_up_to(LIMIT)
    records: list[tuple[int, int, float]] = []

    for n in range(QUOTIENT, LIMIT + 1):
        direct, strehl = certificates(n, QUOTIENT, franel)
        common = gcd(direct, strehl)
        cutoff = (n - QUOTIENT) // (2 * QUOTIENT + 1)
        boundary_carrier = (
            comb(n, cutoff + 1) * comb(n + cutoff + 1, cutoff + 1)
        )
        scaled_difference = factorial(cutoff) ** 2 * (direct - strehl)
        assert scaled_difference % boundary_carrier == 0
        for prime in primes:
            quotient, residue = divmod(n, prime)
            if quotient < QUOTIENT:
                continue
            if quotient > QUOTIENT:
                if prime <= n // (QUOTIENT + 1):
                    continue
                break
            folded = min(residue, prime - 1 - residue)
            assert folded <= (n - QUOTIENT) // (2 * QUOTIENT + 1)
            assert direct % prime == apery[folded] % prime
            assert strehl % prime == apery[folded] % prime
            assert boundary_carrier % prime == 0
            if prime >= 7:
                expected_carrier_valuation = (
                    2
                    if 2 * n + 1 == (2 * QUOTIENT + 1) * prime
                    else 1
                )
                assert (
                    valuation(boundary_carrier, prime)
                    == expected_carrier_valuation
                )
            assert (common % prime == 0) == (apery[folded] % prime == 0)

        rate = log(common) / n if common > 1 else 0.0
        records.append((n, common, rate))

    lower = 10
    while lower < LIMIT:
        upper = min(2 * lower, LIMIT)
        winner = max(
            (record for record in records if lower < record[0] <= upper),
            key=lambda record: record[2],
        )
        print(
            f"({lower},{upper}] max_log_gcd_over_n={winner[2]:.9f} "
            f"at_n={winner[0]} gcd={winner[1]}"
        )
        lower *= 2


if __name__ == "__main__":
    main()
