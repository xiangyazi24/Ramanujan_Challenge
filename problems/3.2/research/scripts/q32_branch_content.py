#!/usr/bin/env python3
"""Measure the strongest q=1 branch ideals in the coefficient family.

Let Gamma_n be the content of the full Legendre--Euler coefficient vector,
and put m=floor((n-1)/3)+1.  The exact branch ideals are

    Gamma_n^- = gcd(Gamma_n, binom(n,m)),
    Gamma_n^+ = gcd(Gamma_n, binom(n+m,m)).

Their top-half prime supports are respectively the direct and reflected q=1
bad primes.  This script verifies the support and reports dyadic maxima.
"""

from __future__ import annotations

from math import comb, gcd, log

from q32_fixed_q_content import truncation_content
from q32_newton import apery_numbers
from q32_strehl_gcd import franel_numbers, primes_up_to


LIMIT = 400


def main() -> None:
    franel = franel_numbers(LIMIT)
    apery = apery_numbers(LIMIT)
    primes = primes_up_to(LIMIT)
    records: list[list[tuple[int, int, float]]] = [[], [], []]

    for n in range(3, LIMIT + 1):
        boundary = (n - 1) // 3 + 1
        content = truncation_content(n, 1, franel)
        values = (
            gcd(content, comb(n, boundary)),
            gcd(content, comb(n + boundary, boundary)),
            content,
        )
        for branch, value in enumerate(values):
            records[branch].append(
                (
                    n,
                    value,
                    log(value) / n if value > 1 else 0.0,
                )
            )

        for prime in primes:
            if not n / 2 < prime <= n:
                continue
            residue = n - prime
            folded = min(residue, prime - 1 - residue)
            bad = apery[folded] % prime == 0
            assert (values[0] % prime == 0) == (
                bad and residue < boundary
            )
            assert (values[1] % prime == 0) == (
                bad and prime - residue <= boundary
            )

    lower = 10
    names = ("direct", "reflected", "content")
    while lower < LIMIT:
        upper = min(2 * lower, LIMIT)
        print(f"range=({lower},{upper}]")
        for name, branch_records in zip(names, records):
            winner = max(
                (
                    record
                    for record in branch_records
                    if lower < record[0] <= upper
                ),
                key=lambda record: record[2],
            )
            print(
                f"{name}: rate={winner[2]:.9f} "
                f"n={winner[0]} value={winner[1]}"
            )
        lower *= 2


if __name__ == "__main__":
    main()
