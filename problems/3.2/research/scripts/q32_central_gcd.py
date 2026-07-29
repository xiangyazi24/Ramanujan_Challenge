#!/usr/bin/env python3
"""Exact scan of gcd(A_n, binom(n,floor(n/2))).

Every prime n/2 < p <= n divides the central binomial coefficient except at
the odd boundary n=2p-1.  That boundary cannot be an Apéry zero for p>=7,
so, apart from the finite p=5 case, the top-half support of this gcd is
exactly the q=1 support.
"""

from __future__ import annotations

from math import comb, gcd

from q32_newton import apery_numbers, log_abs


LIMIT = 3000


def main() -> None:
    apery = apery_numbers(LIMIT)
    records: list[tuple[int, int, float]] = []
    for n in range(2, LIMIT + 1):
        common = gcd(apery[n], comb(n, n // 2))
        ratio = log_abs(common) / n if common > 1 else 0.0
        records.append((n, common, ratio))

    print(f"nontrivial={sum(common > 1 for _, common, _ in records)}")
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
