#!/usr/bin/env python3
"""Scan the folded Newton gcd proposed in Q451.

This is evidence only.  The q=1 radical divides

    gcd(A_n, F_J(n) F_J(-n-1)),  J=floor((n-1)/3),

but the gcd may also contain small primes and primes larger than n.
"""

from __future__ import annotations

import math
from math import gcd

from q32_newton import (
    apery_numbers,
    evaluate_newton,
    forward_differences,
    log_abs,
)


LIMIT = 1200


def main() -> None:
    apery = apery_numbers(LIMIT)
    records: list[tuple[int, int, float]] = []
    for n in range(6, LIMIT + 1):
        j = (n - 1) // 3
        coefficients = forward_differences(apery[: j + 1])
        positive = evaluate_newton(coefficients, n)
        negative = evaluate_newton(coefficients, -n - 1)
        common = gcd(apery[n], positive * negative)
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
