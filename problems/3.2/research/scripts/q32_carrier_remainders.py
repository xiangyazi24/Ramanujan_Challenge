#!/usr/bin/env python3
"""Test Euclidean compression against the two q=1 factorial carriers.

For J=floor((n-1)/3) and m=J+1, put

    S_n = sum_{k<=J} binom(n,k) binom(n+k,k) F_k,
    B_n^- = binom(n,m),
    B_n^+ = binom(n+m,m).

The direct q=1 candidates are carried by B_n^- and the reflected candidates
by B_n^+.  If a centered remainder S_n modulo either carrier had
subexponential height, it would prove the corresponding branch.  This script
measures those remainders; it is a negative diagnostic, not a proof.
"""

from __future__ import annotations

from math import comb, log

from q32_strehl_gcd import franel_numbers


LIMIT = 1200


def centered_remainder(value: int, modulus: int) -> int:
    residue = value % modulus
    return min(residue, modulus - residue)


def main() -> None:
    franel = franel_numbers(LIMIT // 3 + 2)
    records: list[tuple[int, float, float, float, float]] = []

    for n in range(3, LIMIT + 1):
        cutoff = (n - 1) // 3
        boundary = cutoff + 1
        strehl = sum(
            comb(n, k) * comb(n + k, k) * franel[k]
            for k in range(cutoff + 1)
        )
        direct_carrier = comb(n, boundary)
        reflected_carrier = comb(n + boundary, boundary)
        direct_remainder = centered_remainder(strehl, direct_carrier)
        reflected_remainder = centered_remainder(strehl, reflected_carrier)
        full_carrier = direct_carrier * reflected_carrier
        full_remainder = centered_remainder(strehl, full_carrier)
        records.append(
            (
                n,
                log(direct_remainder) / n if direct_remainder else 0.0,
                log(reflected_remainder) / n
                if reflected_remainder
                else 0.0,
                log(full_remainder) / n if full_remainder else 0.0,
                log(full_carrier) / n,
            )
        )

    lower = 10
    while lower < LIMIT:
        upper = min(2 * lower, LIMIT)
        slab = [
            record for record in records if lower < record[0] <= upper
        ]
        last = slab[-1]
        print(
            f"({lower},{upper}] last_n={last[0]} "
            f"direct={last[1]:.9f} reflected={last[2]:.9f} "
            f"full={last[3]:.9f} carrier={last[4]:.9f}"
        )
        lower *= 2


if __name__ == "__main__":
    main()
