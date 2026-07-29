#!/usr/bin/env python3
"""Compare full q=1 branch content with the simple carrier gcd.

The full Legendre--Euler coefficient family can remove nuisance factors from
gcd(S_n,B_n^-/+), but it cannot remove a q=1 bad top-half prime.  This exact
scan records where the two ideals differ.
"""

from __future__ import annotations

from math import comb, gcd

from q32_fixed_q_content import truncation_content
from q32_strehl_gcd import franel_numbers


LIMIT = 400


def main() -> None:
    franel = franel_numbers(LIMIT)
    differences = []
    for n in range(3, LIMIT + 1):
        cutoff = (n - 1) // 3
        boundary = cutoff + 1
        strehl = sum(
            comb(n, index)
            * comb(n + index, index)
            * franel[index]
            for index in range(cutoff + 1)
        )
        content = truncation_content(n, 1, franel)
        carriers = (
            ("direct", comb(n, boundary)),
            ("reflected", comb(n + boundary, boundary)),
        )
        for branch, carrier in carriers:
            full_ideal = gcd(content, carrier)
            simple_ideal = gcd(strehl, carrier)
            assert simple_ideal % full_ideal == 0
            if full_ideal != simple_ideal:
                differences.append(
                    (
                        n,
                        branch,
                        full_ideal,
                        simple_ideal,
                        simple_ideal // full_ideal,
                    )
                )

    print(f"differences={len(differences)} through n={LIMIT}")
    for record in differences[:20]:
        print(record)


if __name__ == "__main__":
    main()
