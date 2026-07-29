#!/usr/bin/env python3
"""Expose the quotient alias in the reflected cutoff-jump construction.

For a reflected zero r=p-1-j and n=q*p+r, the low nested-cutoff jump and
the high nested-cutoff jump do not coincide only for q=1.  With

    H=floor((n-1)/3),
    sigma(h)=ceil((n+1+3h)/2),

the same-index construction aliases every quotient

    q == 1 (mod 3).

The first counterexample uses the Apéry zero 17|A_3:

    q=4, n=81:  low jump = high jump = 20,
    q=7, n=132: low jump = high jump = 37.

The q=1 reflected branch can still be purified by intersecting with the
explicit lcm-interval carrier

    lcm(1,...,floor(2n/3)) / lcm(1,...,floor(n/2)).

For primes above sqrt(n), this carrier has exactly the prime support
(n/2,2n/3], so it removes q=4,7,... aliases at an unavoidable positive
PNT height n/6+o(n).

There is also a cleaner quotient filter using every cutoff transition.
If e=min(q,p-1-q) is the effective degree of Q_q modulo p, the bad prime is
lost exactly e times as J runs from 0 to n.  Consequently primes lost only
once have q=1 or q=p-2.  The latter is a near-square-root nuisance at a
fixed outer index.
"""

from __future__ import annotations

from math import isqrt, lcm

from q32_legendre_content import franel_numbers
from q32_q1_all_cutoff_profile import (
    content_is_zero,
    transform_coefficients_mod,
)


PRIME = 17
FOLDED = 3
QUOTIENTS = range(1, 8)


def ceiling_half(value: int) -> int:
    return (value + 1) // 2


def cutoff_jumps(
    zero_cutoffs: set[int], n: int
) -> tuple[list[int], list[int]]:
    low_limit = (n - 1) // 3
    high_limit = (n - 3) // 3

    low_tails = [
        all(cutoff in zero_cutoffs for cutoff in range(h, low_limit + 1))
        for h in range(low_limit + 1)
    ]

    def sigma(h: int) -> int:
        return ceiling_half(n + 1 + 3 * h)

    high_tails = [
        all(cutoff in zero_cutoffs for cutoff in range(sigma(h), n))
        for h in range(high_limit + 1)
    ]
    low_jumps = [
        h
        for h in range(1, len(low_tails))
        if low_tails[h] and not low_tails[h - 1]
    ]
    high_jumps = [
        h
        for h in range(1, len(high_tails))
        if high_tails[h] and not high_tails[h - 1]
    ]
    return low_jumps, high_jumps


def lcm_through(limit: int) -> int:
    result = 1
    for value in range(1, limit + 1):
        result = lcm(result, value)
    return result


def main() -> None:
    residue = PRIME - 1 - FOLDED
    limit = max(QUOTIENTS) * PRIME + residue
    franel = franel_numbers(limit)
    records = []

    for quotient in QUOTIENTS:
        n = quotient * PRIME + residue
        zero_cutoffs = {
            cutoff
            for cutoff in range(n + 1)
            if content_is_zero(
                transform_coefficients_mod(
                    n, cutoff, PRIME, franel
                ),
                PRIME,
            )
        }
        low_jumps, high_jumps = cutoff_jumps(zero_cutoffs, n)
        common = sorted(set(low_jumps) & set(high_jumps))
        losses = [
            cutoff
            for cutoff in range(1, n + 1)
            if cutoff - 1 in zero_cutoffs
            and cutoff not in zero_cutoffs
        ]
        effective_degree = min(quotient, PRIME - 1 - quotient)
        assert len(losses) == effective_degree

        upper = (2 * n) // 3
        middle = n // 2
        interval_carrier = (
            lcm_through(upper) // lcm_through(middle)
        )
        carrier_contains = interval_carrier % PRIME == 0
        if PRIME > isqrt(n):
            assert carrier_contains == (middle < PRIME <= upper)

        records.append(
            (
                quotient,
                n,
                low_jumps,
                high_jumps,
                common,
                carrier_contains,
                losses,
            )
        )

    for quotient, n, low, high, common, carrier, losses in records:
        expected_common = quotient % 3 == 1
        assert bool(common) == expected_common
        assert carrier == (quotient == 1)
        print(
            f"q={quotient} n={n} low={low} high={high} "
            f"same={common} losses={losses} "
            f"interval_carrier={carrier}"
        )

    # The loss count is independent of the direct/reflected low digit.
    for quotient in (1, 4):
        n = quotient * PRIME + FOLDED
        zero_cutoffs = {
            cutoff
            for cutoff in range(n + 1)
            if content_is_zero(
                transform_coefficients_mod(
                    n, cutoff, PRIME, franel
                ),
                PRIME,
            )
        }
        losses = [
            cutoff
            for cutoff in range(1, n + 1)
            if cutoff - 1 in zero_cutoffs
            and cutoff not in zero_cutoffs
        ]
        assert losses == [
            multiple * PRIME
            for multiple in range(1, quotient + 1)
        ]
        print(f"direct q={quotient} n={n} losses={losses}")


if __name__ == "__main__":
    main()
