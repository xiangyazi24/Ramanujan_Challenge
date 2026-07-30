#!/usr/bin/env python3
"""Exact pair-preserving Newton-margin gcd histories.

For prime endpoints q=d+1 and ell=d+L+1, define

    H_{s,t}=G_{d-s,L+s+t}.

Every stencil contains the fixed core [d,d+L], so all endpoint/core
targets survive while left and right margins are enlarged.  This script
checks the two difference laws and the primitive Pascal cell law, then
computes cumulative rectangle gcds for the five hostile rows.
"""

from argparse import ArgumentParser
from math import comb, gcd, prod

from q32_cartier_packet_audit import shell_batch
from q32_newton_gcd_audit import (
    carrier_from_values,
    forward_difference,
)
from q32_translated_stack_audit import factor_small


INTERVALS = (
    (200, 139, 181, (139, 181)),
    (272, 191, 233, (191, 233)),
    (300, 191, 227, (191, 227)),
    (321, 179, 193, (179, 193)),
    (321, 193, 211, (193, 211)),
    (321, 179, 211, (179, 193, 211)),
    (755, 593, 733, (593, 733)),
)


def available_margins(moment, d, length, requested):
    """Stay in the first quotient cell and keep every stencil length <q."""

    left = min(requested, d - (moment // 2 + 1))
    right = min(requested, moment - (d + length))
    while left >= 0 and length + left + right >= d + 1:
        if left >= right:
            left -= 1
        else:
            right -= 1
    assert left >= 0 and right >= 0
    return left, right


def cumulative_gcd(table):
    rows = len(table)
    columns = len(table[0])
    out = [[0] * columns for _ in range(rows)]
    for s in range(rows):
        for t in range(columns):
            value = abs(table[s][t])
            if s:
                value = gcd(value, out[s - 1][t])
            if t:
                value = gcd(value, out[s][t - 1])
            out[s][t] = value
    return out


def factor_bounded(number, bound=10000):
    """Factor small primes without trial-dividing a huge carrier gcd."""

    number = abs(number)
    factors = {}
    divisor = 2
    while divisor <= bound and divisor * divisor <= number:
        while number % divisor == 0:
            factors[divisor] = factors.get(divisor, 0) + 1
            number //= divisor
        divisor = 3 if divisor == 2 else divisor + 2
    if number > 1:
        if number.bit_length() <= 40:
            factors.update(factor_small(number))
        else:
            factors[f"cofactor<{number.bit_length()} bits>"] = 1
    return factors


def analyze_interval(n, q, ell, preserved_targets, requested):
    moment = n - 1
    d = q - 1
    length = ell - q
    S, T = available_margins(moment, d, length, requested)
    values = shell_batch(
        moment,
        range(d - S, d + length + T + 1),
    )
    assert all(moment // node == 1 for node in values)

    table = [
        [
            carrier_from_values(
                values,
                d - s,
                length + s + t,
            )
            for t in range(T + 1)
        ]
        for s in range(S + 1)
    ]

    for s in range(S + 1):
        for t in range(T + 1):
            order = length + s + t
            node = d - s
            if s:
                difference = forward_difference(
                    [
                        values[index]
                        for index in range(node, node + order + 1)
                    ],
                    order,
                )[0]
                assert table[s][t] - table[s - 1][t] == (
                    (-1) ** order
                    * comb(ell + t, order)
                    * difference
                )
            if t:
                difference = forward_difference(
                    [
                        values[index]
                        for index in range(node, node + order + 1)
                    ],
                    order,
                )[0]
                assert table[s][t] - table[s][t - 1] == (
                    (-1) ** order
                    * comb(ell + t - 1, order)
                    * difference
                )
            if s and t:
                common = gcd(q - s, order)
                assert (
                    order // common * table[s][t]
                    + (q - s) // common * table[s - 1][t]
                    == (ell + t) // common * table[s][t - 1]
                )

    target_product = prod(preserved_targets)
    assert all(
        value % target_product == 0
        for row in table
        for value in row
    )
    gcds = cumulative_gcd(table)

    first_clean = None
    clean_points = []
    for total in range(S + T + 1):
        for s in range(S + 1):
            t = total - s
            if 0 <= t <= T and gcds[s][t] == target_product:
                clean_points.append((s, t))
                if first_clean is None:
                    first_clean = (s, t)

    diagonal = []
    for margin in range(min(S, T) + 1):
        nuisance = gcds[margin][margin] // target_product
        diagonal.append((margin, nuisance, factor_bounded(nuisance)))

    final_nuisance = gcds[S][T] // target_product
    print(
        "MARGIN",
        n,
        (q, ell),
        "targets",
        preserved_targets,
        "shape",
        (S, T),
        "first_clean",
        first_clean,
        "final",
        (final_nuisance, factor_bounded(final_nuisance)),
    )
    print("  diagonal", diagonal)
    return first_clean, final_nuisance


def main():
    parser = ArgumentParser()
    parser.add_argument("--margin", type=int, default=10)
    args = parser.parse_args()

    for interval in INTERVALS:
        analyze_interval(*interval, args.margin)
    print("PASS: pair-preserving Newton-margin audit")


if __name__ == "__main__":
    main()
