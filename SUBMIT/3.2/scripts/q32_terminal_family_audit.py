#!/usr/bin/env python3
"""Exact audit of the terminal Newton family and its Pascal saturation."""

from math import comb, gcd
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from q32_cartier_packet_audit import (  # noqa: E402
    coefficient,
    polytope_points,
    shell_batch,
)


def primes_up_to(limit):
    return [
        p
        for p in range(2, limit + 1)
        if all(p % q for q in range(2, int(p**0.5) + 1))
    ]


def product(values):
    answer = 1
    for value in values:
        answer *= value
    return answer


def newton(values, start, order):
    return sum(
        (-1) ** i
        * comb(start + i, i)
        * comb(start + order + 1, order - i)
        * values[start + i]
        for i in range(order + 1)
    )


def finite_difference(values, start, order):
    return sum(
        (-1) ** (order - i)
        * comb(order, i)
        * values[start + i]
        for i in range(order + 1)
    )


def direct_terminal_ct(moment, order):
    """Return the coefficient expansion of the terminal packet B_L."""

    return sum(
        (-1) ** residue
        * comb(order, residue)
        * coefficient(
            moment,
            (moment - residue) * kappa[0],
            (moment - residue) * kappa[1],
            (moment - residue) * kappa[2],
        )
        for kappa in polytope_points(1)
        for residue in range(order + 1)
    )


def candidate_primorial(moment, width):
    start = moment // 2 + 1
    return product(
        p
        for p in primes_up_to(moment + 1)
        if start + width < p <= moment + 1
    )


def audit_small():
    checks = 0
    for moment in range(6, 23):
        start = moment // 2 + 1
        maximum = moment - start
        values = shell_batch(moment, range(start, moment + 1))
        terminal = {}
        boundary = {}
        for order in range(maximum + 1):
            node = moment - order
            terminal[order] = newton(values, node, order)
            boundary[order] = finite_difference(values, node, order)
            assert boundary[order] == direct_terminal_ct(moment, order)
            checks += 1
            if order:
                assert terminal[order] - terminal[order - 1] == (
                    (-1) ** order
                    * comb(moment + 1, order)
                    * boundary[order]
                )
                checks += 1

        for width in range(1, maximum + 1):
            pascal_gcd = 0
            for shift in range(width):
                pascal_gcd = gcd(
                    pascal_gcd,
                    comb(moment + 1, maximum - shift),
                )
            expected = candidate_primorial(moment, width)
            actual = product(
                p
                for p in primes_up_to(moment + 1)
                if p > (moment + 1) // 2 and pascal_gcd % p == 0
            )
            assert actual == expected
            checks += 1
    return checks


def hostile_records():
    records = []
    for moment in (199, 271, 299, 320):
        start = moment // 2 + 1
        maximum = moment - start
        width = min(38, maximum)
        values = shell_batch(moment, range(start, moment + 1))
        family = [
            newton(values, start + shift, maximum - shift)
            for shift in range(width + 1)
        ]
        family_gcd = 0
        for value in family:
            family_gcd = gcd(family_gcd, value)
        high = tuple(
            p
            for p in primes_up_to(moment + 1)
            if start + width < p <= moment + 1
            and family_gcd % p == 0
        )
        records.append((moment, width, family_gcd, high))
    return records


if __name__ == "__main__":
    for moment, width, family_gcd, high in hostile_records():
        print(
            "HOSTILE",
            moment,
            "WIDTH",
            width,
            "GCD",
            family_gcd,
            "HIGH_PRIMES",
            high,
        )
    print("TERMINAL_FAMILY_CHECKS", audit_small())
    print("Q32_TERMINAL_FAMILY_AUDIT=PASS")
