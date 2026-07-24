#!/usr/bin/env python3
"""Audit the fixed-n Laurent coefficient-shell pair certificate from Q607.

For

    Lambda = (1+x)(1+y)(1+z)((1+y)(1+z)+xyz)/(xyz)

write C_N(e) for the coefficient of x^e1 y^e2 z^e3 in Lambda^N and

    E_N(m) = sum_{v in {-1,0,1}^3} C_N(m v).

Five ambient directions vanish identically in the top-half shell, leaving
the 22 directions in the support of Lambda.  If p is a q=1 hit at row n,
then p divides E_(n-1)(p-1).  Two hits p,q therefore give the elementary
cross-weighted certificate

    q E_(n-1)(p-1) - p E_(n-1)(q-1),

which is divisible by pq.  This script checks the n=321 reflected triple and
shows that the evident certificate retains full linear exponential height.
"""

from __future__ import annotations

from itertools import product
from math import comb, log


def coefficient(n: int, exponent: tuple[int, int, int]) -> int:
    """Exact one-sum formula for [x^a y^b z^c] Lambda^n."""

    first, second, third = exponent
    total = 0
    for index in range(n + 1):
        arguments = (
            index,
            index - first,
            n - second,
            n - third,
        )
        tops = (n, n, 2 * n - index, 2 * n - index)
        if any(
            argument < 0 or argument > top
            for argument, top in zip(arguments, tops)
        ):
            continue
        total += (
            comb(n, index)
            * comb(n, index - first)
            * comb(2 * n - index, n - second)
            * comb(2 * n - index, n - third)
        )
    return total


def shell(n: int, spacing: int) -> int:
    return sum(
        coefficient(
            n,
            tuple(spacing * coordinate for coordinate in direction),
        )
        for direction in product((-1, 0, 1), repeat=3)
    )


def main() -> None:
    row = 321
    power = row - 1
    hits = (179, 193, 211)
    values = {prime: shell(power, prime - 1) for prime in hits}

    for prime in hits:
        value = values[prime]
        assert value % prime == 0
        print(
            f"p={prime} rate={log(value) / row:.12f} "
            f"bits={value.bit_length()}"
        )

    for larger, smaller in ((193, 179), (211, 193), (211, 179)):
        certificate = (
            smaller * values[larger] - larger * values[smaller]
        )
        assert certificate % (larger * smaller) == 0
        print(
            f"p={larger} q={smaller} h={larger-smaller} "
            f"rate={log(abs(certificate)) / row:.12f} "
            f"bits={abs(certificate).bit_length()}"
        )


if __name__ == "__main__":
    main()
