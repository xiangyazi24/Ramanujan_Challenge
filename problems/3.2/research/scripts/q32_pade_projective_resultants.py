#!/usr/bin/env python3
"""Audit resultants of content-stripped Padé numerator polynomials.

For p>H, the full content theorem gives

    v_p(content(P_{H,a})) = the (a+1)-st largest v_p(A_s).

Thus, once a is at least the number of prefix zeros, a direct target
prime remains a root of the content-stripped numerator.  This script asks
whether gcds of adjacent resultants across a growing degree tail compress
those common-root primes.  Exact small cases show quadratic-size
resultant gcds and candidate-window pollution, so the naive global
resultant is not a height compression.
"""

from __future__ import annotations

import argparse
from functools import reduce
from math import factorial, gcd

from sympy import Poly, resultant, symbols

from q32_adjacent_pade_kappa import apery_values, is_prime, newton_coefficients
from q32_pade_family_gcd import (
    candidate_data,
    ceil_two_thirds,
    primitive_pair_fast,
)


VARIABLE = symbols("x")


def primitive_power_polynomial(
    binomial_coefficients: list[int],
) -> Poly:
    """Convert a primitive binomial-basis polynomial to primitive ZZ[x]."""

    degree = len(binomial_coefficients) - 1
    expression = 0
    falling_factorial = 1
    for index, coefficient in enumerate(binomial_coefficients):
        if index:
            falling_factorial *= VARIABLE - (index - 1)
        expression += (
            coefficient
            * (factorial(degree) // factorial(index))
            * falling_factorial
        )
    _, primitive = Poly(expression, VARIABLE).primitive()
    return primitive


def audit_height(
    height: int,
    apery: list[int],
    differences: list[int],
) -> None:
    data, maximum_zero_count = candidate_data(height, apery)
    lower_degree = max(1, maximum_zero_count)
    upper_degree = min(height, ceil_two_thirds(height) + 1)
    if lower_degree >= upper_degree:
        return

    polynomials: list[Poly] = []
    for numerator_degree in range(lower_degree, upper_degree + 1):
        numerator, _ = primitive_pair_fast(
            height,
            numerator_degree,
            differences,
        )
        content = abs(gcd(*numerator))
        numerator = [coefficient // content for coefficient in numerator]
        polynomials.append(primitive_power_polynomial(numerator))

    resultants = [
        abs(int(resultant(left, right, VARIABLE)))
        for left, right in zip(polynomials, polynomials[1:])
    ]
    assert all(resultants)
    common = reduce(gcd, resultants)

    n = 3 * height + 1
    support: list[str] = []
    target_primes: list[int] = []
    for prime, node, _ in data:
        if apery[node] % prime == 0:
            target_primes.append(prime)
        if common % prime == 0:
            tag = "T" if apery[node] % prime == 0 else "G"
            support.append(f"{prime}@{node}:{tag}")
    assert all(common % prime == 0 for prime in target_primes)

    print(
        f"H={height} degrees={lower_degree}..{upper_degree} "
        f"gcd_resultant_bits={common.bit_length()} "
        f"resultant_bits={min(value.bit_length() for value in resultants)}"
        f"..{max(value.bit_length() for value in resultants)} "
        f"candidate_support={','.join(support) if support else '-'}"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("maximum_height", nargs="?", type=int, default=16)
    args = parser.parse_args()
    if args.maximum_height < 5:
        raise SystemExit("maximum_height must be at least 5")

    apery = apery_values(3 * args.maximum_height + 3)
    differences = newton_coefficients(apery)
    for height in range(5, args.maximum_height + 1):
        audit_height(height, apery, differences)


if __name__ == "__main__":
    main()
