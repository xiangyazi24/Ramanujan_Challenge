#!/usr/bin/env python3
"""Audit the numerator-degree 1/3 anchored Padé cross determinant.

The primitive Padé pairs use the integral binomial basis, not ``ZZ[x]``.
Consequently

    P_1 Q_3 - P_3 Q_1 = Phi_H K_H

has a rational linear quotient K_H in the ordinary power basis.  Clearing
its denominator gives an integer polynomial without changing reductions at
candidate primes p>2H, because the denominator has no prime factor above H.

For a good candidate prime common to P_1(3H+1) and P_3(3H+1), differentiating
the cross identity at the moving node proves p divides the cleared value
K_H(3H+1).  This script checks the exact identity and that support implication.
It is an audit, not an asymptotic estimate.
"""

from __future__ import annotations

import argparse
from math import factorial

from sympy import Poly, QQ, symbols

from q32_adjacent_pade_kappa import (
    apery_values,
    evaluate_newton,
    is_prime,
    newton_coefficients,
)
from q32_pade_family_gcd import primitive_pair_fast


X = symbols("x")


def power_polynomial(binomial_coefficients: list[int]) -> Poly:
    """Convert an integral binomial-basis polynomial to ``QQ[x]``."""

    falling_factorial = 1
    expression = 0
    for degree, coefficient in enumerate(binomial_coefficients):
        if degree:
            falling_factorial *= X - (degree - 1)
        expression += QQ(coefficient, factorial(degree)) * falling_factorial
    return Poly(expression, X, domain=QQ)


def node_polynomial(height: int) -> Poly:
    result = Poly(1, X, domain=QQ)
    for node in range(height + 1):
        result *= Poly(X - node, X, domain=QQ)
    return result


def largest_prime_factor(value: int) -> int:
    value = abs(value)
    largest = 1
    divisor = 2
    while divisor * divisor <= value:
        while value % divisor == 0:
            largest = divisor
            value //= divisor
        divisor += 1 if divisor == 2 else 2
    return max(largest, value)


def audit_height(
    height: int,
    apery: list[int],
    differences: list[int],
) -> tuple[int, int, int]:
    numerator_1, denominator_1 = primitive_pair_fast(height, 1, differences)
    numerator_3, denominator_3 = primitive_pair_fast(height, 3, differences)

    cross = (
        power_polynomial(numerator_1) * power_polynomial(denominator_3)
        - power_polynomial(numerator_3) * power_polynomial(denominator_1)
    )
    quotient, remainder = divmod(cross, node_polynomial(height))
    assert remainder.is_zero
    assert not quotient.is_zero
    assert quotient.degree() <= 1

    clearing_denominator, integer_quotient = quotient.clear_denoms(convert=True)
    clearing_denominator = int(clearing_denominator)
    assert largest_prime_factor(clearing_denominator) <= height

    outer_index = 3 * height + 1
    target_count = 0
    good_count = 0
    for node in range(height + 1):
        prime = outer_index - node
        if not is_prime(prime):
            continue
        assert clearing_denominator % prime
        common = (
            evaluate_newton(numerator_1, outer_index) % prime == 0
            and evaluate_newton(numerator_3, outer_index) % prime == 0
        )
        if not common:
            continue
        if apery[node] % prime == 0:
            target_count += 1
        else:
            good_count += 1
            assert int(integer_quotient.eval(outer_index)) % prime == 0

    coefficient_bits = max(
        abs(int(coefficient)).bit_length()
        for coefficient in integer_quotient.all_coeffs()
    )
    print(
        f"H={height} denominator={clearing_denominator} "
        f"K_bits={coefficient_bits} targets={target_count} good={good_count}"
    )
    return target_count, good_count, coefficient_bits


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("maximum_height", nargs="?", type=int, default=32)
    args = parser.parse_args()
    if args.maximum_height < 3:
        raise SystemExit("maximum_height must be at least 3")

    apery = apery_values(3 * args.maximum_height + 4)
    differences = newton_coefficients(apery)
    totals = [0, 0]
    for height in range(3, args.maximum_height + 1):
        target_count, good_count, _ = audit_height(
            height,
            apery,
            differences,
        )
        totals[0] += target_count
        totals[1] += good_count
    print(
        f"heights=3..{args.maximum_height} "
        f"target_checks={totals[0]} good_checks={totals[1]}"
    )


if __name__ == "__main__":
    main()
