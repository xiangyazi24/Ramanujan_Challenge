#!/usr/bin/env python3
"""Exact audit of the adjacent signed Newton--Padé invariant.

For a fixed interpolation height H, let (P_b,Q_b) be the primitive pair
with

    deg P_b <= H-b,  deg Q_b <= b,
    P_b(j) = A_j Q_b(j),  0 <= j <= H.

The cross product of adjacent pairs vanishes at all interpolation nodes:

    P_b Q_(b+1) - P_(b+1) Q_b
        = kappa_b * binom(x,H+1).

Because the left side is integer-valued, kappa_b is the exact integer
obtained by evaluating it at x=H+1.  This script measures its size, the
gcd of adjacent numerator certificates at n=3H+1, and the direct-prime
support of those gcds.  It is an experiment, not an asymptotic proof.
"""

from __future__ import annotations

from functools import reduce
from math import comb, gcd, log

from q32_pade_total_positivity import (
    apery_values,
    entry,
    newton_coefficients,
    primitive_pade_kernel,
)


def evaluate_newton(coefficients: list[int], x: int) -> int:
    return sum(value * comb(x, index) for index, value in enumerate(coefficients))


def primitive_pair(
    height: int,
    denominator_degree: int,
    differences: list[int],
) -> tuple[list[int], list[int]]:
    numerator_degree = height - denominator_degree
    denominator = primitive_pade_kernel(
        height,
        denominator_degree,
        differences,
    )
    numerator = [
        sum(
            denominator[ell] * entry(k, ell, differences)
            for ell in range(min(denominator_degree, k) + 1)
        )
        for k in range(numerator_degree + 1)
    ]
    return numerator, denominator


def is_prime(value: int) -> bool:
    if value < 2:
        return False
    if value % 2 == 0:
        return value == 2
    divisor = 3
    while divisor * divisor <= value:
        if value % divisor == 0:
            return False
        divisor += 2
    return True


def direct_prime_support(
    height: int,
    n: int,
    apery: list[int],
    values: list[int],
) -> tuple[list[int], list[int]]:
    bad: list[int] = []
    good: list[int] = []
    common = reduce(gcd, (abs(value) for value in values), 0)
    for j in range(height + 1):
        prime = n - j
        if not is_prime(prime) or common % prime:
            continue
        if apery[j] % prime == 0:
            bad.append(prime)
        else:
            good.append(prime)
    return bad, good


def log_integer(value: int) -> float:
    """Stable natural logarithm for arbitrarily large positive integers."""

    value = abs(value)
    if value <= 1:
        return 0.0
    shift = max(0, value.bit_length() - 53)
    return log(value >> shift) + shift * log(2)


def audit_height(height: int, apery: list[int], differences: list[int]) -> None:
    n = 3 * height + 1
    pairs = [
        primitive_pair(height, denominator_degree, differences)
        for denominator_degree in range(height + 1)
    ]
    numerator_values = [
        evaluate_newton(numerator, n)
        for numerator, _ in pairs
    ]

    kappas: list[int] = []
    for denominator_degree in range(height):
        numerator, denominator = pairs[denominator_degree]
        next_numerator, next_denominator = pairs[denominator_degree + 1]
        kappa = (
            evaluate_newton(numerator, height + 1)
            * evaluate_newton(next_denominator, height + 1)
            - evaluate_newton(next_numerator, height + 1)
            * evaluate_newton(denominator, height + 1)
        )
        assert kappa
        # Only P_b Q_(b+1) reaches monomial degree H+1.  Comparing
        # leading coefficients in the binomial basis gives this exact
        # primitive-normalized formula.
        assert kappa == (
            comb(height + 1, denominator_degree + 1)
            * numerator[-1]
            * next_denominator[-1]
        )
        kappas.append(kappa)

        # The identity is checked at two points beyond the interpolation grid.
        for x in (height + 2, n):
            cross = (
                evaluate_newton(numerator, x)
                * evaluate_newton(next_denominator, x)
                - evaluate_newton(next_numerator, x)
                * evaluate_newton(denominator, x)
            )
            assert cross == kappa * comb(x, height + 1)

    adjacent = [
        gcd(abs(numerator_values[b]), abs(numerator_values[b + 1]))
        for b in range(height)
    ]
    best_b = min(range(height), key=lambda b: adjacent[b])
    bad, good = direct_prime_support(
        height,
        n,
        apery,
        numerator_values[best_b : best_b + 2],
    )
    assert all(kappas[best_b] % prime == 0 for prime in good)

    triple = [
        reduce(
            gcd,
            (
                abs(numerator_values[b]),
                abs(numerator_values[b + 1]),
                abs(numerator_values[b + 2]),
            ),
        )
        for b in range(max(0, height - 1))
    ]
    if triple:
        best_triple_b = min(range(len(triple)), key=lambda b: triple[b])
        triple_rate = log_integer(triple[best_triple_b]) / n
    else:
        best_triple_b = -1
        triple_rate = 0.0

    print(
        f"H={height:2d} n={n:3d} "
        f"min_adj_b={best_b:2d} "
        f"log_gcd2/n={log_integer(adjacent[best_b]) / n:.6f} "
        f"log_kappa/n={log_integer(kappas[best_b]) / n:.6f} "
        f"min_tri_b={best_triple_b:2d} "
        f"log_gcd3/n={triple_rate:.6f} "
        f"bad={bad} good_pollution={good}"
    )


def main() -> None:
    maximum_height = 30
    apery = apery_values(3 * maximum_height + 3)
    differences = newton_coefficients(apery)
    for height in range(2, maximum_height + 1):
        audit_height(height, apery, differences)


if __name__ == "__main__":
    main()
