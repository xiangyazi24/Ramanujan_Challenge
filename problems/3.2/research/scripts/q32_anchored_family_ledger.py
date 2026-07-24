#!/usr/bin/env python3
"""Audit the candidate support of the anchored Newton--Padé family.

The canonical primitive pairs have integral coefficients in the binomial
basis.  Their ordinary power-basis cross quotient is generally rational:

    P_1 Q_b - P_b Q_1 = Phi_H K_b,  deg K_b <= b-2.

All denominators of K_b have prime factors at most H, so reduction at a
candidate prime p>2H is well defined.  This script checks:

* the exact anchored support split into targets, z=0 pollution, z=2
  pollution, and the content absorber z>B;
* the derivative-defect identity controlling the converse from K_b(j)=0
  to P_b(j)=0;
* all interpolation, quotient, and denominator assertions.

This is a finite exact audit, not an asymptotic estimate.
"""

from __future__ import annotations

import argparse
from collections import Counter
from math import factorial

from sympy import Poly, QQ, Rational, symbols

from q32_adjacent_pade_kappa import (
    apery_values,
    evaluate_newton,
    is_prime,
    newton_coefficients,
)
from q32_anchored_pade_pollution import largest_prime_factor, power_polynomial
from q32_pade_family_gcd import primitive_pair_fast


X = symbols("x")


def rational_mod(value: Rational, prime: int) -> int:
    """Reduce a rational whose denominator is a unit modulo ``prime``."""

    value = Rational(value)
    denominator = int(value.q) % prime
    assert denominator
    return int(value.p) * pow(denominator, -1, prime) % prime


def node_polynomial(height: int) -> Poly:
    result = Poly(1, X, domain=QQ)
    for node in range(height + 1):
        result *= Poly(X - node, X, domain=QQ)
    return result


def audit_height(
    height: int,
    maximum_degree: int,
    apery: list[int],
    differences: list[int],
) -> Counter[str]:
    n = 3 * height + 1
    selected = [1, *range(3, maximum_degree + 1)]
    pairs = {
        degree: primitive_pair_fast(height, degree, differences)
        for degree in selected
    }
    power_pairs = {
        degree: (
            power_polynomial(numerator),
            power_polynomial(denominator),
        )
        for degree, (numerator, denominator) in pairs.items()
    }

    for numerator, denominator in pairs.values():
        for node in range(height + 1):
            assert evaluate_newton(numerator, node) == (
                apery[node] * evaluate_newton(denominator, node)
            )

    phi = node_polynomial(height)
    quotients: dict[int, Poly] = {}
    quotient_denominators: dict[int, int] = {}
    numerator_1, denominator_1 = power_pairs[1]
    for degree in range(3, maximum_degree + 1):
        numerator, denominator = power_pairs[degree]
        cross = numerator_1 * denominator - numerator * denominator_1
        quotient, remainder = divmod(cross, phi)
        assert remainder.is_zero
        assert quotient.is_zero or quotient.degree() <= degree - 2
        clearing_denominator, _ = quotient.clear_denoms(convert=True)
        clearing_denominator = int(clearing_denominator)
        assert largest_prime_factor(clearing_denominator) <= height
        quotients[degree] = quotient
        quotient_denominators[degree] = clearing_denominator

    derivative_1 = numerator_1.diff()
    denominator_derivative_1 = denominator_1.diff()
    counts: Counter[str] = Counter()

    for node in range(height + 1):
        prime = n - node
        if not is_prime(prime):
            continue
        counts["candidates"] += 1
        zero_count = sum(apery[index] % prime == 0 for index in range(height + 1))
        target = apery[node] % prime == 0
        numerator_zeros = {
            degree: evaluate_newton(pairs[degree][0], n) % prime == 0
            for degree in selected
        }
        common = all(numerator_zeros.values())

        if target:
            counts["target"] += 1
            assert common
        elif common:
            if zero_count == 0:
                counts["E0"] += 1
            elif zero_count == 2:
                counts["E2"] += 1
            elif zero_count > maximum_degree:
                counts["Egt"] += 1
            else:
                raise AssertionError(
                    f"unclassified pollution H={height} p={prime} z={zero_count}"
                )

        if not target and zero_count > maximum_degree:
            assert common

        anchor_zero = numerator_zeros[1]
        if not anchor_zero or target:
            continue

        q1_at_node = rational_mod(denominator_1.eval(node), prime)
        assert q1_at_node == 0
        defect = (
            rational_mod(derivative_1.eval(node), prime)
            - (apery[node] % prime)
            * rational_mod(denominator_derivative_1.eval(node), prime)
        ) % prime
        if defect == 0:
            counts["anchor_derivative_defect"] += 1

        all_quotients_zero = True
        for degree, quotient in quotients.items():
            quotient_zero = rational_mod(quotient.eval(node), prime) == 0
            all_quotients_zero &= quotient_zero
            if numerator_zeros[degree]:
                assert quotient_zero
            if defect:
                assert quotient_zero == numerator_zeros[degree]
            else:
                assert quotient_zero

        if all_quotients_zero and not common:
            counts["K_false_positive"] += 1

    denominator_max = max(quotient_denominators.values(), default=1)
    fields = " ".join(
        f"{key}={counts[key]}"
        for key in (
            "candidates",
            "target",
            "E0",
            "E2",
            "Egt",
            "anchor_derivative_defect",
            "K_false_positive",
        )
    )
    print(
        f"H={height} B={maximum_degree} max_K_den={denominator_max} {fields}"
    )
    return counts


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--maximum-degree", "-B", type=int, default=6)
    parser.add_argument("heights", nargs="*", type=int)
    args = parser.parse_args()
    if not args.heights:
        args.heights = list(range(args.maximum_degree, 33))
    if args.maximum_degree < 3:
        raise SystemExit("maximum degree must be at least 3")
    if min(args.heights) < args.maximum_degree:
        raise SystemExit("each height must be at least the maximum degree")

    maximum_height = max(args.heights)
    apery = apery_values(3 * maximum_height + 3)
    differences = newton_coefficients(apery)
    totals: Counter[str] = Counter()
    for height in args.heights:
        totals.update(
            audit_height(
                height,
                args.maximum_degree,
                apery,
                differences,
            )
        )
    print(
        "TOTAL "
        + " ".join(f"{key}={value}" for key, value in sorted(totals.items()))
    )


if __name__ == "__main__":
    main()
