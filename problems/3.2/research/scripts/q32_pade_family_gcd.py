#!/usr/bin/env python3
"""Measure the exact gcd of a growing family of primitive Padé numerators.

For height H and numerator degree a, let (P_{H,a}, Q_{H,a}) be the
primitive integral Newton-basis interpolation pair

    deg P_{H,a} <= a,  deg Q_{H,a} <= H-a,
    P_{H,a}(s) = A_s Q_{H,a}(s),  0 <= s <= H.

This script computes

    G_{H,A}(3H+1) = gcd_{0 <= a <= A} P_{H,a}(3H+1)

for the zero-absorber scale A=ceil(H^(2/3)).  It reports both the total
gcd height and its radical in the candidate window 2H<p<=3H+1.  The
calculation is exploratory evidence, not an asymptotic estimate.
"""

from __future__ import annotations

import argparse
from math import gcd, log

from q32_adjacent_pade_kappa import (
    apery_values,
    evaluate_newton,
    is_prime,
    log_integer,
    newton_coefficients,
    primitive_pair,
)


def ceil_two_thirds(height: int) -> int:
    """Return the least A with A^3 >= H^2."""

    cutoff = 0
    target = height * height
    while cutoff**3 < target:
        cutoff += 1
    return cutoff


def candidate_data(
    height: int,
    apery: list[int],
) -> tuple[list[tuple[int, int, int]], int]:
    """Return (prime, moving node, prefix-zero count) and max zero count."""

    n = 3 * height + 1
    data: list[tuple[int, int, int]] = []
    maximum = 0
    for node in range(height + 1):
        prime = n - node
        if not is_prime(prime):
            continue
        zero_count = sum(
            apery[index] % prime == 0
            for index in range(height + 1)
        )
        data.append((prime, node, zero_count))
        maximum = max(maximum, zero_count)
    return data, maximum


def factor_small(value: int) -> str:
    """Format the exact factorization when the residual gcd is small."""

    value = abs(value)
    if value == 0 or value.bit_length() > 128:
        return "-"
    factors: list[str] = []
    divisor = 2
    while divisor * divisor <= value:
        exponent = 0
        while value % divisor == 0:
            value //= divisor
            exponent += 1
        if exponent:
            factors.append(
                str(divisor) if exponent == 1 else f"{divisor}^{exponent}"
            )
        divisor = 3 if divisor == 2 else divisor + 2
    if value > 1:
        factors.append(str(value))
    return "*".join(factors) if factors else "1"


def audit_height(
    height: int,
    apery: list[int],
    differences: list[int],
) -> None:
    n = 3 * height + 1
    cutoff = min(height, ceil_two_thirds(height))
    data, maximum_zero_count = candidate_data(height, apery)

    common = 0
    milestones = {
        0,
        min(1, cutoff),
        min(2, cutoff),
        cutoff // 2,
        cutoff,
    }
    milestone_rows: list[tuple[int, int, float]] = []

    for numerator_degree in range(cutoff + 1):
        numerator, denominator = primitive_pair(
            height,
            height - numerator_degree,
            differences,
        )
        for node in range(height + 1):
            assert evaluate_newton(numerator, node) == (
                apery[node] * evaluate_newton(denominator, node)
            )
        value = evaluate_newton(numerator, n)
        common = gcd(common, abs(value))
        if numerator_degree in milestones:
            milestone_rows.append(
                (
                    numerator_degree,
                    common.bit_length(),
                    log_integer(common) / height,
                )
            )

    candidate_primes = [
        prime for prime, _, _ in data if common % prime == 0
    ]
    target_primes = [
        prime
        for prime, node, _ in data
        if apery[node] % prime == 0
    ]
    expected_primes = [
        prime
        for prime, node, zero_count in data
        if apery[node] % prime == 0 or zero_count > cutoff
    ]
    assert candidate_primes == expected_primes
    if cutoff >= maximum_zero_count:
        assert candidate_primes == target_primes

    candidate_radical = 1
    for prime in candidate_primes:
        candidate_radical *= prime

    milestones_text = ",".join(
        f"a={degree}:bits={bits}:log/H={rate:.6f}"
        for degree, bits, rate in milestone_rows
    )
    print(
        f"H={height} A={cutoff} max_z={maximum_zero_count} "
        f"G_bits={common.bit_length()} "
        f"G_factor={factor_small(common)} "
        f"log_G/H={log_integer(common) / height:.6f} "
        f"window_count={len(candidate_primes)} "
        f"log_window/H={log_integer(candidate_radical) / height:.6f} "
        f"target_count={len(target_primes)} "
        f"milestones=[{milestones_text}]"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "heights",
        nargs="*",
        type=int,
        default=[10, 15, 20, 25, 30],
    )
    args = parser.parse_args()
    if not args.heights or min(args.heights) < 2:
        raise SystemExit("all heights must be at least 2")

    maximum_height = max(args.heights)
    apery = apery_values(3 * maximum_height + 3)
    differences = newton_coefficients(apery)
    for height in args.heights:
        audit_height(height, apery, differences)


if __name__ == "__main__":
    main()
