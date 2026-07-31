#!/usr/bin/env python3
"""Audit the first p-adic quotient of the adjacent Pascal multipliers.

For 0 <= s < L < p, the binomial coefficient ``binom(p+s, L)``
contains exactly one copy of p.  The quotient has the closed expansion

    binom(p+s, L) / p
      = (-1)^(t-1)/(t binom(L,t))
        prod_{j=1}^s (1+p/j)
        prod_{j=1}^{t-1} (1-p/j),   t=L-s.

The script checks this identity modulo p^2, specializes it to the two
Pascal multipliers in Section 68.18 of Q32_SEPARATION_ANALYSIS.md, and
audits the boundary zero-segment example (n,p)=(147,73).
"""

from fractions import Fraction
from math import comb
from pathlib import Path
import sys


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from q32_cartier_packet_audit import apery, shell_batch


HOSTILE_CASES = (
    (200, 128, 63),
    (272, 180, 63),
    (300, 180, 57),
    (321, 168, 53),
)


def primes_below(limit):
    answer = []
    for candidate in range(2, limit):
        if all(candidate % divisor for divisor in range(2, int(candidate**0.5) + 1)):
            answer.append(candidate)
    return answer


def harmonic(index):
    return sum(
        (Fraction(1, denominator) for denominator in range(1, index + 1)),
        Fraction(0),
    )


def rational_mod(value, modulus):
    assert value.denominator % modulus != 0
    return (
        value.numerator
        * pow(value.denominator, -1, modulus)
        % modulus
    )


def first_quotient_prediction(prime, surplus, length):
    """Return binom(p+surplus,length)/p modulo p^2."""

    assert 0 <= surplus < length < prime
    complement = length - surplus
    leading = Fraction(
        (-1) ** (complement - 1),
        complement * comb(length, complement),
    )
    return rational_mod(
        leading
        * (
            1
            + prime
            * (
                harmonic(surplus)
                - harmonic(complement - 1)
            )
        ),
        prime * prime,
    )


def check_quotient(total, length, prime):
    surplus = total - prime
    assert 0 <= surplus < length < prime and total < 2 * prime
    actual = (comb(total, length) // prime) % (prime * prime)
    predicted = first_quotient_prediction(prime, surplus, length)
    assert actual == predicted


def finite_difference(values, start, order, modulus):
    return sum(
        (-1) ** (order - shift)
        * comb(order, shift)
        * values[start + shift]
        for shift in range(order + 1)
    ) % modulus


def audit_universal_formula():
    checks = 0
    for prime in primes_below(102):
        if prime == 2:
            continue
        for length in range(1, prime):
            for surplus in range(length):
                check_quotient(prime + surplus, length, prime)
                checks += 1
    return checks


def audit_hostile_blocks():
    multiplier_checks = 0
    target_records = []
    for index, D, N in HOSTILE_CASES:
        moment = index - 1
        margin = min(
            D - moment // 2,
            moment - D - N + 2,
        )
        d = D - margin + 1
        length = N + margin - 2
        nodes = range(d, D + length + 1)
        candidates = [
            prime
            for prime in primes_below(D + N + 1)
            if D < prime <= D + N
        ]

        for prime in candidates:
            check_quotient(D + N, length, prime)
            check_quotient(D + length, length, prime)
            multiplier_checks += 2

        targets = [
            prime
            for prime in candidates
            if apery(index - prime) % prime == 0
        ]
        for prime in targets:
            values = shell_batch(moment, nodes, prime)
            left_packet = finite_difference(values, d, length, prime)
            right_packet = finite_difference(values, D, length, prime)
            left_quotient = comb(D + N, length) // prime
            right_quotient = comb(D + length, length) // prime
            sign = -1 if length % 2 else 1
            left_increment = sign * left_quotient * left_packet % prime
            right_increment = sign * right_quotient * right_packet % prime

            # These nonvanishing assertions are calibration for the four
            # hostile blocks, not part of the universal independence proof.
            assert left_increment != 0
            assert right_increment != 0
            target_records.append(
                (
                    index,
                    prime,
                    left_increment,
                    right_increment,
                )
            )
    return multiplier_checks, tuple(target_records)


def audit_boundary_zero_segment():
    index, prime = 147, 73
    moment = index - 1
    assert apery(2) % prime == 0
    values = shell_batch(moment, range(74, 146), prime)
    assert all(value % prime == 0 for value in values.values())

    D, N = (13 * index) // 20, index // 5
    margin = min(
        D - moment // 2,
        moment - D - N + 2,
    )
    assert (D, N, moment, margin) == (95, 29, 146, 22)
    supports = []
    for current_margin in (margin - 1, margin):
        d = D - current_margin + 1
        length = N + current_margin - 2
        assert length < prime
        supports.append(
            (
                current_margin,
                (d, d + length),
                (D, D + length),
            )
        )
    return tuple(supports)


def main():
    universal_checks = audit_universal_formula()
    multiplier_checks, target_records = audit_hostile_blocks()
    boundary_supports = audit_boundary_zero_segment()
    print("UNIVERSAL_QUOTIENT_CHECKS", universal_checks)
    print("HOSTILE_MULTIPLIER_CHECKS", multiplier_checks)
    print("HOSTILE_TARGET_INCREMENTS", target_records)
    print("BOUNDARY_ZERO_SEGMENT_SUPPORTS", boundary_supports)
    print("Q32_PASCAL_FIRST_QUOTIENT_AUDIT=PASS")


if __name__ == "__main__":
    main()
