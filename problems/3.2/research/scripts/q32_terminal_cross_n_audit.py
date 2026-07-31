#!/usr/bin/env python3
"""Exact audit of the terminal n versus n-2 seam.

The script verifies the lower doubled-period seam law

    C_{p+r-1}(p+i) = J_{r-1}(i)              (mod p),
    F_{s+1}^{(p+r)} = S_r                    (mod p),

and checks that direct cross-Casoratians carry only the one universal
candidate-prime copy at the hostile target rows.
"""

from math import comb, gcd
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from q32_cartier_packet_audit import (  # noqa: E402
    LAMBDA,
    apery,
    coefficient,
    shell_fast,
)
from q32_terminal_family_audit import newton, primes_up_to  # noqa: E402
from q32_terminal_turan_hankel_audit import (  # noqa: E402
    terminal_data,
    valuation,
)


def doubled_period(moment, dilation, modulus=None):
    """J_m(i) = CT(Lambda(X)^m Lambda(X^{-i}))."""

    value = sum(
        lambda_coefficient
        * coefficient(
            moment,
            *(dilation * coordinate for coordinate in ray),
            modulus=modulus,
        )
        for ray, lambda_coefficient in LAMBDA.items()
    )
    return value if modulus is None else value % modulus


def seam_scalar(residue, modulus=None):
    """The lower scalar S_r in the first post-selector terminal value."""

    value = sum(
        (-1) ** index
        * comb(residue, index + 1)
        * doubled_period(residue - 1, index, modulus)
        for index in range(residue)
    )
    return value if modulus is None else value % modulus


def dense_seam_audit():
    coefficient_checks = 0
    seam_checks = 0
    for prime in primes_up_to(31):
        if prime < 5:
            continue
        for residue in range(1, prime):
            moment = prime + residue - 1
            values = {}
            for index in range(residue):
                node = prime + index
                actual = shell_fast(moment, node, modulus=prime)
                expected = doubled_period(
                    residue - 1, index, modulus=prime
                )
                assert actual == expected
                values[node] = actual
                coefficient_checks += 1
            actual_seam = newton(values, prime, residue - 1) % prime
            expected_seam = seam_scalar(residue, modulus=prime)
            assert actual_seam == expected_seam
            seam_checks += 1
    return coefficient_checks, seam_checks


def hostile_cross_audit():
    rows = []
    hostile = (
        (168, (97,)),
        (200, (139, 181)),
        (272, (191, 233)),
        (300, (191, 227)),
        (321, (179, 193, 211)),
    )
    width = 5
    for n, targets in hostile:
        start, _, upper, _ = terminal_data(n)
        _, _, lower, _ = terminal_data(n - 2)
        upper_aligned = upper[1 : width + 2]
        lower_aligned = lower[: width + 1]
        cross = [
            upper_aligned[index] * lower_aligned[index + 1]
            - upper_aligned[index + 1] * lower_aligned[index]
            for index in range(width)
        ]
        cross_gcd = 0
        for value in cross:
            cross_gcd = gcd(cross_gcd, value)

        safe_candidates = tuple(
            prime
            for prime in primes_up_to(n - 2)
            if start + width + 1 < prime
        )
        for prime in safe_candidates:
            assert all(value % prime == 0 for value in cross)
        for prime in targets:
            if prime not in safe_candidates:
                continue
            assert [valuation(value, prime) for value in cross] == [
                1
            ] * width

            residue = n - prime
            seam = prime - 1 - start
            assert upper[seam] % prime == apery(residue) % prime == 0
            assert (
                upper[seam + 1] - seam_scalar(residue, modulus=prime)
            ) % prime == 0

        rows.append(
            (
                n,
                targets,
                abs(cross_gcd).bit_length(),
                tuple(
                    (
                        prime,
                        tuple(valuation(value, prime) for value in cross),
                    )
                    for prime in targets
                    if prime in safe_candidates
                ),
            )
        )
    return rows


if __name__ == "__main__":
    coefficient_checks, seam_checks = dense_seam_audit()
    print(
        "DENSE_SEAM_COEFFICIENT_CHECKS",
        coefficient_checks,
        "DENSE_NEWTON_SEAM_CHECKS",
        seam_checks,
    )
    for row in hostile_cross_audit():
        print("HOSTILE_CROSS", row)
    print("Q32_TERMINAL_CROSS_N_AUDIT=PASS")
