#!/usr/bin/env python3
"""Audit the one-dimensional weight-seven finite-MHS endpoint quotient.

For p >= 11, write

    xi  = H(6)/p mod p,
    eta = H(2,4)/p mod p,
    A   = H(2,2,3) mod p,
    B   = H(2,5) mod p.

The exact finite multiple-harmonic-sum identities proved in the
companion research note give

    3 eta = 2 xi,
    3 A   = 14 xi,
    2 B   = -7 xi                         (mod p).

Thus the conservative four-coordinate endpoint presentation at
effective weight seven has actual image of dimension at most one.

This script is only an exact modular audit of those symbolic identities.
It also checks the Bernoulli normalizations and the lifted reversal
identity from which the eta relation follows.
"""

from __future__ import annotations

import argparse
from collections import Counter
from math import isqrt


def primes_at_most(limit: int) -> list[int]:
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[:2] = b"\x00\x00"
    for prime in range(2, isqrt(limit) + 1):
        if sieve[prime]:
            sieve[prime * prime : limit + 1 : prime] = b"\x00" * (
                (limit - prime * prime) // prime + 1
            )
    return [
        candidate
        for candidate in range(2, limit + 1)
        if sieve[candidate]
    ]


def multiple_harmonic_sum(
    exponents: tuple[int, ...], prime: int, modulus: int
) -> int:
    """Return the strict sum H(exponents) modulo ``modulus``."""

    depth = len(exponents)
    partials = [0] * (depth + 1)
    partials[0] = 1
    for value in range(1, prime):
        inverse = pow(value, -1, modulus)
        powers = [
            pow(inverse, exponent, modulus)
            for exponent in exponents
        ]
        for position in range(depth, 0, -1):
            partials[position] = (
                partials[position]
                + partials[position - 1] * powers[position - 1]
            ) % modulus
    return partials[depth]


def bernoulli_mod(index: int, prime: int) -> int:
    """Return B_index modulo p, using B_1=-1/2.

    The routine is used only with 0 <= index < p-1, so all denominators
    in the triangular Bernoulli recurrence are p-adic units.
    """

    inverses = [0, 1] + [
        pow(value, -1, prime) for value in range(2, index + 2)
    ]
    values = [0] * (index + 1)
    values[0] = 1
    for degree in range(1, index + 1):
        binomial = 1
        total = values[0]
        for lower in range(1, degree):
            binomial = (
                binomial
                * (degree + 2 - lower)
                * inverses[lower]
            ) % prime
            total = (total + binomial * values[lower]) % prime
        values[degree] = (
            -total * inverses[degree + 1]
        ) % prime
    return values[index]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prime-limit", type=int, default=1000)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    checks: Counter[str] = Counter()

    for prime in primes_at_most(args.prime_limit):
        if prime < 11:
            continue
        modulus2 = prime**2

        s2 = multiple_harmonic_sum((2,), prime, modulus2)
        s4 = multiple_harmonic_sum((4,), prime, modulus2)
        s6 = multiple_harmonic_sum((6,), prime, modulus2)
        h24 = multiple_harmonic_sum((2, 4), prime, modulus2)
        h42 = multiple_harmonic_sum((4, 2), prime, modulus2)

        h23 = multiple_harmonic_sum((2, 3), prime, prime)
        h223 = multiple_harmonic_sum((2, 2, 3), prime, prime)
        h232 = multiple_harmonic_sum((2, 3, 2), prime, prime)
        h25 = multiple_harmonic_sum((2, 5), prime, prime)
        h43 = multiple_harmonic_sum((4, 3), prime, prime)
        h52 = multiple_harmonic_sum((5, 2), prime, prime)

        assert s2 % prime == s4 % prime == s6 % prime == 0
        assert h24 % prime == h42 % prime == 0
        checks["integrality"] += 5

        # Exact stuffle identities at the required precisions.
        assert (s2 * s4 - h24 - h42 - s6) % modulus2 == 0
        assert (
            (s2 % prime) * h23
            - 2 * h223
            - h232
            - h43
            - h25
        ) % prime == 0
        assert 2 * h232 % prime == 0
        checks["stuffle_and_reversal"] += 3

        # The first-order correction to even-weight reversal.
        assert (
            h24
            - h42
            - prime * (2 * h43 + 4 * h52)
        ) % modulus2 == 0
        checks["lifted_reversal"] += 1

        beta = bernoulli_mod(prime - 7, prime)
        xi = s6 // prime % prime
        eta = h24 // prime % prime

        assert (7 * xi - 6 * beta) % prime == 0
        assert (h25 + 3 * beta) % prime == 0
        assert (h43 + 5 * beta) % prime == 0
        assert (h52 - 3 * beta) % prime == 0
        checks["bernoulli_normalization"] += 4

        assert (3 * eta - 2 * xi) % prime == 0
        assert (3 * h223 - 14 * xi) % prime == 0
        assert (2 * h25 + 7 * xi) % prime == 0
        checks["rank_one_relation"] += 3

    print(f"prime_limit={args.prime_limit}")
    print(
        "checked_primes="
        f"{len([p for p in primes_at_most(args.prime_limit) if p >= 11])}"
    )
    for name in sorted(checks):
        print(f"{name}_checks={checks[name]}")
    print("relations=3*eta-2*xi,3*A-14*xi,2*B+7*xi")
    print("status=exact_audit_of_symbolic_finite_mhs_proof")
    print("failures=0")


if __name__ == "__main__":
    main()
