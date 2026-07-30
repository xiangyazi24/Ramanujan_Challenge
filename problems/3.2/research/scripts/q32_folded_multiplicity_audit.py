#!/usr/bin/env python3
"""Audit the branch/CRT formulas behind folded-node multiplicity pruning.

This is a dependency-free regression checker for Section 64.  It checks
candidate prime assignments only; the Selberg upper-bound sieve used in
the proof is analytic and is not inferred from the finite scan.
"""

from __future__ import annotations

import argparse
from math import gcd, isqrt


def primes_through(limit):
    sieve = bytearray(b"\x01") * (limit + 1)
    if limit >= 0:
        sieve[0] = 0
    if limit >= 1:
        sieve[1] = 0
    for prime in range(2, isqrt(limit) + 1):
        if sieve[prime]:
            start = prime * prime
            sieve[start : limit + 1 : prime] = b"\x00" * (
                (limit - start) // prime + 1
            )
    return [value for value in range(2, limit + 1) if sieve[value]]


def audit(maximum_n, density_limit):
    primes = primes_through(maximum_n)
    small_odd_primes = [prime for prime in primes if 3 <= prime <= 47]

    assignments = 0
    branch_pairs = 0
    central_pairs = 0
    density_pairs = 0
    density_prime_checks = 0
    high_cutoff_checks = 0

    for n in range(5, maximum_n + 1):
        direct = {}
        reflected = {}

        for prime in primes:
            if prime <= isqrt(n):
                continue
            if prime > n:
                break

            quotient, residue = divmod(n, prime)
            folded = min(residue, prime - 1 - residue)

            if residue == folded:
                assert (n - folded) // quotient == prime
                assert (n - folded) % quotient == 0
                assert (2 * quotient + 1) * folded <= n - quotient
                direct.setdefault(folded, []).append((prime, quotient))

            if prime - 1 - residue == folded:
                reflected_label = quotient + 1
                assert (
                    (n + 1 + folded) // reflected_label == prime
                )
                assert (n + 1 + folded) % reflected_label == 0
                assert (
                    (2 * reflected_label - 1) * folded
                    <= n - reflected_label + 1
                )
                reflected.setdefault(folded, []).append(
                    (prime, reflected_label)
                )

            assignments += 1

        assert all(len(values) <= 1 for values in direct.values())
        assert all(len(values) <= 1 for values in reflected.values())

        cutoff = max(2, isqrt(isqrt(n)))
        for folded in direct.keys() & reflected.keys():
            (direct_prime, a), = direct[folded]
            (reflected_prime, b_label), = reflected[folded]

            common = gcd(a, b_label)
            assert (2 * n + 1) % common == 0
            assert (a + b_label) % 2 == 1
            assert (
                a * direct_prime + b_label * reflected_prime
                == 2 * n + 1
            )
            assert folded % (2 * a) == (n - a) % (2 * a)
            assert folded % (2 * b_label) == (
                b_label - n - 1
            ) % (2 * b_label)

            u = a // common
            v = b_label // common
            total = (2 * n + 1) // common
            assert (
                u * direct_prime + v * reflected_prime == total
            )

            if max(a, b_label) > cutoff:
                assert (2 * cutoff - 1) * folded < n
                high_cutoff_checks += 1

            if (
                density_pairs < density_limit
                and gcd(total, u * v) == 1
            ):
                for ell in small_odd_primes:
                    roots = sum(
                        (
                            (direct_prime - 2 * v * parameter)
                            * (
                                reflected_prime
                                + 2 * u * parameter
                            )
                        )
                        % ell
                        == 0
                        for parameter in range(ell)
                    )
                    expected = 1 if (u * v * total) % ell == 0 else 2
                    assert roots == expected
                    density_prime_checks += 1
                density_pairs += 1

            branch_pairs += 1
            central_pairs += direct_prime == reflected_prime

    return {
        "minimum_n": 5,
        "maximum_n": maximum_n,
        "mesoscopic_prime_assignments": assignments,
        "direct_reflected_branch_pairs": branch_pairs,
        "central_branch_pairs": central_pairs,
        "primitive_density_pairs": density_pairs,
        "density_prime_checks": density_prime_checks,
        "high_label_cutoff_checks": high_cutoff_checks,
        "assertions": "passed",
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--maximum-n", type=int, default=5000)
    parser.add_argument("--density-limit", type=int, default=5000)
    args = parser.parse_args()
    if args.maximum_n < 5:
        raise SystemExit("--maximum-n must be at least 5")
    if args.density_limit < 0:
        raise SystemExit("--density-limit must be nonnegative")

    for key, value in audit(
        args.maximum_n, args.density_limit
    ).items():
        print(f"{key:35s} {value}")


if __name__ == "__main__":
    main()
