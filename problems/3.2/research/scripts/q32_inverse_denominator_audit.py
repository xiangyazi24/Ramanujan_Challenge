#!/usr/bin/env python3
"""Exact audits for the inverse-denominator formulation in Section 16.

The default bound keeps the full exact-integer run short.  A bound of 10000
also works, but constructing and reducing all central binomial coefficients is
substantially slower.
"""

from __future__ import annotations

import argparse
import math


def primes_upto(limit: int) -> list[int]:
    sieve = bytearray(b"\x01") * (limit + 1)
    if limit >= 0:
        sieve[0] = 0
    if limit >= 1:
        sieve[1] = 0
    for p in range(2, math.isqrt(limit) + 1):
        if sieve[p]:
            sieve[p * p : limit + 1 : p] = b"\x00" * (
                (limit - p * p) // p + 1
            )
    return [p for p in range(2, limit + 1) if sieve[p]]


def apery_numbers(limit: int) -> list[int]:
    values = [1, 5]
    for n in range(1, limit):
        polynomial = 34 * n**3 + 51 * n**2 + 27 * n + 5
        numerator = polynomial * values[n] - n**3 * values[n - 1]
        denominator = (n + 1) ** 3
        assert numerator % denominator == 0
        values.append(numerator // denominator)
    return values[: limit + 1]


def section_anchor_audit(values: list[int], prime_bound: int = 500) -> int:
    checked = 0
    for d in range(2, 9):
        for a in range(d):
            threshold = max(50, 5 * max(values[:d]) + 1)
            for p in primes_upto(prime_bound):
                if p < threshold:
                    continue
                c = (a - p) % d
                n = p + c
                if n >= len(values):
                    continue
                m = (n - a) // d
                carrier = math.comb(n, m)
                assert c < m < p
                assert carrier % p == 0
                assert (carrier // p) % p != 0
                assert values[n] % p != 0
                checked += 1
    return checked


def coefficient_ratio_audit() -> int:
    checked = 0
    for d in range(2, 9):
        for a in range(d):
            for m in range(40):
                old_binomial = math.comb(d * m + a, m)
                new_binomial = math.comb(d * (m + 1) + a, m + 1)
                numerator = m + 1
                for i in range(1, d):
                    numerator *= (d - 1) * m + a + i
                denominator = 1
                for i in range(1, d + 1):
                    denominator *= d * m + a + i
                assert old_binomial * denominator == new_binomial * numerator
                checked += 1
    return checked


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--limit", type=int, default=2500)
    args = parser.parse_args()
    limit = args.limit
    if limit < 12:
        raise SystemExit("--limit must be at least 12")

    values = apery_numbers(max(limit + 1, 510))
    central = [math.comb(n, n // 2) for n in range(limit + 2)]
    denominators = [
        central[n] // math.gcd(central[n], values[n])
        for n in range(limit + 2)
    ]

    defect = [1] * (limit + 1)
    max_bits = (0, 0)
    max_rate = (0.0, 0)
    comparison_count = 0
    for n in range(1, limit + 1):
        neighbor_gcd = math.gcd(denominators[n - 1], denominators[n + 1])
        defect[n] = neighbor_gcd // math.gcd(neighbor_gcd, denominators[n])
        central_neighbor_gcd = math.gcd(central[n - 1], central[n + 1])
        central_triple_gcd = math.gcd(central_neighbor_gcd, central[n])
        universal_factor = central_neighbor_gcd // central_triple_gcd
        expected_factor = 1 if n % 2 == 0 else math.gcd(n // 2 + 1, 2)
        assert universal_factor == expected_factor
        apery_central_gcd = math.gcd(values[n], central[n])
        assert (apery_central_gcd * universal_factor) % defect[n] == 0
        comparison_count += 1
        if n >= 11:
            max_bits = max(max_bits, (defect[n].bit_length(), n))
            rate = math.log(defect[n]) / n if defect[n] > 1 else 0.0
            max_rate = max(max_rate, (rate, n))

    target_count = 0
    incidence_count = 0
    primes = [p for p in primes_upto(limit) if p >= 7]
    for n in range(11, limit + 1):
        for p in primes:
            if p <= n // 2:
                continue
            if p > n:
                break
            incidence_count += 1
            target = values[n] % p == 0
            in_defect = defect[n] % p == 0
            assert in_defect == target, (n, p, target, defect[n] % p)
            if target:
                assert denominators[n - 1] % p == 0
                assert denominators[n] % p != 0
                assert denominators[n + 1] % p == 0
                target_count += 1

    anchor_count = section_anchor_audit(values)
    ratio_count = coefficient_ratio_audit()

    print(f"PASS: {ratio_count} hypergeometric coefficient ratios")
    print(f"PASS: {anchor_count} residue-section prime anchors")
    print(
        "PASS: "
        f"{comparison_count} exact carrier-to-Apery-gcd comparisons"
    )
    print(
        "PASS: "
        f"{incidence_count} top-half incidences, {target_count} exact holes"
    )
    print(f"max defect bits: {max_bits[0]} at n={max_bits[1]}")
    print(f"max log(defect)/n: {max_rate[0]:.12f} at n={max_rate[1]}")


if __name__ == "__main__":
    main()
