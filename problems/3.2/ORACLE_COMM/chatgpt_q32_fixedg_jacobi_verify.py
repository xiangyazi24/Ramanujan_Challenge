#!/usr/bin/env python3
"""Exact standard-library audit for the Q8377 fixed-g Apéry packet.

The script checks only finite algebra used in
``chatgpt_q32_fixedg_jacobi_report.md``.  It does not claim to verify the
conditional horizontal statement FGJ16.
"""

from __future__ import annotations

import argparse
import hashlib
from dataclasses import dataclass
from math import comb, gcd, prod


@dataclass(frozen=True)
class Label:
    branch: str
    a: int
    s: int
    p: int
    order: int


def prime_flags(limit: int) -> bytearray:
    flags = bytearray(b"\x01") * (limit + 1)
    if limit >= 0:
        flags[0] = 0
    if limit >= 1:
        flags[1] = 0
    q = 2
    while q * q <= limit:
        if flags[q]:
            flags[q * q : limit + 1 : q] = b"\x00" * (
                (limit - q * q) // q + 1
            )
        q += 1
    return flags


def divisors(n: int) -> list[int]:
    out: list[int] = []
    d = 1
    while d * d <= n:
        if n % d == 0:
            out.append(d)
            if d * d != n:
                out.append(n // d)
        d += 1
    return sorted(out)


def apery_values(limit: int) -> list[int]:
    if limit == 0:
        return [1]
    values = [1, 5]
    for n in range(1, limit):
        middle = 34 * n**3 + 51 * n**2 + 27 * n + 5
        numerator = middle * values[n] - n**3 * values[n - 1]
        denominator = (n + 1) ** 3
        value, remainder = divmod(numerator, denominator)
        assert remainder == 0, (n, numerator, denominator)
        values.append(value)
    return values[: limit + 1]


def valuation(value: int, prime: int) -> int:
    out = 0
    while value and value % prime == 0:
        value //= prime
        out += 1
    return out


def packet_labels(m: int, g: int, is_prime: bytearray) -> list[Label]:
    assert m >= 2 and g >= 1 and (m - 1) % g == 0
    e = (m - 1) // g
    exponent = m - 1
    labels: list[Label] = []

    # The midpoint 3a=E belongs to the direct branch.
    for a in range(1, e // 3 + 1):
        if 3 * a > e:
            continue
        s = g * a
        p = m - s
        if not is_prime[p]:
            continue
        order = (p - 1) // gcd(p - 1, exponent)
        labels.append(Label("D", a, s, p, order))

    # The reflected branch is strict at the midpoint.
    for a in range(1, (e - 1) // 3 + 1):
        if 3 * a >= e:
            continue
        numerator = m + 1 + g * a
        if numerator % 2:
            continue
        s = g * a
        p = numerator // 2
        if not is_prime[p]:
            continue
        order = (p - 1) // gcd(p - 1, exponent)
        labels.append(Label("R", a, s, p, order))

    return labels


def audit(maximum_m: int) -> None:
    assert maximum_m >= 143
    is_prime = prime_flags(maximum_m)
    primes = [p for p in range(2, maximum_m + 1) if is_prime[p]]
    b = apery_values(maximum_m)

    packet_checks = 0
    label_checks = 0
    lucas_checks = 0
    order_checks = 0
    product_checks = 0
    height_checks = 0
    selected_occurrences = 0
    p5_exceptions = 0
    digest_rows: list[str] = []

    for m in range(6, maximum_m + 1):
        # Exact height-slack carrier for the strict top strip.
        top_primes = [p for p in primes if p <= m and 2 * p > m + 1]
        top_product = prod(top_primes, start=1)
        central = comb(m, m // 2)
        assert central % top_product == 0, (m, top_primes)
        assert central < 2**m
        assert b[m] >= comb(2 * m, m) ** 2
        assert comb(2 * m, m) * (2 * m + 1) >= 4**m
        height_checks += 1

        for g in divisors(m - 1):
            e = (m - 1) // g
            exponent = m - 1
            labels = packet_labels(m, g, is_prime)

            # Joint injectivity is exact, not merely empirical.
            candidate_primes = [label.p for label in labels]
            assert len(candidate_primes) == len(set(candidate_primes)), (
                m,
                g,
                labels,
            )

            selected: list[Label] = []
            for label in labels:
                a, s, p = label.a, label.s, label.p
                assert p * p > m
                assert 0 < s < p

                if label.branch == "D":
                    assert 3 * a <= e
                    assert p == m - g * a == g * (e - a) + 1
                    assert 3 * p >= 2 * m + 1
                    assert p < m
                    assert 2 * s <= p - 1
                    assert (exponent - s) == p - 1
                    assert exponent % (p - 1) == s % (p - 1)
                    expected_order = (e - a) // gcd(e, a)
                else:
                    assert 3 * a < e
                    assert g * (e + a) % 2 == 0
                    assert 2 * p == m + 1 + g * a
                    assert 2 * p > m + 1
                    assert 3 * p < 2 * m + 1
                    assert 2 * s < p - 1
                    assert exponent + s == 2 * (p - 1)
                    assert exponent % (p - 1) == (-s) % (p - 1)
                    expected_order = (e + a) // gcd(e + a, 2 * e)

                assert label.order == expected_order
                assert label.order == (p - 1) // gcd(p - 1, exponent)
                order_checks += 1

                quotient, residue = divmod(m, p)
                assert quotient == 1
                if label.branch == "D":
                    assert residue == s
                else:
                    assert residue == p - 1 - s

                # Exact two-digit Lucas and reflection checks on integer data.
                assert b[m] % p == (b[quotient] * b[residue]) % p
                assert (b[p - 1 - s] - b[s]) % p == 0
                lucas_checks += 2

                if p >= 7:
                    assert (b[m] % p == 0) == (b[s] % p == 0), (
                        m,
                        g,
                        label,
                    )
                else:
                    p5_exceptions += 1

                if b[s] % p == 0 and p >= 7:
                    selected.append(label)
                    selected_occurrences += 1
                    # A selected recurrence zero has nonzero adjacent states.
                    assert s >= 1
                    assert b[s - 1] % p != 0
                    assert b[s + 1] % p != 0

                digest_rows.append(
                    f"{m}:{g}:{label.branch}:{a}:{s}:{p}:"
                    f"{label.order}:{int(b[s] % p == 0)}"
                )
                label_checks += 1

            selected_product = prod((label.p for label in selected), start=1)
            assert b[m] % selected_product == 0

            selected_by_a: dict[int, list[int]] = {}
            for label in selected:
                selected_by_a.setdefault(label.a, []).append(label.p)

            lower_product = 1
            for a, attached_primes in selected_by_a.items():
                attached_product = prod(attached_primes, start=1)
                s = g * a
                assert b[s] % attached_product == 0
                lower_product *= b[s]
            assert lower_product % selected_product == 0
            product_checks += 1

            packet_checks += 1

    # Actual fixed-g simple-primary and modular-coordinate regressions.
    assert b[5] == 819005
    assert b[5] == 11 * 74455
    assert valuation(b[5], 11) == 1

    labels_16_5 = packet_labels(16, 5, is_prime)
    assert Label("D", 1, 5, 11, 1) in labels_16_5
    assert not any(label.branch == "R" and label.a == 1 for label in labels_16_5)
    assert b[16] % 11 == 0
    assert 16 * 17 - 5 * 6 == 22 * 11

    # The coefficient zero b_5=0 mod 11 is not an evaluation zero of A_11.
    x = (5 * 6) % 11
    hasse_value = sum((b[r] % 11) * pow(x, r, 11) for r in range(11)) % 11
    assert x == 8
    assert hasse_value == 5

    # The same lower zero gives real direct and reflected hits in two packets.
    labels_75_2 = packet_labels(75, 2, is_prime)
    assert any(
        label.branch == "D" and label.a == 1 and label.p == 73
        for label in labels_75_2
    )
    assert b[2] == 73 and b[75] % 73 == 0

    labels_143_2 = packet_labels(143, 2, is_prime)
    assert any(
        label.branch == "R" and label.a == 1 and label.p == 73
        for label in labels_143_2
    )
    assert b[143] % 73 == 0

    # E=37 is prime: all admissible orders are genuinely large.
    for label in labels_75_2:
        if label.branch == "D":
            assert label.order == 37 - label.a
            assert 3 * label.order > 2 * 37
        else:
            assert 2 * label.order >= 37 + label.a
            assert 2 * label.order > 37

    digest = hashlib.sha256("\n".join(digest_rows).encode()).hexdigest()
    print("CHATGPT_Q32_FIXEDG_JACOBI_VERIFY=PASS")
    print("MAXIMUM_M", maximum_m)
    print("PACKET_CHECKS", packet_checks)
    print("LABEL_CHECKS", label_checks)
    print("LUCAS_REFLECTION_CHECKS", lucas_checks)
    print("ORDER_CHECKS", order_checks)
    print("PRODUCT_CHECKS", product_checks)
    print("HEIGHT_SLACK_CHECKS", height_checks)
    print("SELECTED_OCCURRENCES", selected_occurrences)
    print("P5_FINITE_EXCEPTIONS", p5_exceptions)
    print("DIGEST", digest)
    print("FGJ16", "CONDITIONAL_NOT_VERIFIED")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--maximum-m", type=int, default=220)
    arguments = parser.parse_args()
    audit(arguments.maximum_m)
