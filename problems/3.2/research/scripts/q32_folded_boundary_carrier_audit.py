#!/usr/bin/env python3
"""Audit the folded-index and folded-slack carrier identities."""

from __future__ import annotations

import argparse
from math import isqrt


def primes_upto(limit: int) -> list[int]:
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[:2] = b"\x00\x00"
    for prime in range(2, isqrt(limit) + 1):
        if sieve[prime]:
            sieve[prime * prime : limit + 1 : prime] = (
                b"\x00" * ((limit - prime * prime) // prime + 1)
            )
    return [prime for prime in range(2, limit + 1) if sieve[prime]]


def apery_values(limit: int) -> list[int]:
    values = [1, 5]
    for index in range(1, limit):
        polynomial = 34 * index**3 + 51 * index**2 + 27 * index + 5
        numerator = (
            polynomial * values[index] - index**3 * values[index - 1]
        )
        denominator = (index + 1) ** 3
        assert numerator % denominator == 0
        values.append(numerator // denominator)
    return values[: limit + 1]


def audit(maximum_n: int) -> None:
    primes = primes_upto(maximum_n)
    apery = apery_values(maximum_n)

    target_checks = 0
    direct_checks = 0
    reflected_checks = 0
    slack_products = 0
    index_products = 0

    for n in range(10, maximum_n + 1):
        slack_radical: dict[int, int] = {}
        index_radical: dict[int, int] = {}

        for prime in primes:
            if not isqrt(n) < prime <= n:
                continue
            remainder = n % prime
            if apery[remainder] % prime:
                continue

            folded = min(remainder, prime - 1 - remainder)
            slack = prime - 2 * folded
            assert 1 <= slack <= prime
            assert slack % 2 == 1
            assert apery[folded] % prime == 0

            if remainder == folded:
                assert (n - folded) % prime == 0
                assert (2 * n + slack) % prime == 0
                direct_checks += 1
            else:
                assert remainder == prime - 1 - folded
                assert (n + folded + 1) % prime == 0
                assert (2 * n + 2 - slack) % prime == 0
                reflected_checks += 1

            slack_radical[slack] = slack_radical.get(slack, 1) * prime
            index_radical[folded] = index_radical.get(folded, 1) * prime
            target_checks += 1

        for slack, radical in slack_radical.items():
            carrier = (2 * n + slack) * (2 * n + 2 - slack)
            assert carrier % radical == 0
            slack_products += 1

        for folded, radical in index_radical.items():
            carrier = (n - folded) * (n + folded + 1)
            assert carrier % radical == 0
            index_products += 1

    print(
        {
            "maximum_n": maximum_n,
            "target_checks": target_checks,
            "direct_checks": direct_checks,
            "reflected_checks": reflected_checks,
            "slack_products": slack_products,
            "index_products": index_products,
            "failures": 0,
        }
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--maximum-n", type=int, default=800)
    audit(parser.parse_args().maximum_n)
