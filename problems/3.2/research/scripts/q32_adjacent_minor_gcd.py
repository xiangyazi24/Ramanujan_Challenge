#!/usr/bin/env python3
"""Audit the long-prefix adjacent-minor gcd for Apéry top targets.

For

    t_k(n) = binom(n,k)^2 binom(n+k,k)^2,
    S_k(n) = sum_{j=0}^k t_j(n),

write the reduced ratio t_(k+1)/t_k = a_k/b_k and define

    C_(n,k) = b_k S_k - a_k S_(k-1).

The normalized Hankel identity shows that C_(n,k) is a positive integer
divisible by every top-half target prime.  This script checks that fact,
computes the gcd of the full adjacent-minor family, and compares it with
the older low/upper truncation gcd.  It supplies finite evidence only.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from math import gcd, isqrt, log, prod


@dataclass(frozen=True)
class Record:
    index: int
    target_radical: int
    adjacent_gcd: int
    adjacent_top_part: int
    truncation_gcd: int
    combined_gcd: int


def primes_at_most(limit: int) -> list[int]:
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[:2] = b"\x00\x00"
    for prime in range(2, isqrt(limit) + 1):
        if sieve[prime]:
            sieve[prime * prime : limit + 1 : prime] = b"\x00" * (
                (limit - prime * prime) // prime + 1
            )
    return [candidate for candidate in range(2, limit + 1) if sieve[candidate]]


def apery_numbers(limit: int) -> list[int]:
    values = [1, 5]
    for index in range(1, limit):
        numerator = (
            (34 * index**3 + 51 * index**2 + 27 * index + 5)
            * values[index]
            - index**3 * values[index - 1]
        )
        denominator = (index + 1) ** 3
        assert numerator % denominator == 0
        values.append(numerator // denominator)
    return values


def summands_and_prefixes(index: int) -> tuple[list[int], list[int]]:
    summands = [1]
    prefixes = [1]
    for cutoff in range(index):
        numerator = (
            summands[-1]
            * (index - cutoff) ** 2
            * (index + cutoff + 1) ** 2
        )
        denominator = (cutoff + 1) ** 4
        assert numerator % denominator == 0
        summands.append(numerator // denominator)
        prefixes.append(prefixes[-1] + summands[-1])
    return summands, prefixes


def adjacent_carrier(
    index: int,
    cutoff: int,
    summands: list[int],
    prefixes: list[int],
) -> int:
    raw_a = (index - cutoff) ** 2 * (index + cutoff + 1) ** 2
    raw_b = (cutoff + 1) ** 4
    common = gcd(raw_a, raw_b)
    coefficient_a = raw_a // common
    coefficient_b = raw_b // common
    carrier = (
        coefficient_b * prefixes[cutoff]
        - coefficient_a * prefixes[cutoff - 1]
    )

    term_gcd = gcd(summands[cutoff], summands[cutoff + 1])
    hankel = (
        prefixes[cutoff - 1] * prefixes[cutoff + 1]
        - prefixes[cutoff] ** 2
    )
    assert -hankel == term_gcd * carrier
    assert carrier > 0
    return carrier


def audit_index(
    index: int,
    primes: list[int],
    apery: list[int],
) -> Record | None:
    height = index // 3
    midpoint = index // 2
    if height + 1 >= midpoint:
        return None

    summands, prefixes = summands_and_prefixes(index)
    candidate_primes = [
        prime for prime in primes if index // 2 < prime <= index
    ]
    target_primes = [
        prime for prime in candidate_primes if apery[index] % prime == 0
    ]
    target_radical = prod(target_primes)

    adjacent_gcd = 0
    for cutoff in range(height + 1, midpoint):
        carrier = adjacent_carrier(index, cutoff, summands, prefixes)
        assert carrier % target_radical == 0
        adjacent_gcd = gcd(adjacent_gcd, carrier)

    adjacent_top_part = prod(
        prime ** valuation(adjacent_gcd, prime) for prime in candidate_primes
    )
    lower_candidates = prod(
        prime for prime in candidate_primes if prime <= 2 * height + 1
    )
    smith_cap = top_primorial = prod(candidate_primes)
    smith_cap = smith_cap**2 * lower_candidates**2
    if index >= 60:
        assert adjacent_top_part == gcd(prefixes[height], smith_cap)

    sharp_height = (index - 1) // 3
    low_prefix = prefixes[sharp_height]
    upper_suffix = sum(summands[(index + 1) // 2 :])
    truncation_gcd = gcd(low_prefix, upper_suffix)
    combined_gcd = gcd(adjacent_gcd, truncation_gcd)

    assert (upper_suffix - 4 * low_prefix) % top_primorial**2 == 0
    assert (prefixes[-1] - 5 * low_prefix) % top_primorial**2 == 0
    assert truncation_gcd % target_radical == 0
    assert combined_gcd % target_radical == 0

    return Record(
        index=index,
        target_radical=target_radical,
        adjacent_gcd=adjacent_gcd,
        adjacent_top_part=adjacent_top_part,
        truncation_gcd=truncation_gcd,
        combined_gcd=combined_gcd,
    )


def rate(value: int, index: int) -> float:
    return log(value) / index if value > 1 else 0.0


def valuation(value: int, prime: int) -> int:
    exponent = 0
    while value % prime == 0:
        value //= prime
        exponent += 1
    return exponent


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--limit", type=int, default=500)
    parser.add_argument("--block-start", type=int, default=100)
    args = parser.parse_args()
    if args.limit < 16:
        raise SystemExit("--limit must be at least 16")

    primes = primes_at_most(args.limit)
    apery = apery_numbers(args.limit)
    records = [
        record
        for index in range(16, args.limit + 1)
        if (record := audit_index(index, primes, apery)) is not None
    ]

    exact_adjacent = sum(
        record.adjacent_gcd == record.target_radical for record in records
    )
    exact_combined = sum(
        record.combined_gcd == record.target_radical for record in records
    )
    print(f"records={len(records)}")
    print(f"adjacent_exact_target_radical={exact_adjacent}")
    print(f"combined_exact_target_radical={exact_combined}")

    lower = args.block_start
    while lower < args.limit:
        upper = min(2 * lower, args.limit)
        block = [
            record for record in records if lower < record.index <= upper
        ]
        if block:
            adjacent_winner = max(
                block,
                key=lambda record: rate(record.adjacent_gcd, record.index),
            )
            combined_winner = max(
                block,
                key=lambda record: rate(record.combined_gcd, record.index),
            )
            print(
                f"({lower},{upper}] "
                f"adjacent_max={rate(adjacent_winner.adjacent_gcd, adjacent_winner.index):.9f} "
                f"at_n={adjacent_winner.index} "
                f"adjacent_extra={adjacent_winner.adjacent_gcd // adjacent_winner.target_radical} "
                f"combined_max={rate(combined_winner.combined_gcd, combined_winner.index):.9f} "
                f"combined_at_n={combined_winner.index} "
                f"combined_extra={combined_winner.combined_gcd // combined_winner.target_radical}"
            )
        lower *= 2


if __name__ == "__main__":
    main()
