#!/usr/bin/env python3
"""Index-first exact search for short Apéry zero quadruples.

Prime-first scans bound ``p`` and iterate the recurrence modulo every prime.
This complementary scan bounds the integer indices but not the prime.  It
computes exact Apéry numbers, factors every pair gcd in the requested index
window, and then reconstructs the complete zero prefix for every prime factor.
Consequently the absence of a reported short off-center quadruple certifies
all primes, including primes above the prime-first cutoff, provided the four
indices and their span lie in the requested box.

The default box is deliberately modest; larger values are useful as a long
background computation.
"""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from itertools import combinations
from math import comb, gcd, isqrt


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-index", type=int, default=3_000)
    parser.add_argument("--max-span", type=int, default=300)
    parser.add_argument("--progress-step", type=int, default=500)
    return parser.parse_args()


def apery_numbers(max_index: int) -> list[int]:
    assert max_index >= 1
    values = [1, 5]
    for index in range(1, max_index):
        numerator = (
            (34 * index**3 + 51 * index**2 + 27 * index + 5)
            * values[index]
            - index**3 * values[index - 1]
        )
        denominator = (index + 1) ** 3
        assert numerator % denominator == 0
        values.append(numerator // denominator)
    return values


def primes_up_to(limit: int) -> list[int]:
    if limit < 2:
        return []
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[:2] = b"\x00\x00"
    for prime in range(2, isqrt(limit) + 1):
        if sieve[prime]:
            start = prime * prime
            sieve[start : limit + 1 : prime] = b"\x00" * (
                (limit - start) // prime + 1
            )
    return [value for value in range(2, limit + 1) if sieve[value]]


def factor_integer(value: int, trial_primes: list[int]) -> tuple[int, ...]:
    """Return distinct prime factors, with an exact reconstruction audit."""

    original = value
    factors = []
    exponents = []
    for prime in trial_primes:
        if prime * prime > value:
            break
        if value % prime:
            continue
        exponent = 0
        while value % prime == 0:
            value //= prime
            exponent += 1
        factors.append(prime)
        exponents.append(exponent)
    if value > 1:
        factors.append(value)
        exponents.append(1)
    reconstructed = 1
    for prime, exponent in zip(factors, exponents):
        reconstructed *= prime**exponent
    assert reconstructed == original
    return tuple(factors)


def main() -> None:
    args = parse_args()
    assert args.max_index >= 4
    assert 3 <= args.max_span <= args.max_index
    values = apery_numbers(args.max_index)
    assert values[:5] == [1, 5, 73, 1445, 33001]

    gcd_values: set[int] = set()
    nontrivial_pairs = 0
    total_pairs = 0
    for left in range(args.max_index + 1):
        right_stop = min(args.max_index, left + args.max_span)
        for right in range(left + 2, right_stop + 1):
            total_pairs += 1
            common = gcd(values[left], values[right])
            if common > 1:
                nontrivial_pairs += 1
                gcd_values.add(common)
        if args.progress_step and left and left % args.progress_step == 0:
            print(
                "progress "
                f"left={left} pairs={total_pairs} "
                f"unique_nontrivial_gcds={len(gcd_values)}",
                flush=True,
            )

    maximum_gcd = max(gcd_values, default=1)
    trial_primes = primes_up_to(isqrt(maximum_gcd))
    candidate_primes: set[int] = set()
    for common in gcd_values:
        candidate_primes.update(factor_integer(common, trial_primes))
    candidate_primes.difference_update({2, 3, 5})

    active = {}
    all_windows = []
    off_center = []
    short_off_center = []
    four_zero_subsets = 0
    block_four_subsets = 0
    block_excess = 0
    for prime in sorted(candidate_primes):
        top = min(args.max_index, prime - 1)
        zeros = tuple(
            index for index in range(top + 1) if values[index] % prime == 0
        )
        if len(zeros) >= 2 and any(
            right - left <= args.max_span
            for left, right in zip(zeros, zeros[1:])
        ):
            active[prime] = zeros
        four_zero_subsets += sum(
            chain[-1] - chain[0] <= args.max_span
            for chain in combinations(zeros, 4)
        )
        block_counts = Counter(
            index // (args.max_span + 1) for index in zeros
        )
        block_four_subsets += sum(
            comb(count, 4) for count in block_counts.values()
            if count >= 4
        )
        block_excess += sum(
            max(count - 3, 0) for count in block_counts.values()
        )
        for offset in range(len(zeros) - 3):
            chain = zeros[offset : offset + 4]
            span = chain[-1] - chain[0]
            if span > args.max_span:
                continue
            centered = tuple(
                index
                for index in range(3)
                if chain[index] + chain[index + 1] == prime - 1
            )
            row = (prime, chain, span, centered)
            all_windows.append(row)
            if not centered:
                off_center.append(row)
                if span * span <= prime:
                    short_off_center.append(row)

    # Regression anchor: whenever the requested box contains it, the short
    # p=1049 windows are all explained by their unique centered adjacent pair.
    if args.max_index >= 1048 and args.max_span >= 36:
        assert active.get(1049) == (494, 504, 508, 540, 544, 554)
        assert any(row[0] == 1049 for row in all_windows)
        assert all(row[3] for row in all_windows if row[0] == 1049)

    closest = None
    if off_center:
        closest = min(
            off_center,
            key=lambda row: (
                Fraction(row[2] ** 2, row[0]), row[0], row[1]
            ),
        )

    print(f"MAX_INDEX={args.max_index}")
    print(f"MAX_SPAN={args.max_span}")
    print(f"TOTAL_PAIRS={total_pairs}")
    print(f"NONTRIVIAL_PAIR_GCDS={nontrivial_pairs}")
    print(f"UNIQUE_NONTRIVIAL_GCDS={len(gcd_values)}")
    print(f"MAXIMUM_GCD_BITS={maximum_gcd.bit_length()}")
    print(f"CANDIDATE_PRIMES={len(candidate_primes)}")
    print(f"ACTIVE_PAIR_PRIMES={len(active)}")
    print(
        "ACTIVE_ZERO_COUNT_DISTRIBUTION="
        f"{sorted(Counter(map(len, active.values())).items())}"
    )
    print(f"CONSECUTIVE_QUADRUPLES={len(all_windows)}")
    print(f"FOUR_ZERO_SUBSETS={four_zero_subsets}")
    print(f"BLOCK_FOUR_SUBSETS={block_four_subsets}")
    print(f"BLOCK_EXCESS={block_excess}")
    print(f"OFF_CENTER_QUADRUPLES={len(off_center)}")
    print(f"SHORT_OFF_CENTER_QUADRUPLES={len(short_off_center)}")
    if closest is not None:
        ratio = Fraction(closest[2] ** 2, closest[0])
        print(f"CLOSEST_OFF_CENTER={closest}")
        print(
            "CLOSEST_SPAN_SQUARED_OVER_P="
            f"{ratio.numerator}/{ratio.denominator}"
        )
    if short_off_center:
        for row in short_off_center[:20]:
            print(f"SHORT_COUNTEREXAMPLE={row}")
        raise AssertionError("short off-center actual Apéry quadruple found")
    print("ACTUAL_QUADRUPLE_INDEX_SCAN PASS")


if __name__ == "__main__":
    main()
