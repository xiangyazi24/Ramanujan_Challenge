#!/usr/bin/env python3
"""Independent small-scale verifier for ``long_bridge_incidence_scan.cpp``.

The production scanner counts from sorted occurrence lists.  This verifier
instead loops over every admissible pair and triple of indices.  It also
evaluates every relevant gap continuant pointwise and checks

    N_h(x) == 0  <=>  pi_p(x) == pi_p(x+h).

The default dyadic interval is deliberately small enough for the quadratic
reference algorithm.  Its exact aggregate is retained as a regression.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import hashlib
import json
import math


@dataclass
class Bin:
    lower: int
    upper: int
    all_short: int = 0
    four_gap_valid: int = 0
    selected_first: int = 0


def primes_in_dyadic_interval(scale: int) -> list[int]:
    sieve = bytearray(b"\x01") * (2 * scale + 1)
    sieve[:2] = b"\x00\x00"
    for prime in range(2, math.isqrt(2 * scale) + 1):
        if not sieve[prime]:
            continue
        start = prime * prime
        sieve[start : 2 * scale + 1 : prime] = b"\x00" * (
            (2 * scale - start) // prime + 1
        )
    return [
        value
        for value in range(max(7, scale + 1), 2 * scale + 1)
        if sieve[value]
    ]


def make_bins(scale: int) -> list[Bin]:
    height = math.isqrt(2 * scale)
    near_cutoff = max(2, int(height / math.sqrt(math.log(scale))))
    bins = [Bin(2, near_cutoff)]
    if near_cutoff < height:
        bins.append(Bin(near_cutoff + 1, height))
    lower = height + 1
    while lower <= 2 * scale:
        upper = min(2 * scale, 2 * lower - 1)
        bins.append(Bin(lower, upper))
        lower = upper + 1
    return bins


def solution_pair(prime: int) -> tuple[list[int], list[int]]:
    apery = [1, 5 % prime]
    companion = [0, 1]
    for index in range(1, prime - 1):
        coefficient = (34 * index**3 + 51 * index**2 + 27 * index + 5) % prime
        denominator_inverse = pow(index + 1, -3, prime)
        apery.append(
            (
                coefficient * apery[index]
                - index**3 * apery[index - 1]
            )
            * denominator_inverse
            % prime
        )
        companion.append(
            (
                coefficient * companion[index]
                - index**3 * companion[index - 1]
            )
            * denominator_inverse
            % prime
        )
    return apery, companion


def projective_states(
    prime: int, apery: list[int], companion: list[int]
) -> list[int]:
    states = []
    assert len(apery) == len(companion)
    for left, right in zip(apery, companion):
        assert left or right
        states.append(prime if right == 0 else left * pow(right, -1, prime) % prime)
    assert states == list(reversed(states))
    return states


def check_gap_continuants(prime: int, states: list[int]) -> None:
    """Evaluate N_h(x) numerically, without constructing any polynomials."""

    for start in range(prime):
        previous = 0  # N_0(start)
        current = 1  # N_1(start)
        for height in range(1, prime - start):
            assert (current == 0) == (states[start] == states[start + height])
            argument = start + height
            coefficient = (
                (2 * argument + 1)
                * (17 * argument**2 + 17 * argument + 5)
            ) % prime
            following = (
                coefficient * current - pow(argument, 6, prime) * previous
            ) % prime
            previous, current = current, following


def carrier_prefix(prime: int, apery: list[int], height: int) -> list[bool]:
    bad = [False] * (height + 1)
    previous, current = 0, 1
    for span in range(1, height + 1):
        following = (34 * current - previous) % prime
        previous, current = current, following
        bad[span] = bad[span - 1] or apery[span] == 0 or current == 0
    return bad


def locate_bin(bins: list[Bin], gap: int) -> Bin | None:
    for bridge_bin in bins:
        if bridge_bin.lower <= gap <= bridge_bin.upper:
            return bridge_bin
    return None


def scan_prime_direct(prime: int, template: list[Bin]) -> tuple[int, int, int, list[Bin]]:
    apery, companion = solution_pair(prime)
    states = projective_states(prime, apery, companion)
    check_gap_continuants(prime, states)
    height = math.isqrt(prime)
    carrier_bad = carrier_prefix(prime, apery, height)
    bins = [Bin(item.lower, item.upper) for item in template]
    short_pairs = 0
    four_chains = 0
    selected_chains = 0
    prior_occurrences: dict[int, list[int]] = {}

    for middle in range(prime):
        state = states[middle]
        prior = prior_occurrences.setdefault(state, [])
        short_starts = [
            start
            for start in prior
            if 2 <= middle - start <= height
        ]
        short_pairs += len(short_starts)

        four_gap_valid = False
        selected_first = False
        if len(prior) >= 3:
            x0, x1, x2 = prior[-3:]
            a, b, c = x1 - x0, x2 - x1, middle - x2
            span = middle - x0
            four_gap_valid = span <= height and min(a, b, c) >= 2
            selected_first = (
                four_gap_valid
                and not (a == b == c)
                and x0 + x1 != prime - 1
                and x1 + x2 != prime - 1
                and x2 + middle != prime - 1
                and not carrier_bad[span]
            )
        four_chains += int(four_gap_valid)
        selected_chains += int(selected_first)

        for later in range(middle + 2, prime):
            if states[later] != state:
                continue
            bridge_bin = locate_bin(bins, later - middle)
            assert bridge_bin is not None
            bridge_bin.all_short += len(short_starts)
            bridge_bin.four_gap_valid += int(four_gap_valid)
            bridge_bin.selected_first += int(selected_first)
        prior.append(middle)

    return short_pairs, four_chains, selected_chains, bins


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--x", type=int, default=70)
    args = parser.parse_args()
    if args.x < 7:
        raise SystemExit("X must be at least 7")

    template = make_bins(args.x)
    totals = {
        "short_pairs": 0,
        "four_chains": 0,
        "selected_chains": 0,
    }
    bins = [Bin(item.lower, item.upper) for item in template]
    primes = primes_in_dyadic_interval(args.x)
    for prime in primes:
        short_pairs, four_chains, selected_chains, prime_bins = scan_prime_direct(
            prime, template
        )
        totals["short_pairs"] += short_pairs
        totals["four_chains"] += four_chains
        totals["selected_chains"] += selected_chains
        assert len(bins) == len(prime_bins)
        for total_bin, prime_bin in zip(bins, prime_bins):
            total_bin.all_short += prime_bin.all_short
            total_bin.four_gap_valid += prime_bin.four_gap_valid
            total_bin.selected_first += prime_bin.selected_first

    result = {
        "X": args.x,
        "primes": len(primes),
        **totals,
        "bins": [vars(item) for item in bins],
    }
    encoded = json.dumps(result, sort_keys=True, separators=(",", ":")).encode()
    print(json.dumps(result, sort_keys=True))
    print(f"sha256={hashlib.sha256(encoded).hexdigest()}")

    if args.x == 70:
        assert result == {
            "X": 70,
            "primes": 15,
            "short_pairs": 199,
            "four_chains": 1,
            "selected_chains": 0,
            "bins": [
                {
                    "lower": 2,
                    "upper": 5,
                    "all_short": 17,
                    "four_gap_valid": 0,
                    "selected_first": 0,
                },
                {
                    "lower": 6,
                    "upper": 11,
                    "all_short": 21,
                    "four_gap_valid": 0,
                    "selected_first": 0,
                },
                {
                    "lower": 12,
                    "upper": 23,
                    "all_short": 41,
                    "four_gap_valid": 0,
                    "selected_first": 0,
                },
                {
                    "lower": 24,
                    "upper": 47,
                    "all_short": 73,
                    "four_gap_valid": 1,
                    "selected_first": 0,
                },
                {
                    "lower": 48,
                    "upper": 95,
                    "all_short": 60,
                    "four_gap_valid": 0,
                    "selected_first": 0,
                },
                {
                    "lower": 96,
                    "upper": 140,
                    "all_short": 10,
                    "four_gap_valid": 0,
                    "selected_first": 0,
                },
            ],
        }
        print("regression=PASS")


if __name__ == "__main__":
    main()
