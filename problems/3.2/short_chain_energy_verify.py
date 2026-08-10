#!/usr/bin/env python3
"""Audit the short-chain edge charge and its level-set energy bound."""

from __future__ import annotations

from collections import defaultdict
import hashlib
import json
import math


LIMIT = 499
EXPECTED_SHA256 = "0796e0f34777aba4a631767076dbedb71503eef4422e19553ad09e702be6566c"


def primes_up_to(limit: int) -> list[int]:
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[:2] = b"\x00\x00"
    for prime in range(2, math.isqrt(limit) + 1):
        if sieve[prime]:
            sieve[prime * prime :: prime] = b"\x00" * (
                (limit - prime * prime) // prime + 1
            )
    return [prime for prime in range(7, limit + 1) if sieve[prime]]


def coefficient(index: int, prime: int) -> int:
    return (
        34 * index**3 + 51 * index**2 + 27 * index + 5
    ) % prime


def solution(
    prime: int, initial_zero: int, initial_one: int
) -> list[int]:
    values = [initial_zero % prime, initial_one % prime]
    for index in range(1, prime - 1):
        values.append(
            (
                coefficient(index, prime) * values[index]
                - index**3 * values[index - 1]
            )
            * pow(index + 1, -3, prime)
            % prime
        )
    return values


def projective_states(prime: int) -> list[int]:
    apery = solution(prime, 1, 5)
    companion = solution(prime, 0, 1)
    states = [
        prime if right == 0 else left * pow(right, -1, prime) % prime
        for left, right in zip(apery, companion)
    ]
    assert states == list(reversed(states))
    return states


def audit_prime(prime: int) -> list[int]:
    height = math.isqrt(prime)
    states = projective_states(prime)
    fibers: dict[int, list[int]] = defaultdict(list)
    for index, state in enumerate(states):
        fibers[state].append(index)

    global_consecutive = [0] * (height + 1)
    edge_counts: list[int] = []
    total_windows = 0
    edge_energy = 0
    window_energy = 0

    for occurrences in fibers.values():
        short_edges: list[tuple[int, int]] = []
        for edge_index, (left, right) in enumerate(
            zip(occurrences, occurrences[1:])
        ):
            gap = right - left
            if gap <= height:
                short_edges.append((edge_index, gap))
                global_consecutive[gap] += 1

        edge_count = len(short_edges)
        edge_counts.append(edge_count)
        assert sum(gap for _, gap in short_edges) <= prime - 1

        incidence = defaultdict(int)
        window_count = 0
        for start in range(max(0, len(occurrences) - 3)):
            if occurrences[start + 3] - occurrences[start] > height:
                continue
            window_count += 1
            for edge_index in range(start, start + 3):
                gap = occurrences[edge_index + 1] - occurrences[edge_index]
                assert gap <= height
                incidence[edge_index] += 1

        assert sum(incidence.values()) == 3 * window_count
        assert all(multiplicity <= 3 for multiplicity in incidence.values())
        assert 3 * window_count <= 3 * edge_count
        total_windows += window_count
        edge_energy += edge_count * edge_count
        window_energy += window_count * window_count

    for gap in range(1, height + 1):
        all_pairs = sum(
            states[start] == states[start + gap]
            for start in range(prime - gap)
        )
        assert global_consecutive[gap] <= all_pairs
        assert all_pairs <= 3 * (gap - 1)

    assert sum(global_consecutive) <= 3 * height * (height - 1) // 2

    level = 0
    while 4 * (1 << level) <= height:
        scale = 1 << level
        population = sum(
            edge_count * height >= 4 * prime * scale
            for edge_count in edge_counts
        )
        assert 16 * prime * scale**3 * population <= 3 * height**3
        level += 1

    assert window_energy <= edge_energy
    assert edge_energy <= 30 * prime * height
    return [
        prime,
        height,
        total_windows,
        edge_energy,
        window_energy,
        max(edge_counts, default=0),
    ]


def main() -> None:
    rows = [audit_prime(prime) for prime in primes_up_to(LIMIT)]
    encoded = json.dumps(rows, separators=(",", ":")).encode()
    digest = hashlib.sha256(encoded).hexdigest()
    assert digest == EXPECTED_SHA256
    max_edge_ratio = max(row[3] / (row[0] * row[1]) for row in rows)
    max_window_ratio = max(row[4] / (row[0] * row[1]) for row in rows)
    print(
        "SHORT_CHAIN_ENERGY_VERIFY"
        f" primes={len(rows)} total_windows={sum(row[2] for row in rows)}"
        f" max_edge_ratio={max_edge_ratio:.9f}"
        f" max_window_ratio={max_window_ratio:.9f}"
    )
    print(f"sha256={digest}")


if __name__ == "__main__":
    main()
