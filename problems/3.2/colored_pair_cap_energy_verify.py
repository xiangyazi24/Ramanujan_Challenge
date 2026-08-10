#!/usr/bin/env python3
"""Verify the colored row-shift model that saturates the pH energy bound."""

from __future__ import annotations

from collections import defaultdict
import hashlib
import json
import math


Q_LIMIT = 61
EXPECTED_SHA256 = "707f006bbafe1da79d69fb859f76b78b885d238c9bacda97e47043981600e39e"


def is_prime(value: int) -> bool:
    if value < 2:
        return False
    if value % 2 == 0:
        return value == 2
    divisor = 3
    while divisor * divisor <= value:
        if value % divisor == 0:
            return False
        divisor += 2
    return True


def next_prime(value: int) -> int:
    candidate = value + 1
    if candidate % 2 == 0:
        candidate += 1
    while not is_prime(candidate):
        candidate += 2
    return candidate


def primes_between(lower: int, upper: int) -> list[int]:
    return [value for value in range(lower, upper + 1) if is_prime(value)]


def audit_parameter(row_prime: int) -> list[int]:
    rows = row_prime - 2
    height = 6 * row_prime
    packet_length = rows * row_prime
    prime = next_prime(40 * row_prime**2)
    orbit_height = math.isqrt(prime)
    assert prime < 80 * row_prime**2
    assert height * height < prime
    assert prime - 1 > 4 * (1 + packet_length)

    shifts = [0]
    for row in range(rows - 1):
        shifts.append((shifts[-1] + row + 3) % row_prime)

    colors = [-1] * prime
    left_by_color: list[list[int]] = [[] for _ in range(row_prime)]
    for row in range(rows):
        occupied_columns = set()
        for color in range(row_prime):
            column = (color + shifts[row]) % row_prime
            occupied_columns.add(column)
            left = 1 + row * row_prime + column
            right = prime - 1 - left
            assert colors[left] == colors[right] == -1
            colors[left] = colors[right] = color
            left_by_color[color].append(left)
        assert len(occupied_columns) == row_prime

    next_color = row_prime
    for index in range((prime + 1) // 2):
        reflected = prime - 1 - index
        if colors[index] != -1:
            assert colors[index] == colors[reflected]
            continue
        colors[index] = colors[reflected] = next_color
        next_color += 1

    assert all(color >= 0 for color in colors)
    assert colors == list(reversed(colors))
    fibers: dict[int, list[int]] = defaultdict(list)
    for index, color in enumerate(colors):
        fibers[color].append(index)
    fiber_energy = sum(len(occurrences) ** 2 for occurrences in fibers.values())
    assert fiber_energy**2 <= prime**3

    pair_counts = [0] * prime
    for occurrences in fibers.values():
        for left_index, left in enumerate(occurrences):
            for right in occurrences[left_index + 1 :]:
                pair_counts[right - left] += 1
    for gap in range(1, prime):
        assert pair_counts[gap] <= 3 * (gap - 1)

    quotient_windows = row_prime - 5
    separated_energy = 0
    conditioned_incidence = 0
    far_cutoff = int(orbit_height / math.sqrt(math.log(prime)))
    for color, left_occurrences in enumerate(left_by_color):
        assert len(left_occurrences) == rows
        assert fibers[color][:rows] == left_occurrences
        windows = []
        for start in range(quotient_windows):
            window = tuple(left_occurrences[start : start + 4])
            gaps = tuple(
                window[index + 1] - window[index] for index in range(3)
            )
            assert min(gaps) >= 3
            assert len(set(gaps)) == 3
            assert window[-1] - window[0] < height
            assert all(
                window[index] + window[index + 1] != prime - 1
                for index in range(3)
            )
            reflected = tuple(prime - 1 - item for item in reversed(window))
            assert all(colors[item] == color for item in window + reflected)
            windows.append((window[0], window[-1], reflected[0], reflected[-1]))

        for first_start in range(quotient_windows):
            first = tuple(left_occurrences[first_start : first_start + 4])
            reflected_first = {
                prime - 1 - item for item in first
            }
            support = set(first) | reflected_first
            for second_start in range(first_start + 1, quotient_windows):
                second = tuple(
                    left_occurrences[second_start : second_start + 4]
                )
                external = second[0]
                if external <= first[-1] + far_cutoff:
                    continue
                assert external not in support
                assert tuple(
                    fibers[color][second_start : second_start + 4]
                ) == second
                assert second[-1] - second[0] < height <= orbit_height
                assert all(colors[item] == color for item in second)
                conditioned_incidence += 1

        color_energy = 0
        for first_index, first in enumerate(windows):
            for second_index, second in enumerate(windows):
                if first_index == second_index:
                    continue
                overlaps = not (
                    first[1] < second[0] or second[1] < first[0]
                ) or not (
                    first[3] < second[2] or second[3] < first[2]
                )
                if not overlaps:
                    color_energy += 1
        assert color_energy == (row_prime - 8) * (row_prime - 9)
        separated_energy += color_energy

    assert separated_energy == row_prime * (row_prime - 8) * (row_prime - 9)
    assert 1920 * separated_energy >= prime * height
    assert 3840 * conditioned_incidence >= prime * height
    return [
        row_prime,
        prime,
        height,
        max(pair_counts),
        separated_energy,
        conditioned_incidence,
        fiber_energy,
        next_color,
    ]


def main() -> None:
    rows = [
        audit_parameter(row_prime)
        for row_prime in primes_between(19, Q_LIMIT)
    ]
    encoded = json.dumps(rows, separators=(",", ":")).encode()
    digest = hashlib.sha256(encoded).hexdigest()
    assert digest == EXPECTED_SHA256
    min_energy_ratio = min(row[4] / (row[1] * row[2]) for row in rows)
    min_conditioned_ratio = min(row[5] / (row[1] * row[2]) for row in rows)
    print(
        "COLORED_PAIR_CAP_ENERGY_VERIFY"
        f" parameters={len(rows)} min_energy_ratio={min_energy_ratio:.9f}"
        f" min_conditioned_ratio={min_conditioned_ratio:.9f}"
    )
    print(f"sha256={digest}")


if __name__ == "__main__":
    main()
