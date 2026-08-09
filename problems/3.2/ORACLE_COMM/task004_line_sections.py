#!/usr/bin/env python3
"""Line-section statistics for the Apéry zero-pair incidence cloud.

For each p, the cloud contains the point (r,h) for every ordered pair
(r,r') in Z_p^2, with h=r'-r mod p represented by 1,...,p.  Since 0 is
not an Apéry zero, every cloud point has r != 0 and hence lies on exactly
one nonvertical line h=c r through the origin.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable, Sequence


HERE = Path(__file__).resolve().parent


def primes_upto(limit: int) -> list[int]:
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[0:2] = b"\x00\x00"
    for prime in range(2, math.isqrt(limit) + 1):
        if sieve[prime]:
            start = prime * prime
            sieve[start : limit + 1 : prime] = b"\x00" * (
                (limit - start) // prime + 1
            )
    return [n for n in range(2, limit + 1) if sieve[n]]


def apery_zeros(prime: int) -> tuple[int, ...]:
    """Zeros via beta_n=(n!)^3 b_n and its division-free recurrence."""

    previous, current = 1 % prime, 5 % prime
    zeros: list[int] = []
    if previous == 0:
        zeros.append(0)
    if current == 0:
        zeros.append(1)
    for index in range(1, prime - 1):
        polynomial = (
            34 * index**3 + 51 * index**2 + 27 * index + 5
        ) % prime
        following = (
            polynomial * current - pow(index, 6, prime) * previous
        ) % prime
        previous, current = current, following
        if current == 0:
            zeros.append(index + 1)
    return tuple(zeros)


def apery_zeros_divided(prime: int) -> tuple[int, ...]:
    """Independent check using the original Apéry recurrence."""

    previous, current = 1 % prime, 5 % prime
    zeros: list[int] = []
    if previous == 0:
        zeros.append(0)
    if current == 0:
        zeros.append(1)
    for index in range(1, prime - 1):
        polynomial = (
            34 * index**3 + 51 * index**2 + 27 * index + 5
        ) % prime
        numerator = (
            polynomial * current - pow(index, 3, prime) * previous
        ) % prime
        denominator = pow(index + 1, 3, prime)
        following = numerator * pow(denominator, -1, prime) % prime
        previous, current = current, following
        if current == 0:
            zeros.append(index + 1)
    return tuple(zeros)


def line_sections(prime: int, zeros: Sequence[int]) -> list[int]:
    """Counts on h=c r for c in F_p; the vertical line has count zero."""

    counts = [0] * prime
    for r in zeros:
        if r == 0:
            raise AssertionError("b_0=1, so r=0 cannot occur")
        inverse = pow(r, -1, prime)
        for rp in zeros:
            h = (rp - r) % prime
            slope = h * inverse % prime
            counts[slope] += 1
    return counts


def row(prime: int, zeros: Sequence[int]) -> dict[str, object]:
    z = len(zeros)
    counts = line_sections(prime, zeros)
    cloud_size = sum(counts)
    if cloud_size != z * z or counts[0] != z:
        raise AssertionError((prime, z, cloud_size, counts[0]))

    energy = sum(value * value for value in counts)
    nonzero_energy = sum(value * value for value in counts[1:])
    off_cloud_size = z * (z - 1)
    random_raw = cloud_size + cloud_size * (cloud_size - 1) / prime
    # For an off-diagonal pair, r'/r is neither 0 nor 1.  Hence its slope
    # c=r'/r-1 ranges over p-2 values: c=0 is diagonal and c=-1 would force
    # r'=0.  The occupancy benchmark must use these p-2 admissible lines.
    random_nonzero = (
        off_cloud_size
        + off_cloud_size * (off_cloud_size - 1) / (prime - 2)
    )
    random_forced_diagonal = z * z + random_nonzero
    distribution = Counter(counts)
    pencil_distribution = distribution.copy()
    pencil_distribution[0] += 1  # the vertical line r=0
    nonzero_distribution = Counter(counts[1:])
    return {
        "p": prime,
        "Z": z,
        "cloud_size": cloud_size,
        "pencil_line_count": prime + 1,
        "vertical_line_count": 0,
        "slope_zero_count": counts[0],
        "max_line": max(counts),
        "max_nonzero_line": max(counts[1:], default=0),
        "line_energy": energy,
        "nonzero_line_energy": nonzero_energy,
        "random_raw_energy": f"{random_raw:.12g}",
        "raw_energy_ratio": f"{energy / random_raw:.12g}",
        "forced_diagonal_prediction": f"{random_forced_diagonal:.12g}",
        "forced_diagonal_ratio": f"{energy / random_forced_diagonal:.12g}",
        "random_nonzero_energy": f"{random_nonzero:.12g}",
        "nonzero_energy_ratio": f"{nonzero_energy / random_nonzero:.12g}",
        "line_count_distribution": json.dumps(dict(sorted(distribution.items()))),
        "pencil_line_count_distribution": json.dumps(
            dict(sorted(pencil_distribution.items()))
        ),
        "nonzero_line_distribution": json.dumps(
            dict(sorted(nonzero_distribution.items()))
        ),
    }


def aggregate(rows: Sequence[dict[str, object]]) -> dict[str, object]:
    total_energy = sum(int(item["line_energy"]) for item in rows)
    total_forced = sum(float(item["forced_diagonal_prediction"]) for item in rows)
    total_nonzero = sum(int(item["nonzero_line_energy"]) for item in rows)
    total_nonzero_prediction = sum(float(item["random_nonzero_energy"]) for item in rows)
    by_z: dict[int, list[dict[str, object]]] = defaultdict(list)
    finite_slope_distribution: Counter[int] = Counter()
    pencil_distribution: Counter[int] = Counter()
    for item in rows:
        by_z[int(item["Z"])].append(item)
        finite_slope_distribution.update(
            {int(k): int(v) for k, v in json.loads(
                str(item["line_count_distribution"])
            ).items()}
        )
        pencil_distribution.update(
            {int(k): int(v) for k, v in json.loads(
                str(item["pencil_line_count_distribution"])
            ).items()}
        )
    groups = []
    for z, group in sorted(by_z.items()):
        observed = sum(int(item["line_energy"]) for item in group)
        predicted = sum(float(item["forced_diagonal_prediction"]) for item in group)
        observed_nonzero = sum(int(item["nonzero_line_energy"]) for item in group)
        predicted_nonzero = sum(float(item["random_nonzero_energy"]) for item in group)
        groups.append(
            {
                "Z": z,
                "primes": len(group),
                "energy": observed,
                "forced_prediction": predicted,
                "ratio": observed / predicted,
                "nonzero_energy": observed_nonzero,
                "nonzero_prediction": predicted_nonzero,
                "nonzero_ratio": observed_nonzero / predicted_nonzero,
                "max_nonzero_line": max(int(item["max_nonzero_line"]) for item in group),
            }
        )
    return {
        "active_primes": len(rows),
        "total_line_energy": total_energy,
        "total_forced_diagonal_prediction": total_forced,
        "aggregate_forced_diagonal_ratio": total_energy / total_forced,
        "total_nonzero_line_energy": total_nonzero,
        "total_nonzero_prediction": total_nonzero_prediction,
        "aggregate_nonzero_ratio": total_nonzero / total_nonzero_prediction,
        "max_nonzero_line": max(int(item["max_nonzero_line"]) for item in rows),
        "finite_slope_distribution": dict(sorted(finite_slope_distribution.items())),
        "full_pencil_distribution": dict(sorted(pencil_distribution.items())),
        "by_Z": groups,
    }


def write_csv(path: Path, rows: Iterable[dict[str, object]]) -> None:
    rows = list(rows)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--limit", type=int, default=2_000)
    parser.add_argument(
        "--csv", type=Path, default=HERE / "task004_line_sections.csv"
    )
    parser.add_argument(
        "--summary", type=Path, default=HERE / "task004_line_summary.json"
    )
    args = parser.parse_args()

    rows = []
    for prime in primes_upto(args.limit):
        if prime < 7:
            continue
        zeros = apery_zeros(prime)
        divided_zeros = apery_zeros_divided(prime)
        if zeros != divided_zeros:
            raise AssertionError(("recurrence mismatch", prime, zeros, divided_zeros))
        if len(zeros) >= 2:
            rows.append(row(prime, zeros))

    summary = {"limit": args.limit, **aggregate(rows)}
    write_csv(args.csv, rows)
    args.summary.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    print(f"CSV: {args.csv}")
    print(f"summary: {args.summary}")


if __name__ == "__main__":
    main()
