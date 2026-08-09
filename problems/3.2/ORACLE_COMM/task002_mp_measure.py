#!/usr/bin/env python3
"""Codex Task 002(C): multiplicative-collision measurements.

The task asks for signed frequencies 1 <= |k| <= K and excludes k=k', but
its proposed benchmark includes the diagonal K Z.  We therefore record:

* the literal ordered off-diagonal statistic;
* the total statistic including k=k';
* the benchmark printed in the task;
* cardinality-correct versions using L=2K signed frequencies.

All collision counts are exact integers.  The script uses only stdlib.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter
from pathlib import Path
from statistics import median
from typing import Iterable, Sequence


HERE = Path(__file__).resolve().parent


def primes_upto(limit: int) -> list[int]:
    sieve = bytearray(b"\x01") * (limit + 1)
    if limit >= 0:
        sieve[0] = 0
    if limit >= 1:
        sieve[1] = 0
    for prime in range(2, math.isqrt(limit) + 1):
        if sieve[prime]:
            start = prime * prime
            sieve[start : limit + 1 : prime] = b"\x00" * (
                (limit - start) // prime + 1
            )
    return [n for n in range(2, limit + 1) if sieve[n]]


def apery_zeros(prime: int) -> tuple[int, ...]:
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


def integer_cuberoot(value: int) -> int:
    root = int(round(value ** (1.0 / 3.0)))
    while (root + 1) ** 3 <= value:
        root += 1
    while root**3 > value:
        root -= 1
    return root


def frequency_set(prime: int, half_width: int, signed: bool) -> tuple[int, ...]:
    half_width = min(half_width, (prime - 1) // 2)
    positive = tuple(range(1, half_width + 1))
    if not signed:
        return positive
    return positive + tuple(range(prime - half_width, prime))


def product_collision(
    prime: int, zeros: Sequence[int], half_width: int, signed: bool
) -> dict[str, int | float]:
    frequencies = frequency_set(prime, half_width, signed)
    histogram: Counter[int] = Counter()
    for k in frequencies:
        for residue in zeros:
            histogram[(k * residue) % prime] += 1
    total = sum(count * count for count in histogram.values())
    diagonal = len(frequencies) * len(zeros)
    off_diagonal = total - diagonal
    if off_diagonal < 0:
        raise AssertionError((prime, zeros, half_width, total, diagonal))
    return {
        "frequency_count": len(frequencies),
        "total": total,
        "diagonal": diagonal,
        "off_diagonal": off_diagonal,
    }


def row(prime: int, zeros: Sequence[int], scale: str, half_width: int) -> dict[str, object]:
    signed = product_collision(prime, zeros, half_width, signed=True)
    positive = product_collision(prime, zeros, half_width, signed=False)
    z = len(zeros)
    signed_count = int(signed["frequency_count"])

    # This is exactly the expression requested, despite its mismatch with
    # both the exclusion k != k' and the 2K signed frequency cardinality.
    task_prediction = half_width * half_width * z * z / prime + half_width * z
    signed_total_prediction = signed_count * z + (signed_count * z) ** 2 / prime
    signed_off_prediction = signed_count * (signed_count - 1) * z * z / prime
    fixed_margin_off_prediction = (
        signed_count * (signed_count - 1) * z * (z - 1) / (prime - 2)
        if prime > 2
        else 0.0
    )
    fixed_margin_total_prediction = (
        int(signed["diagonal"]) + fixed_margin_off_prediction
    )
    return {
        "p": prime,
        "Z": z,
        "scale": scale,
        "K0": half_width,
        "signed_frequency_count": signed_count,
        "literal_signed_offdiag": signed["off_diagonal"],
        "signed_total_including_diag": signed["total"],
        "signed_diagonal": signed["diagonal"],
        "task_prediction": f"{task_prediction:.12g}",
        "literal_over_task_prediction": (
            f"{int(signed['off_diagonal']) / task_prediction:.12g}"
            if task_prediction
            else "nan"
        ),
        "signed_total_prediction": f"{signed_total_prediction:.12g}",
        "signed_total_ratio": (
            f"{int(signed['total']) / signed_total_prediction:.12g}"
            if signed_total_prediction
            else "nan"
        ),
        "signed_off_poisson_prediction": f"{signed_off_prediction:.12g}",
        "signed_off_poisson_ratio": (
            f"{int(signed['off_diagonal']) / signed_off_prediction:.12g}"
            if signed_off_prediction
            else "nan"
        ),
        "signed_off_fixed_margin_prediction": f"{fixed_margin_off_prediction:.12g}",
        "signed_off_fixed_margin_ratio": (
            f"{int(signed['off_diagonal']) / fixed_margin_off_prediction:.12g}"
            if fixed_margin_off_prediction
            else "nan"
        ),
        "signed_total_fixed_margin_prediction": (
            f"{fixed_margin_total_prediction:.12g}"
        ),
        "signed_total_fixed_margin_ratio": (
            f"{int(signed['total']) / fixed_margin_total_prediction:.12g}"
            if fixed_margin_total_prediction
            else "nan"
        ),
        "positive_offdiag": positive["off_diagonal"],
        "positive_total": positive["total"],
    }


def aggregate(rows: Sequence[dict[str, object]], scale: str) -> dict[str, object]:
    selected = [entry for entry in rows if entry["scale"] == scale and int(entry["Z"]) >= 2]

    def total(field: str) -> float:
        return sum(float(entry[field]) for entry in selected)

    observed_off = total("literal_signed_offdiag")
    observed_total = total("signed_total_including_diag")
    task_prediction = total("task_prediction")
    corrected_off = total("signed_off_fixed_margin_prediction")
    poisson_total = total("signed_total_prediction")
    fixed_margin_total = total("signed_total_fixed_margin_prediction")
    literal_ratios = [
        float(entry["literal_over_task_prediction"])
        for entry in selected
        if entry["literal_over_task_prediction"] != "nan"
    ]
    return {
        "scale": scale,
        "active_primes": len(selected),
        "sum_literal_offdiag": int(observed_off),
        "sum_task_prediction": task_prediction,
        "aggregate_literal_over_task": observed_off / task_prediction,
        "median_literal_over_task": median(literal_ratios),
        "sum_signed_total": int(observed_total),
        "sum_signed_total_poisson_prediction": poisson_total,
        "aggregate_signed_total_poisson_ratio": observed_total / poisson_total,
        "sum_signed_total_fixed_margin_prediction": fixed_margin_total,
        "aggregate_signed_total_fixed_margin_ratio": (
            observed_total / fixed_margin_total
        ),
        "sum_fixed_margin_off_prediction": corrected_off,
        "aggregate_off_fixed_margin_ratio": (
            observed_off / corrected_off if corrected_off else float("nan")
        ),
    }


def write_csv(path: Path, rows: Iterable[dict[str, object]]) -> None:
    rows = list(rows)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--limit", type=int, default=10_000)
    parser.add_argument(
        "--csv", type=Path, default=HERE / "task002_mp_measurements.csv"
    )
    parser.add_argument(
        "--summary", type=Path, default=HERE / "task002_mp_summary.json"
    )
    args = parser.parse_args()

    rows: list[dict[str, object]] = []
    prime_count = 0
    active_count = 0
    for prime in primes_upto(args.limit):
        if prime < 7:
            continue
        prime_count += 1
        zeros = apery_zeros(prime)
        if len(zeros) >= 2:
            active_count += 1
        widths = {
            "sqrt": math.isqrt(prime),
            "cuberoot": integer_cuberoot(prime),
        }
        for scale, width in widths.items():
            rows.append(row(prime, zeros, scale, width))

    summaries = [aggregate(rows, scale) for scale in ("sqrt", "cuberoot")]
    payload = {
        "limit": args.limit,
        "prime_count_p_ge_7": prime_count,
        "primes_with_Z_ge_2": active_count,
        "summaries": summaries,
    }
    write_csv(args.csv, rows)
    args.summary.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))
    print(f"CSV: {args.csv}")
    print(f"summary: {args.summary}")


if __name__ == "__main__":
    main()
