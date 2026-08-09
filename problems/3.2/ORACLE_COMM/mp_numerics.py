#!/usr/bin/env python3
"""Codex Task 003: positive-frequency M_p(k,k') and ratio energy.

All counts are exact integers.  Besides the literal task normalizations, the
output includes diagonal-consistent and fixed-cardinality random benchmarks;
the distinctions are mathematically substantial when Z(p) is small.
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
    """Independent check using the original divided Apéry recurrence."""

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


def falling_factorial(value: int, length: int) -> int:
    return math.prod(range(value - length + 1, value + 1))


def fixed_subset_mult_energy_expectation(prime: int, size: int) -> float:
    """Exact mean over uniform size-Z subsets of F_p^times.

    Ratio-energy quadruples have respectively 1, 2, 3, or 4 distinct
    entries. For N=p-1, their counts are
      N, 2N^2-N, 2N(N-2), N(N-2)^2.
    A quadruple with t distinct entries is present with probability
    (Z)_t/(N)_t.
    """

    group_size = prime - 1
    counts = (
        group_size,
        2 * group_size * group_size - group_size,
        2 * group_size * (group_size - 2),
        group_size * (group_size - 2) ** 2,
    )
    expectation = 0.0
    for distinct, count in enumerate(counts, start=1):
        if size >= distinct:
            expectation += (
                count
                * falling_factorial(size, distinct)
                / falling_factorial(group_size, distinct)
            )
    return expectation


def integer_cuberoot(value: int) -> int:
    root = int(round(value ** (1.0 / 3.0)))
    while (root + 1) ** 3 <= value:
        root += 1
    while root**3 > value:
        root -= 1
    return root


def small_frequency_energy(
    prime: int, zeros: Sequence[int], width: int
) -> tuple[int, int, int]:
    """Return (off-diagonal, total, exact diagonal) for k,k'=1,...,K."""

    histogram: Counter[int] = Counter()
    for k in range(1, width + 1):
        for residue in zeros:
            histogram[k * residue % prime] += 1
    total = sum(count * count for count in histogram.values())
    diagonal = width * len(zeros)
    return total - diagonal, total, diagonal


def multiplicative_energy(prime: int, zeros: Sequence[int]) -> tuple[int, int, int]:
    """Return energy, maximum multiplicity, and maximum away from k=1."""

    histogram: Counter[int] = Counter()
    for residue in zeros:
        inverse = pow(residue, -1, prime)
        for target in zeros:
            histogram[target * inverse % prime] += 1
    return (
        sum(count * count for count in histogram.values()),
        max(histogram.values(), default=0),
        max((count for ratio, count in histogram.items() if ratio != 1), default=0),
    )


def make_row(prime: int, zeros: Sequence[int]) -> dict[str, object]:
    z = len(zeros)
    output: dict[str, object] = {"p": prime, "Z": z, "zeros": " ".join(map(str, zeros))}
    for label, width in (
        ("sqrt", math.isqrt(prime)),
        ("cuberoot", integer_cuberoot(prime)),
    ):
        off, total, diagonal = small_frequency_energy(prime, zeros, width)
        task_prediction = width * width * z * z / prime + width * z
        fixed_off_prediction = width * (width - 1) * z * (z - 1) / (prime - 2)
        fixed_total_prediction = diagonal + fixed_off_prediction
        output.update(
            {
                f"K_{label}": width,
                f"off_{label}": off,
                f"total_{label}": total,
                f"diag_{label}": diagonal,
                f"task_prediction_{label}": f"{task_prediction:.12g}",
                f"off_over_task_{label}": f"{off / task_prediction:.12g}",
                f"fixed_off_prediction_{label}": f"{fixed_off_prediction:.12g}",
                f"off_over_fixed_{label}": (
                    f"{off / fixed_off_prediction:.12g}"
                    if fixed_off_prediction
                    else "nan"
                ),
                f"fixed_total_prediction_{label}": f"{fixed_total_prediction:.12g}",
                f"total_over_fixed_{label}": f"{total / fixed_total_prediction:.12g}",
            }
        )

    energy, max_ratio_multiplicity, max_nonidentity_multiplicity = multiplicative_energy(
        prime, zeros
    )
    off_pairs = z * (z - 1)
    occupancy_prediction = (
        z * z
        + off_pairs
        + off_pairs * (off_pairs - 1) / (prime - 2)
    )
    fixed_subset_prediction = fixed_subset_mult_energy_expectation(prime, z)
    task_mult_prediction = prime * z * z
    output.update(
        {
            "E_mult": energy,
            "max_M_k1": max_ratio_multiplicity,
            "max_M_k1_nonidentity": max_nonidentity_multiplicity,
            "task_mult_prediction_pZ2": task_mult_prediction,
            "E_mult_over_pZ2": f"{energy / task_mult_prediction:.12g}",
            "occupancy_mult_prediction": f"{occupancy_prediction:.12g}",
            "E_mult_over_occupancy": f"{energy / occupancy_prediction:.12g}",
            "fixed_subset_mult_prediction": f"{fixed_subset_prediction:.12g}",
            "E_mult_over_fixed_subset": f"{energy / fixed_subset_prediction:.12g}",
            "minimum_mult_energy": 2 * z * z - z,
            "excess_over_minimum": energy - (2 * z * z - z),
        }
    )
    return output


def aggregate(rows: Sequence[dict[str, object]]) -> dict[str, object]:
    result: dict[str, object] = {"active_primes": len(rows)}
    for label in ("sqrt", "cuberoot"):
        observed_off = sum(int(row[f"off_{label}"]) for row in rows)
        task_prediction = sum(float(row[f"task_prediction_{label}"]) for row in rows)
        fixed_off = sum(float(row[f"fixed_off_prediction_{label}"]) for row in rows)
        observed_total = sum(int(row[f"total_{label}"]) for row in rows)
        fixed_total = sum(float(row[f"fixed_total_prediction_{label}"]) for row in rows)
        result[label] = {
            "observed_off": observed_off,
            "task_prediction": task_prediction,
            "off_over_task": observed_off / task_prediction,
            "fixed_off_prediction": fixed_off,
            "off_over_fixed": observed_off / fixed_off,
            "observed_total": observed_total,
            "fixed_total_prediction": fixed_total,
            "total_over_fixed": observed_total / fixed_total,
            "zero_off_primes": sum(int(row[f"off_{label}"]) == 0 for row in rows),
            "median_off_over_task": median(float(row[f"off_over_task_{label}"]) for row in rows),
        }

    observed_mult = sum(int(row["E_mult"]) for row in rows)
    task_mult = sum(float(row["task_mult_prediction_pZ2"]) for row in rows)
    occupancy_mult = sum(float(row["occupancy_mult_prediction"]) for row in rows)
    fixed_subset_mult = sum(
        float(row["fixed_subset_mult_prediction"]) for row in rows
    )
    exceptional = [row for row in rows if int(row["excess_over_minimum"]) > 0]
    result["multiplicative_energy"] = {
        "observed": observed_mult,
        "task_pZ2_prediction": task_mult,
        "observed_over_task": observed_mult / task_mult,
        "occupancy_prediction": occupancy_mult,
        "observed_over_occupancy": observed_mult / occupancy_mult,
        "fixed_subset_prediction": fixed_subset_mult,
        "observed_over_fixed_subset": observed_mult / fixed_subset_mult,
        "sidon_primes": len(rows) - len(exceptional),
        "non_sidon_primes": len(exceptional),
        "max_M_k1": max(int(row["max_M_k1"]) for row in rows),
        "max_M_k1_nonidentity": max(
            int(row["max_M_k1_nonidentity"]) for row in rows
        ),
        "exceptions": [
            {
                "p": int(row["p"]),
                "Z": int(row["Z"]),
                "E_mult": int(row["E_mult"]),
                "minimum": int(row["minimum_mult_energy"]),
                "excess": int(row["excess_over_minimum"]),
                "max_M_k1": int(row["max_M_k1"]),
                "max_M_k1_nonidentity": int(row["max_M_k1_nonidentity"]),
            }
            for row in exceptional
        ],
    }
    return result


def write_csv(path: Path, rows: Iterable[dict[str, object]]) -> None:
    rows = list(rows)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--limit", type=int, default=5_000)
    parser.add_argument("--csv", type=Path, default=HERE / "mp_numerics.csv")
    parser.add_argument("--summary", type=Path, default=HERE / "mp_numerics_summary.json")
    args = parser.parse_args()

    rows = []
    primes_checked = 0
    for prime in primes_upto(args.limit):
        if prime < 7:
            continue
        primes_checked += 1
        zeros = apery_zeros(prime)
        divided_zeros = apery_zeros_divided(prime)
        if zeros != divided_zeros:
            raise AssertionError(("recurrence mismatch", prime, zeros, divided_zeros))
        if len(zeros) >= 2:
            rows.append(make_row(prime, zeros))

    summary = {
        "limit": args.limit,
        "primes_checked_p_ge_7": primes_checked,
        **aggregate(rows),
    }
    write_csv(args.csv, rows)
    args.summary.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    print(f"CSV: {args.csv}")
    print(f"summary: {args.summary}")


if __name__ == "__main__":
    main()
