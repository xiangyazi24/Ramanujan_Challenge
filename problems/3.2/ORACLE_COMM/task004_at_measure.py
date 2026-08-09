#!/usr/bin/env python3
"""Numerical test of the anti-tail (AT) scale for dyadic prime blocks."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


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


def measure(x: int, zero_sets: dict[int, tuple[int, ...]]) -> dict[str, object]:
    length = x * x
    loads = bytearray(length)
    primes = [prime for prime in zero_sets if x < prime <= 2 * x]
    lam = 0.0
    incidence_mass = 0
    for prime in primes:
        zeros = zero_sets[prime]
        lam += len(zeros) / prime
        incidence_mass += len(zeros)
        for residue in zeros:
            for m in range(residue, length, prime):
                if loads[m] == 255:
                    raise OverflowError("K_X(m) exceeded byte storage")
                loads[m] += 1

    maximum = max(loads, default=0)
    argmax = loads.index(maximum) if loads else 0
    max_count = loads.count(maximum) if loads else 0
    target = x ** (2.0 / 3.0) * lam
    hm2_scale = x * lam
    return {
        "X": x,
        "prime_count": len(primes),
        "active_prime_count": sum(bool(zero_sets[p]) for p in primes),
        "sum_Z": incidence_mass,
        "lambda_X": f"{lam:.15g}",
        "max_K": maximum,
        "first_argmax": argmax,
        "number_of_argmax": max_count,
        "X_2over3_lambda": f"{target:.15g}",
        "max_over_X_2over3_lambda": f"{maximum / target:.15g}" if target else "nan",
        "X_lambda": f"{hm2_scale:.15g}",
        "max_over_X_lambda": f"{maximum / hm2_scale:.15g}" if hm2_scale else "nan",
        "mean_K": f"{sum(loads) / length:.15g}",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--j-min", type=int, default=4)
    parser.add_argument("--j-max", type=int, default=12)
    parser.add_argument(
        "--csv", type=Path, default=HERE / "task004_at_measurements.csv"
    )
    parser.add_argument(
        "--summary", type=Path, default=HERE / "task004_at_summary.json"
    )
    args = parser.parse_args()

    xs = [1 << exponent for exponent in range(args.j_min, args.j_max + 1)]
    primes = [p for p in primes_upto(2 * max(xs)) if p >= 7]
    zero_sets = {prime: apery_zeros(prime) for prime in primes}
    rows = [measure(x, zero_sets) for x in xs]

    with args.csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    payload = {"j_min": args.j_min, "j_max": args.j_max, "rows": rows}
    args.summary.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))
    print(f"CSV: {args.csv}")
    print(f"summary: {args.summary}")


if __name__ == "__main__":
    main()
