#!/usr/bin/env python3
"""Measure fixed prefix-zero tails in the direct candidate window.

For n=3H+1 and 2H<p<=3H+1, let

    z_p(H) = #{0<=s<=H : p divides A_s}.

Fixed-degree anchored Padé certificates retain every good candidate prime
with z_p(H) above a fixed threshold through numerator content.  This script
uses the canonical complete (p,s) zero-pair bank to measure that obstruction.
It supplies finite evidence only; it does not infer an asymptotic.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from math import log
from pathlib import Path
from struct import Struct


PAIR = Struct("<II")
DEFAULT_DATA = Path(__file__).parents[2] / "data_zp_pairs.bin"
DEFAULT_HEIGHTS = (20_000, 50_000, 100_000, 200_000, 400_000, 600_000)


def read_pairs(path: Path) -> dict[int, list[int]]:
    raw = path.read_bytes()
    if len(raw) % PAIR.size:
        raise ValueError("zero-pair bank has a partial record")
    zero_sets: dict[int, list[int]] = defaultdict(list)
    for offset in range(0, len(raw), PAIR.size):
        prime, node = PAIR.unpack_from(raw, offset)
        zero_sets[prime].append(node)
    return zero_sets


def audit_height(
    height: int,
    threshold: int,
    zero_sets: dict[int, list[int]],
) -> None:
    outer_index = 3 * height + 1
    rows: list[tuple[int, int, bool]] = []
    for prime, nodes in zero_sets.items():
        if not 2 * height < prime <= outer_index:
            continue
        zero_count = sum(node <= height for node in nodes)
        if zero_count < threshold:
            continue
        moving_node = outer_index - prime
        rows.append((prime, zero_count, moving_node in nodes))

    mass = sum(log(prime) for prime, _, _ in rows)
    good_mass = sum(log(prime) for prime, _, target in rows if not target)
    histogram = {
        zero_count: sum(row_zero_count == zero_count for _, row_zero_count, _ in rows)
        for zero_count in sorted({zero_count for _, zero_count, _ in rows})
    }
    print(
        f"H={height} threshold={threshold} count={len(rows)} "
        f"mass/H={mass / height:.9f} "
        f"good_mass/H={good_mass / height:.9f} "
        f"targets={sum(target for _, _, target in rows)} "
        f"histogram={histogram}"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", type=Path, default=DEFAULT_DATA)
    parser.add_argument("--threshold", type=int, default=4)
    parser.add_argument("heights", nargs="*", type=int)
    args = parser.parse_args()
    if args.threshold < 1:
        raise SystemExit("threshold must be positive")

    zero_sets = read_pairs(args.data)
    heights = tuple(args.heights) or DEFAULT_HEIGHTS
    for height in heights:
        audit_height(height, args.threshold, zero_sets)


if __name__ == "__main__":
    main()
