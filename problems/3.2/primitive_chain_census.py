#!/usr/bin/env python3
"""Audit short primitive four-zero chains in the saved Apéry census.

``data_zp_pairs.bin`` is the repository's headerless little-endian sequence
of strictly increasing ``uint32`` pairs ``(p,r)`` with ``b_r == 0 (mod p)``.
For each stored prime this script enumerates sliding quadruples in the *full*
ordered zero set, so every reported gap is genuinely consecutive.  It checks
whether a quadruple contains the unique reflection-centered adjacent pair and
whether its total span is at most ``sqrt(p)``.

This is a finite audit, not a proof of the short-chain statement.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from fractions import Fraction
import hashlib
from pathlib import Path
import struct
from typing import Optional


HERE = Path(__file__).resolve().parent
DATA = HERE / "data_zp_pairs.bin"
EXPECTED_SHA256 = "8746d0b400c1b669b001eae955c602908a10c9ee4cb3cac62c6676ea2ddd874d"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=DATA)
    parser.add_argument(
        "--expected-sha256",
        help="optional expected digest for a noncanonical pair file",
    )
    return parser.parse_args()


def load(
    data: Path, expected_sha256: Optional[str]
) -> tuple[dict[int, tuple[int, ...]], str, int]:
    raw = data.read_bytes()
    if len(raw) % 8:
        raise AssertionError("partial uint32 pair in zero-set data")
    digest = hashlib.sha256(raw).hexdigest()
    if expected_sha256 is not None:
        assert digest == expected_sha256

    grouped: dict[int, list[int]] = defaultdict(list)
    previous = (-1, -1)
    for pair in struct.iter_unpack("<II", raw):
        prime, residue = pair
        assert pair > previous
        assert residue < prime
        grouped[prime].append(residue)
        previous = pair
    return (
        {prime: tuple(residues) for prime, residues in grouped.items()},
        digest,
        len(raw) // 8,
    )


def main() -> None:
    args = parse_args()
    canonical = args.data.resolve() == DATA.resolve()
    expected = args.expected_sha256 or (EXPECTED_SHA256 if canonical else None)
    zero_sets, digest, records = load(args.data, expected)
    quadruples = []
    off_center = []
    short_off_center = []

    for prime, zeros in zero_sets.items():
        assert tuple(sorted(set(zeros))) == zeros
        zero_set = set(zeros)
        assert all(prime - 1 - residue in zero_set for residue in zeros)
        assert all(residue + 1 not in zero_set for residue in zeros)

        for index in range(len(zeros) - 3):
            chain = zeros[index : index + 4]
            gaps = tuple(chain[j + 1] - chain[j] for j in range(3))
            centered = tuple(
                j for j in range(3) if chain[j] + chain[j + 1] == prime - 1
            )
            row = (prime, chain, gaps, centered)
            quadruples.append(row)
            if not centered:
                off_center.append(row)
                span = chain[-1] - chain[0]
                if span * span <= prime:
                    short_off_center.append(row)

    first = min(
        off_center,
        key=lambda row: (row[1][-1] - row[1][0], row[0], row[1]),
    )
    closest = min(
        off_center,
        key=lambda row: (
            Fraction((row[1][-1] - row[1][0]) ** 2, row[0]),
            row[0],
            row[1],
        ),
    )
    closest_ratio_squared = Fraction(
        (closest[1][-1] - closest[1][0]) ** 2, closest[0]
    )
    maximum = max(zero_sets.items(), key=lambda item: (len(item[1]), -item[0]))

    print(f"data={args.data}")
    print(f"sha256={digest}")
    print(
        f"records={records} active_primes={len(zero_sets)} "
        f"maximum_prime={max(zero_sets)}"
    )
    print(f"zero_count_distribution={sorted(Counter(map(len, zero_sets.values())).items())}")
    print(f"maximum_zero_set=p{maximum[0]}:{maximum[1]}")
    print(
        f"consecutive_quadruples={len(quadruples)} "
        f"off_center={len(off_center)} short_off_center={len(short_off_center)}"
    )
    print(f"first_off_center_by_span={first}")
    print(f"closest_off_center_to_sqrt_scale={closest}")
    print(
        "closest_span_squared_over_p="
        f"{closest_ratio_squared.numerator}/{closest_ratio_squared.denominator} "
        f"approx={float(closest_ratio_squared):.12f}"
    )
    if short_off_center:
        first_short = min(
            short_off_center,
            key=lambda row: (row[0], row[1]),
        )
        print(f"first_short_off_center={first_short}")
        raise AssertionError("short off-center consecutive chain found")
    print("PRIMITIVE_CHAIN_CENSUS PASS")


if __name__ == "__main__":
    main()
