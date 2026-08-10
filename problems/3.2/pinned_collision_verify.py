#!/usr/bin/env python3
"""Exact census for the half-integer pinning mechanism.

For an even total span ``s=a+b+c``, reflection forces
``N_s(-(s+1)/2)=0``.  At this pinned starting point the current aligned
pencil has first and middle values

    N_a(-(s+1)/2),
    N_b(-(s+1)/2+a)

(with the paper's central deflation on ``c=a``).  This script computes the
primitive integer gcd of their rational numerators, removes the structural
carrier, and measures its total mass.  Up to a separately selectable height
it also computes the exact aligned content and checks whether every surviving
prime is supported on the pinned gcd.

The calculation is exact.  Ordinary Python automatically re-executes the
file under SageMath.
"""

from __future__ import annotations

import argparse
import math
import os
from pathlib import Path
import shutil
import sys


SAGE_CACHE = "/tmp/ramanujan-pinned-collision-sage"
os.environ.setdefault("DOT_SAGE", SAGE_CACHE)
Path(os.environ["DOT_SAGE"]).mkdir(parents=True, exist_ok=True)

try:
    from sage.all import QQ, RR, ZZ, gcd  # type: ignore[import-not-found]
except ModuleNotFoundError:
    sage = shutil.which("sage")
    if sage is None:
        raise SystemExit("SageMath is required (the `sage` executable was not found)")
    environment = os.environ.copy()
    environment["DOT_SAGE"] = SAGE_CACHE
    os.execvpe(sage, [sage, "-python", __file__, *sys.argv[1:]], environment)

sys.path.insert(0, str(Path(__file__).resolve().parent))
import aligned_corank_verify as acv  # noqa: E402
import meso_explore as me  # noqa: E402


def numerator_at(polynomial, argument) -> int:
    """Return the signed numerator of an exact rational evaluation."""

    return int(QQ(polynomial(argument)).numerator())


def pinned_values(polynomials: list, a: int, b: int, c: int) -> tuple[int, int]:
    """Return the two integer numerators defining the pinned carrier."""

    span = a + b + c
    assert span % 2 == 0
    start = -QQ(span + 1) / 2
    first = polynomials[a]
    if c == a:
        middle = acv.central_deflation(polynomials[b], b)
    else:
        middle = polynomials[b]
    return (
        numerator_at(first, start),
        numerator_at(middle, start + a),
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-height", type=int, default=40)
    parser.add_argument(
        "--aligned-check-height",
        type=int,
        default=20,
        help="compute expensive aligned resultants only through this height",
    )
    parser.add_argument("--show-records", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    maximum = args.max_height
    check_height = min(args.aligned_check_height, maximum)
    if maximum < 8 or check_height < 0:
        raise SystemExit("invalid height")

    polynomials = me.build_gap_polynomials(maximum)
    apery = me.apery_values(maximum)
    pell = me.pell_values(maximum + 1)
    triples = [triple for triple in acv.gap_triples(maximum) if len(set(triple)) > 1]

    pinned_raw: dict[tuple[int, int, int], int] = {}
    for a, b, c in triples:
        if (a + b + c) % 2:
            pinned_raw[(a, b, c)] = 1
            continue
        first_value, middle_value = pinned_values(polynomials, a, b, c)
        pinned_raw[(a, b, c)] = int(gcd(abs(first_value), abs(middle_value)))

    print("H triples even_span pinned_nontrivial pinned_mass/H^3 pal_mass/H^3 nonpal_mass/H^3")
    maximum_records = []
    for height in range(8, maximum + 1):
        carrier = acv.structural_carrier(height, apery, pell)
        height_triples = [triple for triple in triples if sum(triple) <= height]
        mass = RR.zero()
        pal_mass = RR.zero()
        nonpal_mass = RR.zero()
        nontrivial = 0
        records = []
        for triple in height_triples:
            reduced = acv.remove_supported_part(pinned_raw[triple], carrier)
            if reduced > 1:
                contribution = RR(reduced).log()
                mass += contribution
                if triple[0] == triple[2]:
                    pal_mass += contribution
                else:
                    nonpal_mass += contribution
                nontrivial += 1
                records.append((*triple, int(reduced)))
        print(
            f"{height:2d} {len(height_triples):6d} "
            f"{sum(sum(triple) % 2 == 0 for triple in height_triples):6d} "
            f"{nontrivial:6d} {float(mass / height**3):.9f} "
            f"{float(pal_mass / height**3):.9f} "
            f"{float(nonpal_mass / height**3):.9f}"
        )
        if height == maximum:
            maximum_records = records

    if check_height >= 8:
        carrier = acv.structural_carrier(check_height, apery, pell)
        unpinned_records = []
        checked_nontrivial = 0
        for triple in triples:
            if sum(triple) > check_height:
                continue
            datum = acv.aligned_datum(polynomials, *triple)
            reduced_content = acv.remove_supported_part(datum.content, carrier)
            if reduced_content == 1:
                continue
            checked_nontrivial += 1
            unpinned = acv.remove_supported_part(
                reduced_content, pinned_raw[triple]
            )
            if unpinned > 1:
                unpinned_records.append((*triple, int(reduced_content), int(unpinned)))
        print(
            f"aligned_check_height={check_height} "
            f"aligned_nontrivial={checked_nontrivial} "
            f"unpinned_records={unpinned_records}"
        )

    if args.show_records:
        print(f"height={maximum} pinned_records={maximum_records}")


if __name__ == "__main__":
    main()
