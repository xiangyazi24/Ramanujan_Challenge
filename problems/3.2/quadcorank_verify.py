#!/usr/bin/env python3
"""Exact audit of the two quadruple-corank gcd masses.

The adjacent certificate pair

    gcd(S^*_(a,b), S^*_(b,c))

has a self-gcd diagonal when ``a == c``.  This script separates that
diagonal from the off-diagonal mass and compares it with the skipped-triple
pair

    gcd(S^*_(a,b), S^*_(a,b+c)).

Here ``S_(d,r) = |Res_x(N_d(x), N_r(x+d))|`` and ``S^*`` is obtained by
removing every prime supported on

    U_H = product_(j<=H) j! b_j V_(j+1).

All resultants, gcds, and saturations are exact integers.  Floating-point
logarithms are printed only as reconnaissance; the canonical digest is over
the exact reduced gcds.  Ordinary Python automatically re-executes the file
under SageMath.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
from pathlib import Path
import os
import shutil
import sys


SAGE_CACHE = "/tmp/ramanujan-quadcorank-sage"
os.environ.setdefault("DOT_SAGE", SAGE_CACHE)
Path(os.environ["DOT_SAGE"]).mkdir(parents=True, exist_ok=True)

try:
    from sage.all import RR, ZZ, factorial, gcd  # type: ignore[import-not-found]
except ModuleNotFoundError:
    sage = shutil.which("sage")
    if sage is None:
        raise SystemExit("SageMath is required (the `sage` executable was not found)")
    environment = os.environ.copy()
    environment["DOT_SAGE"] = SAGE_CACHE
    os.execvpe(sage, [sage, "-python", __file__, *sys.argv[1:]], environment)

sys.path.insert(0, str(Path(__file__).resolve().parent))
import meso_explore as me  # noqa: E402


DEFAULT_MAX_HEIGHT = 20
EXPECTED_HEIGHT_20_DIGEST = (
    "e1b61c3e46dc326e8a214af08d53a1fea0ec24fae2bfc552bc8f42472e8c1a93"
)


def gap_triples(height: int) -> list[tuple[int, int, int]]:
    """Return all ``a,b,c >= 2`` with ``a+b+c <= height``."""

    return [
        (a, b, c)
        for a in range(2, height - 3)
        for b in range(2, height - a - 1)
        for c in range(2, height - a - b + 1)
    ]


def structural_carrier(height: int, apery: list[int], pell: list[int]):
    carrier = ZZ.one()
    for index in range(1, height + 1):
        carrier *= factorial(index) * apery[index] * pell[index + 1]
    return carrier


def remove_supported_part(value, carrier):
    """Remove the full prime-power part of ``value`` supported on ``carrier``."""

    reduced = ZZ(value)
    while True:
        common = gcd(reduced, carrier)
        if common == 1:
            return reduced
        reduced //= common


def exact_log(value) -> RR:
    return RR(value).log()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-height", type=int, default=DEFAULT_MAX_HEIGHT)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    maximum = args.max_height
    if maximum < 8:
        raise SystemExit("max-height must be at least 8")

    polynomials = me.build_gap_polynomials(maximum)
    resultants = me.exact_resultants(polynomials, maximum)
    apery = me.apery_values(maximum)
    pell = me.pell_values(maximum + 1)

    # The exact absolute symmetry is essential to the old self-gcd diagonal.
    for d in range(2, maximum - 1):
        for r in range(2, maximum - d + 1):
            forward = me.separated_resultant(polynomials, d, r)
            reverse = me.separated_resultant(polynomials, r, d)
            assert forward == reverse == me.lookup_resultant(resultants, d, r)

    print(
        "H triples diag old_diag/H^4 old_off/H^3 "
        "old_all/H^3 skip/H^3 skip_avg skip_nontrivial"
    )
    maximum_payload: list[str] = []
    for height in range(8, maximum + 1):
        triples = gap_triples(height)
        carrier = structural_carrier(height, apery, pell)
        needed_pairs = {
            pair
            for a, b, c in triples
            for pair in ((a, b), (b, c), (a, b + c))
        }
        reduced = {
            pair: remove_supported_part(
                me.lookup_resultant(resultants, *pair), carrier
            )
            for pair in needed_pairs
        }

        old_diagonal = RR.zero()
        old_off_diagonal = RR.zero()
        skipped_mass = RR.zero()
        skipped_nontrivial = 0
        payload: list[str] = []
        diagonal_count = 0
        for a, b, c in triples:
            first = reduced[a, b]
            adjacent = reduced[b, c]
            skipped = reduced[a, b + c]
            old_gcd = gcd(first, adjacent)
            skipped_gcd = gcd(first, skipped)

            # Since c >= 2, the skipped pair can never be the same unordered
            # parameter pair as (a,b).  This removes the automatic self-gcd.
            assert tuple(sorted((a, b))) != tuple(sorted((a, b + c)))
            if a == c:
                assert first == adjacent and old_gcd == first
                old_diagonal += exact_log(old_gcd)
                diagonal_count += 1
            else:
                old_off_diagonal += exact_log(old_gcd)
            skipped_mass += exact_log(skipped_gcd)
            skipped_nontrivial += int(skipped_gcd > 1)
            payload.append(f"{a},{b},{c},{old_gcd},{skipped_gcd}")

        if height == maximum:
            maximum_payload = payload
        old_all = old_diagonal + old_off_diagonal
        print(
            f"{height:2d} {len(triples):7d} {diagonal_count:4d} "
            f"{float(old_diagonal / height**4):.9f} "
            f"{float(old_off_diagonal / height**3):.9f} "
            f"{float(old_all / height**3):.9f} "
            f"{float(skipped_mass / height**3):.9f} "
            f"{float(skipped_mass / len(triples)):.9f} "
            f"{skipped_nontrivial:5d}"
        )

    digest = sha256("\n".join(maximum_payload).encode("ascii")).hexdigest()
    print(f"height={maximum} exact_gcd_sha256={digest}")
    if maximum == DEFAULT_MAX_HEIGHT:
        assert digest == EXPECTED_HEIGHT_20_DIGEST
        print("QUADCORANK_VERIFY PASS")


if __name__ == "__main__":
    main()
