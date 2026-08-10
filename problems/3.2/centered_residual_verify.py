#!/usr/bin/env python3
"""Classify the first fully-deflated residual mechanisms.

The canonical lists below are produced independently by exact runs of
``fully_deflated_corank_verify.py --show-records`` at ambient heights 28,
32, and 36.  For every listed prime/triple this script reconstructs the common roots of the fully
deflated direct pencil and then enumerates all return times inside the selected
span.  It verifies two distinct mechanisms: skipped residual chains whose
refinement contains a reflection-centered adjacent pair, and primitive
off-center phantoms whose common root is not an Apéry-number zero.

All arithmetic is exact.  Ordinary Python automatically re-executes the file
under SageMath.
"""

from __future__ import annotations

import os
from pathlib import Path
import shutil
import sys


SAGE_CACHE = "/tmp/ramanujan-centered-residual-sage"
os.environ.setdefault("DOT_SAGE", SAGE_CACHE)
Path(os.environ["DOT_SAGE"]).mkdir(parents=True, exist_ok=True)

try:
    from sage.all import GF, PolynomialRing, gcd  # type: ignore[import-not-found]
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
import pairpalm_verify as ppv  # noqa: E402


RECORDS_FD_H28 = (
    (2, 4, 18, 653),
    (2, 4, 20, 653),
    (2, 10, 12, 131),
    (2, 18, 4, 653),
    (2, 18, 6, 653),
    (4, 18, 2, 653),
    (6, 18, 2, 653),
    (10, 2, 12, 131),
    (12, 2, 10, 131),
    (12, 10, 2, 131),
    (18, 4, 2, 653),
    (20, 4, 2, 653),
)

# New skipped records visible at ambient height 32.  The middle reduced
# content is 24649 = 157^2; the last coordinate below is deliberately the
# characteristic prime rather than the full content.
SKIPPED_FD_H32_NEW = (
    (4, 14, 14, 157),
    (14, 4, 14, 157),
    (14, 14, 4, 157),
)

# (a,b,c,p,root,expected actual Apéry zero set).  The first two records occur
# in the height-32 content census.  The remaining entries include every new
# primitive record in the height-36 full-content audit and the height-76
# endpoint phantom pair.
PRIMITIVE_PHANTOMS = (
    (13, 3, 15, 431, 279, (56, 374)),
    (15, 3, 13, 431, 120, (56, 374)),
    (6, 14, 13, 499, 417, (67, 170, 203, 295, 328, 431)),
    (13, 14, 6, 499, 48, (67, 170, 203, 295, 328, 431)),
    (4, 14, 17, 157, 83, ()),
    (17, 14, 4, 157, 38, ()),
    (5, 20, 10, 1297, 360, (459, 530, 766, 837)),
    (10, 20, 5, 1297, 901, (459, 530, 766, 837)),
    (3, 40, 31, 7411, 4681, ()),
    (31, 40, 3, 7411, 2655, ()),
)


def reduction(polynomial, ring):
    return ring(polynomial)


def centered_pair(root, left: int, right: int, prime: int) -> bool:
    """Whether the two absolute zero positions sum to ``p-1``."""

    return (2 * int(root) + left + right + 1) % prime == 0


def main() -> None:
    polynomials = me.build_gap_polynomials(74)
    root_count = 0
    print("a b c p root zero_offsets intermediate centered_pairs")
    skipped_records = RECORDS_FD_H28 + SKIPPED_FD_H32_NEW
    for a, b, c, prime in skipped_records:
        field = GF(prime)
        ring = PolynomialRing(field, "z")

        direct_zz = (
            acv.central_deflation(polynomials[a], a),
            me.shifted(acv.central_deflation(polynomials[b], b), a),
            me.shifted(acv.central_deflation(polynomials[c], c), a + b),
        )
        direct = [reduction(polynomial, ring) for polynomial in direct_zz]
        common = gcd(gcd(direct[0], direct[1]), direct[2]).monic()
        roots = common.roots(multiplicities=False)
        assert roots, (a, b, c, prime, common)
        assert len(roots) == common.degree(), (
            "non-rational common roots require a separate audit",
            a,
            b,
            c,
            prime,
            common,
        )

        span = a + b + c
        selected = (0, a, a + b, span)
        for root in roots:
            root_count += 1
            zero_offsets = [0]
            for h in range(1, span + 1):
                if reduction(polynomials[h], ring)(root) == 0:
                    zero_offsets.append(h)
            assert all(offset in zero_offsets for offset in selected)

            intermediate = [offset for offset in zero_offsets if offset not in selected]
            selected_centered = [
                (left, right)
                for left, right in zip(selected, selected[1:])
                if centered_pair(root, left, right, prime)
            ]
            centered_pairs = [
                (left, right)
                for left, right in zip(zero_offsets, zero_offsets[1:])
                if centered_pair(root, left, right, prime)
            ]

            # Full deflation removes a simple center root from each selected
            # gap.  Every surviving root in the census is instead a skipped
            # subchain of a longer reflection-symmetric return orbit.
            assert not selected_centered, (
                a,
                b,
                c,
                prime,
                int(root),
                selected_centered,
            )
            assert intermediate, (a, b, c, prime, int(root), zero_offsets)
            assert centered_pairs, (a, b, c, prime, int(root), zero_offsets)
            print(
                a,
                b,
                c,
                prime,
                int(root),
                zero_offsets,
                intermediate,
                centered_pairs,
            )

    print(
        "CENTERED_RESIDUAL_VERIFY PASS "
        f"records={len(skipped_records)} roots={root_count}"
    )

    for a, b, c, prime, expected_root, expected_apery_zeros in PRIMITIVE_PHANTOMS:
        field = GF(prime)
        ring = PolynomialRing(field, "z")
        direct_zz = (
            acv.central_deflation(polynomials[a], a),
            me.shifted(acv.central_deflation(polynomials[b], b), a),
            me.shifted(acv.central_deflation(polynomials[c], c), a + b),
        )
        direct = [reduction(polynomial, ring) for polynomial in direct_zz]
        common = gcd(gcd(direct[0], direct[1]), direct[2]).monic()
        roots = common.roots(multiplicities=False)
        assert roots == [field(expected_root)], (a, b, c, prime, common, roots)
        assert common == ring.gen() - field(expected_root), (
            a,
            b,
            c,
            prime,
            common,
        )

        span = a + b + c
        root = roots[0]
        selected = [0, a, a + b, span]
        zero_offsets = [0] + [
            h
            for h in range(1, span + 1)
            if reduction(polynomials[h], ring)(root) == 0
        ]
        assert zero_offsets == selected, (a, b, c, prime, zero_offsets)
        assert not any(
            centered_pair(root, left, right, prime)
            for left, right in zip(zero_offsets, zero_offsets[1:])
        )

        carrier = acv.structural_carrier(
            span,
            me.apery_values(span),
            me.pell_values(span + 1),
        )
        # The displayed common F_p-root makes every coefficient of
        # Res(F,G+T*J) divisible by p.  Checking p outside the structural
        # carrier is therefore the exact local condition needed here; the
        # separate full-content audit computes the integer resultant itself.
        assert carrier % prime != 0
        assert ppv.apery_zeros_division(prime) == expected_apery_zeros
        assert ppv.apery_zeros_cleared(prime) == expected_apery_zeros
        assert expected_root not in expected_apery_zeros
        print(
            "PRIMITIVE_PHANTOM PASS "
            f"triple={(a, b, c)} p={prime} root={int(root)} "
            f"zero_offsets={zero_offsets} "
            f"apery_zero_set={expected_apery_zeros}"
        )


if __name__ == "__main__":
    main()
