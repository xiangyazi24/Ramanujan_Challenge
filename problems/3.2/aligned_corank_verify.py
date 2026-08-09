#!/usr/bin/env python3
"""Exact audit of the aligned quadruple-corank content.

For a non-arithmetic-progression gap triple ``(a,b,c)``, put

    F = N_a(x),  H = N_(b+c)(x+a),

and take ``G = N_b(x+a)`` off the palindromic slice.  On ``c == a``
with ``a != b``, take ``G = N^o_b(x+a)`` instead, where the central
linear factor is removed from ``N_b`` when ``b`` is even.  This script
computes the exact integer

    C_(a,b,c) = content_T Res_x(F, G + T*H).

The reduced content ``C*`` is obtained by removing every prime supported
on

    U_H = product_(j<=H) j! b_j V_(j+1).

All polynomial arithmetic, contents, resultants, gcds, and saturations are
exact.  Floating-point logarithms are printed only as reconnaissance.  The
script also verifies, term by term, that ``C*`` divides the older saturated
scalar gcd.  Ordinary Python automatically re-executes the file under
SageMath.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from hashlib import sha256
import os
from pathlib import Path
import shutil
import sys


SAGE_CACHE = "/tmp/ramanujan-aligned-corank-sage"
os.environ.setdefault("DOT_SAGE", SAGE_CACHE)
Path(os.environ["DOT_SAGE"]).mkdir(parents=True, exist_ok=True)

try:
    from sage.all import (  # type: ignore[import-not-found]
        PolynomialRing,
        RR,
        ZZ,
        factorial,
        gcd,
    )
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
    "c7174f329e889bbb1f999a235015042987e14575a2a9bff111fe2ed4e180af4f"
)

PARAMETER_RING = PolynomialRing(ZZ, "T")
T = PARAMETER_RING.gen()
RESULTANT_RING = PolynomialRing(PARAMETER_RING, "y")


@dataclass(frozen=True)
class AlignedDatum:
    """Exact data independent of the saturation height."""

    content: int
    first_resultant: int
    second_resultant: int


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

    reduced = abs(ZZ(value))
    while True:
        common = gcd(reduced, carrier)
        if common == 1:
            return reduced
        reduced //= common


def lift_polynomial(polynomial):
    """Coerce a polynomial in ``ZZ[x]`` into ``ZZ[T][y]``."""

    return RESULTANT_RING([PARAMETER_RING(coefficient) for coefficient in polynomial])


def central_deflation(polynomial, index: int):
    if index % 2:
        return polynomial
    quotient, remainder = polynomial.quo_rem(2 * me.X + index + 1)
    assert remainder == 0
    return quotient


def integer_content(polynomial) -> int:
    content = ZZ.zero()
    for coefficient in polynomial:
        content = gcd(content, ZZ(coefficient))
    return abs(int(content))


def aligned_datum(polynomials: list, a: int, b: int, c: int) -> AlignedDatum:
    """Compute one exact aligned content and audit its endpoint coefficients."""

    assert not (a == b == c)
    first_zz = polynomials[a]
    if c == a:
        second_zz = me.shifted(central_deflation(polynomials[b], b), a)
    else:
        second_zz = me.shifted(polynomials[b], a)
    third_zz = me.shifted(polynomials[b + c], a)
    assert third_zz.degree() > second_zz.degree()

    first = lift_polynomial(first_zz)
    second = lift_polynomial(second_zz)
    third = lift_polynomial(third_zz)
    resultant = first.resultant(second + T * third)
    assert resultant != 0

    first_resultant = abs(int(first_zz.resultant(second_zz)))
    second_resultant = abs(int(first_zz.resultant(third_zz)))
    degree_padding = third_zz.degree() - second_zz.degree()
    padded_constant = int(first_zz.leading_coefficient()) ** degree_padding

    # These identities certify that the content is an aligned refinement of
    # the old scalar gcd, including the generic-degree padding at T = 0.
    assert abs(int(resultant[0])) == padded_constant * first_resultant
    assert resultant.degree() == first_zz.degree()
    assert abs(int(resultant[resultant.degree()])) == second_resultant

    content = integer_content(resultant)
    assert content > 0
    assert (padded_constant * first_resultant) % content == 0
    assert second_resultant % content == 0
    return AlignedDatum(content, first_resultant, second_resultant)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-height", type=int, default=DEFAULT_MAX_HEIGHT)
    parser.add_argument(
        "--show-records",
        action="store_true",
        help="print every nontrivial reduced content at the maximum height",
    )
    parser.add_argument(
        "--skip-digest-check",
        action="store_true",
        help="do not enforce the canonical digest at the default height",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    maximum = args.max_height
    if maximum < 8:
        raise SystemExit("max-height must be at least 8")

    polynomials = me.build_gap_polynomials(maximum)
    apery = me.apery_values(maximum)
    pell = me.pell_values(maximum + 1)
    triples = [
        triple for triple in gap_triples(maximum) if len(set(triple)) > 1
    ]

    data: dict[tuple[int, int, int], AlignedDatum] = {}
    for index, triple in enumerate(triples, start=1):
        data[triple] = aligned_datum(polynomials, *triple)
        if index % 250 == 0:
            print(f"computed {index}/{len(triples)} aligned resultants", flush=True)

    print(
        "H non_ap aligned/H^3 scalar/H^3 progression/H^3 "
        "aligned_nontrivial"
    )
    maximum_payload: list[str] = []
    maximum_records: list[tuple[int, int, int, int]] = []
    for height in range(8, maximum + 1):
        carrier = structural_carrier(height, apery, pell)
        aligned_mass = RR.zero()
        scalar_mass = RR.zero()
        progression_mass = RR.zero()
        nontrivial = 0
        height_triples = [triple for triple in triples if sum(triple) <= height]
        payload: list[str] = []
        records: list[tuple[int, int, int, int]] = []

        for triple in height_triples:
            datum = data[triple]
            reduced_content = remove_supported_part(datum.content, carrier)
            reduced_first = remove_supported_part(datum.first_resultant, carrier)
            reduced_second = remove_supported_part(datum.second_resultant, carrier)
            reduced_scalar = gcd(reduced_first, reduced_second)

            assert reduced_scalar % reduced_content == 0
            aligned_mass += RR(reduced_content).log()
            scalar_mass += RR(reduced_scalar).log()
            nontrivial += int(reduced_content > 1)
            payload.append(f"{triple[0]},{triple[1]},{triple[2]},{reduced_content}")
            if reduced_content > 1:
                records.append((*triple, int(reduced_content)))

        for a in range(2, height // 3 + 1):
            if 3 * a > height:
                continue
            progression = me.separated_resultant(polynomials, a, a)
            progression_mass += RR(
                remove_supported_part(progression, carrier)
            ).log()

        if height == maximum:
            maximum_payload = payload
            maximum_records = records
        print(
            f"{height:2d} {len(height_triples):6d} "
            f"{float(aligned_mass / height**3):.9f} "
            f"{float(scalar_mass / height**3):.9f} "
            f"{float(progression_mass / height**3):.9f} "
            f"{nontrivial:5d}"
        )

    digest = sha256("\n".join(maximum_payload).encode("ascii")).hexdigest()
    print(f"height={maximum} exact_content_sha256={digest}")
    if args.show_records:
        print(f"height={maximum} nontrivial_records={maximum_records}")
    if maximum == DEFAULT_MAX_HEIGHT and not args.skip_digest_check:
        assert digest == EXPECTED_HEIGHT_20_DIGEST
        print("ALIGNED_CORANK_VERIFY PASS")


if __name__ == "__main__":
    main()
