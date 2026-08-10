#!/usr/bin/env python3
"""Exact audit of the fully deflated consecutive-segment pencil.

For a gap triple ``(a,b,c)`` put

    F = N^o_a(x),
    G = N^o_b(x+a),
    J = N^o_c(x+a+b),

where ``N^o_h`` is obtained by removing the reflection-center linear factor
when ``h`` is even.  Thus every polynomial represents one of the three
adjacent gaps of a four-zero chain, and every centered adjacent pair is
removed symmetrically.

The script computes

    C^fd_(a,b,c) = content_T Res_x(F, G + T*J)

and removes the same structural carrier as the paper.  It also verifies the
two endpoint coefficients and their nonvanishing.  All arithmetic is exact;
ordinary Python automatically re-executes the file under SageMath.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from hashlib import sha256
import os
from pathlib import Path
import shutil
import sys


SAGE_CACHE = "/tmp/ramanujan-fully-deflated-sage"
os.environ.setdefault("DOT_SAGE", SAGE_CACHE)
Path(os.environ["DOT_SAGE"]).mkdir(parents=True, exist_ok=True)

try:
    from sage.all import PolynomialRing, RR, ZZ, gcd  # type: ignore[import-not-found]
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


PARAMETER_RING = PolynomialRing(ZZ, "T")
T = PARAMETER_RING.gen()
RESULTANT_RING = PolynomialRing(PARAMETER_RING, "y")

EXPECTED_DIGESTS = {
    20: "4e4c87ab9b630c3590c8d5a788051aa1dc91617a3107b034e64d2cc5be217872",
    28: "26472e0d3c81f8e099a0d8c51881c735b67a6e6b0c757e3f330f0d6a04f3e399",
    32: "437f7328a4ea7d78a62b3e781165a76a31f277f7ec6983f9f6c3a4e67ee05efa",
    36: "196327788a7ea9adbbb141a5c6161ae96125ea5b8242bc6603352080b32aac76",
}


@dataclass(frozen=True)
class Datum:
    content: int
    first_resultant: int
    second_resultant: int


def lift(polynomial):
    return RESULTANT_RING([PARAMETER_RING(coefficient) for coefficient in polynomial])


def datum(polynomials: list, a: int, b: int, c: int) -> Datum:
    first_zz = acv.central_deflation(polynomials[a], a)
    second_zz = me.shifted(acv.central_deflation(polynomials[b], b), a)
    third_zz = me.shifted(acv.central_deflation(polynomials[c], c), a + b)

    first = lift(first_zz)
    second = lift(second_zz)
    third = lift(third_zz)
    resultant = first.resultant(second + T * third)
    assert resultant != 0

    first_resultant = abs(int(first_zz.resultant(second_zz)))
    second_resultant = abs(int(first_zz.resultant(third_zz)))
    assert first_resultant > 0 and second_resultant > 0

    generic_degree = max(second_zz.degree(), third_zz.degree())
    constant_padding = generic_degree - second_zz.degree()
    top_padding = generic_degree - third_zz.degree()
    expected_constant = (
        abs(int(first_zz.leading_coefficient())) ** constant_padding
        * first_resultant
    )
    expected_top = (
        abs(int(first_zz.leading_coefficient())) ** top_padding
        * second_resultant
    )
    assert abs(int(resultant[0])) == expected_constant
    assert abs(int(resultant[resultant.degree()])) == expected_top

    content = acv.integer_content(resultant)
    assert content > 0
    assert expected_constant % content == 0
    assert expected_top % content == 0
    return Datum(content, first_resultant, second_resultant)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-height", type=int, default=20)
    parser.add_argument("--show-records", action="store_true")
    parser.add_argument(
        "--skip-digest-check",
        action="store_true",
        help="do not enforce a recorded digest at canonical audited heights",
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
    triples = [triple for triple in acv.gap_triples(maximum) if len(set(triple)) > 1]
    data = {}
    for index, triple in enumerate(triples, start=1):
        data[triple] = datum(polynomials, *triple)
        if index % 250 == 0:
            print(f"computed {index}/{len(triples)} fully deflated pencils", flush=True)

    print("H non_ap mass/H^3 nontrivial")
    maximum_payload = []
    maximum_records = []
    for height in range(8, maximum + 1):
        carrier = acv.structural_carrier(height, apery, pell)
        mass = RR.zero()
        records = []
        payload = []
        height_triples = [triple for triple in triples if sum(triple) <= height]
        for triple in height_triples:
            reduced = acv.remove_supported_part(data[triple].content, carrier)
            mass += RR(reduced).log()
            payload.append(f"{triple[0]},{triple[1]},{triple[2]},{int(reduced)}")
            if reduced > 1:
                records.append((*triple, int(reduced)))
        print(
            f"{height:2d} {len(height_triples):6d} "
            f"{float(mass / height**3):.9f} {len(records):5d}"
        )
        if height == maximum:
            maximum_payload = payload
            maximum_records = records

    digest = sha256("\n".join(maximum_payload).encode("ascii")).hexdigest()
    print(f"height={maximum} exact_content_sha256={digest}")
    if not args.skip_digest_check and maximum in EXPECTED_DIGESTS:
        assert digest == EXPECTED_DIGESTS[maximum]
    if args.show_records:
        print(f"height={maximum} nontrivial_records={maximum_records}")


if __name__ == "__main__":
    main()
