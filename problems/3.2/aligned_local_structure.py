#!/usr/bin/env python3
"""Audit the local structure behind the aligned pencil contents.

For every gap triple that is not all equal, up to the requested height,
this script computes the exact aligned content used in
``aligned_corank_verify.py``.  For each prime surviving the structural
saturation it then computes, over ``GF(p)``, the common gcd of ``F,G,J`` and
its squarefree degree.  Thus the output distinguishes

* coefficient-content valuation ``v_p(C)``;
* common-factor degree (with modular multiplicity); and
* the number of distinct common geometric roots.

The calculation is exact.  Ordinary Python automatically re-executes the
file under SageMath.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import shutil
import sys


SAGE_CACHE = "/tmp/ramanujan-aligned-local-sage"
os.environ.setdefault("DOT_SAGE", SAGE_CACHE)
Path(os.environ["DOT_SAGE"]).mkdir(parents=True, exist_ok=True)

try:
    from sage.all import GF, PolynomialRing, ZZ, factorial, gcd  # type: ignore[import-not-found]
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


def valuation(value, prime: int) -> int:
    """Return the exponent of ``prime`` in a nonzero integer ``value``."""

    exponent = 0
    remaining = abs(ZZ(value))
    while remaining % prime == 0:
        remaining //= prime
        exponent += 1
    return exponent


def aligned_polynomials(polynomials: list, a: int, b: int, c: int):
    """Return the integer polynomials ``F,G,J`` from the paper."""

    first = polynomials[a]
    if c == a:
        second = me.shifted(acv.central_deflation(polynomials[b], b), a)
    else:
        second = me.shifted(polynomials[b], a)
    third = me.shifted(polynomials[b + c], a)
    return first, second, third


def common_factor_data(first, second, third, prime: int):
    """Return modular gcd degree, squarefree degree, and factor degrees."""

    field = GF(prime)
    ring = PolynomialRing(field, "z")
    reductions = [ring(polynomial) for polynomial in (first, second, third)]
    common = gcd(gcd(reductions[0], reductions[1]), reductions[2]).monic()
    radical = common.radical().monic()
    factors = tuple(
        sorted((int(factor.degree()), int(multiplicity)) for factor, multiplicity in common.factor())
    )
    return int(common.degree()), int(radical.degree()), factors


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-height", type=int, default=20)
    parser.add_argument(
        "--all-supported",
        action="store_true",
        help="also print structurally supported prime factors of the raw content",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    height = args.max_height
    if height < 8:
        raise SystemExit("max-height must be at least 8")

    polynomials = me.build_gap_polynomials(height)
    apery = me.apery_values(height)
    pell = me.pell_values(height + 1)
    carrier = acv.structural_carrier(height, apery, pell)

    records = []
    triples = [triple for triple in acv.gap_triples(height) if len(set(triple)) > 1]
    for index, triple in enumerate(triples, start=1):
        a, b, c = triple
        datum = acv.aligned_datum(polynomials, a, b, c)
        reduced = acv.remove_supported_part(datum.content, carrier)
        factor_source = ZZ(datum.content if args.all_supported else reduced)
        if factor_source == 1:
            continue
        first, second, third = aligned_polynomials(polynomials, a, b, c)
        for prime, _ in factor_source.factor():
            p = int(prime)
            raw_v = valuation(datum.content, p)
            reduced_v = valuation(reduced, p)
            gcd_degree, distinct_degree, factors = common_factor_data(
                first, second, third, p
            )
            assert raw_v >= distinct_degree
            if reduced_v:
                assert gcd_degree > 0
            records.append(
                (
                    a,
                    b,
                    c,
                    p,
                    raw_v,
                    reduced_v,
                    gcd_degree,
                    distinct_degree,
                    factors,
                )
            )
        if index % 250 == 0:
            print(f"computed {index}/{len(triples)} aligned pencils", flush=True)

    print(
        "a b c p v_p(C) v_p(C*) deg_gcd deg_radical factor_degrees_and_multiplicities"
    )
    for record in records:
        print(*record)
    print(
        f"height={height} triples={len(triples)} "
        f"reported_prime_triples={len(records)}"
    )


if __name__ == "__main__":
    main()
