#!/usr/bin/env python3
"""Exact endpoint-resultant scan for primitive fully-deflated returns.

For ``a,b,c >= 2`` put ``s=a+b+c`` and

    F = N^o_a(x),  G = N^o_b(x+a),  J = N^o_c(x+a+b).

A common root modulo ``p`` forces ``p`` to divide both endpoint resultants
``Res(F,G)`` and ``Res(F,J)``.  Their integer gcd is only a candidate
carrier: the two resultants may vanish at different roots.  This script
uses that inexpensive superset, removes the structural carrier ``U_s``,
and then checks every candidate in the actual finite field.  It filters
center roots and intermediate projective returns and, for manageable
primes, compares the surviving root with the distinguished Apéry zero set.

The endpoint-gcd method is deliberately not identified with the content of
``Res(F,G+T*J)``.  All polynomial and integer arithmetic is exact.  Ordinary
Python automatically re-executes the file under SageMath.

This also covers chains defined only by the raw return set
``{h : N_h(x)=0}``: the continuant congruence for ``N_(u+v)`` modulo
``N_u`` forces the next shifted gap polynomial to vanish whenever the exact
return set makes the intervening multiplier ``N_(u+1)(x)`` nonzero.
"""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from hashlib import sha256
import os
from pathlib import Path
import shutil
import sys


SAGE_CACHE = "/tmp/ramanujan-primitive-fd-candidate-sage"
os.environ.setdefault("DOT_SAGE", SAGE_CACHE)
Path(os.environ["DOT_SAGE"]).mkdir(parents=True, exist_ok=True)

try:
    from sage.all import GF, PolynomialRing, ZZ, gcd  # type: ignore[import-not-found]
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


EXPECTED_DIGESTS = {
    # Exact short-range scan: only primes p > (a+b+c)^2 are classified.
    (36, 6, False, 2_000_000):
        "7279c7675994a3709a2039a0b04073f88e3fca346d08fe94a39f5fa4a3ee76e9",
    (40, 37, False, 2_000_000):
        "d7b607d7b08252a0aecb665bb3faec0ea379be177bb334a10536e8632c57f1b7",
    (44, 41, False, 2_000_000):
        "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855",
    (52, 45, False, 2_000_000):
        "9d4a63bfe29b388a01c5b440daf5bb734184cebb5c1cc7a4e996d33b12f89ec3",
    (60, 53, False, 2_000_000):
        "3ecf3218d6a71c660ea9e182963caf109a968f19ed3b102e11187d16e9cc2eaf",
    (68, 61, False, 2_000_000):
        "7591cd1cd5772a33484cf1cf94d145e4bb9d14ae72e090a0bc58337a8a9ba690",
    (76, 69, False, 2_000_000):
        "ea92e1052ac8eb33c823f1721d46a54731b8dab1765e2948fffb92dd45fb6f8d",
    (84, 77, False, 2_000_000):
        "784a819756e125cf60d09614139cc037fa237d70dbeb1514e694cc7e71aa6afb",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-height", type=int, default=32)
    parser.add_argument(
        "--min-span",
        type=int,
        default=6,
        help="only scan triples with a+b+c at least this value",
    )
    parser.add_argument(
        "--include-small-primes",
        action="store_true",
        help="classify candidate primes p <= (a+b+c)^2 as well",
    )
    parser.add_argument(
        "--apery-limit",
        type=int,
        default=2_000_000,
        help="compute the distinguished Apéry zero set only up to this prime",
    )
    parser.add_argument("--show-records", action="store_true")
    return parser.parse_args()


def centered(root: int, left: int, right: int, prime: int) -> bool:
    return (2 * root + left + right + 1) % prime == 0


def main() -> None:
    args = parse_args()
    maximum = args.max_height
    if maximum < 8:
        raise SystemExit("max-height must be at least 8")
    if not 6 <= args.min_span <= maximum:
        raise SystemExit("min-span must lie between 6 and max-height")

    polynomials = me.build_gap_polynomials(maximum)
    deflated = [None] + [
        acv.central_deflation(polynomials[index], index)
        for index in range(1, maximum + 1)
    ]
    apery = me.apery_values(maximum)
    pell = me.pell_values(maximum + 1)
    carriers = {
        height: acv.structural_carrier(height, apery, pell)
        for height in range(6, maximum + 1)
    }

    first_endpoint_cache = {}
    second_endpoint_cache = {}
    zero_set_cache = {}
    category_counts: Counter[str] = Counter()
    payload = []
    records = []
    candidate_prime_triples = 0
    tested_prime_triples = 0
    maximum_primitive_ratio = None
    maximum_primitive_record = None

    triples = [
        triple
        for triple in acv.gap_triples(maximum)
        if len(set(triple)) > 1 and sum(triple) >= args.min_span
    ]
    for index, (a, b, c) in enumerate(triples, start=1):
        span = a + b + c
        first = deflated[a]

        first_key = (a, b)
        if first_key not in first_endpoint_cache:
            second = me.shifted(deflated[b], a)
            resultant = abs(ZZ(first.resultant(second)))
            assert resultant > 0
            first_endpoint_cache[first_key] = resultant
        first_resultant = first_endpoint_cache[first_key]

        second_key = (a, c, span)
        if second_key not in second_endpoint_cache:
            third = me.shifted(deflated[c], a + b)
            resultant = abs(ZZ(first.resultant(third)))
            assert resultant > 0
            second_endpoint_cache[second_key] = resultant
        second_resultant = second_endpoint_cache[second_key]

        candidate = acv.remove_supported_part(
            gcd(first_resultant, second_resultant), carriers[span]
        )
        if index % 250 == 0:
            print(f"computed {index}/{len(triples)} endpoint-gcd triples", flush=True)
        if candidate == 1:
            continue

        for prime_sage, exponent in ZZ(candidate).factor():
            prime = int(prime_sage)
            candidate_prime_triples += 1
            short_span = prime > span * span
            if not args.include_small_primes and not short_span:
                continue
            tested_prime_triples += 1

            field = GF(prime)
            ring = PolynomialRing(field, "z")
            local = (
                ring(first),
                ring(me.shifted(deflated[b], a)),
                ring(me.shifted(deflated[c], a + b)),
            )
            assert all(polynomial != 0 for polynomial in local)
            common = gcd(gcd(local[0], local[1]), local[2]).monic()
            roots = common.roots(multiplicities=False)
            if not roots:
                category = (
                    "endpoint_false_positive"
                    if common.degree() == 0
                    else "triple_extension_only"
                )
                category_counts[category] += 1
                payload.append(
                    f"{a},{b},{c},{prime},{int(exponent)},none,"
                    f"{category}"
                )
                continue

            for root_field in roots:
                root = int(root_field)
                center_pairs = (
                    centered(root, 0, a, prime),
                    centered(root, a, a + b, prime),
                    centered(root, a + b, span, prime),
                )
                raw_ring = ring
                return_offsets = [0] + [
                    level
                    for level in range(1, span + 1)
                    if raw_ring(polynomials[level])(root_field) == 0
                ]
                selected = [0, a, a + b, span]
                nonwrapping = root + span < prime

                if any(center_pairs):
                    category = "center"
                elif return_offsets != selected:
                    category = "nonprimitive"
                else:
                    if prime <= args.apery_limit:
                        if prime not in zero_set_cache:
                            divided = ppv.apery_zeros_division(prime)
                            cleared = ppv.apery_zeros_cleared(prime)
                            assert divided == cleared
                            zero_set_cache[prime] = divided
                        actual = root in zero_set_cache[prime]
                        if actual and nonwrapping:
                            category = "primitive_actual"
                        elif actual:
                            category = "primitive_actual_wrapping"
                        else:
                            category = "primitive_phantom"
                    else:
                        category = "primitive_actual_unknown"

                category_counts[category] += 1
                record = (
                    a,
                    b,
                    c,
                    prime,
                    int(exponent),
                    root,
                    int(short_span),
                    int(nonwrapping),
                    tuple(return_offsets),
                    category,
                )
                records.append(record)
                if category.startswith("primitive_") and nonwrapping:
                    ratio = Fraction(prime, span * span)
                    if (
                        maximum_primitive_ratio is None
                        or ratio > maximum_primitive_ratio
                    ):
                        maximum_primitive_ratio = ratio
                        maximum_primitive_record = record
                payload.append(
                    ",".join(
                        map(
                            str,
                            (
                                a,
                                b,
                                c,
                                prime,
                                int(exponent),
                                root,
                                int(short_span),
                                int(nonwrapping),
                                ":".join(map(str, return_offsets)),
                                category,
                            ),
                        )
                    )
                )

    digest = sha256("\n".join(payload).encode("ascii")).hexdigest()
    print(
        f"height={maximum} min_span={args.min_span} triples={len(triples)} "
        f"candidate_prime_triples={candidate_prime_triples} "
        f"tested_prime_triples={tested_prime_triples}"
    )
    print(f"categories={sorted(category_counts.items())}")
    print(f"payload_sha256={digest}")
    if maximum_primitive_ratio is None:
        print("max_primitive_p_over_span_squared=none")
    else:
        print(
            "max_primitive_p_over_span_squared="
            f"{maximum_primitive_ratio.numerator}/"
            f"{maximum_primitive_ratio.denominator} "
            f"approx={float(maximum_primitive_ratio):.12f} "
            f"record={maximum_primitive_record}"
        )
    expected = EXPECTED_DIGESTS.get(
        (maximum, args.min_span, args.include_small_primes, args.apery_limit)
    )
    if expected is not None:
        assert digest == expected
    if args.show_records:
        for record in records:
            print(record)


if __name__ == "__main__":
    main()
