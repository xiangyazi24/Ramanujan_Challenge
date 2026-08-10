#!/usr/bin/env python3
"""Measure exact carrier deflation in separated Apéry resultants.

This is reconnaissance, not an asymptotic proof.  For 2 <= s <= H and
2 <= G <= Gmax it computes both

    R_full = Res_x(N_s(x), N_G(x+s))

and

    R_ctr = Res_x(N_s^circ(x), N_G^circ(x+s)),

where the forced central linear factor is removed for every even index.
The multiplicativity of the resultant is checked through
``R_ctr | R_full``.  The script then removes from ``R_ctr`` every prime
power supported on the structural carrier

    U_H = product_{j<=H} j! b_j V_{j+1}.

Only bit lengths and an integrity digest are printed.
"""

from __future__ import annotations

import argparse
import hashlib
import math
import os
from pathlib import Path
import shutil
import sys


SAGE_CACHE = "/tmp/ramanujan-separated-deflation-sage"
os.environ.setdefault("DOT_SAGE", SAGE_CACHE)
Path(os.environ["DOT_SAGE"]).mkdir(parents=True, exist_ok=True)

try:
    from sage.all import ZZ  # type: ignore[import-not-found]
except ModuleNotFoundError:
    sage = shutil.which("sage")
    if sage is None:
        raise SystemExit("SageMath is required")
    environment = os.environ.copy()
    environment["DOT_SAGE"] = SAGE_CACHE
    os.execvpe(sage, [sage, "-python", __file__, *sys.argv[1:]], environment)


sys.path.insert(0, str(Path(__file__).resolve().parent))
import meso_explore as me  # noqa: E402


def central_deflate(poly, index: int):
    if index % 2:
        return poly
    quotient, remainder = poly.quo_rem(2 * me.X + index + 1)
    assert remainder == 0
    return quotient


def prime_support(value: int) -> set[int]:
    return {int(prime) for prime, _ in ZZ(abs(value)).factor()}


def structural_support(height: int) -> set[int]:
    support = set(me.prime_sieve(height))
    for value in me.apery_values(height)[1:]:
        support.update(prime_support(value))
    for value in me.pell_values(height + 1)[1:]:
        support.update(prime_support(value))
    return support


def remove_supported_part(value: int, support: set[int]) -> int:
    residual = abs(value)
    for prime in support:
        while residual % prime == 0:
            residual //= prime
    return residual


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--height", type=int, default=12)
    parser.add_argument("--bridge-max", type=int, default=36)
    parser.add_argument(
        "--show-diagonal",
        action="store_true",
        help="also print the h=h logarithmic decomposition",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.height < 2 or args.bridge_max < 2:
        raise SystemExit("height and bridge-max must be at least 2")

    polynomials = me.build_gap_polynomials(max(args.height, args.bridge_max))
    support = structural_support(args.height)
    rows: list[tuple[int, int, int, int, int]] = []
    digest_rows: list[str] = []
    diagonal_values: dict[int, tuple[int, int]] = {}

    for span in range(2, args.height + 1):
        first_full = polynomials[span]
        first = central_deflate(polynomials[span], span)
        for bridge in range(2, args.bridge_max + 1):
            second_full = polynomials[bridge]
            second = central_deflate(polynomials[bridge], bridge)
            full_value = abs(
                int(first_full.resultant(me.shifted(second_full, span)))
            )
            center_value = abs(int(first.resultant(me.shifted(second, span))))
            assert full_value and center_value
            assert full_value % center_value == 0
            residual = remove_supported_part(center_value, support)
            full_bits = full_value.bit_length()
            center_bits = center_value.bit_length()
            residual_bits = residual.bit_length()
            rows.append((span, bridge, full_bits, center_bits, residual_bits))
            if span == bridge:
                diagonal_values[span] = (full_value, center_value)
            digest_rows.append(
                f"{span},{bridge},{full_value},{center_value},{residual}"
            )

    total_full = sum(row[2] for row in rows)
    total_center = sum(row[3] for row in rows)
    total_residual = sum(row[4] for row in rows)
    fully_removed = sum(row[4] == 1 for row in rows)
    worst = max(rows, key=lambda row: row[4] / row[3])
    best_center = max(rows, key=lambda row: row[2] - row[3])
    best_carrier = max(rows, key=lambda row: row[3] - row[4])
    digest = hashlib.sha256("\n".join(digest_rows).encode("ascii")).hexdigest()

    print(
        "SEPARATED_DEFLATION"
        f" H={args.height} Gmax={args.bridge_max}"
        f" pairs={len(rows)} carrier_primes={len(support)}"
    )
    print(f"sha256={digest}")
    print(
        f"aggregate full_bits={total_full} center_bits={total_center}"
        f" residual_bits={total_residual}"
        f" center_fraction={total_center / total_full:.9f}"
        f" carrier_fraction={total_residual / total_center:.9f}"
        f" total_fraction={total_residual / total_full:.9f}"
        f" fully_removed={fully_removed}"
    )
    print(
        "largest_residual_fraction"
        f" s={worst[0]} G={worst[1]}"
        f" center_bits={worst[3]} residual_bits={worst[4]}"
        f" fraction={worst[4] / worst[3]:.9f}"
    )
    print(
        "largest_center_saving"
        f" s={best_center[0]} G={best_center[1]}"
        f" full_bits={best_center[2]} center_bits={best_center[3]}"
        f" saved_bits={best_center[2] - best_center[3]}"
    )
    print(
        "largest_carrier_saving"
        f" s={best_carrier[0]} G={best_carrier[1]}"
        f" center_bits={best_carrier[3]} residual_bits={best_carrier[4]}"
        f" saved_bits={best_carrier[3] - best_carrier[4]}"
    )

    for lower in range(2, args.bridge_max + 1, args.height):
        upper = min(args.bridge_max, lower + args.height - 1)
        bucket = [row for row in rows if lower <= row[1] <= upper]
        bucket_full = sum(row[2] for row in bucket)
        bucket_center = sum(row[3] for row in bucket)
        bucket_residual = sum(row[4] for row in bucket)
        print(
            f"bridge_bucket=[{lower},{upper}]"
            f" full_bits={bucket_full} center_bits={bucket_center}"
            f" residual_bits={bucket_residual}"
            f" center_fraction={bucket_center / bucket_full:.9f}"
            f" carrier_fraction={bucket_residual / bucket_center:.9f}"
            f" total_fraction={bucket_residual / bucket_full:.9f}"
        )

    if args.show_diagonal:
        print(
            "diagonal_columns="
            "h,full_ln,formal_lc_ln,center_ln,residual_ln,"
            "residual_over_h2logh,formal_lc_integer_divides"
        )
        for index, (full_value, center_value) in diagonal_values.items():
            polynomial = polynomials[index]
            degree = int(polynomial.degree())
            leading = int(polynomial.leading_coefficient())
            formal_leading_part = leading ** (2 * degree)
            residual = remove_supported_part(
                center_value, structural_support(index)
            )
            print(
                f"diagonal h={index}"
                f" full_ln={math.log(full_value):.9f}"
                f" formal_lc_ln={math.log(formal_leading_part):.9f}"
                f" center_ln={math.log(center_value):.9f}"
                f" residual_ln={math.log(residual):.9f}"
                " residual_over_h2logh="
                f"{math.log(residual) / (index * index * math.log(index)):.9f}"
                " formal_lc_integer_divides="
                f"{int(full_value % formal_leading_part == 0)}"
            )


if __name__ == "__main__":
    main()
