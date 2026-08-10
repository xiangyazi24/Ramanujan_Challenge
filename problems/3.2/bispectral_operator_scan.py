#!/usr/bin/env python3
"""Exclude bounded polynomial-coefficient difference eigenoperators.

For the positional-gap continuants ``N_h``, form the exact linear system for

    L = sum_{j=-r}^r A_j(x) T^j,       deg A_j <= d,
    L N_h = lambda_h N_h               (1 <= h <= m).

The system has integer coefficients.  If its reduction modulo one prime has
one-dimensional kernel, then its rational kernel is also one-dimensional;
the scalar identity operator already spans that kernel.  Thus the computation
is a rigorous exclusion for the displayed finite ansatz, not evidence about
arbitrary rational-coefficient or higher-order operators.
"""

from __future__ import annotations

import argparse
import hashlib
import os
from pathlib import Path
import shutil
import sys


SAGE_CACHE = "/tmp/ramanujan-bispectral-sage"
os.environ.setdefault("DOT_SAGE", SAGE_CACHE)
Path(os.environ["DOT_SAGE"]).mkdir(parents=True, exist_ok=True)

try:
    from sage.all import GF, matrix, vector  # type: ignore[import-not-found]
except ModuleNotFoundError:
    sage = shutil.which("sage")
    if sage is None:
        raise SystemExit("SageMath is required")
    environment = os.environ.copy()
    environment["DOT_SAGE"] = SAGE_CACHE
    os.execvpe(sage, [sage, "-python", __file__, *sys.argv[1:]], environment)


sys.path.insert(0, str(Path(__file__).resolve().parent))
import meso_explore as me  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--radius", type=int, default=8)
    parser.add_argument("--degree", type=int, default=30)
    parser.add_argument("--height", type=int, default=20)
    parser.add_argument("--prime", type=int, default=1_000_003)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.radius < 0 or args.degree < 0 or args.height < 1:
        raise SystemExit(
            "radius and degree must be nonnegative; height must be positive"
        )

    field = GF(args.prime)
    polynomial_ring = field["x"]
    variable = polynomial_ring.gen()
    polynomials = [
        polynomial_ring([field(coefficient) for coefficient in polynomial.list()])
        for polynomial in me.build_gap_polynomials(args.height)
    ]

    shifts = list(range(-args.radius, args.radius + 1))
    operator_terms = [
        (shift, degree)
        for shift in shifts
        for degree in range(args.degree + 1)
    ]
    operator_column = {
        term: column for column, term in enumerate(operator_terms)
    }
    lambda_start = len(operator_terms)
    column_count = lambda_start + args.height
    rows = []

    for height in range(1, args.height + 1):
        polynomial = polynomials[height]
        translated = {
            shift: polynomial(variable + shift) for shift in shifts
        }
        for power in range(polynomial.degree() + args.degree + 1):
            row = [field.zero()] * column_count
            for shift, shifted in translated.items():
                lower = max(0, power - shifted.degree())
                upper = min(args.degree, power)
                for degree in range(lower, upper + 1):
                    row[operator_column[shift, degree]] = shifted[power - degree]
            if power <= polynomial.degree():
                row[lambda_start + height - 1] = -polynomial[power]
            rows.append(row)

    system = matrix(field, rows)
    identity = vector(field, column_count)
    identity[operator_column[0, 0]] = 1
    for height in range(args.height):
        identity[lambda_start + height] = 1
    assert system * identity == 0

    digest = hashlib.sha256()
    digest.update(
        f"{args.prime},{args.radius},{args.degree},{args.height}\n".encode("ascii")
    )
    for row in rows:
        digest.update(",".join(str(int(entry)) for entry in row).encode("ascii"))
        digest.update(b"\n")

    rank = int(system.rank())
    nullity = int(system.ncols() - rank)
    print(
        "BISPECTRAL_OPERATOR_SCAN"
        f" prime={args.prime} radius={args.radius} degree={args.degree}"
        f" height={args.height} rows={system.nrows()} cols={system.ncols()}"
        f" rank={rank} nullity={nullity}"
    )
    print(f"sha256={digest.hexdigest()}")
    if nullity == 1:
        print("verdict=SCALAR_IDENTITY_ONLY_OVER_Q")
    else:
        print("verdict=INCONCLUSIVE_NONTRIVIAL_MODULAR_KERNEL")


if __name__ == "__main__":
    main()
