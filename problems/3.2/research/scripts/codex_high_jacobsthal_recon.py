#!/usr/bin/env python3
"""Empirical Jacobsthal/apparition scan for the Apéry square-root branches.

The natural finite windows are

    tau:   0 <= j <= (p - 1) / 2,  p mod 24 in {1, 5, 7, 11},
    sigma: 0 <= j <= (p - 3) / 2,  p mod 24 in {13, 17, 19, 23}.

Run without arguments for the summary and with ``--tables`` for the complete
zero-set tables in the two quarter-zero classes.
"""

from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import dataclass
from math import gcd, isqrt


CLASSES = (1, 5, 7, 11, 13, 17, 19, 23)


@dataclass(frozen=True)
class Row:
    p: int
    residue: int
    branch: str
    degree: int
    quarter: int
    zeros: tuple[int, ...]
    rep_23: tuple[tuple[int, int], ...]
    rep_16: tuple[tuple[int, int], ...]


def primes_below(limit: int) -> list[int]:
    sieve = bytearray(b"\1") * limit
    sieve[:2] = b"\0\0"
    for d in range(2, isqrt(limit - 1) + 1):
        if sieve[d]:
            sieve[d * d : limit : d] = b"\0" * (((limit - 1 - d * d) // d) + 1)
    return [p for p in range(5, limit) if sieve[p]]


def branch_for(residue: int) -> str:
    return "tau" if residue in (1, 5, 7, 11) else "sigma"


def coefficients_mod(p: int, branch: str, degree: int) -> list[int]:
    """Return tau or sigma through ``degree`` using its order-two recurrence."""
    if degree == 0:
        return [1]
    linear = 170 if branch == "tau" else 238
    constant = 107 if branch == "tau" else 209
    low_shift = 1 if branch == "tau" else 3
    values = [1, (5 if branch == "tau" else 39) * pow(2, -1, p) % p]
    for n in range(degree - 1):
        middle = 68 * n * n + linear * n + constant
        low = (2 * n + low_shift) ** 2
        numerator = 2 * middle * values[-1] - low * values[-2]
        values.append(numerator * pow(4 * (n + 2) ** 2, -1, p) % p)
    return values


def representations(p: int, a: int, b: int) -> tuple[tuple[int, int], ...]:
    """Enumerate positive solutions of p = a*x^2 + b*y^2."""
    answers = []
    for x in range(1, isqrt((p - 1) // a) + 1):
        remainder = p - a * x * x
        if remainder <= 0 or remainder % b:
            continue
        y = isqrt(remainder // b)
        if y > 0 and b * y * y == remainder:
            answers.append((x, y))
    return tuple(answers)


def scan(limit: int = 3000) -> list[Row]:
    rows = []
    for p in primes_below(limit):
        residue = p % 24
        branch = branch_for(residue)
        degree = (p - 1) // 2 if branch == "tau" else (p - 3) // 2
        quarter = (p - 1) // 4 if p % 4 == 1 else (p - 3) // 4
        values = coefficients_mod(p, branch, degree)
        rows.append(
            Row(
                p,
                residue,
                branch,
                degree,
                quarter,
                tuple(j for j, value in enumerate(values) if value == 0),
                representations(p, 2, 3),
                representations(p, 1, 6),
            )
        )
    return rows


def apery_mod(p: int, degree: int) -> list[int]:
    values = [1, 5]
    for n in range(1, degree):
        numerator = (2 * n + 1) * (17 * n * n + 17 * n + 5) * values[n]
        numerator -= n**3 * values[n - 1]
        quotient, remainder = divmod(numerator, (n + 1) ** 3)
        assert remainder == 0
        values.append(quotient)
    return [value % p for value in values[: degree + 1]]


def direct_sqrt(coefficients: list[int], p: int) -> list[int]:
    root = [1]
    inverse_two = pow(2, -1, p)
    for n in range(1, len(coefficients)):
        convolution = sum(root[j] * root[n - j] for j in range(1, n))
        root.append((coefficients[n] - convolution) * inverse_two % p)
    return root


def direct_div_q(coefficients: list[int], p: int) -> list[int]:
    quotient = []
    for n, value in enumerate(coefficients):
        if n >= 1:
            value += 34 * quotient[n - 1]
        if n >= 2:
            value -= quotient[n - 2]
        quotient.append(value % p)
    return quotient


def cross_check() -> None:
    """Check both fast recurrences against the defining series for small p."""
    for p in primes_below(200):
        for branch in ("tau", "sigma"):
            degree = (p - 3) // 2
            f_coefficients = apery_mod(p, degree)
            radicand = (
                f_coefficients
                if branch == "tau"
                else direct_div_q(f_coefficients, p)
            )
            assert coefficients_mod(p, branch, degree) == direct_sqrt(radicand, p)


def rational_grid_max(rows: list[Row], max_denominator: int = 24) -> int:
    """Largest hit count for floor/ceil/nearest(k*degree/q), excluding endpoints."""
    best = 0
    for q in range(2, max_denominator + 1):
        for k in range(1, q):
            for rounding in ("floor", "ceil", "nearest"):
                hits = 0
                for row in rows:
                    if rounding == "floor":
                        position = k * row.degree // q
                    elif rounding == "ceil":
                        position = (k * row.degree + q - 1) // q
                    else:
                        position = (2 * k * row.degree + q) // (2 * q)
                    hits += position in row.zeros
                best = max(best, hits)
    return best


def best_bounded_rational_fit(
    samples: list[tuple[int, int]], bound: int = 48
) -> int:
    """Best support of ratio == a/b (mod p), |a|,b <= bound."""
    best = 0
    for denominator in range(1, bound + 1):
        for numerator in range(-bound, bound + 1):
            if gcd(abs(numerator), denominator) != 1:
                continue
            hits = sum(
                (ratio * denominator - numerator) % p == 0
                for p, ratio in samples
            )
            best = max(best, hits)
    return best


def print_value_hunt(rows: list[Row]) -> None:
    features = (
        lambda x, y, p: x,
        lambda x, y, p: y,
        lambda x, y, p: x * y,
        lambda x, y, p: x * pow(y, -1, p) % p,
        lambda x, y, p: 4 * x * x,
        lambda x, y, p: 3 * y * y,
    )

    class_5 = [row for row in rows if row.residue == 5 and row.p >= 29]
    best_near_quarter = {}
    for offset in (-2, -1, 1, 2):
        best = 0
        for feature in features:
            samples = []
            for row in class_5:
                x, y = row.rep_23[0]
                value = coefficients_mod(row.p, "tau", row.quarter + 2)[
                    row.quarter + offset
                ]
                scale = feature(x, y, row.p) % row.p
                samples.append((row.p, value * pow(scale, -1, row.p) % row.p))
            best = max(best, best_bounded_rational_fit(samples))
        best_near_quarter[offset] = best
    print(
        "class 5 near-quarter candidate fits "
        "(x,y,xy,x/y,4x^2,3y^2; |num|,den <= 48): "
        + ", ".join(
            f"offset {offset:+}: {hits}/{len(class_5)}"
            for offset, hits in best_near_quarter.items()
        )
    )

    class_1 = [row for row in rows if row.residue == 1]
    best = 0
    for feature in features:
        samples = []
        for row in class_1:
            x, y = row.rep_16[0]
            value = coefficients_mod(row.p, "tau", row.quarter)[row.quarter]
            scale = feature(x, y, row.p) % row.p
            samples.append((row.p, value * pow(scale, -1, row.p) % row.p))
        best = max(best, best_bounded_rational_fit(samples))
    print(
        "class 1 quarter-value candidate fits with the same grid: "
        f"best {best}/{len(class_1)}"
    )


def representation_linear_position_max(
    rows: list[Row], center: str, coefficient_bound: int = 64
) -> int:
    """Best support of center +/- (a*x+b*y) among actual zeros."""
    best = 0
    for a in range(-coefficient_bound, coefficient_bound + 1):
        for b in range(-coefficient_bound, coefficient_bound + 1):
            if a == b == 0:
                continue
            hits = 0
            for row in rows:
                x, y = row.rep_23[0]
                base = row.quarter if center == "quarter" else 0
                displacement = a * x + b * y
                hits += displacement != 0 and (
                    base + displacement in row.zeros or base - displacement in row.zeros
                )
            best = max(best, hits)
    return best


def print_summary(rows: list[Row]) -> None:
    print("class  branch  primes  quarter zeros  primes with zeros  |Y_p| histogram")
    for residue in CLASSES:
        class_rows = [row for row in rows if row.residue == residue]
        histogram = dict(sorted(Counter(len(row.zeros) for row in class_rows).items()))
        quarter_zeros = sum(row.quarter in row.zeros for row in class_rows)
        nonempty = sum(bool(row.zeros) for row in class_rows)
        print(
            f"{residue:>5}  {class_rows[0].branch:>6}  {len(class_rows):>6}"
            f"  {quarter_zeros:>13}  {nonempty:>16}  {histogram}"
        )

    class_13 = [row for row in rows if row.residue == 13]
    class_23 = [row for row in rows if row.residue == 23]
    assert all(
        (row.p - 5) // 8 in row.zeros and (3 * row.p - 7) // 8 in row.zeros
        for row in class_13
    )
    assert all(
        (row.p - 7) // 8 in row.zeros
        and row.quarter in row.zeros
        and (3 * row.p - 5) // 8 in row.zeros
        for row in class_23
    )
    print("\nUniversal sigma octant laws:")
    print(f"  p = 13 mod 24: 2/2 positions for all {len(class_13)} primes")
    print(f"  p = 23 mod 24: 3/3 positions for all {len(class_23)} primes")

    for residue in (11, 17):
        class_rows = [row for row in rows if row.residue == residue]
        print(
            f"class {residue} rational-position grid (denominator <= 24): "
            f"best hit count {rational_grid_max(class_rows)}/{len(class_rows)}"
        )

    class_11 = [row for row in rows if row.residue == 11]
    print(
        "class 11 quarter +/- (a*x+b*y), |a|,|b| <= 64: "
        f"best hit count {representation_linear_position_max(class_11, 'quarter')}"
        f"/{len(class_11)}"
    )
    class_5_extra = [row for row in rows if row.residue == 5 and len(row.zeros) > 1]
    print(
        "class 5 extra-pair displacement a*x+b*y, |a|,|b| <= 64: "
        f"best hit count {representation_linear_position_max(class_5_extra, 'quarter')}"
        f"/{len(class_5_extra)}"
    )
    print_value_hunt(rows)


def print_tables(rows: list[Row]) -> None:
    for residue in (5, 23):
        print(f"\n### p = {residue} (mod 24)")
        print("| p | (x,y) for 2x^2+3y^2 | Y_p |")
        print("|---:|:---:|:---|")
        for row in rows:
            if row.residue != residue:
                continue
            representation = ", ".join(map(str, row.rep_23)) or "none"
            print(f"| {row.p} | {representation} | {list(row.zeros)} |")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tables", action="store_true")
    arguments = parser.parse_args()
    cross_check()
    rows = scan()
    assert len(rows) == 428
    print("Direct-series recurrence cross-check: PASS (all primes 5 <= p < 200)")
    print_summary(rows)
    if arguments.tables:
        print_tables(rows)


if __name__ == "__main__":
    main()
