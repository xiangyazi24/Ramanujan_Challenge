#!/usr/bin/env python3
"""Factor and classify output from q32_pade_family_gcd_fast.cpp.

The C++ scanner prints one whitespace-separated RESULT line per height.  This
postprocessor factors the exact family gcd, separates factors into the ranges
requested in Q889, and records basic tests against A_{3H+1} and the binomial
carrier.
"""

from __future__ import annotations

import argparse
import subprocess
from dataclasses import dataclass
from pathlib import Path

from sympy import factorint


@dataclass(frozen=True)
class Record:
    fields: dict[str, str]

    @property
    def height(self) -> int:
        return int(self.fields["H"])

    @property
    def cutoff(self) -> int:
        return int(self.fields["A"])

    @property
    def gcd(self) -> int:
        return int(self.fields["G"])


def parse_result(line: str) -> Record:
    pieces = line.strip().split()
    assert pieces and pieces[0] == "RESULT", line
    fields: dict[str, str] = {}
    for piece in pieces[1:]:
        key, value = piece.split("=", 1)
        fields[key] = value
    return Record(fields)


def format_factorization(factors: dict[int, int]) -> str:
    if not factors:
        return "1"
    return "*".join(
        str(prime) if exponent == 1 else f"{prime}^{exponent}"
        for prime, exponent in sorted(factors.items())
    )


def range_factorization(
    factors: dict[int, int], lower_exclusive: int | None, upper_inclusive: int | None
) -> str:
    selected = {
        prime: exponent
        for prime, exponent in factors.items()
        if (lower_exclusive is None or lower_exclusive < prime)
        and (upper_inclusive is None or prime <= upper_inclusive)
    }
    return format_factorization(selected)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("binary", type=Path)
    parser.add_argument("heights", nargs="+", type=int)
    args = parser.parse_args()

    completed = subprocess.run(
        [str(args.binary), *(str(height) for height in args.heights)],
        check=True,
        text=True,
        stdout=subprocess.PIPE,
    )

    print(completed.stdout, end="")
    print("\n# FACTORED SUMMARY")
    for line in completed.stdout.splitlines():
        if not line.startswith("RESULT "):
            continue
        record = parse_result(line)
        height = record.height
        n = 3 * height + 1
        value = record.gcd
        factors = {int(p): int(e) for p, e in factorint(value).items()}
        factor_text = format_factorization(factors)
        low = range_factorization(factors, None, height)
        middle = range_factorization(factors, height, 2 * height)
        candidate = range_factorization(factors, 2 * height, n)
        high = range_factorization(factors, n, None)
        gcd_an = int(record.fields["gcd_G_An"])
        gcd_b = int(record.fields["gcd_G_B"])
        core = int(record.fields["direct_core"])
        targets = record.fields["targets"]
        support = record.fields["candidate_support"]
        assert support == record.fields["expected_support"]

        print(
            f"H={height} A={record.cutoff} "
            f"G={value} factor={factor_text} "
            f"range_le_H={low} range_H_to_2H={middle} "
            f"range_2H_to_n={candidate} range_gt_n={high} "
            f"targets={targets} candidate_support={support} "
            f"G_divides_A_n={gcd_an == value} "
            f"G_divides_binomial={gcd_b == value} "
            f"G_equals_core={core == value} "
            f"core_divides_G={value % core == 0}"
        )


if __name__ == "__main__":
    main()
