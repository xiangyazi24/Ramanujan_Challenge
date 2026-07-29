#!/usr/bin/env python3
"""Audit the first Euclidean-remainder certificates for the central carrier.

If D divides both the Apéry number A_n and the central binomial B_n, then D
divides the least signed residue of n^d A_n modulo B_n.  Subexponential
residues for a fixed d would therefore prove the desired gcd bound.  This
script measures their exact Archimedean size.
"""

from __future__ import annotations

from math import comb, log

from q32_newton import apery_numbers


SAMPLES = (60, 120, 240, 480, 720, 1000)
MAX_DEGREE = 8


def signed_residue(value: int, modulus: int) -> int:
    residue = value % modulus
    return min(residue, modulus - residue)


def main() -> None:
    apery = apery_numbers(max(SAMPLES))
    for n in SAMPLES:
        central = comb(n, n // 2)
        rates = []
        for degree in range(MAX_DEGREE + 1):
            residue = signed_residue(pow(n, degree, central) * apery[n], central)
            rates.append(log(residue) / n if residue else float("-inf"))
        rendered = " ".join(f"{rate:.9f}" for rate in rates)
        print(
            f"n={n} log_B_over_n={log(central) / n:.9f} "
            f"residue_rates_d0_to_d{MAX_DEGREE}={rendered}"
        )


if __name__ == "__main__":
    main()
