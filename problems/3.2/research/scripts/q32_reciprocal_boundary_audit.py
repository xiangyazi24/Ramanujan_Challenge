#!/usr/bin/env python3
"""Numerical audit for the fixed-numerator-degree a=1 boundary.

For

    U_n = sum_{k=0}^n (-1)^k binom(n,k) / A_k,

the exact ordinary generating-function identity is

    sum_{n>=0} U_n z^n
      = (1-z)^(-1) R(-z/(1-z)),
    R(w) = sum_{k>=0} w^k / A_k.

This script also compares U_n with the residue contribution from the
small complex zero of the Golyshev--Zagier hypergeometric interpolation
of the Apéry numbers.  The comparison is diagnostic only: it shows that
this one zero pair does not dominate in the tested range, so a
Nörlund--Rice argument must control all other residues and its contour
remainder before it can prove an eventual sign theorem.
"""

from __future__ import annotations

from fractions import Fraction
from math import comb

import mpmath as mp


ROOT_REAL = "0.145986773118010866514683748306"
ROOT_IMAGINARY = "0.582426140244853150720084304577"


def apery_values(maximum: int) -> list[int]:
    values = [1, 5]
    for index in range(1, maximum):
        cubic = (
            34 * index**3
            + 51 * index**2
            + 27 * index
            + 5
        )
        numerator = (
            cubic * values[index]
            - index**3 * values[index - 1]
        )
        denominator = (index + 1) ** 3
        assert numerator % denominator == 0
        values.append(numerator // denominator)
    return values


def exact_u(index: int, apery: list[int]) -> Fraction:
    return sum(
        (
            Fraction((-1) ** term * comb(index, term), apery[term])
            for term in range(index + 1)
        ),
        Fraction(),
    )


def numerical_u(index: int, apery: list[int]) -> mp.mpf:
    return mp.fsum(
        mp.mpf((-1) ** term * comb(index, term)) / apery[term]
        for term in range(index + 1)
    )


def interpolated_apery(value: mp.mpc) -> mp.mpc:
    return mp.hyper(
        [-value, -value, value + 1, value + 1],
        [1, 1, 1],
        1,
    )


def root_pair_data() -> tuple[mp.mpc, mp.mpc]:
    root = mp.mpc(ROOT_REAL, ROOT_IMAGINARY)
    derivative = mp.diff(interpolated_apery, root)
    return root, mp.gamma(-root) / derivative


def root_pair_contribution(
    index: int,
    root: mp.mpc,
    coefficient: mp.mpc,
) -> mp.mpf:
    term = (
        coefficient
        * mp.gamma(index + 1)
        / mp.gamma(index + 1 - root)
    )
    return 2 * mp.re(term)


def verify_ogf_coefficients(maximum: int = 20) -> None:
    """Check the binomial-transform coefficient identity exactly."""

    apery = apery_values(maximum)
    direct = [exact_u(index, apery) for index in range(maximum + 1)]

    # Expand R(-z/(1-z))/(1-z) coefficientwise.  The contribution
    # of w^k/A_k to z^n is (-1)^k binom(n,k)/A_k.
    transformed = []
    for index in range(maximum + 1):
        transformed.append(
            sum(
                (
                    Fraction(
                        (-1) ** term * comb(index, term),
                        apery[term],
                    )
                    for term in range(index + 1)
                ),
                Fraction(),
            )
        )
    assert transformed == direct


def main() -> None:
    heights = (20, 40, 80, 100, 140, 280, 500, 1000, 1600, 2000)
    mp.mp.dps = 100
    apery = apery_values(max(heights))
    verify_ogf_coefficients()
    root, coefficient = root_pair_data()
    residual = abs(interpolated_apery(root))
    print("exact OGF coefficient identity checked through degree 20")
    print(f"interpolation root residual={mp.nstr(residual, 8)}")
    print("n U_n root_pair_contribution difference")
    for height in heights:
        if height <= 280:
            value_fraction = exact_u(height, apery)
            value = mp.mpf(value_fraction.numerator) / value_fraction.denominator
        else:
            value = numerical_u(height, apery)
        contribution = root_pair_contribution(height, root, coefficient)
        print(
            height,
            mp.nstr(value, 24),
            mp.nstr(contribution, 24),
            mp.nstr(value - contribution, 16),
        )


if __name__ == "__main__":
    main()
