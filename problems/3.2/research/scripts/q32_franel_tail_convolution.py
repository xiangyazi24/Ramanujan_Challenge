#!/usr/bin/env python3
"""Verify the exact convolution form of the Franel tail inverse.

For a cutoff J, let

    P_J(x) = sum_{i=0}^J (-1)^i binom(J,i) F_i x^i

and let z_d be the recursively defined tail-left-inverse moments.  The
apparently banded, nonconstant-coefficient system is equivalent to

    P_J(x) * sum_{d>J} z_d x^d/d
      = (-1)^J sum_{k>J} F_k x^k/(k binom(k-1,J)).

This finite exact audit checks the coefficient identity both from the
recursion and by direct rational-series division.  The identity itself is
proved algebraically in the working notes; this script guards signs and
index shifts.
"""

from __future__ import annotations

from fractions import Fraction
from math import comb

from q32_franel_tail_lattice import integral_left_inverse_moments
from q32_strehl_gcd import franel_numbers


LIMIT = 160


def main() -> None:
    franel = franel_numbers(LIMIT)
    for cutoff in range(LIMIT):
        moments = integral_left_inverse_moments(cutoff, LIMIT, franel)
        polynomial = [
            (-1) ** index
            * comb(cutoff, index)
            * franel[index]
            for index in range(cutoff + 1)
        ]
        quotient = [Fraction(0) for _ in range(LIMIT + 1)]
        for degree in range(cutoff + 1, LIMIT + 1):
            right = Fraction(
                (-1) ** cutoff * franel[degree],
                degree * comb(degree - 1, cutoff),
            )
            quotient[degree] = right - sum(
                polynomial[index] * quotient[degree - index]
                for index in range(1, cutoff + 1)
            )
            assert degree * quotient[degree] == moments[degree], (
                cutoff,
                degree,
                degree * quotient[degree],
                moments[degree],
            )
    print(
        "verified Franel tail convolution identity for "
        f"0 <= J < k <= {LIMIT}"
    )


if __name__ == "__main__":
    main()
