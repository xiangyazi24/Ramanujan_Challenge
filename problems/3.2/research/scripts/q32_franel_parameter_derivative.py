#!/usr/bin/env python3
"""Audit the parameter-derivative reformulation of the Franel tail inverse.

Define

    Phi(alpha,x) = sum_{m>=0} (-1)^m binom(alpha,m) F_m x^m

and evaluate its alpha derivative at an integer J.  If D_J denotes that
derivative and D_J^{<=J} its degree-J truncation, then

    W_J = -D_J/P_J + D_J^{<=J}/P_J,
    P_J = Phi(J,x).

Consequently tail integrality follows if the two Euler derivatives on the
right are integral.  Both are integral for all tested Franel cutoffs.  The
same assertion fails for arbitrary integer sequences, so the audit does not
mistake a formal unit-series fact for the missing theorem.
"""

from __future__ import annotations

from fractions import Fraction
from math import comb

from q32_strehl_gcd import franel_numbers


LIMIT = 180


def derivative_binomial(integer: int, degree: int) -> Fraction:
    """Return d/dalpha binom(alpha,degree) at alpha=integer."""

    if degree == 0:
        return Fraction(0)
    if degree <= integer:
        return Fraction(comb(integer, degree)) * sum(
            Fraction(1, integer - offset) for offset in range(degree)
        )
    return Fraction(
        (-1) ** (degree - integer - 1),
        degree * comb(degree - 1, integer),
    )


def divide_by_unit(
    numerator: list[Fraction], denominator: list[int]
) -> list[Fraction]:
    quotient = [Fraction(0) for _ in numerator]
    for degree in range(len(numerator)):
        quotient[degree] = numerator[degree] - sum(
            denominator[index] * quotient[degree - index]
            for index in range(1, min(degree, len(denominator) - 1) + 1)
        )
    return quotient


def main() -> None:
    franel = franel_numbers(LIMIT)
    for cutoff in range(LIMIT):
        polynomial = [
            (-1) ** degree
            * comb(cutoff, degree)
            * franel[degree]
            for degree in range(cutoff + 1)
        ]
        derivative = [
            Fraction((-1) ** degree * franel[degree])
            * derivative_binomial(cutoff, degree)
            for degree in range(LIMIT + 1)
        ]
        truncated = [
            value if degree <= cutoff else Fraction(0)
            for degree, value in enumerate(derivative)
        ]
        full_quotient = divide_by_unit(derivative, polynomial)
        truncated_quotient = divide_by_unit(truncated, polynomial)
        for degree in range(1, LIMIT + 1):
            assert (degree * full_quotient[degree]).denominator == 1, (
                "full",
                cutoff,
                degree,
                degree * full_quotient[degree],
            )
            assert (
                degree * truncated_quotient[degree]
            ).denominator == 1, (
                "truncated",
                cutoff,
                degree,
                degree * truncated_quotient[degree],
            )
    print(
        "verified integral Franel parameter logarithmic derivatives for "
        f"0 <= J < k <= {LIMIT}"
    )


if __name__ == "__main__":
    main()
