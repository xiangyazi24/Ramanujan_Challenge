#!/usr/bin/env python3
"""Verify the exact polynomial resultant behind the factorial carrier.

For fixed J, B_J(X)=L_(J+1)(X) has roots
0,...,J,-1,...,-J-1, and the truncated Strehl polynomial S_J takes the
values A_0,...,A_J twice on those roots.  Clearing denominators therefore
makes the resultant an explicit product of squared Apéry numbers.  Its
quadratic logarithmic height is the exact obstruction to the naive global
resultant route.
"""

from __future__ import annotations

from math import factorial, prod

from sympy import Poly, QQ, resultant, symbols

from q32_newton import apery_numbers
from q32_strehl_gcd import franel_numbers


MAX_CUTOFF = 8
X = symbols("X")


def legendre_kernel(index: int):
    numerator = prod(X - root for root in range(index))
    numerator *= prod(X + root for root in range(1, index + 1))
    return numerator / factorial(index) ** 2


def main() -> None:
    apery = apery_numbers(MAX_CUTOFF)
    franel = franel_numbers(MAX_CUTOFF)
    for cutoff in range(MAX_CUTOFF + 1):
        strehl = sum(
            legendre_kernel(index) * franel[index]
            for index in range(cutoff + 1)
        )
        carrier = factorial(cutoff + 1) ** 2 * legendre_kernel(cutoff + 1)
        scaled_strehl = factorial(cutoff) ** 2 * strehl
        carrier_poly = Poly(carrier, X, domain=QQ)
        strehl_poly = Poly(scaled_strehl, X, domain=QQ)
        assert carrier_poly.monic() == carrier_poly
        assert all(coefficient.q == 1 for coefficient in strehl_poly.all_coeffs())

        expected = factorial(cutoff) ** (4 * (cutoff + 1))
        expected *= prod(apery[index] ** 2 for index in range(cutoff + 1))
        actual = int(resultant(carrier_poly, strehl_poly))
        assert actual == expected
        print(
            f"J={cutoff} resultant_digits={len(str(abs(actual)))}"
        )


if __name__ == "__main__":
    main()
