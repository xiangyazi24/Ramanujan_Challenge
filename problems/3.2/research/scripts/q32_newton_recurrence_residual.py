#!/usr/bin/env python3
"""Test the quartic recurrence residual of local Newton interpolants.

An interpolating polynomial through consecutive Apéry values satisfies the
Apéry recurrence at every interior node.  Its polynomial recurrence residual
therefore has the complete interior-node product as a factor and a quotient
of degree at most four.  This script checks whether a q=1 hit supplies an
extra factor in that quartic quotient.  It does not in the three reflected
hits at n=321.
"""

from __future__ import annotations

import sympy as sp

from q32_newton import apery_numbers


def main() -> None:
    n = 321
    apery = apery_numbers(n)
    variable = sp.symbols("x")
    recurrence_coefficient = (
        34 * variable**3
        + 51 * variable**2
        + 27 * variable
        + 5
    )
    blocks = (
        (0, 40, 36, 179),
        (41, 40, 64, 193),
        (82, 24, 100, 211),
    )

    for start, degree, hit_index, prime in blocks:
        interpolant = sp.interpolate(
            [
                (start + offset, apery[start + offset])
                for offset in range(degree + 1)
            ],
            variable,
        )
        residual = sp.expand(
            (variable + 1) ** 3
            * interpolant.subs(variable, variable + 1)
            - recurrence_coefficient * interpolant
            + variable**3
            * interpolant.subs(variable, variable - 1)
        )
        universal_factor = sp.prod(
            variable - (start + offset)
            for offset in range(1, degree)
        )
        quotient, remainder = sp.div(
            sp.Poly(residual, variable),
            sp.Poly(universal_factor, variable),
        )
        assert remainder.is_zero
        assert quotient.degree() <= 4

        value = sp.cancel(quotient.as_expr().subs(variable, n))
        numerator, denominator = sp.fraction(value)
        assert denominator % prime
        assert numerator % prime
        print(
            f"start={start} degree={degree} hit={hit_index} "
            f"prime={prime} quotient_degree={quotient.degree()} "
            f"quotient_mod_prime="
            f"{int(numerator % prime) * pow(int(denominator % prime), -1, prime) % prime}"
        )


if __name__ == "__main__":
    main()
