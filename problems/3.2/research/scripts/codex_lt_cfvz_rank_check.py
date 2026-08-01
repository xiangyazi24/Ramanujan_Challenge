#!/usr/bin/env python3
"""Mechanical guards for the CFVZ hypergeometric-rank adjudication.

The source-level facts (CFVZ, arXiv:2510.23298v1, pp. 2--4 and p. 9) are:
the Franel series h satisfies a second-order equation, is a rational pullback
of 2F1(1/3,2/3;1), and the Apéry series pulls back to (1+x)h(x)^2.

This script checks the two elementary hypergeometric identities used to keep
the ranks straight.  It is a regression/normalization check, not a substitute
for the cited source or for Clausen's identity as a theorem.
"""

from __future__ import annotations

from math import comb

import sympy as sp


Y, X, T = sp.symbols("y x t")


def truncate(poly: sp.Expr, variable: sp.Symbol, order: int) -> sp.Expr:
    """Expand at zero and discard terms of degree >= order."""
    return sp.series(poly, variable, 0, order).removeO().expand()


def hypergeom_coefficient(parameters: tuple[sp.Rational, ...], n: int) -> sp.Rational:
    """Coefficient of z^n in pFq(parameters; all lower parameters 1)."""
    numerator = sp.prod(sp.rf(parameter, n) for parameter in parameters)
    # There are len(parameters)-1 lower parameters, all equal to 1, plus n!.
    denominator = sp.factorial(n) ** len(parameters)
    return sp.cancel(numerator / denominator)


def check_literal_cancellation(order: int = 30) -> None:
    """Check 3F2(1/3,2/3,1;1,1)=2F1(1/3,2/3;1)."""
    for n in range(order):
        three_f_two = hypergeom_coefficient(
            (sp.Rational(1, 3), sp.Rational(2, 3), sp.Rational(1)), n
        )
        two_f_one = hypergeom_coefficient(
            (sp.Rational(1, 3), sp.Rational(2, 3)), n
        )
        assert three_f_two == two_f_one
    print(
        "VERIFIED literal 3F2(1/3,2/3,1;1,1) coefficient cancellation "
        f"through degree {order - 1}: it is the rank-2 2F1 series"
    )


def check_clausen_square(order: int = 24) -> None:
    """Check the Clausen form of the genuine symmetric-square period."""
    two_f_one = sum(
        hypergeom_coefficient(
            (sp.Rational(1, 3), sp.Rational(2, 3)), n
        )
        * Y**n
        for n in range(order)
    )
    clausen_argument = 4 * Y * (1 - Y)
    three_f_two = sum(
        hypergeom_coefficient(
            (sp.Rational(1, 3), sp.Rational(2, 3), sp.Rational(1, 2)), n
        )
        * clausen_argument**n
        for n in range(order)
    )
    assert truncate(two_f_one**2 - three_f_two, Y, order) == 0
    print(
        "VERIFIED Clausen square through degree "
        f"{order - 1}: 2F1(1/3,2/3;1;y)^2 = "
        "3F2(1/3,2/3,1/2;1,1;4y(1-y))"
    )


def apery(n: int) -> int:
    return sum(comb(n, k) ** 2 * comb(n + k, n) ** 2 for k in range(n + 1))


def franel(n: int) -> int:
    return sum(comb(n, k) ** 3 for k in range(n + 1))


def check_cfvz_pullback(order: int = 14) -> None:
    """Check f_A(phi(x))=(1+x)h(x)^2 as a formal series."""
    phi = X * (1 - 8 * X) / (1 + X)
    apery_series = sum(apery(n) * T**n for n in range(order))
    franel_series = sum(franel(n) * X**n for n in range(order))
    difference = apery_series.subs(T, phi) - (1 + X) * franel_series**2
    assert truncate(difference, X, order) == 0
    print(
        "VERIFIED CFVZ Apéry/Franel pullback f_A(phi(x))=(1+x)h(x)^2 "
        f"through degree {order - 1}"
    )


def main() -> None:
    check_literal_cancellation()
    check_clausen_square()
    check_cfvz_pullback()


if __name__ == "__main__":
    main()
