#!/usr/bin/env python3
"""Certify the Franel tail-lattice inverse for cutoffs J=1 and J=2.

Let F(x) be the ordinary generating function of the Franel numbers.  It
satisfies

  2(1+4x)F + (-1+14x+24x^2)F'
    + x(1+x)(-1+8x)F'' = 0.

For a fixed cutoff J, the candidate inverse moments z_k satisfy an ODE
L_J Z_J = F - sum_{k=0}^J F_k x^k.  This script verifies, symbolically modulo
the displayed Franel ODE, explicit rational expressions for Z_1 and Z_2.
Their denominators have constant term one, and the only scalar denominator
in the J=2 numerator is removed by the elementary evenness of F_k for k>0.
Consequently both inverse-moment sequences are integral.
"""

from __future__ import annotations

import sympy as sp


x = sp.symbols("x")
franel_initial = (1, 2, 10)
ode_denominator = x * (1 + x) * (8 * x - 1)
second_f = -2 * (1 + 4 * x) / ode_denominator
second_fp = -(-1 + 14 * x + 24 * x**2) / ode_denominator


Triple = tuple[sp.Expr, sp.Expr, sp.Expr]


def add(left: Triple, right: Triple) -> Triple:
    return tuple(
        sp.factor(left[index] + right[index]) for index in range(3)
    )  # type: ignore[return-value]


def scale(scalar: sp.Expr, value: Triple) -> Triple:
    return tuple(
        sp.factor(scalar * coordinate) for coordinate in value
    )  # type: ignore[return-value]


def theta(value: Triple) -> Triple:
    """Apply x*d/dx modulo the Franel differential equation."""

    coefficient_f, coefficient_fp, rational = value
    return (
        sp.factor(
            x
            * (
                sp.diff(coefficient_f, x)
                + coefficient_fp * second_f
            )
        ),
        sp.factor(
            x
            * (
                coefficient_f
                + sp.diff(coefficient_fp, x)
                + coefficient_fp * second_fp
            )
        ),
        sp.factor(x * sp.diff(rational, x)),
    )


def shifted_theta(shift: int, value: Triple) -> Triple:
    return add(theta(value), scale(sp.Integer(shift), value))


def compose_shifts(shifts: tuple[int, ...], value: Triple) -> Triple:
    for shift in shifts:
        value = shifted_theta(shift, value)
    return value


def assert_equal(left: Triple, right: Triple) -> None:
    assert all(
        sp.cancel(left[index] - right[index]) == 0 for index in range(3)
    )


def main() -> None:
    # J=1.  P_1=1-2x and Z_1=(A_1 F+B_1 F'+C_1)/P_1^2.
    polynomial_1 = 1 - 2 * x
    solution_1 = (
        (1 + x) * (1 - 8 * x) / polynomial_1**2,
        x * (1 + x) * (1 - 8 * x) / polynomial_1**2,
        (-1 + 3 * x - 4 * x**2) / polynomial_1**2,
    )
    operator_1 = add(
        scale(-1, shifted_theta(-1, solution_1)),
        scale(2 * x, shifted_theta(1, solution_1)),
    )
    assert_equal(operator_1, (sp.Integer(1), sp.Integer(0), -1 - 2 * x))

    # J=2.  P_2=1-4x+10x^2 and
    # Z_2=(A_2 F+B_2 F'+C_2)/P_2^2.
    polynomial_2 = 1 - 4 * x + 10 * x**2
    solution_2 = (
        (
            (1 + x)
            * (8 * x - 1)
            * (4 * x**2 + 6 * x - 1)
            / polynomial_2**2
        ),
        (
            x
            * (1 + x)
            * (8 * x - 1)
            * (8 * x**2 + 14 * x - 3)
            / (2 * polynomial_2**2)
        ),
        (
            -100 * x**4
            + 26 * x**3
            - 9 * x**2
            + 8 * x
            - 1
        )
        / polynomial_2**2,
    )
    operator_2 = add(
        scale(
            sp.Rational(1, 2),
            compose_shifts((-2, -1), solution_2),
        ),
        add(
            scale(
                -2 * x,
                compose_shifts((-1, 1), solution_2),
            ),
            scale(
                5 * x**2,
                compose_shifts((1, 2), solution_2),
            ),
        ),
    )
    assert_equal(
        operator_2,
        (sp.Integer(1), sp.Integer(0), -1 - 2 * x - 10 * x**2),
    )
    print("certified Franel tail inverse OGFs for J=1 and J=2")


if __name__ == "__main__":
    main()
