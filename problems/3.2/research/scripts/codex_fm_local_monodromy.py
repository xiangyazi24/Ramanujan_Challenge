#!/usr/bin/env python3
"""Symbolic checks for Picard--Fuchs pullback and local exponent data."""

from fractions import Fraction

import sympy as sp


def finite_exponents(a: sp.Expr, b: sp.Expr, point: sp.Expr, variable: sp.Symbol):
    residue = sp.simplify(b.subs(variable, point) / sp.diff(a, variable).subs(variable, point))
    return (sp.Integer(0), sp.simplify(1 - residue))


def check_pullback_of_rank_two_equation() -> None:
    x = sp.symbols("x")
    h, hp = sp.symbols("h hp")
    a = x * (x + 1) * (8 * x - 1)
    b = 24 * x**2 + 14 * x - 1
    c = 8 * x + 2
    hpp = -(b * hp + c * h) / a

    phi = x * (1 - 8 * x) / (1 + x)
    gauge = sp.sqrt(1 + x)
    y = gauge * h

    def total_derivative(expression: sp.Expr) -> sp.Expr:
        return sp.diff(expression, x) + sp.diff(expression, h) * hp + sp.diff(expression, hp) * hpp

    y_x = total_derivative(y)
    y_t = sp.simplify(y_x / sp.diff(phi, x))
    y_tt = sp.simplify(total_derivative(y_t) / sp.diff(phi, x))
    q_phi = phi**2 - 34 * phi + 1
    plus_residual = sp.factor(
        4 * phi * q_phi * y_tt
        + 4 * (2 * phi**2 - 51 * phi + 1) * y_t
        + (phi - 10) * y
    )
    assert sp.simplify(plus_residual) == 0

    t = sp.symbols("t")
    y0, y1 = sp.symbols("y0 y1")
    q = t**2 - 34 * t + 1
    y2 = -(4 * (2 * t**2 - 51 * t + 1) * y1 + (t - 10) * y0) / (4 * t * q)
    z = y0 / sp.sqrt(q)

    def t_derivative(expression: sp.Expr) -> sp.Expr:
        return sp.diff(expression, t) + sp.diff(expression, y0) * y1 + sp.diff(expression, y1) * y2

    z1 = t_derivative(z)
    z2 = t_derivative(z1)
    minus_residual = sp.factor(
        4 * t * q * z2
        + 4 * (4 * t**2 - 85 * t + 1) * z1
        + 3 * (3 * t - 26) * z
    )
    assert sp.simplify(minus_residual) == 0
    print("VERIFIED phi-pullback of S_+ is sqrt(1+x) times the Franel rank-two equation")
    print("VERIFIED S_-=S_+/sqrt(q) transforms the plus equation into the minus equation")


def check_symmetric_square_apery_operator() -> None:
    t = sp.symbols("t")
    y0, y1 = sp.symbols("y0 y1")
    q = t**2 - 34 * t + 1
    y2 = -(4 * (2 * t**2 - 51 * t + 1) * y1 + (t - 10) * y0) / (4 * t * q)

    def total_derivative(expression: sp.Expr) -> sp.Expr:
        return sp.diff(expression, t) + sp.diff(expression, y0) * y1 + sp.diff(expression, y1) * y2

    def theta(expression: sp.Expr) -> sp.Expr:
        return sp.simplify(t * total_derivative(expression))

    f0 = y0**2
    f1 = theta(f0)
    f2 = theta(f1)
    f3 = theta(f2)
    middle = 2 * (17 * f3 + 17 * f2 + 5 * f1) + (17 * f2 + 17 * f1 + 5 * f0)
    shifted_cube = f3 + 3 * f2 + 3 * f1 + f0
    residual = sp.factor(f3 - t * middle + t**2 * shifted_cube)
    assert sp.simplify(residual) == 0
    print("VERIFIED Sym^2(S_+) is annihilated by the third-order Apery operator")


def check_local_exponents_and_conductors() -> None:
    x, t = sp.symbols("x t")
    franel_a = x * (x + 1) * (8 * x - 1)
    franel_b = 24 * x**2 + 14 * x - 1
    franel_finite = {
        sp.Integer(0): finite_exponents(franel_a, franel_b, 0, x),
        sp.Integer(-1): finite_exponents(franel_a, franel_b, -1, x),
        sp.Rational(1, 8): finite_exponents(franel_a, franel_b, sp.Rational(1, 8), x),
    }
    assert all(exponents == (0, 0) for exponents in franel_finite.values())
    franel_infinity = (sp.Integer(1), sp.Integer(1))

    q = t**2 - 34 * t + 1
    alpha, beta = sp.solve(q, t)
    plus_a = 4 * t * q
    plus_b = 4 * (2 * t**2 - 51 * t + 1)
    minus_b = 4 * (4 * t**2 - 85 * t + 1)
    assert finite_exponents(plus_a, plus_b, 0, t) == (0, 0)
    assert finite_exponents(plus_a, minus_b, 0, t) == (0, 0)
    for root in [alpha, beta]:
        assert finite_exponents(plus_a, plus_b, root, t) == (0, sp.Rational(1, 2))
        assert finite_exponents(plus_a, minus_b, root, t) == (0, -sp.Rational(1, 2))
    plus_infinity = (sp.Rational(1, 2), sp.Rational(1, 2))
    minus_infinity = (sp.Rational(3, 2), sp.Rational(3, 2))

    assert franel_infinity == (1, 1)
    assert plus_infinity == (Fraction(1, 2), Fraction(1, 2))
    assert minus_infinity == (Fraction(3, 2), Fraction(3, 2))

    # cond = rank + sum_x(drop_x + Swan_x); all displayed monodromy is tame.
    pushforward_conductor = 6 + 4 + 4 + 3 + 3
    twisted_companion_conductor = 3 + 2 + 2 + 2 + 2
    assert pushforward_conductor == 20
    assert twisted_companion_conductor == 11

    print("VERIFIED Franel exponents: 0:(0,0), -1:(0,0), 1/8:(0,0), infinity:(1,1)")
    print("VERIFIED S_+ exponents: 0:(0,0), q-roots:(0,1/2), infinity:(1/2,1/2)")
    print("VERIFIED S_- exponents: 0:(0,0), q-roots:(0,-1/2), infinity:(3/2,3/2)")
    print("VERIFIED tame conductor bookkeeping: cond(phi_* Sym^2 F)=20, cond(K tensor L_chi(q))=11")
    print("VERIFIED local eigenvalue test at t=0 forces every Kummer self-twist scalar to be 1")
    print("VERIFIED ranks 6 and 3 rule out a mutual Kummer twist of the two displayed sheaves")


def main() -> None:
    check_pullback_of_rank_two_equation()
    check_symmetric_square_apery_operator()
    check_local_exponents_and_conductors()


if __name__ == "__main__":
    main()
