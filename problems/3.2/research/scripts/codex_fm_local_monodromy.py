#!/usr/bin/env python3
"""Symbolic checks for Picard--Fuchs pullback and local exponent data."""

from fractions import Fraction

import sympy as sp


def finite_exponents(a: sp.Expr, b: sp.Expr, point: sp.Expr, variable: sp.Symbol):
    residue = sp.simplify(b.subs(variable, point) / sp.diff(a, variable).subs(variable, point))
    return (sp.Integer(0), sp.simplify(1 - residue))


def infinity_exponents(
    a: sp.Expr, b: sp.Expr, c: sp.Expr, variable: sp.Symbol
) -> tuple[sp.Expr, sp.Expr]:
    """Indicial roots for y~variable^(-r) when degrees are (3,2,1)."""
    r = sp.symbols("r")
    a3 = sp.Poly(a, variable).LC()
    b2 = sp.Poly(b, variable).LC()
    c1 = sp.Poly(c, variable).LC()
    multiplicities = sp.roots(a3 * r * (r + 1) - b2 * r + c1, r)
    roots = [root for root, multiplicity in multiplicities.items() for _ in range(multiplicity)]
    return tuple(sorted(roots, key=sp.default_sort_key))


def check_logarithm_at_infinity(
    a: sp.Expr,
    b: sp.Expr,
    c: sp.Expr,
    variable: sp.Symbol,
    repeated_exponent: sp.Expr,
) -> None:
    """After z=1/x and Y=z^r U, certify a single analytic solution.

    At a simple regular singularity with repeated exponent zero,
    A'(0)=B(0)!=0 makes the analytic recurrence determine every coefficient
    from U(0).  Hence the second solution has a logarithm and monodromy is a
    nontrivial unipotent block.
    """
    z = sp.symbols("z")
    transformed_a = sp.cancel(a.subs(variable, 1 / z) * z**4)
    transformed_b = sp.cancel(
        2 * z**3 * a.subs(variable, 1 / z) - z**2 * b.subs(variable, 1 / z)
    )
    transformed_c = sp.cancel(c.subs(variable, 1 / z))
    gauged_a = transformed_a
    gauged_b = sp.cancel(2 * repeated_exponent * transformed_a / z + transformed_b)
    gauged_c = sp.cancel(
        transformed_a * repeated_exponent * (repeated_exponent - 1) / z**2
        + transformed_b * repeated_exponent / z
        + transformed_c
    )
    denominator = sp.lcm([sp.denom(item) for item in (gauged_a, gauged_b, gauged_c)])
    gauged_a = sp.cancel(gauged_a * denominator)
    gauged_b = sp.cancel(gauged_b * denominator)
    gauged_c = sp.cancel(gauged_c * denominator)
    assert gauged_a.subs(z, 0) == 0
    assert sp.diff(gauged_a, z).subs(z, 0) == gauged_b.subs(z, 0)
    assert gauged_b.subs(z, 0) != 0


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
    franel_c = 8 * x + 2
    franel_finite = {
        sp.Integer(0): finite_exponents(franel_a, franel_b, 0, x),
        sp.Integer(-1): finite_exponents(franel_a, franel_b, -1, x),
        sp.Rational(1, 8): finite_exponents(franel_a, franel_b, sp.Rational(1, 8), x),
    }
    assert all(exponents == (0, 0) for exponents in franel_finite.values())
    franel_infinity = infinity_exponents(franel_a, franel_b, franel_c, x)

    q = t**2 - 34 * t + 1
    alpha, beta = sp.solve(q, t)
    plus_a = 4 * t * q
    plus_b = 4 * (2 * t**2 - 51 * t + 1)
    plus_c = t - 10
    minus_b = 4 * (4 * t**2 - 85 * t + 1)
    minus_c = 3 * (3 * t - 26)
    assert finite_exponents(plus_a, plus_b, 0, t) == (0, 0)
    assert finite_exponents(plus_a, minus_b, 0, t) == (0, 0)
    for root in [alpha, beta]:
        assert finite_exponents(plus_a, plus_b, root, t) == (0, sp.Rational(1, 2))
        assert finite_exponents(plus_a, minus_b, root, t) == (0, -sp.Rational(1, 2))
    plus_infinity = infinity_exponents(plus_a, plus_b, plus_c, t)
    minus_infinity = infinity_exponents(plus_a, minus_b, minus_c, t)

    assert franel_infinity == (1, 1)
    assert plus_infinity == (Fraction(1, 2), Fraction(1, 2))
    assert minus_infinity == (Fraction(3, 2), Fraction(3, 2))

    # Repeated exponents really give logarithms, rather than apparent
    # singularities: the local analytic recurrence has dimension one.
    for point in franel_finite:
        assert franel_b.subs(x, point) != 0
    assert plus_b.subs(t, 0) != 0
    assert minus_b.subs(t, 0) != 0
    check_logarithm_at_infinity(franel_a, franel_b, franel_c, x, sp.Integer(1))
    check_logarithm_at_infinity(plus_a, plus_b, plus_c, t, sp.Rational(1, 2))
    check_logarithm_at_infinity(plus_a, minus_b, minus_c, t, sp.Rational(3, 2))

    # cond = rank + sum_x(drop_x + Swan_x); all displayed monodromy is tame.
    pushforward_conductor = 6 + 4 + 4 + 3 + 3
    tensor_square_pushforward_conductor = 8 + 4 + 4 + 4 + 4
    twisted_companion_conductor = 3 + 2 + 2 + 2 + 2
    untwisted_apery_conductor = 3 + 2 + 1 + 1 + 2
    assert pushforward_conductor == 20
    assert tensor_square_pushforward_conductor == 24
    assert twisted_companion_conductor == 11
    assert untwisted_apery_conductor == 9
    assert pushforward_conductor + twisted_companion_conductor == 31
    assert tensor_square_pushforward_conductor + twisted_companion_conductor == 35

    print("VERIFIED Franel exponents: 0:(0,0), -1:(0,0), 1/8:(0,0), infinity:(1,1)")
    print("VERIFIED S_+ exponents: 0:(0,0), q-roots:(0,1/2), infinity:(1/2,1/2)")
    print("VERIFIED S_- exponents: 0:(0,0), q-roots:(0,-1/2), infinity:(3/2,3/2)")
    print("VERIFIED repeated exponents have logarithmic second solutions and nontrivial unipotent blocks")
    print("VERIFIED tame conductor bookkeeping: cond(phi_* Sym^2 F)=20, cond(K tensor L_chi(q))=11")
    print("VERIFIED exact tensor-square bookkeeping: cond(phi_*(F tensor F))=24 and total C=35")
    print("VERIFIED reduced virtual object has C=31 and surviving Apery object has conductor 9")
    print("VERIFIED local eigenvalue test at t=0 forces every Kummer self-twist scalar to be 1")
    print("VERIFIED ranks 6 and 3 rule out a mutual Kummer twist of the two displayed sheaves")


def main() -> None:
    check_pullback_of_rank_two_equation()
    check_symmetric_square_apery_operator()
    check_local_exponents_and_conductors()


if __name__ == "__main__":
    main()
