#!/usr/bin/env python3
"""Falsification certificate for ``CODEX_SPEC_laststand_norun.md``.

The specification identifies the projective point ``(b_n:c_n)`` with a
state propagated by the scalar Apéry companion matrix.  Those are different
objects: the companion matrix propagates ``(u_{n-1},u_n)`` for one solution
``u``.  This verifier first reproduces the advertised symbolic polynomial
``F`` for that *state* problem, then checks the claimed action and collapse on
the actual two-solution orbit.  Both checks fail, the latter at an explicit
live two-run modulo 997.

The final gate intentionally exits nonzero.  Continuing to the requested
mixed-branch enumeration and forty-prime survey would certify a theorem about
the wrong dynamical object.
"""

from __future__ import annotations

from fractions import Fraction
import math
import sys

import sympy as sp


EXPECTED_G_TERMS = (
    ((8, 0), 108),
    ((7, 1), 432), ((7, 0), 864),
    ((6, 2), 648), ((6, 1), 3024), ((6, 0), 2763),
    ((5, 3), 432), ((5, 2), 3888), ((5, 1), 8289), ((5, 0), 4482),
    ((4, 4), 108), ((4, 3), 2160), ((4, 2), 8826),
    ((4, 1), 11205), ((4, 0), 3849),
    ((3, 4), 432), ((3, 3), 3837), ((3, 2), 9384),
    ((3, 1), 7698), ((3, 0), 1644),
    ((2, 4), 537), ((2, 3), 2871), ((2, 2), 4719),
    ((2, 1), 2466), ((2, 0), 276),
    ((1, 4), 210), ((1, 3), 870), ((1, 2), 1038), ((1, 1), 276),
    ((0, 4), 32), ((0, 3), 108), ((0, 2), 92),
)


def apery_P(n):
    return 34 * n**3 + 51 * n**2 + 27 * n + 5


def symbolic_state_certificate() -> None:
    """Reproduce the advertised F, but only for a companion-state orbit."""

    r, d = sp.symbols("r d")
    a = lambda x: apery_P(x) / (x + 1) ** 3
    be = lambda x: x**3 / (x + 1) ** 3

    # v=(A:B) is obtained from det(M_r v,M_{r+d} v)=0 when B != 0.
    A = a(r) - a(r + d)
    B = be(r) - be(r + d)
    A1 = a(r + 1) - a(r + d + 1)
    B1 = be(r + 1) - be(r + d + 1)
    state_run_identity = B * B1 - (-be(r) * A + a(r) * B) * A1
    numerator, denominator = sp.cancel(state_run_identity).as_numer_denom()
    G_expr = sp.cancel(numerator / (-24 * d**2))
    G = sp.Poly(G_expr, r, d)
    terms = tuple((monomial, int(coefficient)) for monomial, coefficient in G.terms())

    assert terms == EXPECTED_G_TERMS
    assert sp.Poly(numerator, r, d).degree(r) == 8
    assert sp.Poly(numerator, r, d).degree(d) == 6
    assert sp.Poly(numerator, r, d).total_degree() == 10
    assert sp.factor(denominator) == (
        (r + 1) ** 3 * (r + 2) ** 3
        * (d + r + 1) ** 3 * (d + r + 2) ** 3
    )
    assert numerator != 0

    print("SYMBOLIC STATE CERTIFICATE: PASS")
    print("F numerator = -24*d^2*G; degrees (r,d,total) = (8,6,10)")
    print("G exact coefficient table ((r-degree,d-degree), coefficient):")
    print(terms)
    print("factor(F numerator) =", sp.factor(numerator))


def symbolic_action_counterexample() -> None:
    """Disprove xi_(n+1)=M_n xi_n already over Q at n=1."""

    # b_0=1,b_1=5 and c_0=0,c_1=6.  At n=1 the recurrence gives
    # b_2=73 and c_2=351/4.
    b1, c1 = Fraction(5), Fraction(6)
    b2, c2 = Fraction(73), Fraction(351, 4)
    alpha1, beta1 = Fraction(117, 8), Fraction(1, 8)
    image = (c1, -beta1 * b1 + alpha1 * c1)
    determinant = image[0] * c2 - image[1] * b2

    assert image == (Fraction(6), Fraction(697, 8))
    assert determinant == Fraction(-46669, 8)
    assert determinant != 0
    print("FOUNDATIONAL ACTION GATE: FAIL")
    print("M_1*(b_1,c_1) =", image, "while (b_2,c_2) =", (b2, c2))
    print("projective determinant =", determinant, "!= 0")


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    limit = math.isqrt(n)
    return all(n % q for q in range(3, limit + 1, 2))


def apery_orbit(p: int) -> tuple[list[int], list[int]]:
    """Compute the two solutions with c_1=6 on 0,...,p-2."""

    N = p - 2
    b = [0] * (N + 1)
    c = [0] * (N + 1)
    b[0], b[1] = 1, 5 % p
    c[0], c[1] = 0, 6 % p
    for n in range(1, N):
        inverse_cube = pow(n + 1, -3, p)
        b[n + 1] = (apery_P(n) * b[n] - n**3 * b[n - 1]) * inverse_cube % p
        c[n + 1] = (apery_P(n) * c[n] - n**3 * c[n - 1]) * inverse_cube % p
    return b, c


def projective_key(x: int, y: int, p: int) -> tuple[int, int]:
    if y % p:
        return (x * pow(y, -1, p) % p, 1)
    assert x % p
    return (1, 0)


def mod_alpha(n: int, p: int) -> int:
    return apery_P(n) * pow(n + 1, -3, p) % p


def mod_beta(n: int, p: int) -> int:
    return n**3 * pow(n + 1, -3, p) % p


def casoratian_and_degeneracy_certificates() -> None:
    """Check the two independent algebraic claims that remain valid."""

    r, d, z = sp.symbols("r d z")
    beta_cross = sp.expand(r**3 * (r + d + 1) ** 3 - (r + d) ** 3 * (r + 1) ** 3)
    linear_cross = r * (r + d + 1) - z * (r + d) * (r + 1)
    assert sp.expand(linear_cross.subs(z, 1)) == -d

    # Use the exact polynomial identity X^3-Y^3=(X-Y)(X^2+XY+Y^2),
    # which is the product over z^3=1 without representing the roots.
    X = r * (r + d + 1)
    Y = (r + d) * (r + 1)
    assert sp.expand(beta_cross - (X - Y) * (X**2 + X * Y + Y**2)) == 0
    assert sp.Poly(linear_cross, r).degree() == 2

    alpha_difference_numerator = sp.cancel(
        apery_P(r) / (r + 1) ** 3
        - apery_P(r + d) / (r + d + 1) ** 3
    ).as_numer_denom()[0]
    alpha_polynomial = sp.Poly(alpha_difference_numerator, r)
    assert alpha_polynomial.degree(r) == 4
    assert alpha_polynomial.LC() == -51 * d

    # W_n=det((b_n,c_n),(b_(n+1),c_(n+1))) obeys
    # W_n=beta_n W_(n-1), W_0=6, hence W_n=6/(n+1)^3.
    for n in range(1, 21):
        product = sp.prod(sp.Rational(j, j + 1) ** 3 for j in range(1, n + 1))
        assert sp.cancel(6 * product) == sp.Rational(6, (n + 1) ** 3)
    print("GAP-ONE CASORATIAN: PASS; W_r = 6/(r+1)^3")
    print("BETA DEGENERACY: PASS; z=1 gives d=0, each nontrivial z^3=1 gives <=2 r")
    print("ALPHA DEGENERACY: PASS; characteristic-zero r-degree 4, leading coefficient -51*d")
    print("uniform simultaneous alpha/beta degeneracy bound: <=4 r per d for p>17")


def live_collapse_counterexample() -> None:
    """Verify an actual adjacent pair in Z_182 that violates xi=v."""

    p, d = 997, 182
    assert is_prime(p)
    b, c = apery_orbit(p)

    # Independently check the exact gap-one Casoratian at every legal index.
    for n in range(p - 2):
        W = (b[n] * c[n + 1] - b[n + 1] * c[n]) % p
        assert W == 6 * pow(n + 1, -3, p) % p

    witnesses = []
    for r in (248, 565):
        assert r + d + 1 <= p - 2
        for j in (0, 1):
            assert projective_key(b[r + j], c[r + j], p) == projective_key(
                b[r + d + j], c[r + d + j], p
            )
        assert c[r] != 0

        A = (mod_alpha(r, p) - mod_alpha(r + d, p)) % p
        B = (mod_beta(r, p) - mod_beta(r + d, p)) % p
        xi = projective_key(b[r], c[r], p)
        predicted = projective_key(A, B, p)
        assert xi != predicted
        witnesses.append((r, b[r], c[r], xi, A, B, predicted))

    expected = [
        (248, 717, 994, (758, 1), 384, 960, (798, 1)),
        (565, 763, 153, (409, 1), 164, 714, (165, 1)),
    ]
    assert witnesses == expected
    print("LIVE COLLAPSE GATE: FAIL")
    for witness in witnesses:
        print("p=997 d=182 witness (r,b_r,c_r,xi,A,B,v) =", witness)


def main() -> int:
    symbolic_state_certificate()
    symbolic_action_counterexample()
    casoratian_and_degeneracy_certificates()
    live_collapse_counterexample()
    print("MIXED-BRANCH ENUMERATION: SKIPPED -- it would use the false action gate")
    print("40-PRIME SURVEY: SKIPPED -- the mandated 2-run prediction is already false")
    print("FINAL GATE: FAIL -- CODEX_SPEC_laststand_norun.md conflates two projective orbits")
    return 1


if __name__ == "__main__":
    sys.exit(main())
