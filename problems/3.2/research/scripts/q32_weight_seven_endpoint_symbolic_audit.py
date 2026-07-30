#!/usr/bin/env python3
"""Symbolic audit for the all-m seventh endpoint rank-one proof."""

from __future__ import annotations

import sympy as sp


def main() -> None:
    m, j = sp.symbols("m j")
    a = m - j
    c = m + j

    anchor_x = -2 * j * m**2
    anchor_y = -j * m**2 * (2 * j**2 + m**2 - 2)
    anchor_u = 2 * j**2 * m**4

    kernel = -sp.Rational(4, 5) * j * m**2 * (
        2 * j**2 + m**2 + 8
    )

    direct_x = anchor_x + a**2
    direct_y = anchor_y + a**2 * (3 * a**2 - 3 * a - m**2)
    direct_u = anchor_u - a**2 * m**2 * (1 + 2 * j)
    direct_z = -4 * m**3 * a**2
    direct_c = (
        4
        - sp.Rational(12, 5) * a
        + sp.Rational(12, 5) * a**2
        - sp.Rational(4, 5) * m**2
        + 8 * m**3
    )
    direct_ell = kernel + a**2 * direct_c

    reflected_x = anchor_x + c**2
    reflected_y = anchor_y + c**2 * (3 * c**2 + 3 * c - m**2)
    reflected_u = anchor_u - c**2 * m**2 * (1 + 2 * j)
    reflected_z = 4 * m**3 * c**2
    reflected_c = (
        4
        + sp.Rational(12, 5) * c
        + sp.Rational(12, 5) * c**2
        - sp.Rational(4, 5) * m**2
        - 8 * m**3
    )
    reflected_ell = kernel + c**2 * reflected_c

    assert sp.factor(
        10 * direct_z
        - 4 * direct_y
        - 20 * direct_x
        + 5 * direct_ell
    ) == 0
    assert sp.factor(
        10 * reflected_z
        - 4 * reflected_y
        - 20 * reflected_x
        + 5 * reflected_ell
    ) == 0

    direct_q = (
        direct_y + 5 * direct_u + sp.Rational(5, 4) * direct_ell
    )
    direct_density = sp.Rational(1, 12) * m**3 * (
        60 * m**3
        - 14 * m**2
        - 51
        + (22 * m**2 + 3) * (m - j) ** 2 / (m + j) ** 2
    )
    direct_ratio = (
        (m - j) * (m + j + 1) / (j + 1) ** 2
    ) ** 2
    direct_certificate_polynomial = (
        6 * j**3
        - 15 * j**2 * m**2
        - 18 * j**2 * m
        - 30 * j * m**3
        + 12 * j * m**2
        + 9 * j
        - 15 * m**4
        + 58 * m**3
        + 12 * m
    )
    direct_certificate = (
        j**4
        * direct_certificate_polynomial
        / (3 * (m + j) ** 2)
    )
    assert sp.factor(
        direct_ratio * direct_certificate.subs(j, j + 1)
        - direct_certificate
        - direct_q
        + direct_density
    ) == 0

    reflected_q = (
        reflected_y
        + 5 * reflected_u
        + sp.Rational(5, 4) * reflected_ell
    )
    reflected_density = sp.Rational(1, 12) * m**3 * (
        60 * m**3
        + 14 * m**2
        + 51
        - (22 * m**2 + 3) * (m + j) ** 2 / (m - j) ** 2
    )
    reflected_ratio = (
        (m - 1 - j) * (m + j) / (j + 1) ** 2
    ) ** 2
    reflected_certificate_polynomial = (
        6 * j**3
        - 15 * j**2 * m**2
        + 18 * j**2 * m
        + 30 * j * m**3
        + 12 * j * m**2
        + 9 * j
        - 15 * m**4
        - 58 * m**3
        - 12 * m
    )
    reflected_certificate = (
        j**4
        * reflected_certificate_polynomial
        / (3 * (m - j) ** 2)
    )
    assert sp.factor(
        reflected_ratio * reflected_certificate.subs(j, j + 1)
        - reflected_certificate
        - reflected_q
        + reflected_density
    ) == 0

    reflected_terminal = sp.factor(
        (2 * m - 1) ** 2
        * reflected_certificate_polynomial.subs(j, m)
        / 3
    )
    expected_terminal = sp.factor(
        -m * (2 * m - 1) ** 2 * (22 * m**2 + 3) / 3
    )
    assert reflected_terminal == expected_terminal

    print("termwise_plane_identities=2")
    print("gosper_certificates=2")
    print("terminal_boundary_checks=1")
    print("failures=0")


if __name__ == "__main__":
    main()
