#!/usr/bin/env python3
"""Exact audit of the per-ray split of the seam scalar S_r.

The seam theorem (S1)--(S3) of `/tmp/P32_TERMINAL_CROSS_N_FINAL.md` reduces the
post-selector boundary datum to the pair ``gcd(b_r, S_r)`` with

    S_r = b_r - sum_kappa lambda_kappa CT[ Lambda^{r-1} X^{-(r-1)kappa} (X^kappa-1)^r ].

This script proves and checks the further factorisation

    S_r = b_r - T_r,
    T_r = sum_kappa lambda_kappa U_kappa(r),
    U_kappa(r) = CT[ G_kappa^{r-1} (X^kappa - 1) ],     G_kappa = Lambda * (1 - X^{-kappa}),

which follows from the regrouping

    Lambda^{r-1} X^{-(r-1)kappa} (X^kappa-1)^r
      = (Lambda X^{-kappa})^{r-1} (X^kappa-1)^{r-1} (X^kappa-1)
      = (Lambda (1 - X^{-kappa}))^{r-1} (X^kappa-1).

Consequence: each U_kappa is the constant term of the powers of ONE fixed Laurent
polynomial (times a fixed Laurent polynomial), hence D-finite in r; therefore S_r is
D-finite as a finite Z-linear combination of 21 such pieces plus b_r.  Its own operator
is the LCLM of the pieces, which is why direct guessing on S_r finds nothing (checked:
no operator of order <= 10 and degree <= 20 with data through r = 260), while the
individual pieces have operators of order 4-6 and degree 10-13.

The binomial expansion of G_kappa^{r-1} gives the computationally fast form

    U_kappa(r) = sum_{j=0}^{r-1} (-1)^j C(r-1,j) ( c_{r-1}((j-1)kappa) - c_{r-1}(j kappa) ),
    c_m(eta) = [X^eta] Lambda^m,

which the script also cross-checks against honest Laurent-polynomial arithmetic.

Checks performed
----------------
1. polynomial form vs. binomial form of U_kappa, all rays, 2 <= r <= 7;
2. sum_kappa lambda_kappa U_kappa = b_r - S_r against the reference seam scalar,
   exactly for 2 <= r <= 7 and modulo a 61-bit prime for 2 <= r <= 40.
"""

from collections import defaultdict
from math import comb
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from q32_cartier_packet_audit import LAMBDA, apery, coefficient  # noqa: E402
from q32_terminal_cross_n_audit import seam_scalar  # noqa: E402

PRIME = (1 << 61) - 1


def multiply(left, right):
    out = defaultdict(int)
    for u, a in left.items():
        for v, b in right.items():
            out[(u[0] + v[0], u[1] + v[1], u[2] + v[2])] += a * b
    return {k: v for k, v in out.items() if v}


def power(base, exponent):
    out = {(0, 0, 0): 1}
    for _ in range(exponent):
        out = multiply(out, base)
    return out


def ray_polynomial(ray):
    """G_kappa = Lambda * (1 - X^{-kappa})."""

    out = defaultdict(int)
    for point, weight in LAMBDA.items():
        out[point] += weight
        shifted = (point[0] - ray[0], point[1] - ray[1], point[2] - ray[2])
        out[shifted] -= weight
    return {k: v for k, v in out.items() if v}


def ray_term_polynomial(ray, residue):
    """U_kappa(r) by honest Laurent-polynomial arithmetic."""

    tail = {ray: 1, (0, 0, 0): -1}
    base = ray_polynomial(ray)
    return multiply(power(base, residue - 1), tail).get((0, 0, 0), 0)


def ray_term_binomial(ray, residue, modulus=None):
    """U_kappa(r) by the binomial expansion of G_kappa^{r-1}."""

    moment = residue - 1
    total = 0
    for index in range(moment + 1):
        sign = 1 if index % 2 == 0 else -1
        shifted = tuple((index - 1) * c for c in ray)
        plain = tuple(index * c for c in ray)
        total += sign * comb(moment, index) * (
            coefficient(moment, *shifted, modulus=modulus)
            - coefficient(moment, *plain, modulus=modulus)
        )
        if modulus:
            total %= modulus
    return total if modulus is None else total % modulus


def rays():
    return tuple(point for point in sorted(LAMBDA) if any(point))


def audit_forms(rmax=7):
    checks = 0
    for residue in range(2, rmax + 1):
        for ray in rays():
            assert ray_term_polynomial(ray, residue) == ray_term_binomial(
                ray, residue
            )
            checks += 1
    return checks


def audit_split_exact(rmax=7):
    checks = 0
    for residue in range(2, rmax + 1):
        total = sum(
            LAMBDA[ray] * ray_term_polynomial(ray, residue) for ray in rays()
        )
        assert apery(residue) - total == seam_scalar(residue)
        checks += 1
    return checks


def audit_split_modular(rmax=40, modulus=PRIME):
    checks = 0
    for residue in range(2, rmax + 1):
        total = 0
        for ray in rays():
            total += LAMBDA[ray] * ray_term_binomial(
                ray, residue, modulus=modulus
            )
            total %= modulus
        left = (apery(residue) - total) % modulus
        assert left == seam_scalar(residue) % modulus
        checks += 1
    return checks


if __name__ == "__main__":
    print("RAYS", len(rays()))
    print("FORM_CHECKS", audit_forms())
    print("EXACT_SPLIT_CHECKS", audit_split_exact())
    print("MODULAR_SPLIT_CHECKS", audit_split_modular())
    print("Q32_SEAM_RAY_SPLIT_AUDIT=PASS")
