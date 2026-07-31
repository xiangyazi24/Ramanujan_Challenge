#!/usr/bin/env python3
"""Numerically test the corrected smooth-quotient L operator on CTs.

The stable quotient operator is q(E)^3, where q is the product of the
eleven non-nilpotent degree-two characteristic polynomials.  This probe
checks whether it already annihilates the terminal CT sequence (it need
not: ideal lifting predicts a moment-(M+1) source).
"""

from fractions import Fraction
from math import comb

from q32_cartier_packet_audit import shell_batch


def multiply(left, right):
    out = [Fraction(0)] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            out[i + j] += a * b
    return out


# Coefficients are in increasing powers of E.
QUADRATICS = (
    (1, -6, 1),
    (-1, 2, 1),
    (2, -4, 1),
    (Fraction(-1, 2), 0, 1),
    (Fraction(7, 4), -3, 1),
    (1, 2, 1),
    (-1, -2, 1),
    (Fraction(1, 2), -2, 1),
    (Fraction(1, 4), -1, 1),
    (-2, 0, 1),
    (Fraction(-1, 4), -1, 1),
)

q = [Fraction(1)]
for factor in QUADRATICS:
    q = multiply(q, tuple(Fraction(c) for c in factor))
operator = multiply(multiply(q, q), q)
assert len(q) - 1 == 22
assert len(operator) - 1 == 66


def terminal_values(moment, max_order):
    shell = shell_batch(
        moment, range(moment - max_order, moment + 1)
    )
    return [
        sum(
            (-1) ** residue
            * comb(order, residue)
            * shell[moment - residue]
            for residue in range(order + 1)
        )
        for order in range(max_order + 1)
    ]


for moment in (140, 150):
    values = terminal_values(moment, 70)
    for start in range(5):
        residual = sum(
            coefficient * values[start + shift]
            for shift, coefficient in enumerate(operator)
        )
        print(
            "moment", moment,
            "start", start,
            "zero", residual == 0,
            "bits", abs(residual.numerator).bit_length(),
            "denominator", residual.denominator,
        )
