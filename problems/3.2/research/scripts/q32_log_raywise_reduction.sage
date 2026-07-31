#!/usr/bin/env sage
"""Raywise symbolic same-moment reduction for the corrected log route.

The computation never expands ``(1-X^kappa)^L``.  A term is represented
as

    A(X;M,L) X^{-(M-j)kappa}(1-X^kappa)^{L-k}.

The two CT-preserving projectors produce only 0 <= j,k <= 2.  Applying
the stable order-66 L operator can therefore be factored raywise.  With
base order L-4, the remaining fixed Laurent decoration lies in the
smooth ideal for every ray, including the two nilpotent ray classes.
It is divided by

    y-z, x-2z-1, z^2-1/2,

and reduced by logarithmic integration by parts at the SAME moment M.
"""

from collections import defaultdict
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from q32_cartier_packet_audit import polytope_points

S.<MM, LL> = PolynomialRing(QQ)
K = S.fraction_field()
LR.<x, y, z> = LaurentPolynomialRing(K, 3)


def theta(polynomial, variable):
    return variable * polynomial.derivative(variable)


def vector_field(H, polynomial):
    return sum(
        H[index] * theta(polynomial, variable)
        for index, variable in enumerate((x, y, z))
    )


def divergence(H):
    return sum(
        theta(H[index], variable)
        for index, variable in enumerate((x, y, z))
    )


def dot_kappa(H, kappa):
    return sum(K(kappa[index]) * H[index] for index in range(3))


G_ABC_H = (-x - 1, LR.zero(), (x + 1) * (z + 1))
G_ABC_h = 2 * x * z - x + 2 * z
G_BC_H = (
    -(x + 1)
    * (2 * y * z**2 + y * z + 2 * z**2 - 4 * y + 6 * z + 1)
    / x,
    (y + 1)
    * (2 * x * z**2 + 5 * x * z + 2 * z**2 - 2 * x + 5 * z)
    / x,
    -(z + 1)
    * (
        4 * x * y - x * z + 2 * y * z + 2 * z**2
        - 3 * x + 5 * y + 2 * z
    )
    / x,
)
G_BC_h = -(4 * z**3 + 4 * z**2 + x) / x

# Minimum-degree polynomial fields of the three triangular smooth
# generators.  Their weights are exactly g_y,g_x,g_z.
g_y = y - z
g_x = x - 2 * z - 1
g_z = z**2 - K(1) / 2
smooth_generators = (g_y, g_x, g_z)

H_y = (LR.zero(), y + 1, -z - 1)
H_x = (
    x * y * z + x * z + y * z + x + z + 1,
    -x * y * z - x * z + y + 1,
    -x * z + y * z - x + y - z - 1,
)
H_z = (
    -K(1) / 2 * x * y**2 * z - K(1) / 2 * y**2 * z,
    K(1) / 2 * x * y**2 * z
    + K(1) / 2 * x * y * z
    - K(1) / 2 * y**2
    - K(1) / 2 * y * z
    - K(1) / 2 * y
    - K(1) / 2 * z,
    -K(1) / 2 * y**2 * z
    - K(1) / 2 * y**2
    + K(1) / 2 * y * z
    + z**2
    + K(1) / 2 * y
    + K(3) / 2 * z
    + K(1) / 2,
)
weight_fields = (H_y, H_x, H_z)
weight_betas = tuple(divergence(H) for H in weight_fields)
F = (
    (1 + x) * (1 + y) * (1 + z)
    * ((1 + y) * (1 + z) + x * y * z)
)
for generator, H in zip(smooth_generators, weight_fields):
    h = generator + sum(H)
    assert vector_field(H, F) == h * F


def twisted_terms(terms, H, h, kappa):
    """Apply T_{H,M} to a decorated ray-term dictionary."""

    alpha = h - sum(H)
    beta = divergence(H)
    kappa_H = dot_kappa(H, kappa)
    out = defaultdict(LR.zero)
    for (endpoint_drop, order_drop), decoration in terms.items():
        same = (
            vector_field(H, decoration)
            + (
                beta
                + MM * alpha
                - (MM - endpoint_drop) * kappa_H
            )
            * decoration
        )
        out[(endpoint_drop, order_drop)] += same
        out[(endpoint_drop + 1, order_drop + 1)] += (
            -(LL - order_drop) * kappa_H * decoration
        )
    return dict(out)


def apply_projector(terms, H, h, sign, kappa):
    twisted = twisted_terms(terms, H, h, kappa)
    out = defaultdict(LR.zero)
    for key, value in terms.items():
        out[key] += value
    for key, value in twisted.items():
        out[key] += sign * value / (MM + 1)
    return {key: value for key, value in out.items() if value}


PE.<E> = PolynomialRing(QQ)
quadratics = (
    E**2 - 6 * E + 1,
    E**2 + 2 * E - 1,
    E**2 - 4 * E + 2,
    E**2 - QQ(1) / 2,
    E**2 - 3 * E + QQ(7) / 4,
    (E + 1)**2,
    E**2 - 2 * E - 1,
    E**2 - 2 * E + QQ(1) / 2,
    (E - QQ(1) / 2)**2,
    E**2 - 2,
    E**2 - E - QQ(1) / 4,
)
q = prod(quadratics)
stable = q**3
assert stable.degree() == 66


def shift_L(polynomial, amount):
    return polynomial.map_coefficients(
        lambda coefficient: K(coefficient(LL=LL + amount))
    )


PR.<X, Y, Z> = PolynomialRing(K, order="degrevlex")


def divide_smooth_laurent(polynomial):
    """Canonical triangular division after clearing Laurent exponents."""

    exponent_list = list(polynomial.dict())
    clear = tuple(
        max(0, -min(exponent[index] for exponent in exponent_list))
        for index in range(3)
    )
    cleared = sum(
        coefficient
        * X**(exponent[0] + clear[0])
        * Y**(exponent[1] + clear[1])
        * Z**(exponent[2] + clear[2])
        for exponent, coefficient in polynomial.dict().items()
    )
    after_y = PR(cleared(Y=Z))
    quotient_y = PR((cleared - after_y) // (Y - Z))
    after_x = PR(after_y(X=2 * Z + 1))
    quotient_x = PR((after_y - after_x) // (X - 2 * Z - 1))
    quotient_z, remainder = after_x.quo_rem(Z**2 - K(1) / 2)

    clearing = x**clear[0] * y**clear[1] * z**clear[2]

    def back(value):
        return LR(value(X=x, Y=y, Z=z)) / clearing

    quotients = (back(quotient_y), back(quotient_x), back(quotient_z))
    assert polynomial == sum(
        smooth_generators[index] * quotients[index]
        for index in range(3)
    ) + back(remainder)
    return quotients, back(remainder)


BASE_DROP = 4
records = []
total_same = LR.zero()
total_lower = LR.zero()
for kappa in polytope_points(1):
    if kappa == (0, 0, 0):
        continue
    monomial = x**kappa[0] * y**kappa[1] * z**kappa[2]
    m = 1 - monomial
    terms = {(0, 0): LR.one()}
    terms = apply_projector(
        terms, G_ABC_H, G_ABC_h, -1, kappa
    )
    terms = apply_projector(
        terms, G_BC_H, G_BC_h, +1, kappa
    )
    assert max(key[0] for key in terms) <= 2
    assert max(key[1] for key in terms) <= 2

    decoration = LR.zero()
    for shift in range(67):
        coefficient = K(stable[shift])
        for (endpoint_drop, order_drop), value in terms.items():
            decoration += (
                coefficient
                * shift_L(value, shift)
                * monomial**endpoint_drop
                * m**(shift + BASE_DROP - order_drop)
            )
    quotients, remainder = divide_smooth_laurent(decoration)
    assert remainder == 0

    same = LR.zero()
    lower = LR.zero()
    for H, beta, quotient in zip(
        weight_fields, weight_betas, quotients
    ):
        kappa_H = dot_kappa(H, kappa)
        same += (
            vector_field(H, quotient)
            + (beta - MM * kappa_H) * quotient
        )
        lower += -(LL - BASE_DROP) * kappa_H * quotient * monomial
    same = -same / MM
    lower = -lower / MM
    total_same += same
    total_lower += lower
    records.append(
        (
            kappa,
            len(decoration.dict()),
            tuple(len(value.dict()) for value in quotients),
            len(same.dict()),
            len(lower.dict()),
        )
    )

print("Q32_LOG_RAYWISE_REDUCTION=PASS")
print("BASE_ORDERS", "L-4", "L-5")
print("RAY_RECORDS", records)
print("MAX_DECORATION_SUPPORT", max(record[1] for record in records))
print("MAX_QUOTIENT_SUPPORT", max(max(record[2]) for record in records))
print("MAX_SOURCE_SUPPORTS",
      max(record[3] for record in records),
      max(record[4] for record in records))
