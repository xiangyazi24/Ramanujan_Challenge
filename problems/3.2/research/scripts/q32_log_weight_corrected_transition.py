#!/usr/bin/env python3
"""Exact transition graph after cancelling the three diagonal log loops.

The three prescribed logarithmic fields have weights

    g_y=y-z,  g_x=x-2z-1,  g_z=z^2-1/2.

Weight-zero Koszul fields may be added without changing those weights.
The explicit gauges discovered in ``q32_log_weight_pair_gauge.py``
cancel the constant self coefficient in every channel.  This script
computes the *full* canonical quotient transition after those changes,
including Laurent z-shifts, and prints the directed support graph used
for a termination/cycle audit.
"""

from collections import defaultdict

import sympy as sp

x, y, z = sp.symbols("x y z")
a, b, c = sp.symbols("a b c")

g_y = y - z
g_x = x - 2 * z - 1
g_z = z**2 - sp.Rational(1, 2)
weights = (g_y, g_x, g_z)

H_y = (0, y + 1, -z - 1)
H_x = (
    x * y * z + x * z + y * z + x + z + 1,
    -x * y * z - x * z + y + 1,
    -x * z + y * z - x + y - z - 1,
)
H_z = (
    -sp.Rational(1, 2) * x * y**2 * z
    - sp.Rational(1, 2) * y**2 * z,
    sp.Rational(1, 2) * x * y**2 * z
    + sp.Rational(1, 2) * x * y * z
    - sp.Rational(1, 2) * y**2
    - sp.Rational(1, 2) * y * z
    - sp.Rational(1, 2) * y
    - sp.Rational(1, 2) * z,
    -sp.Rational(1, 2) * y**2 * z
    - sp.Rational(1, 2) * y**2
    + sp.Rational(1, 2) * y * z
    + z**2
    + sp.Rational(1, 2) * y
    + sp.Rational(3, 2) * z
    + sp.Rational(1, 2),
)
fields = (H_y, H_x, H_z)


def theta(polynomial, variable):
    return sp.expand(variable * sp.diff(polynomial, variable))


def scale_field(scalar, field):
    return tuple(sp.expand(scalar * entry) for entry in field)


def add_fields(*summands):
    return tuple(
        sp.expand(sum(field[index] for field in summands))
        for index in range(3)
    )


def subtract_fields(left, right):
    return tuple(sp.expand(left[index] - right[index]) for index in range(3))


pair_indices = ((0, 1), (0, 2), (1, 2))
pair_fields = []
for left, right in pair_indices:
    pair_fields.append(
        subtract_fields(
            scale_field(weights[right], fields[left]),
            scale_field(weights[left], fields[right]),
        )
    )
K_yx, K_yz, K_xz = pair_fields

corrected_fields = (
    add_fields(H_y, K_xz),
    add_fields(
        H_x,
        K_yx,
        scale_field(-2 + 3 * x - x**2, K_xz),
    ),
    add_fields(
        H_z,
        scale_field(-sp.Rational(1, 2), K_yz),
        scale_field(
            -sp.Rational(19, 4) * z**-3
            + 3 * z**-2
            - sp.Rational(1, 2) * z**-1,
            K_xz,
        ),
    ),
)


def laurent_terms(polynomial):
    """Return {(i,j,k): coeff} for a Laurent polynomial."""

    polynomial = sp.expand(polynomial)
    out = defaultdict(lambda: sp.Integer(0))
    for term in sp.Add.make_args(polynomial):
        powers = term.as_powers_dict()
        exponent = (
            int(powers.get(x, 0)),
            int(powers.get(y, 0)),
            int(powers.get(z, 0)),
        )
        coefficient = sp.cancel(
            term / (x**exponent[0] * y**exponent[1] * z**exponent[2])
        )
        assert not coefficient.free_symbols.intersection({x, y, z})
        out[exponent] += coefficient
    return {exponent: sp.factor(coefficient)
            for exponent, coefficient in out.items() if coefficient != 0}


def z_remainder(polynomial):
    """Remainder a+b*z modulo z^2-1/2 for a Laurent polynomial."""

    out = sp.Integer(0)
    for (ix, iy, iz), coefficient in laurent_terms(polynomial).items():
        assert ix == iy == 0
        quotient, parity = divmod(iz, 2)
        out += coefficient * sp.Rational(1, 2) ** quotient * z**parity
    return sp.expand(out)


def canonical_quotients(polynomial):
    after_y = sp.expand(polynomial.subs(y, z))
    quotient_y = sp.cancel((polynomial - after_y) / g_y)
    assert sp.expand(polynomial - after_y - g_y * quotient_y) == 0

    after_x = sp.expand(after_y.subs(x, 2 * z + 1))
    quotient_x = sp.cancel((after_y - after_x) / g_x)
    assert sp.expand(after_y - after_x - g_x * quotient_x) == 0

    remainder = z_remainder(after_x)
    quotient_z = sp.cancel((after_x - remainder) / g_z)
    assert sp.expand(after_x - remainder - g_z * quotient_z) == 0
    return tuple(map(sp.expand, (quotient_y, quotient_x, quotient_z))), remainder


# The weight identities for the three original fields are checked against
# Lambda in q32_log_weight_pair_gauge.py.  Each correction here is a
# Koszul combination g_j H_i-g_i H_j and therefore has weight zero.


all_edges = []
for source, field in enumerate(corrected_fields):
    beta = sum(
        theta(field[index], (x, y, z)[index])
        for index in range(3)
    )
    coefficient = sp.expand(a * field[0] + b * field[1] + c * field[2] + beta)
    quotients, remainder = canonical_quotients(coefficient)
    print("SOURCE", ("y", "x", "z")[source])
    for target, quotient in enumerate(quotients):
        terms = laurent_terms(quotient)
        constant = sp.factor(terms.pop((0, 0, 0), 0))
        print(
            " TARGET",
            ("y", "x", "z")[target],
            "CONSTANT",
            constant,
            "SHIFTS",
            sorted(terms.items()),
        )
        assert constant == 0 if target == source else True
        for shift, coefficient_value in terms.items():
            all_edges.append((source, target, shift, coefficient_value))
    print(" REMAINDER", sp.factor(remainder))


print("EDGE_COUNT", len(all_edges))
for source, target, shift, coefficient in all_edges:
    print(
        "EDGE",
        ("y", "x", "z")[source],
        "->",
        ("y", "x", "z")[target],
        shift,
        coefficient,
    )
print("Q32_LOG_WEIGHT_CORRECTED_TRANSITION=PASS")
