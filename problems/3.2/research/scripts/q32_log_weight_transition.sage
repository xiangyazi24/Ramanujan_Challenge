#!/usr/bin/env sage
"""Exact monomial transition for same-moment logarithmic reduction.

Let I_i(a,b,c)=CT[Lambda^M g_i x^a y^b z^c] for

    g_y=y-z, g_x=x-2z-1, g_z=z^2-1/2.

For a logarithmic field of weight g_i,

    M I_i = -CT[Lambda^M (V_i+beta_i)(x^a y^b z^c)].

This script canonically divides the small coefficient

    a H_x + b H_y + c H_z + beta

by (g_y,g_x,g_z), exposing every shifted ideal state and the
two-dimensional remainder in span(1,z).
"""

S.<aa, bb, cc, MM> = PolynomialRing(QQ)
K = S.fraction_field()
R.<x, y, z> = PolynomialRing(K, order="degrevlex")

g_y = y - z
g_x = x - 2 * z - 1
g_z = z**2 - K(1) / 2
generators = (g_y, g_x, g_z)

H_y = (
    R.zero(),
    y + 1,
    -z - 1,
)
beta_y = y - z

H_x = (
    x * y * z + x * z + y * z + x + z + 1,
    -x * y * z - x * z + y + 1,
    -x * z + y * z - x + y - z - 1,
)
beta_x = y * z + x + y - z

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
beta_z = (
    K(1) / 2 * x * y**2 * z
    + K(1) / 2 * x * y * z
    - K(1) / 2 * y**2 * z
    - y**2
    + 2 * z**2
    - K(1) / 2 * y
    + K(3) / 2 * z
)

fields = (H_y, H_x, H_z)
betas = (beta_y, beta_x, beta_z)


def canonical_division(polynomial):
    after_y = R(polynomial(y=z))
    quotient_y = R((polynomial - after_y) // g_y)
    after_x = R(after_y(x=2 * z + 1))
    quotient_x = R((after_y - after_x) // g_x)
    quotient_z, remainder = after_x.quo_rem(g_z)
    assert polynomial == (
        g_y * quotient_y
        + g_x * quotient_x
        + g_z * quotient_z
        + remainder
    )
    assert remainder.degree(x) <= 0 and remainder.degree(y) <= 0
    assert remainder.degree(z) <= 1
    return (quotient_y, quotient_x, quotient_z), remainder


def shift_ledger(polynomial):
    return sorted(
        (
            tuple(int(entry) for entry in exponent),
            factor(coefficient),
        )
        for exponent, coefficient in polynomial.dict().items()
    )


all_transitions = []
for field_index, (H, beta) in enumerate(zip(fields, betas)):
    coefficient = aa * H[0] + bb * H[1] + cc * H[2] + beta
    quotients, remainder = canonical_division(coefficient)
    print("FIELD", field_index, ("y", "x", "z")[field_index])
    print("COEFFICIENT", coefficient)
    for target_index, quotient in enumerate(quotients):
        ledger = shift_ledger(quotient)
        print(
            " TO", target_index, ("y", "x", "z")[target_index],
            ledger,
        )
        all_transitions.extend(
            (field_index, target_index, shift, value)
            for shift, value in ledger
        )
    print(" REMAINDER", shift_ledger(remainder))

# A scalar monomial order cannot be strictly decreased: each diagonal
# transition has a zero shift with a generically nonzero coefficient.
for index in range(3):
    zero_shifts = [
        coefficient
        for source, target, shift, coefficient in all_transitions
        if source == index and target == index and shift == (0, 0, 0)
    ]
    print("DIAGONAL_ZERO_SHIFT", index, zero_shifts)
    assert zero_shifts

print("Q32_LOG_WEIGHT_TRANSITION=PASS")
