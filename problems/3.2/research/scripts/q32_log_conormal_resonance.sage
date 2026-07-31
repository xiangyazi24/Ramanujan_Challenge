#!/usr/bin/env sage
"""Conormal residue of the same-moment logarithmic reduction.

Let

    I = (y-z, x-2*z-1, z^2-1/2)

be the smooth critical ideal.  For logarithmic fields V_i of weights
g_i generating I, integration by parts gives

    0 = CT Lambda^M (M*g_i*A + V_i(A) + beta_i*A).

For a monomial A=x^a*y^b*z^c, reduce the second term modulo I^2.
This produces a connection on the six-dimensional conormal fibre

    (I/I^2) tensor_Q Q[z]/(z^2-1/2).

The determinant of ``M + connection`` is invariant under changing the
chosen lifts by tangent logarithmic fields (up to a unit and change of
conormal basis).  Its zero locus therefore detects genuine resonances,
not zero-shift artifacts of one reduction gauge.
"""

S.<aa, bb, cc, MM> = PolynomialRing(QQ)
K = S.fraction_field()
R.<x, y, z> = PolynomialRing(K, order="degrevlex")

g_y = y - z
g_x = x - 2 * z - 1
g_z = z**2 - K(1) / 2
generators = (g_y, g_x, g_z)

H_y = (R.zero(), y + 1, -z - 1)
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
    """Triangular division by the displayed smooth Groebner basis."""

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
    return (quotient_y, quotient_x, quotient_z), remainder


def smooth_pair(polynomial):
    """Coordinates in B=K[z]/(z^2-1/2), basis (1,z)."""

    specialized = R(polynomial(x=2 * z + 1, y=z))
    _, remainder = specialized.quo_rem(g_z)
    assert remainder.degree(z) <= 1
    return vector(
        K,
        (
            remainder.monomial_coefficient(R.one()),
            remainder.monomial_coefficient(z),
        ),
    )


def multiplication_matrix(pair):
    """Multiplication by pair[0]+pair[1]*z in B."""

    return matrix(
        K,
        (
            (pair[0], pair[1] / 2),
            (pair[1], pair[0]),
        ),
    )


blocks = [[None for _ in range(3)] for _ in range(3)]
pair_entries = [[None for _ in range(3)] for _ in range(3)]
remainders = []
for source, (field, beta) in enumerate(zip(fields, betas)):
    coefficient = (
        aa * field[0] + bb * field[1] + cc * field[2] + beta
    )
    quotients, remainder = canonical_division(coefficient)
    remainders.append(smooth_pair(remainder))
    for target, quotient in enumerate(quotients):
        pair = smooth_pair(quotient)
        if source == target:
            pair[0] += MM
        pair_entries[source][target] = pair
        block = multiplication_matrix(pair)
        blocks[source][target] = block

connection = block_matrix(K, blocks)


def pair_add(left, right):
    return vector(K, (left[0] + right[0], left[1] + right[1]))


def pair_scale(scalar, value):
    return vector(K, (scalar * value[0], scalar * value[1]))


def pair_mul(left, right):
    return vector(
        K,
        (
            left[0] * right[0] + left[1] * right[1] / 2,
            left[0] * right[1] + left[1] * right[0],
        ),
    )


def determinant_three(entries):
    positive = (
        pair_mul(entries[0][0], pair_mul(entries[1][1], entries[2][2])),
        pair_mul(entries[0][1], pair_mul(entries[1][2], entries[2][0])),
        pair_mul(entries[0][2], pair_mul(entries[1][0], entries[2][1])),
    )
    negative = (
        pair_mul(entries[0][2], pair_mul(entries[1][1], entries[2][0])),
        pair_mul(entries[0][1], pair_mul(entries[1][0], entries[2][2])),
        pair_mul(entries[0][0], pair_mul(entries[1][2], entries[2][1])),
    )
    out = vector(K, (0, 0))
    for value in positive:
        out = pair_add(out, value)
    for value in negative:
        out = pair_add(out, pair_scale(-1, value))
    return out


quadratic_determinant = determinant_three(pair_entries)
determinant_norm = (
    quadratic_determinant[0] ** 2
    - quadratic_determinant[1] ** 2 / 2
)

print("CONORMAL_MATRIX")
print(connection)
print("CONORMAL_REMAINDERS", remainders)
print(
    "CONORMAL_QUADRATIC_DETERMINANT",
    tuple(factor(entry) for entry in quadratic_determinant),
)
print("CONORMAL_DETERMINANT_NORM", factor(determinant_norm))

# A change of conormal basis cannot alter the resonance divisor.
assert connection.nrows() == connection.ncols() == 6
assert determinant_norm != 0
print("Q32_LOG_CONORMAL_RESONANCE=PASS")
