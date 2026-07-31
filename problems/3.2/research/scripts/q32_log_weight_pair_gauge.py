#!/usr/bin/env python3
"""Lightweight exact audit of tangent gauges built from weight pairs.

If H_i has logarithmic weight g_i, then

    K_ij = g_j H_i - g_i H_j

has weight zero.  We compute the induced diagonal Groebner-quotient
functional exactly and solve for constant linear combinations of the
three pair gauges which cancel the old self-loops.
"""

import sympy as sp

x, y, z = sp.symbols("x y z")
a, b, c = sp.symbols("a b c")

g_y = y - z
g_x = x - 2 * z - 1
g_z = z**2 - sp.Rational(1, 2)
weights = (g_y, g_x, g_z)

H_y = (
    0,
    y + 1,
    -z - 1,
)
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


def subtract_fields(left, right):
    return tuple(sp.expand(left[i] - right[i]) for i in range(3))


pair_names = ("yx", "yz", "xz")
pair_indices = ((0, 1), (0, 2), (1, 2))
pair_fields = []
for left, right in pair_indices:
    pair_fields.append(
        subtract_fields(
            scale_field(weights[right], fields[left]),
            scale_field(weights[left], fields[right]),
        )
    )


def canonical_quotients(polynomial):
    polynomial = sp.expand(polynomial)
    after_y = sp.expand(polynomial.subs(y, z))
    quotient_y = sp.cancel((polynomial - after_y) / g_y)
    assert sp.expand(polynomial - after_y - g_y * quotient_y) == 0

    after_x = sp.expand(after_y.subs(x, 2 * z + 1))
    quotient_x = sp.cancel((after_y - after_x) / g_x)
    assert sp.expand(after_y - after_x - g_x * quotient_x) == 0

    quotient_z, remainder = sp.div(after_x, g_z, z)
    assert sp.expand(after_x - g_z * quotient_z - remainder) == 0
    return tuple(map(sp.expand, (quotient_y, quotient_x, quotient_z)))


def diagonal_functional(field, target_index):
    beta = sum(
        theta(field[index], (x, y, z)[index])
        for index in range(3)
    )
    coefficient = sp.expand(
        a * field[0] + b * field[1] + c * field[2] + beta
    )
    quotient = canonical_quotients(coefficient)[target_index]
    constant = sp.Poly(quotient, x, y, z).coeff_monomial(1)
    polynomial = sp.Poly(sp.expand(constant), a, b, c)
    return sp.Matrix(
        [
            polynomial.coeff_monomial(a),
            polynomial.coeff_monomial(b),
            polynomial.coeff_monomial(c),
            polynomial.coeff_monomial(1),
        ]
    )


old_diagonals = (
    sp.Matrix([0, 1, 0, 1]),
    sp.Matrix([1, 0, -1, 1]),
    sp.Matrix([-sp.Rational(1, 2), 0, 1, 2]),
)

print("PAIR_GAUGE_DIAGONAL_SIGNATURES")
for target_index, target_name in enumerate(("y", "x", "z")):
    columns = [
        diagonal_functional(field, target_index)
        for field in pair_fields
    ]
    matrix = sp.Matrix.hstack(*columns)
    print("TARGET", target_name)
    for name, column in zip(pair_names, columns):
        print(" ", name, tuple(column))
    print(" RANK", matrix.rank())
    solution = sp.linsolve((matrix, -old_diagonals[target_index]))
    print(" SOLUTION", solution)
    for candidate in solution:
        corrected = old_diagonals[target_index] + matrix * sp.Matrix(candidate)
        assert corrected == sp.zeros(4, 1)
        break


def monomials_up_to(degree):
    return [
        x**i * y**j * z**k
        for i in range(degree + 1)
        for j in range(degree + 1 - i)
        for k in range(degree + 1 - i - j)
    ]


def choose_zero_parameter_solution(solution_set):
    if solution_set is sp.EmptySet:
        return None
    candidate = next(iter(solution_set))
    parameters = set().union(*(entry.free_symbols for entry in candidate))
    parameters.difference_update({x, y, z, a, b, c})
    return tuple(sp.expand(entry.subs({parameter: 0 for parameter in parameters}))
                 for entry in candidate)


print("POLYNOMIAL_MULTIPLIER_SEARCH")
for target_index, target_name in enumerate(("y", "x", "z")):
    print("TARGET", target_name)
    for degree in range(4):
        decorated = [
            (name, monomial, scale_field(monomial, field))
            for name, field in zip(pair_names, pair_fields)
            for monomial in monomials_up_to(degree)
        ]
        columns = [
            diagonal_functional(field, target_index)
            for _, _, field in decorated
        ]
        matrix = sp.Matrix.hstack(*columns)
        solution_set = sp.linsolve(
            (matrix, -old_diagonals[target_index])
        )
        candidate = choose_zero_parameter_solution(solution_set)
        print(
            " DEGREE", degree,
            "COLUMNS", len(columns),
            "RANK", matrix.rank(),
            "EXISTS", candidate is not None,
        )
        if candidate is None:
            if degree == 3:
                left_kernel = matrix.T.nullspace()
                print(
                    " LEFT_INVARIANTS",
                    [tuple(vector) for vector in left_kernel],
                )
                print(
                    " OLD_PAIRINGS",
                    [
                        (vector.T * old_diagonals[target_index])[0]
                        for vector in left_kernel
                    ],
                )
            continue
        corrected = old_diagonals[target_index] + matrix * sp.Matrix(candidate)
        assert corrected == sp.zeros(4, 1)
        nonzero = [
            (decorated[index][0], decorated[index][1], coefficient)
            for index, coefficient in enumerate(candidate)
            if coefficient != 0
        ]
        tangent = tuple(
            sp.expand(
                sum(
                    candidate[index] * decorated[index][2][coordinate]
                    for index in range(len(decorated))
                )
            )
            for coordinate in range(3)
        )
        print(" NONZERO", nonzero)
        print(
            " TANGENT_DEGREES",
            tuple(sp.Poly(entry, x, y, z).total_degree()
                  if entry != 0 else -1 for entry in tangent),
        )
        print(" SIGNATURE", tuple(diagonal_functional(tangent, target_index)))
        break

print("Q32_LOG_WEIGHT_PAIR_GAUGE_AUDIT=PASS")


def laurent_z_quotient_constant(polynomial):
    """Constant of Q in Pbar=(z^2-1/2)Q+R, deg(R)<2.

    Here Pbar is obtained by x=2z+1,y=z and is a finite Laurent
    polynomial in z.  Negative powers contribute no z^0 term to Q;
    z^(2h), h>=1, contributes 2^(1-h).
    """

    restricted = sp.expand(polynomial.subs({x: 2 * z + 1, y: z}))
    out = sp.Rational(0)
    for term in sp.Add.make_args(restricted):
        power = int(term.as_powers_dict().get(z, 0))
        coefficient = sp.cancel(term / z**power)
        assert z not in coefficient.free_symbols
        if power >= 2 and power % 2 == 0:
            out += coefficient * sp.Rational(1, 2)**(power // 2 - 1)
    return sp.simplify(out)


def z_diagonal_laurent_signature(field):
    beta = sum(
        theta(field[index], (x, y, z)[index])
        for index in range(3)
    )
    return sp.Matrix([
        laurent_z_quotient_constant(field[0]),
        laurent_z_quotient_constant(field[1]),
        laurent_z_quotient_constant(field[2]),
        laurent_z_quotient_constant(beta),
    ])


# The polynomial-multiplier invariant is broken by allowed torus Laurent
# gauges.  This explicit combination cancels the z-channel diagonal.
K_yz = pair_fields[1]
K_xz = pair_fields[2]
z_correction = tuple(
    sp.expand(
        -sp.Rational(1, 2) * K_yz[index]
        + (
            -sp.Rational(19, 4) * z**-3
            + 3 * z**-2
            - sp.Rational(1, 2) * z**-1
        ) * K_xz[index]
    )
    for index in range(3)
)
z_signature = z_diagonal_laurent_signature(z_correction)
assert z_signature == sp.Matrix([
    sp.Rational(1, 2), 0, -1, -2
])
assert old_diagonals[2] + z_signature == sp.zeros(4, 1)
print("LAURENT_Z_CORRECTION")
print(
    " COMBINATION",
    "-1/2*K_yz+(-19/4*z^-3+3*z^-2-1/2*z^-1)*K_xz",
)
print(" SIGNATURE", tuple(z_signature))
print("Q32_LOG_WEIGHT_LAURENT_GAUGE=PASS")
