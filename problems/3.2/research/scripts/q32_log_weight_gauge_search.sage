#!/usr/bin/env sage
"""Search polynomial tangent gauges that remove logarithmic self-loops.

For each smooth generator

    g_y = y-z, g_x = x-2z-1, g_z = z^2-1/2,

we already have a logarithmic field of that weight.  Adding a tangent
field K with

    K_x q_x + K_y q_y + K_z q_z = 0,
    q_i = theta_i(F)-F,

does not change the weight.  This script solves the exact coefficient
linear system for a bounded-degree polynomial K and asks whether the
constant monomial in the corresponding diagonal Groebner quotient can
cancel the previously observed self-loop.

The search is affine and exact over QQ; it does not first compute a
possibly expensive syzygy-module basis.
"""

from itertools import product

R.<x, y, z> = PolynomialRing(QQ)
F = (
    (1 + x) * (1 + y) * (1 + z)
    * ((1 + y) * (1 + z) + x * y * z)
)


def theta(polynomial, variable):
    return variable * polynomial.derivative(variable)


q = tuple(theta(F, variable) - F for variable in (x, y, z))

A.<aa, bb, cc> = PolynomialRing(QQ)
K = A.fraction_field()
RK.<X, Y, Z> = PolynomialRing(K)
g_y = Y - Z
g_x = X - 2 * Z - 1
g_z = Z**2 - K(1) / 2


def into_rk(polynomial):
    return RK(polynomial(x=X, y=Y, z=Z))


def canonical_quotients(polynomial):
    after_y = RK(polynomial(Y=Z))
    quotient_y = RK((polynomial - after_y) // g_y)
    after_x = RK(after_y(X=2 * Z + 1))
    quotient_x = RK((after_y - after_x) // g_x)
    quotient_z, remainder = after_x.quo_rem(g_z)
    assert polynomial == (
        g_y * quotient_y
        + g_x * quotient_x
        + g_z * quotient_z
        + remainder
    )
    return quotient_y, quotient_x, quotient_z


def monomials_up_to(degree):
    return [
        x**i * y**j * z**k
        for i in range(degree + 1)
        for j in range(degree + 1 - i)
        for k in range(degree + 1 - i - j)
    ]


def objective_signature(coordinate, monomial, target_index):
    """Coefficients of aa,bb,cc,1 in the target diagonal quotient."""

    variables = (x, y, z)
    exponents = (aa, bb, cc)
    term = (
        exponents[coordinate] * into_rk(monomial)
        + into_rk(theta(monomial, variables[coordinate]))
    )
    quotient = canonical_quotients(term)[target_index]
    constant = A(quotient.dict().get((0, 0, 0), 0))
    return (
        constant.monomial_coefficient(aa),
        constant.monomial_coefficient(bb),
        constant.monomial_coefficient(cc),
        constant.monomial_coefficient(A.one()),
    )


# Old diagonal coefficients D_i in M I_i = -D_i I_i - ...
targets = (
    (QQ(0), QQ(1), QQ(0), QQ(1)),       # bb + 1
    (QQ(1), QQ(0), QQ(-1), QQ(1)),      # aa - cc + 1
    (QQ(-1) / 2, QQ(0), QQ(1), QQ(2)),  # (-aa+2cc+4)/2
)


def solve_degree(degree, target_index):
    monomials = monomials_up_to(degree)
    unknowns = [
        (coordinate, monomial)
        for coordinate in range(3)
        for monomial in monomials
    ]

    syzygy_columns = [
        monomial * q[coordinate]
        for coordinate, monomial in unknowns
    ]
    support = sorted(
        set().union(*(set(column.dict()) for column in syzygy_columns))
    )
    signatures = [
        objective_signature(coordinate, monomial, target_index)
        for coordinate, monomial in unknowns
    ]

    rows = [
        [column.dict().get(exponent, 0) for column in syzygy_columns]
        for exponent in support
    ]
    rhs = [QQ.zero()] * len(rows)
    for signature_index in range(4):
        rows.append([
            signature[signature_index] for signature in signatures
        ])
        rhs.append(-targets[target_index][signature_index])

    matrix_left = matrix(QQ, rows)
    vector_right = vector(QQ, rhs)
    rank = matrix_left.rank()
    augmented_rank = matrix_left.augment(vector_right).rank()
    if rank != augmented_rank:
        return None, len(unknowns), rank

    solution = matrix_left.solve_right(vector_right)
    field = tuple(
        sum(
            solution[index] * monomial
            for index, (entry_coordinate, monomial) in enumerate(unknowns)
            if entry_coordinate == coordinate
        )
        for coordinate in range(3)
    )
    assert sum(field[index] * q[index] for index in range(3)) == 0
    signature = [QQ.zero()] * 4
    for index, (coordinate, monomial) in enumerate(unknowns):
        if solution[index]:
            entry = objective_signature(coordinate, monomial, target_index)
            for signature_index in range(4):
                signature[signature_index] += (
                    solution[index] * entry[signature_index]
                )
    assert tuple(signature[index] + targets[target_index][index]
                 for index in range(4)) == (0, 0, 0, 0)
    return field, len(unknowns), rank


for target_index, name in enumerate(("y", "x", "z")):
    print("TARGET", name)
    found = False
    for degree in range(9):
        field, unknown_count, rank = solve_degree(degree, target_index)
        print(
            " DEGREE", degree,
            "UNKNOWNS", unknown_count,
            "RANK", rank,
            "EXISTS", field is not None,
        )
        if field is not None:
            print(" FIELD", field)
            found = True
            break
    print(" FOUND", found)

print("Q32_LOG_WEIGHT_GAUGE_SEARCH=PASS")
