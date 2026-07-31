#!/usr/bin/env sage
"""Search minimum-degree logarithmic fields with prescribed weights.

A field H has logarithmic weight alpha exactly when

    sum_i H_i (theta_i(F)-F) = alpha*F.

This is a finite linear system for polynomial (or bounded Laurent)
coefficients H_i.
"""

from itertools import product

R.<x, y, z> = PolynomialRing(QQ, order="degrevlex")
F = (
    (1 + x) * (1 + y) * (1 + z)
    * ((1 + y) * (1 + z) + x * y * z)
)


def theta(polynomial, variable):
    return variable * polynomial.derivative(variable)


critical = tuple(
    theta(F, variable) - F for variable in (x, y, z)
)
targets = (y - z, x - 2 * z - 1, z**2 - QQ(1) / 2)


def monomials_up_to(degree):
    return [
        x**a * y**b * z**c
        for a in range(degree + 1)
        for b in range(degree + 1 - a)
        for c in range(degree + 1 - a - b)
    ]


def solve_polynomial_field(target, degree):
    monomials = monomials_up_to(degree)
    columns = [
        monomial * critical[coordinate]
        for coordinate in range(3)
        for monomial in monomials
    ]
    rhs_polynomial = target * F
    support = sorted(
        set(rhs_polynomial.dict()).union(
            *(set(column.dict()) for column in columns)
        )
    )
    coefficient_matrix = matrix(
        QQ,
        [
            [column.dict().get(exponent, 0)
             for column in columns]
            for exponent in support
        ],
    )
    rhs = vector(
        QQ,
        [rhs_polynomial.dict().get(exponent, 0)
         for exponent in support],
    )
    if coefficient_matrix.rank() != coefficient_matrix.augment(rhs).rank():
        return None
    solution = coefficient_matrix.solve_right(rhs)
    H = tuple(
        sum(
            solution[coordinate * len(monomials) + index] * monomial
            for index, monomial in enumerate(monomials)
        )
        for coordinate in range(3)
    )
    assert sum(H[index] * critical[index] for index in range(3)) == target * F
    beta = sum(
        theta(H[index], (x, y, z)[index])
        for index in range(3)
    )
    return H, beta


for target in targets:
    print("TARGET", target)
    for degree in range(5):
        solution = solve_polynomial_field(target, degree)
        print(" DEGREE", degree, "EXISTS", solution is not None)
        if solution is not None:
            H, beta = solution
            print(" H", H)
            print(" beta", beta)
            print(" beta_minus_alpha", beta - target)
            break

print("Q32_LOWDEGREE_WEIGHT_FIELDS=PASS")
