#!/usr/bin/env sage
"""Lift the simple smooth-critical generators to logarithmic fields.

For the correct Apéry polynomial, the logarithmic weights of a Singular
syzygy basis generate the smooth critical ideal after torus saturation.
This script expresses the simple Groebner generators

    y-z, x-2z-1, z^2-1/2

as Laurent combinations of those weights, constructs the corresponding
logarithmic vector fields, and records their divergences.
"""

R.<x, y, z> = PolynomialRing(QQ, order="degrevlex")
F = (
    (1 + x) * (1 + y) * (1 + z)
    * ((1 + y) * (1 + z) + x * y * z)
)


def theta(polynomial, variable):
    return variable * polynomial.derivative(variable)


rows = list(
    R.ideal(
        theta(F, x), theta(F, y), theta(F, z), -F
    ).syzygy_module()
)
fields = [tuple(row[index] for index in range(3)) for row in rows]
multipliers = [row[3] for row in rows]
weights = [
    multipliers[index] - sum(fields[index])
    for index in range(len(rows))
]

R4.<X, Y, Z, U> = PolynomialRing(QQ, order="degrevlex")


def into_r4(polynomial):
    return R4(polynomial(x=X, y=Y, z=Z))


weight4 = [into_r4(weight) for weight in weights]
generators4 = tuple(weight4) + (1 - U * X * Y * Z,)
ideal4 = R4.ideal(generators4)
targets = (y - z, x - 2 * z - 1, z**2 - QQ(1) / 2)

LR.<lx, ly, lz> = LaurentPolynomialRing(QQ, 3)


def inverse_substitution(polynomial):
    return sum(
        QQ(coefficient)
        * lx**(exponent[0] - exponent[3])
        * ly**(exponent[1] - exponent[3])
        * lz**(exponent[2] - exponent[3])
        for exponent, coefficient in polynomial.dict().items()
    )


def into_laurent(polynomial):
    return LR(polynomial(x=lx, y=ly, z=lz))


for target_index, target in enumerate(targets):
    lift_matrix = singular.lift(
        singular(ideal4), singular(into_r4(target))
    ).sage()
    lift4 = [
        R4(lift_matrix[index, 0])
        for index in range(len(generators4))
    ]
    assert into_r4(target) == sum(
        lift4[index] * generators4[index]
        for index in range(len(generators4))
    )
    coefficients = [
        inverse_substitution(lift4[index])
        for index in range(len(rows))
    ]
    H = tuple(
        sum(
            coefficients[row_index]
            * into_laurent(fields[row_index][coordinate])
            for row_index in range(len(rows))
        )
        for coordinate in range(3)
    )
    h = sum(
        coefficients[row_index]
        * into_laurent(multipliers[row_index])
        for row_index in range(len(rows))
    )
    alpha = h - sum(H)
    beta = sum(
        variable * H[index].derivative(variable)
        for index, variable in enumerate((lx, ly, lz))
    )
    assert alpha == into_laurent(target)
    assert sum(
        H[index]
        * variable
        * into_laurent(F).derivative(variable)
        for index, variable in enumerate((lx, ly, lz))
    ) == h * into_laurent(F)
    print("TARGET", target_index, target)
    print("COEFFICIENTS", coefficients)
    print("H", H)
    print("h", h)
    print("alpha", alpha)
    print("beta", beta)
    print("beta_minus_alpha", beta - alpha)

print("Q32_LOG_WEIGHT_FIELDS=PASS")
