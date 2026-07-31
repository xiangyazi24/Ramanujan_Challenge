#!/usr/bin/env sage
"""Scratch probe for the correct Apéry Laurent polynomial.

This file intentionally contains only exact ideal computations.  It is
used to rebuild the logarithmic-reduction claims after correcting

    D = (1+y)(1+z) + xyz.
"""

from itertools import combinations
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from q32_cartier_packet_audit import polytope_points

R.<x, y, z> = PolynomialRing(QQ, order="degrevlex")
A = 1 + x
B = 1 + y
C = 1 + z
D = (1 + y) * (1 + z) + x * y * z
FACTORS = (A, B, C, D)
F = prod(FACTORS)


def torus_saturation(ideal):
    return ideal.saturation(R.ideal(x * y * z))[0]


critical = torus_saturation(
    R.ideal(
        x * F.derivative(x) - F,
        y * F.derivative(y) - F,
        z * F.derivative(z) - F,
    )
)
singular = critical + R.ideal(F)
smooth = critical.saturation(R.ideal(F))[0]
intersection = singular.intersection(smooth)

print("critical dimension", critical.dimension())
print("singular dimension", singular.dimension())
print("smooth dimension", smooth.dimension())
if smooth.dimension() == 0:
    print("smooth vector dimension", smooth.vector_space_dimension())
print("critical == intersection",
      all(intersection.reduce(g) == 0 for g in critical.gens())
      and all(critical.reduce(g) == 0 for g in intersection.gens()))
print("singular radical", singular.radical() == singular)
print("singular minimal associated primes")
for prime in singular.minimal_associated_primes():
    print(" ", prime.dimension(), list(prime.gens()))

pair_intersection = None
for left, right in combinations(range(4), 2):
    pair = torus_saturation(R.ideal(FACTORS[left], FACTORS[right]))
    print("pair", left, right, "dim", pair.dimension(),
          "gens", list(pair.gens()))
    pair_intersection = (
        pair if pair_intersection is None
        else pair_intersection.intersection(pair)
    )
print("singular == pair intersection",
      all(pair_intersection.reduce(g) == 0 for g in singular.gens())
      and all(singular.reduce(g) == 0 for g in pair_intersection.gens()))

# Direct smooth algebra with an inverse of xyzF.
RC.<X, Y, Z, T> = PolynomialRing(QQ, order="degrevlex")
FC = (1 + X) * (1 + Y) * (1 + Z) * (
    (1 + Y) * (1 + Z) + X * Y * Z
)
smooth_inverse = RC.ideal(
    X * FC.derivative(X) - FC,
    Y * FC.derivative(Y) - FC,
    Z * FC.derivative(Z) - FC,
    1 - T * FC * X * Y * Z,
)
print("inverse dimension", smooth_inverse.dimension())
if smooth_inverse.dimension() == 0:
    print("inverse vector dimension", smooth_inverse.vector_space_dimension())
    print("inverse normal basis", smooth_inverse.normal_basis())


def theta(poly, variable):
    return variable * poly.derivative(variable)


syzygies = R.ideal(
    x * F.derivative(x),
    y * F.derivative(y),
    z * F.derivative(z),
    -F,
).syzygy_module()
weight_generators = []
print("syzygies")
for index, row in enumerate(syzygies):
    H = tuple(row[j] for j in range(3))
    h = row[3]
    aa = h - sum(H)
    weight_generators.append(aa)
    bb = sum(theta(H[j], (x, y, z)[j]) for j in range(3))
    print(" G", index, "H", H, "h", h, "a", aa, "b", bb)
    assert sum(H[j] * theta(F, (x, y, z)[j])
               for j in range(3)) == h * F
    for name, ideal in (
        ("AB", torus_saturation(R.ideal(A, B))),
        ("AC", torus_saturation(R.ideal(A, C))),
        ("AD", torus_saturation(R.ideal(A, D))),
        ("BC", torus_saturation(R.ideal(B, C))),
    ):
        print("   ", name,
              "H", [ideal.reduce(entry) for entry in H],
              "a", ideal.reduce(aa),
              "b", ideal.reduce(bb))

weight_ideal = torus_saturation(R.ideal(weight_generators))
print("weight ideal GB", weight_ideal.groebner_basis())
print("weight ideal == smooth",
      all(weight_ideal.reduce(g) == 0 for g in smooth.gens())
      and all(smooth.reduce(g) == 0 for g in weight_ideal.gens()))


normal_basis = smooth_inverse.normal_basis()
inverse_monomials = {
    X: T * FC * Y * Z,
    Y: T * FC * X * Z,
    Z: T * FC * X * Y,
}


def torus_monomial(point):
    out = RC.one()
    for variable, exponent in zip((X, Y, Z), point):
        if exponent == 1:
            out *= variable
        elif exponent == -1:
            out *= inverse_monomials[variable]
        else:
            assert exponent == 0
    return smooth_inverse.reduce(out)


def quotient_coordinates(polynomial):
    remainder = smooth_inverse.reduce(polynomial)
    return vector(
        QQ,
        [remainder.monomial_coefficient(monomial)
         for monomial in normal_basis],
    )


def multiplication_matrix(element):
    return matrix(
        QQ,
        [list(quotient_coordinates(element * monomial))
         for monomial in normal_basis],
    ).transpose()


PE.<E> = PolynomialRing(QQ)
charpolys = {}
for point in polytope_points(1):
    if point == (0, 0, 0):
        continue
    base = smooth_inverse.reduce(1 - torus_monomial(point))
    charpoly = multiplication_matrix(base).charpoly(E)
    charpolys.setdefault(charpoly, []).append(point)

annihilator = PE.one()
for charpoly in charpolys:
    annihilator = lcm(annihilator, charpoly)
print("characteristic classes", len(charpolys))
for charpoly, points in charpolys.items():
    print(" ", charpoly, points)
print("annihilator", annihilator)
print("annihilator degree", annihilator.degree())
print("annihilator valuation", annihilator.valuation())
print("annihilator endpoints", annihilator[0],
      annihilator.leading_coefficient())
