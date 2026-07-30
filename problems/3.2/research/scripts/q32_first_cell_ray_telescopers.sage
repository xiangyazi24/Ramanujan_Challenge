#!/usr/bin/env sage
"""Exact P-recursive decomposition of the complete first-cell shell.

For d=M-r>M/2,

    C_M(d)-b_M = sum_{kappa in P(Z), kappa != 0} c_M(d*kappa).

The 21 nonzero lattice rays of the Newton polytope fall into 14
y/z-symmetry classes.  Each ray coefficient is one proper hypergeometric
sum in t.  This script derives an exact Ore telescoper in r for every
class, verifies it against direct integer sums, and proves the uniform
order bound

    order(C_M(M-r)-b_M) <= 38.

The bound is the sum of the 14 individual orders; computing the large
least common left multiple is not needed for this rigorous upper bound.
"""

from itertools import combinations
from math import comb

from ore_algebra import *
from ore_algebra import nullspace
from ore_algebra.ore_algebra import OreAlgebra_generic


def _associated_commutative_algebra(self):
    try:
        return self._commutative_ring
    except AttributeError:
        self._commutative_ring = PolynomialRing(
            self.base_ring(), self.variable_names()
        )
        return self._commutative_ring


OreAlgebra_generic.associated_commutative_algebra = (
    _associated_commutative_algebra
)

_original_kronecker = nullspace.kronecker


def _qq_safe_kronecker(subsolver, presolver=None):
    return nullspace.clear(_original_kronecker(subsolver, presolver))


nullspace.kronecker = _qq_safe_kronecker

R.<M, r, t> = QQ[]
OA.<Sr, St> = OreAlgebra(R)
d = M - r

# Representative, multiplicity under y/z interchange.  The multiplicities
# sum to the 21 nonzero lattice points of P.
RAY_CLASSES = (
    ((-1, -1, -1), 1),
    ((-1, -1, 0), 2),
    ((-1, -1, 1), 2),
    ((-1, 0, 0), 1),
    ((-1, 0, 1), 2),
    ((-1, 1, 1), 1),
    ((0, -1, -1), 1),
    ((0, -1, 0), 2),
    ((0, -1, 1), 2),
    ((0, 0, 1), 2),
    ((0, 1, 1), 1),
    ((1, 0, 0), 1),
    ((1, 0, 1), 2),
    ((1, 1, 1), 1),
)
assert sum(multiplicity for _, multiplicity in RAY_CLASSES) == 21


def shifted_binomial_ratio(upper, lower, shift):
    """binom(upper,lower+shift)/binom(upper,lower), shift in {-1,0,1}."""

    if shift == 1:
        return (upper - lower) / (lower + 1)
    if shift == -1:
        return lower / (upper - lower + 1)
    assert shift == 0
    return R.one()


def ray_operator(point):
    """Derive an exact r-telescoper for c_M((M-r)*point)."""

    u, v, w = point
    kx = t - u * d
    upper = 2 * M - t
    ky = M - v * d
    kz = M - w * d

    # Ratio under t -> t+1.  The first two binomials have fixed upper M;
    # the final two have upper 2M-t, which decreases by one.
    t_ratio = (
        (M - t)
        / (t + 1)
        * (M - kx)
        / (kx + 1)
        * (upper - ky)
        / upper
        * (upper - kz)
        / upper
    )

    # Under r -> r+1, d decreases by one and the three moving lower
    # arguments change by u,v,w.
    r_ratio = (
        shifted_binomial_ratio(M, kx, u)
        * shifted_binomial_ratio(upper, ky, v)
        * shifted_binomial_ratio(upper, kz, w)
    )

    ideal = OA.ideal(
        [
            t_ratio.denominator() * St - t_ratio.numerator(),
            r_ratio.denominator() * Sr - r_ratio.numerator(),
        ]
    )
    telescopers, certificates = ideal.ct(
        St - 1,
        certificates=True,
        early_termination=True,
        iteration_limit=16,
    )
    assert len(telescopers) == len(certificates) == 1
    return telescopers[0], certificates[0]


def C(n, k):
    return comb(n, k) if 0 <= k <= n else 0


def ray_value(moment, residue, point):
    u, v, w = point
    node = moment - residue
    return sum(
        C(moment, index)
        * C(moment, index - node * u)
        * C(2 * moment - index, moment - node * v)
        * C(2 * moment - index, moment - node * w)
        for index in range(moment + 1)
    )


def apery(moment):
    return sum(
        C(moment, index)^2 * C(moment + index, index)^2
        for index in range(moment + 1)
    )


def shell_fast(moment, node):
    quotient = moment // node
    out = 0
    for index in range(moment + 1):
        x_packet = sum(
            C(moment, moment - index + node * u)
            for u in range(-quotient, quotient + 1)
        )
        yz_packet = sum(
            C(2 * moment - index, moment - index + node * v)
            for v in range(-quotient, quotient + 1)
        )
        out += C(moment, index) * x_packet * yz_packet^2
    return out


def evaluate_coefficient(polynomial, moment, residue):
    return ZZ(polynomial(M=moment, r=residue))


def audit_ray(operator, point):
    for moment in list(range(8, 27)) + [37, 52]:
        values = [
            ray_value(moment, residue, point)
            for residue in range(moment + 1)
        ]
        for residue in range(moment - operator.order()):
            total = sum(
                evaluate_coefficient(operator[shift], moment, residue)
                * values[residue + shift]
                for shift in range(operator.order() + 1)
            )
            assert total == 0, (point, moment, residue, total)


operators = []
for point, multiplicity in RAY_CLASSES:
    operator, certificate = ray_operator(point)
    assert operator.order() <= 3
    audit_ray(operator, point)
    operators.append(operator)
    print(
        "RAY",
        point,
        "MULTIPLICITY",
        multiplicity,
        "ORDER",
        operator.order(),
    )
    print("  TRAILING", factor(operator[0]))
    print("  LEADING", factor(operator[operator.order()]))
    print("  CERTIFICATE_DEGREE", certificate.degree(St))

# The raw scalar modules have no common solution: every pair has unit
# greatest common right divisor.  This rules out a reduction coming from a
# literal shared scalar solution, but not a rational gauge/intertwiner or a
# special relation among the fourteen distinguished ray solutions.
for left, right in combinations(operators, 2):
    assert left.gcrd(right).order() == 0
print("PAIRWISE_GCRD_COUNT", binomial(len(operators), 2))

# Independent check that the 14-class ray sum is the complete nonconstant
# first-cell shell.
for moment in range(3, 22):
    for residue in range((moment - 1) // 2 + 1):
        node = moment - residue
        correction = sum(
            multiplicity * ray_value(moment, residue, point)
            for point, multiplicity in RAY_CLASSES
        )
        assert shell_fast(moment, node) == apery(moment) + correction

order_bound = sum(operator.order() for operator in operators)
assert order_bound == 38
print("FULL_FIRST_CELL_ORDER_BOUND", order_bound)
print("PASS: 14 exact ray telescopers and complete first-cell decomposition")
