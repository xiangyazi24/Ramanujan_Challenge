#!/usr/bin/env sage
"""Certify the order-three recurrence for the doubled period J_n.

Let

    T(n,k) = binom(n,k)^2 binom(2*n-k,n)^2,
    J_n = CT(Lambda(X)^n (Lambda(X^2)+40)).

The shifted-binomial formula for the coefficients of Lambda^n rewrites
J_n as one proper hypergeometric sum

    J_n = sum_k T(n,k) Q(n,k).

This script constructs Q exactly, derives a Zeilberger certificate with
``ore_algebra``, verifies ideal membership, and checks that the resulting
operator is exactly the primitive integer operator stored by
``q32_doubled_period_recurrence_guess.py``.

The certificate multiplier has only the displayed removable boundary
poles.  For n >= 5, the reciprocal-factorial continuation of T kills
the lower tail, is finite at k=n+1,n+2,n+3, and kills the upper tail
from k=n+4 onward.  Thus its forward difference sums to zero over all
integer k.  The cases n=0,...,4 are checked directly.
"""

import ast
from math import comb
from pathlib import Path

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

R = QQ["n,k"]
n, k = R.gens()
K = R.fraction_field()
OA = OreAlgebra(R, names=("Sn", "Sk"))
Sn, Sk = OA.gens()


def multiply(left, right):
    out = {}
    for u, left_coefficient in left.items():
        for v, right_coefficient in right.items():
            exponent = tuple(u[index] + v[index] for index in range(3))
            out[exponent] = (
                out.get(exponent, 0)
                + left_coefficient * right_coefficient
            )

    return out


def lambda_polynomial():
    one = {(0, 0, 0): 1}
    x = {(0, 0, 0): 1, (1, 0, 0): 1}
    y = {(0, 0, 0): 1, (0, 1, 0): 1}
    z = {(0, 0, 0): 1, (0, 0, 1): 1}
    bracket = multiply(y, z)
    bracket[(1, 1, 1)] = bracket.get((1, 1, 1), 0) + 1
    numerator = one
    for factor in (x, y, z, bracket):
        numerator = multiply(numerator, factor)

    return {
        (u[0] - 1, u[1] - 1, u[2] - 1): coefficient
        for u, coefficient in numerator.items()
    }


LAMBDA = lambda_polynomial()
assert len(LAMBDA) == 22
assert sum(LAMBDA.values()) == 40


def shifted_binomial_ratio(upper, lower, shift):
    """Return binom(upper,lower+shift)/binom(upper,lower)."""

    if shift == 0:
        return R.one()
    if shift == 2:
        return (upper - lower) * (upper - lower - 1) / (
            (lower + 1) * (lower + 2)
        )
    if shift == -2:
        return lower * (lower - 1) / (
            (upper - lower + 1) * (upper - lower + 2)
        )
    raise ValueError(shift)


Q = K(40)
for (u, v, w), weight in LAMBDA.items():
    Q += (
        weight
        * shifted_binomial_ratio(n, k, 2 * u)
        * shifted_binomial_ratio(2 * n - k, n, 2 * v)
        * shifted_binomial_ratio(2 * n - k, n, 2 * w)
    )
Q = K(Q)

expected_q_denominator = (
    (k + 1)
    * (k + 2)
    * (n + 1)^2
    * (n + 2)^2
    * (k - n - 2)^3
    * (k - n - 1)^3
)
assert Q.denominator() == expected_q_denominator
assert len(Q.numerator().monomials()) == 91

# T(n,k)=binom(n,k)^2*binom(2*n-k,n)^2 has these two shift ratios.
T_n_ratio = (
    (2 * n + 2 - k) * (2 * n + 1 - k) / (n + 1 - k)^2
)^2
T_k_ratio = (n - k)^4 / ((k + 1)^2 * (2 * n - k)^2)
F_n_ratio = K(T_n_ratio * Q(n=n + 1, k=k) / Q)
F_k_ratio = K(T_k_ratio * Q(n=n, k=k + 1) / Q)

ideal = OA.ideal(
    [
        F_n_ratio.denominator() * Sn - F_n_ratio.numerator(),
        F_k_ratio.denominator() * Sk - F_k_ratio.numerator(),
    ]
)
telescopers, certificates = ideal.ct(
    Sk - 1,
    certificates=True,
    early_termination=True,
    iteration_limit=12,
)
assert len(telescopers) == len(certificates) == 1
operator = telescopers[0]
certificate = certificates[0]

# This is the exact WZ identity P-(Sk-1)C in Ann(F), not a fit to data.
assert OA(operator) - (Sk - 1) * certificate in ideal
assert operator.order() == 3
assert all(
    operator[shift].denominator() == 1
    and operator[shift].numerator().degree() == 21
    for shift in range(4)
)


def stored_candidate():
    candidate_path = Path(__file__).with_name(
        "q32_doubled_period_recurrence_guess.py"
    )
    tree = ast.parse(candidate_path.read_text(encoding="utf-8"))
    for statement in tree.body:
        if not isinstance(statement, ast.Assign):
            continue
        if any(
            isinstance(target, ast.Name) and target.id == "RECURRENCE"
            for target in statement.targets
        ):
            return ast.literal_eval(statement.value)

    raise AssertionError("RECURRENCE was not found")


candidate_polynomials = [
    sum(
        ZZ(coefficient) * n^degree
        for degree, coefficient in enumerate(row)
    )
    for row in stored_candidate()
]
normalizing_ratio = K(operator[0] / candidate_polynomials[0])
assert normalizing_ratio == -4
assert all(
    operator[shift] == -4 * candidate_polynomials[shift]
    for shift in range(4)
)

# The certificate has degree zero in both shift operators.  If C(n,k) is
# its scalar and F=T*Q, the telescoping antidifference is C*T*Q.
assert certificate.degree(Sn) == 0
assert certificate.degree(Sk) == 0
assert len(certificate.coefficients()) == 1
certificate_scalar = certificate.coefficients()[0]
certificate_multiplier = K(certificate_scalar * Q)
expected_certificate_denominator = (
    (k + 1)
    * n
    * (n + 1)
    * (n + 2)
    * (n + 4)
    * (n + 3)^2
    * (k - n - 5)^3
    * (k - n - 4)^3
    * (k - n - 3)^4
    * (k - n - 2)^4
    * (k - n - 1)^4
)
assert (
    certificate_multiplier.denominator()
    == expected_certificate_denominator
)
assert certificate_multiplier.numerator().degree(n) == 45
assert certificate_multiplier.numerator().degree(k) == 24
assert certificate_multiplier.numerator() % k == 0
assert certificate_multiplier.numerator() % (k + 1) != 0


def C(upper, lower):
    return comb(upper, lower) if 0 <= lower <= upper else 0


def direct_summand(index, summation_index):
    base = (
        C(index, summation_index)^2
        * C(2 * index - summation_index, index)^2
    )
    out = 40 * base
    for (u, v, w), weight in LAMBDA.items():
        out += (
            weight
            * C(index, summation_index)
            * C(index, summation_index + 2 * u)
            * C(2 * index - summation_index, index + 2 * v)
            * C(2 * index - summation_index, index + 2 * w)
        )

    return out


def rational_summand(index, summation_index):
    base = (
        C(index, summation_index)^2
        * C(2 * index - summation_index, index)^2
    )

    return ZZ(base * Q(n=index, k=summation_index))


values = []
for index in range(8):
    for summation_index in range(index + 1):
        assert (
            rational_summand(index, summation_index)
            == direct_summand(index, summation_index)
        )
    values.append(
        sum(
            direct_summand(index, summation_index)
            for summation_index in range(index + 1)
        )
    )

# The boundary argument above applies from n=5.  These exact checks close
# the five initial recurrence positions where its generic certificate has
# a parameter pole or too-short upper tail.
for index in range(5):
    assert sum(
        candidate_polynomials[shift](n=index) * values[index + shift]
        for shift in range(4)
    ) == 0

print("Q32_DOUBLED_PERIOD_TELESCOPER=PASS")
print("Q_NUMERATOR_TERMS", len(Q.numerator().monomials()))
print("ORDER", operator.order())
print("COEFFICIENT_DEGREE", 21)
print("NORMALIZING_RATIO_TO_STORED_OPERATOR", normalizing_ratio)
print("CERTIFICATE_NUMERATOR_BIDEGREE", (45, 24))
print("INITIAL_BOUNDARY_CASES", 5)
