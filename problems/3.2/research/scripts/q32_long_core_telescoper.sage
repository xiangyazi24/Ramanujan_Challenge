#!/usr/bin/env sage
"""Exact creative telescopers for the two long first-cell shell rays.

Section 57 writes the full-support part of C_M(M-r)-b_M as

    2 A_M(r) + B_M(r),

where

    A_M(r) = sum_t binom(M,t)^2 binom(2M-t,M) binom(2M-t,r),
    B_M(r) = sum_t binom(M,t)^2 binom(2M-t,r)^2.

Both sums have order-three recurrences in r.  This script derives the
operators (and certificates) over QQ(M,r,t) using exact Ore creative
telescoping, prints their factored coefficients, and checks them against
direct integer sums.  No guessed recurrence is used.

The two small monkey patches work around Sage 10.9 / ore_algebra 0.5
compatibility issues:

* the default associated commutative algebra tries to build an unsupported
  nested Singular ring;
* the default multivariate QQ nullspace path omits the denominator-clearing
  wrapper.
"""

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


def derive_operator(t_ratio_den, t_ratio_num, r_ratio_den, r_ratio_num):
    """Return the first telescoper and its exact certificate."""

    ideal = OA.ideal(
        [
            t_ratio_den * St - t_ratio_num,
            r_ratio_den * Sr - r_ratio_num,
        ]
    )
    telescopers, certificates = ideal.ct(
        St - 1,
        certificates=True,
        early_termination=True,
        iteration_limit=15,
    )
    assert len(telescopers) == len(certificates) == 1
    operator = telescopers[0]
    assert operator.order() == 3
    return operator, certificates[0]


A_operator, A_certificate = derive_operator(
    (t + 1)^2 * (2 * M - t)^2,
    (M - t)^3 * (2 * M - t - r),
    r + 1,
    2 * M - t - r,
)

B_operator, B_certificate = derive_operator(
    (t + 1)^2 * (2 * M - t)^2,
    (M - t)^2 * (2 * M - t - r)^2,
    (r + 1)^2,
    (2 * M - t - r)^2,
)

# The long core is 2*A_M+B_M.  Closure under addition therefore supplies
# an annihilator of order at most 3+3=6 (the least common left multiple).
# Their greatest common right divisor is one, so the operator lclm itself
# has order exactly six.  We deliberately do not construct its enormous
# coefficients here.
assert A_operator.gcrd(B_operator).order() == 0
LONG_CORE_ORDER_BOUND = A_operator.order() + B_operator.order()


def C(n, k):
    return comb(n, k) if 0 <= k <= n else 0


def A_value(moment, residue):
    return sum(
        C(moment, index)^2
        * C(2 * moment - index, moment)
        * C(2 * moment - index, residue)
        for index in range(moment + 1)
    )


def B_value(moment, residue):
    return sum(
        C(moment, index)^2
        * C(2 * moment - index, residue)^2
        for index in range(moment + 1)
    )


def evaluate_coefficient(polynomial, moment, residue):
    return ZZ(polynomial(M=moment, r=residue))


def audit_operator(operator, value):
    for moment in list(range(4, 30)) + [50, 73]:
        values = [value(moment, residue) for residue in range(2 * moment + 1)]
        for residue in range(0, 2 * moment - operator.order()):
            total = sum(
                evaluate_coefficient(operator[shift], moment, residue)
                * values[residue + shift]
                for shift in range(operator.order() + 1)
            )
            assert total == 0, (moment, residue, total)


audit_operator(A_operator, A_value)
audit_operator(B_operator, B_value)

print("A_ORDER", A_operator.order())
for shift in range(A_operator.order() + 1):
    print("A_COEFF", shift, factor(A_operator[shift]))

print("B_ORDER", B_operator.order())
for shift in range(B_operator.order() + 1):
    print("B_COEFF", shift, factor(B_operator[shift]))

print("A_CERTIFICATE_DEGREE", A_certificate.degree(St))
print("B_CERTIFICATE_DEGREE", B_certificate.degree(St))
assert LONG_CORE_ORDER_BOUND == 6
print("LONG_CORE_ORDER_BOUND", LONG_CORE_ORDER_BOUND)
print("PASS: exact order-three long-core telescopers")
