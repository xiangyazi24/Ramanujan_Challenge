#!/usr/bin/env sage
"""Exact, fast audit for Q5938.

Use the explicit certified order-three long-ray operator (60.9) and conjugate
it through the TRUE Newton row.  This is enough to expose the first exact
characteristic-zero obstruction; no numerical recurrence guessing occurs.

For f_d and fixed L,

  G_{d,L}(f)=(-1)^L (d+1) binom(d+L+1,L) Delta_d^L(f_d/(d+1)).

Set g_L(r)=G_{M-r-L,L}(f), so only forward shifts f(r+j) occur.
The companion-state observability calculation below constructs the exact
primitive annihilator of g_L and factors its end coefficients.
"""

from math import factorial
from ore_algebra import *

R.<M, r> = QQ[]
K = R.fraction_field()
OA.<Sr> = OreAlgebra(R)

# Certified operator (60.9) for
# A_M(r)=sum_t binom(M,t)^2 binom(2M-t,M) binom(2M-t,r).
a0 = -(r + 1) * (r - 2*M)^2
a1 = -(r - M + 1) * (
    -M^2 - 6*M*r + 3*r^2 - 10*M + 6*r + 3
)
a2 = -(r + 2) * (
    M^2 - 6*M*r + 3*r^2 - 11*M + 9*r + 7
)
a3 = -(r - M + 2) * (r + 2) * (r + 3)
P = a0 + a1*Sr + a2*Sr^2 + a3*Sr^3
assert P.order() == 3


def sigma(value, amount=1):
    return K(value)(M=M, r=r + amount)


def sigma_matrix(A, amount=1):
    return A.apply_map(lambda x: sigma(x, amount))


def poly_binomial(upper, lower):
    out = K.one()
    for j in range(lower):
        out *= (upper - j) / (j + 1)
    return K(out)


def true_newton_forward_coefficients(length):
    """g_L(r)=sum_j coeff[j] f(r+j), exactly."""
    d0 = M - r - length
    return tuple(
        K(
            (-1)^(length-j)
            * poly_binomial(d0 + length - j, length - j)
            * poly_binomial(d0 + length + 1, j)
        )
        for j in range(length + 1)
    )


def audit_true_newton_identity(max_length=16):
    Rd.<d> = QQ[]
    for length in range(max_length + 1):
        Krow = (d + 1) * poly_binomial(K(d + length + 1), length)
        for i in range(length + 1):
            direct = (
                (-1)^i
                * prod(d + j for j in range(1, i + 1)) / factorial(i)
                * prod(d + length + 2 - j for j in range(1, length-i + 1))
                / factorial(length-i)
            )
            gauged = (
                (-1)^length * Krow
                * (-1)^(length-i) * binomial(length, i)
                / (d + i + 1)
            )
            assert K(direct - gauged) == 0
    print("TRUE_NEWTON_GAUGE_IDENTITY", max_length + 1, "PASS")


def companion(operator):
    rho = operator.order()
    A = matrix(K, rho, rho, 0)
    for i in range(rho - 1):
        A[i, i + 1] = 1
    for j in range(rho):
        A[rho - 1, j] = -K(operator[j]) / K(operator[rho])
    return A


def observation(A, coefficients):
    rho = A.nrows()
    e0 = vector(K, [1] + [0] * (rho - 1))
    out = vector(K, [0] * rho)
    transition = identity_matrix(K, rho)
    for j, coefficient_value in enumerate(coefficients):
        out += coefficient_value * (e0 * transition)
        transition = sigma_matrix(A, j) * transition
    return out


def annihilator_relation(A, c):
    rho = A.nrows()
    rows = []
    transition = identity_matrix(K, rho)
    for j in range(rho + 1):
        shifted = vector(K, [sigma(x, j) for x in c])
        rows.append(shifted * transition)
        transition = sigma_matrix(A, j) * transition
    O = matrix(K, rows)
    ker = O.left_kernel()
    assert ker.dimension() >= 1
    return ker.basis()[0], O


def primitive_relation(relation):
    den = lcm([K(x).denominator() for x in relation if x])
    polys = [R(den*K(x)) for x in relation]
    common = gcd([x for x in polys if x])
    polys = [R(x/common) for x in polys]
    qden = lcm(
        c.denominator()
        for f in polys
        for c in f.coefficients()
    )
    polys = [R(qden*f) for f in polys]
    zcontent = gcd([
        ZZ(c) for f in polys for c in f.coefficients() if c
    ])
    polys = tuple(R(f/zcontent) for f in polys)
    return polys, R(den), R(common)


def alias_factors(f, length):
    rem = R(f)
    out = []
    for c in range(-length - 8, 10):
        lin = R(M-r+c)
        e = 0
        while rem and rem % lin == 0:
            rem = R(rem/lin)
            e += 1
        if e:
            out.append((c, e))
    return tuple(out)


def coefficient_bits(polys):
    return max(
        [1] + [
            abs(ZZ(c.numerator())).nbits()
            for f in polys for c in f.coefficients() if c
        ]
    )


A = companion(P)
audit_true_newton_identity()

for length in range(0, 17):
    c = observation(A, true_newton_forward_coefficients(length))
    relation, O = annihilator_relation(A, c)
    primitive, den, common = primitive_relation(relation)
    assert O.rank() == 3
    assert len(primitive) == 4
    record = {
        "L": length,
        "order": 3,
        "degrees": tuple(f.total_degree() for f in primitive),
        "coefficient_bits": coefficient_bits(primitive),
        "denominator_degree": den.total_degree(),
        "removed_common_degree": common.total_degree(),
        "trailing_alias": alias_factors(primitive[0], length),
        "leading_alias": alias_factors(primitive[3], length),
        "trailing": str(factor(primitive[0])),
        "leading": str(factor(primitive[3])),
    }
    print("TRUE_NEWTON_CONJUGATE", record)

print("Q5938_FAST_AUDIT=PASS")
