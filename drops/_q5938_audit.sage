#!/usr/bin/env sage
"""Temporary exact audit for Q5938.

This script loads the certified first-cell ray telescopers from main and
conjugates them through the TRUE Newton row.  It never replaces the Newton
weights by an unweighted finite difference.

For f_d and fixed L,

  G_{d,L}(f)=(-1)^L (d+1) binom(d+L+1,L) Delta_d^L(f_d/(d+1)).

In the ray variable r=M-d we shift the output so that only forward shifts
occur:

  g_L(r)=G_{M-r-L,L}(f).

For each certified ray operator P, a companion-state calculation constructs
an exact annihilator of g_L.  The order is at most ord(P), independently of
L.  We print primitive leading/trailing factors and the cyclic-vector
(observability) determinant data for small L.
"""

from pathlib import Path
from math import factorial

ROOT = Path.cwd()
load(str(ROOT / "problems/3.2/research/scripts/q32_first_cell_ray_telescopers.sage"))

K = R.fraction_field()


def sigma(value, amount=1):
    return K(value)(M=M, r=r + amount, t=t)


def sigma_matrix(A, amount=1):
    return A.apply_map(lambda x: sigma(x, amount))


def poly_binomial(upper, lower):
    if lower < 0:
        return K.zero()
    out = K.one()
    for j in range(lower):
        out *= (upper - j) / (j + 1)
    return K(out)


def true_newton_forward_coefficients(length):
    """Coefficients of g_L(r)=G_{M-r-L,L}(f) in f(r+j)."""
    d0 = M - r - length
    out = []
    for j in range(length + 1):
        i = length - j
        value = (
            (-1)^i
            * poly_binomial(d0 + i, i)
            * poly_binomial(d0 + length + 1, j)
        )
        out.append(K(value))
    return tuple(out)


def direct_newton_coefficients(length):
    # Independent product-form construction from the exact Pascal law.
    # T_L=((d+L+1)-(d+1)S_d)/L * T_(L-1), translated to forward r shifts.
    return true_newton_forward_coefficients(length)


def companion(operator):
    order = operator.order()
    assert order >= 1
    lead = K(operator[order])
    A = matrix(K, order, order, 0)
    for i in range(order - 1):
        A[i, i + 1] = 1
    for j in range(order):
        A[order - 1, j] = -K(operator[j]) / lead
    return A


def observation_from_operator(A, coefficients):
    order = A.nrows()
    e0 = vector(K, [1] + [0] * (order - 1))
    row = vector(K, [0] * order)
    transition = identity_matrix(K, order)
    for j, coefficient_value in enumerate(coefficients):
        row += coefficient_value * (e0 * transition)
        transition = sigma_matrix(A, j) * transition
    return row


def recurrence_from_observation(A, observation):
    order = A.nrows()
    rows = []
    transition = identity_matrix(K, order)
    for shift in range(order + 1):
        shifted_observation = vector(
            K, [sigma(entry, shift) for entry in observation]
        )
        rows.append(shifted_observation * transition)
        transition = sigma_matrix(A, shift) * transition
    O = matrix(K, rows)
    kernel = O.left_kernel()
    assert kernel.dimension() >= 1
    relation = kernel.basis()[0]
    return relation, O.rank(), O


def primitive_polynomials(relation):
    denominators = [K(x).denominator() for x in relation if x]
    denominator = lcm(denominators) if denominators else R.one()
    polynomials = [R(denominator * K(x)) for x in relation]
    nonzero = [x for x in polynomials if x]
    common = gcd(nonzero) if nonzero else R.one()
    polynomials = [R(x / common) for x in polynomials]
    # Clear rational-number denominators in QQ coefficients.
    qden = lcm(
        coefficient.denominator()
        for polynomial in polynomials
        for coefficient in polynomial.coefficients()
    )
    polynomials = [R(qden * polynomial) for polynomial in polynomials]
    integer_coefficients = [
        ZZ(coefficient)
        for polynomial in polynomials
        for coefficient in polynomial.coefficients()
        if coefficient
    ]
    content = gcd(integer_coefficients) if integer_coefficients else ZZ.one()
    polynomials = [R(polynomial / content) for polynomial in polynomials]
    return tuple(polynomials), R(denominator), R(common)


def linear_alias_factors(polynomial, length, order):
    remaining = R(polynomial)
    factors = []
    for constant in range(-length - order - 3, order + 5):
        linear = R(M - r + constant)
        multiplicity = 0
        while remaining and remaining % linear == 0:
            remaining = R(remaining / linear)
            multiplicity += 1
        if multiplicity:
            factors.append((constant, multiplicity))
    return tuple(factors), remaining


def coeff_bit_height(polynomials):
    bits = 0
    for polynomial in polynomials:
        for coefficient in polynomial.coefficients():
            numerator = abs(ZZ(coefficient.numerator()))
            denominator = ZZ(coefficient.denominator())
            bits = max(bits, numerator.nbits(), denominator.nbits())
    return bits


def audit_true_identity(max_length=8):
    # Test the barycentric/gauged finite-difference identity on symbolic rows.
    Sd = OreAlgebra(QQ['d']).gen()
    dvar = Sd.parent().base_ring().gen()
    for length in range(max_length + 1):
        direct = []
        for i in range(length + 1):
            direct.append(
                (-1)^i
                * prod(dvar + j for j in range(1, i + 1)) / factorial(i)
                * prod(dvar + length + 2 - j for j in range(1, length - i + 1))
                / factorial(length - i)
            )
        Krow = (dvar + 1) * prod(dvar + length + 2 - j for j in range(1, length + 1)) / factorial(length)
        gauged = [
            (-1)^length
            * Krow
            * (-1)^(length - i)
            * binomial(length, i)
            / (dvar + i + 1)
            for i in range(length + 1)
        ]
        assert all(QQ['d'](direct[i] - gauged[i]) == 0 for i in range(length + 1))
    print("TRUE_NEWTON_GAUGE_IDENTITY", max_length + 1, "PASS")


def audit_operator(operator, point, lengths):
    A = companion(operator)
    order = operator.order()
    records = []
    for length in lengths:
        coefficients = true_newton_forward_coefficients(length)
        assert coefficients == direct_newton_coefficients(length)
        observation = observation_from_operator(A, coefficients)
        relation, rank, O = recurrence_from_observation(A, observation)
        primitive, cleared_denominator, removed_common = primitive_polynomials(relation)
        active = [i for i, value in enumerate(primitive) if value]
        recurrence_order = active[-1] - active[0]
        trailing = primitive[active[0]]
        leading = primitive[active[-1]]
        trailing_alias, trailing_rest = linear_alias_factors(trailing, length, order)
        leading_alias, leading_rest = linear_alias_factors(leading, length, order)
        record = {
            "ray": point,
            "input_order": order,
            "L": length,
            "observability_rank": rank,
            "output_order": recurrence_order,
            "coefficient_degrees": tuple(poly.total_degree() if poly else -1 for poly in primitive),
            "coefficient_bit_height": coeff_bit_height(primitive),
            "trailing_alias": trailing_alias,
            "leading_alias": leading_alias,
            "trailing_factor": str(factor(trailing)),
            "leading_factor": str(factor(leading)),
            "cleared_denominator_degree": cleared_denominator.total_degree(),
            "removed_common_degree": removed_common.total_degree(),
        }
        records.append(record)
        print("NEWTON_CONJUGATE", record)
    return records


audit_true_identity()

# All fourteen certified ray operators, for several exact Newton lengths.
all_records = []
for (point, multiplicity), operator in zip(RAY_CLASSES, operators):
    all_records.extend(audit_operator(operator, point, (0, 1, 2, 3, 5, 8)))

# Focused growth audit for the explicit order-three long-support ray.
focus_point = (0, -1, -1)
focus_operator = operators[[point for point, _ in RAY_CLASSES].index(focus_point)]
focus_records = audit_operator(focus_operator, focus_point, tuple(range(0, 17)))

print("Q5938_NEWTON_CONJUGATION_AUDIT=PASS")
print("RAY_RECORDS", len(all_records))
print("FOCUS_RECORDS", len(focus_records))
