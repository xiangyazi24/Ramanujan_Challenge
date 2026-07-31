#!/usr/bin/env sage
"""Canonical finite-field Jacobian lift for the corrected log route.

For one fixed (M,L) this constructs the two projected terminal packet,
applies the stable order-66 smooth-quotient operator, lifts the result
through

    <theta_x(F)-F, theta_y(F)-F, theta_z(F)-F, 1-u*x*y*z>,

and forms the exact moment-(M+1) source ``sum theta_i(C_i)``.  It then
tests a deliberately generous bounded L-window made from the adjacent
packet and the five Lambda powers occurring in the W integrand.

This is a falsifiable modular experiment, not a characteristic-zero
existence proof: a different Jacobian lift can change the displayed
source by a twisted exact term.
"""

from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from q32_cartier_packet_audit import polytope_points


prime = 1009
field = GF(prime)
M = 10
L0 = 5

LR.<x, y, z> = LaurentPolynomialRing(field, 3)
Lambda = (
    (1 + x) * (1 + y) * (1 + z)
    * ((1 + y) * (1 + z) + x * y * z)
    / (x * y * z)
)
F = x * y * z * Lambda


def theta(polynomial, variable):
    return variable * polynomial.derivative(variable)


def vector_field(coefficients, polynomial):
    return sum(
        coefficient * theta(polynomial, variable)
        for coefficient, variable in zip(coefficients, (x, y, z))
    )


def divergence(coefficients):
    return sum(
        theta(coefficient, variable)
        for coefficient, variable in zip(coefficients, (x, y, z))
    )


G_ABC_H = (-x - 1, LR.zero(), (x + 1) * (z + 1))
G_ABC_h = 2 * x * z - x + 2 * z
G_BC_H = (
    -(x + 1)
    * (2 * y * z**2 + y * z + 2 * z**2 - 4 * y + 6 * z + 1)
    / x,
    (y + 1)
    * (2 * x * z**2 + 5 * x * z + 2 * z**2 - 2 * x + 5 * z)
    / x,
    -(z + 1)
    * (
        4 * x * y - x * z + 2 * y * z + 2 * z**2
        - 3 * x + 5 * y + 2 * z
    )
    / x,
)
G_BC_h = -(4 * z**3 + 4 * z**2 + x) / x


def twisted(coefficients, multiplier, polynomial):
    aa = multiplier - sum(coefficients)
    bb = divergence(coefficients)
    return (
        vector_field(coefficients, polynomial)
        + (bb + M * aa) * polynomial
    )


def packet(endpoint, order):
    out = LR.zero()
    for point in polytope_points(1):
        monomial = x**point[0] * y**point[1] * z**point[2]
        out += monomial**(-endpoint) * (1 - monomial)**order
    return out


def projected_packet(order):
    value = packet(M, order)
    value -= twisted(G_ABC_H, G_ABC_h, value) / (M + 1)
    value += twisted(G_BC_H, G_BC_h, value) / (M + 1)
    return value


def convolution(left, right):
    out = [field.zero()] * (len(left) + len(right) - 1)
    for i, aa in enumerate(left):
        for j, bb in enumerate(right):
            out[i + j] += aa * bb
    return out


quadratics = (
    (1, -6, 1),
    (-1, 2, 1),
    (2, -4, 1),
    (field(-1) / 2, 0, 1),
    (field(7) / 4, -3, 1),
    (1, 2, 1),
    (-1, -2, 1),
    (field(1) / 2, -2, 1),
    (field(1) / 4, -1, 1),
    (-2, 0, 1),
    (field(-1) / 4, -1, 1),
)
q = [field.one()]
for factor in quadratics:
    q = convolution(q, [field(entry) for entry in factor])
stable_operator = convolution(convolution(q, q), q)
assert len(stable_operator) == 67

P = sum(
    stable_operator[shift] * projected_packet(L0 + shift)
    for shift in range(67)
)

# Clear all Laurent exponents, then use Singular's canonical lift in a
# polynomial ring carrying an inverse of xyz.
exponents = list(P.dict())
clear = tuple(max(0, -min(exp[index] for exp in exponents))
              for index in range(3))
R4.<X, Y, Z, U> = PolynomialRing(field, order="degrevlex")


def cleared_polynomial(polynomial):
    return sum(
        field(coefficient)
        * X**(exponent[0] + clear[0])
        * Y**(exponent[1] + clear[1])
        * Z**(exponent[2] + clear[2])
        for exponent, coefficient in polynomial.dict().items()
    )


F4 = (
    (1 + X) * (1 + Y) * (1 + Z)
    * ((1 + Y) * (1 + Z) + X * Y * Z)
)
generators = (
    X * F4.derivative(X) - F4,
    Y * F4.derivative(Y) - F4,
    Z * F4.derivative(Z) - F4,
    1 - U * X * Y * Z,
)
ideal4 = R4.ideal(generators)
target = cleared_polynomial(P)
lifted_matrix = singular.lift(
    singular(ideal4), singular(target)
).sage()
lifted = [R4(lifted_matrix[index, 0]) for index in range(4)]
assert target == sum(
    lifted[index] * generators[index] for index in range(4)
)


def substitute_inverse(polynomial):
    return sum(
        field(coefficient)
        * x**(exponent[0] - exponent[3])
        * y**(exponent[1] - exponent[3])
        * z**(exponent[2] - exponent[3])
        for exponent, coefficient in polynomial.dict().items()
    )


clearing_monomial = x**clear[0] * y**clear[1] * z**clear[2]
B_lift = [
    substitute_inverse(lifted[index]) / clearing_monomial
    for index in range(3)
]
critical_generators = (
    theta(F, x) - F,
    theta(F, y) - F,
    theta(F, z) - F,
)
assert P == sum(
    B_lift[index] * critical_generators[index]
    for index in range(3)
)

C_lift = [x * y * z * entry for entry in B_lift]
source = sum(
    theta(C_lift[index], (x, y, z)[index])
    for index in range(3)
)

# A generous exact-polynomial W envelope after multiplying the
# moment-(M+1) source by Lambda.  It contains each of the five
# Lambda-powers separately, rather than only their prescribed U_n
# combination, and also contains the Phi_M term.
window = 6
candidates = []
for shift in range(window + 1):
    for power in range(1, 6):
        candidates.append(
            Lambda**power
            * (
                packet(M - 1, L0 + shift)
                + packet(M + 1, L0 + shift)
            )
        )
    candidates.append(packet(M, L0 + shift))


def span_contains(vectors, target_vector):
    support = sorted(
        set(target_vector.dict()).union(
            *(set(vector_value.dict()) for vector_value in vectors)
        )
    )
    matrix_value = matrix(
        field,
        [[vector_value[exponent] for vector_value in vectors]
         for exponent in support],
    )
    rhs = vector(field, [target_vector[exponent] for exponent in support])
    return (
        matrix_value.rank(),
        matrix_value.augment(rhs).rank(),
        len(support),
    )


span_ranks = span_contains(candidates, Lambda * source)
print("Q32_LOG_CANONICAL_SOURCE_PROBE=PASS")
print("MODULUS", prime, "M", M, "L", L0)
print("PROJECTED_OPERATOR_SUPPORT", len(P.dict()))
print("CLEARING_EXPONENT", clear)
print("LIFT_SUPPORTS", [len(entry.dict()) for entry in B_lift])
print("SOURCE_SUPPORT", len(source.dict()))
print("SOURCE_EXPONENT_BOX",
      tuple(
          (
              min(exp[index] for exp in source.dict()),
              max(exp[index] for exp in source.dict()),
          )
          for index in range(3)
      ))
print("GENEROUS_W_ENVELOPE_RANKS", span_ranks)
