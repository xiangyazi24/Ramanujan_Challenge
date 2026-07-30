#!/usr/bin/env sage
"""Certify the rational Apéry quotient of the doubled-period module.

This script corrects the bounded polynomial-gauge test in
``q32_doubled_period_factor_audit.sage``.  The order-two left factor of
the certified doubled-period operator is rationally gauge equivalent to
the shifted Apéry operator.  For the distinguished solutions,

    G(n)^(-1) (K_n, K_{n+1})^T
        = -10080 (b_n, b_{n+1})^T,

where K is the first-order right-factor transform of J.

The script also constructs a transparent Ore-CRT operator U with

    U(J)_n = M(n) b_n,       U(b)_n = 0,

where U has five primitive integral polynomial coefficients of degree
32 and M is an integral polynomial of degree 30.  The resulting shell
carrier is target-preserving, but it is universally aliased to the old
carrier modulo every candidate prime.  With ``--blocks``, the four
hostile Newton blocks are evaluated exactly.
"""

import ast
import sys
from math import gcd as integer_gcd
from math import prod
from pathlib import Path

from ore_algebra import *


SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

from q32_cartier_packet_audit import apery, shell_batch, shell_fast


R = QQ["n"]
n = R.gen()
K = R.fraction_field()
OA = OreAlgebra(R, names=("Sn",))
Sn = OA.gen()
OAK = OreAlgebra(K, names=("Sn",))


def stored_candidate():
    tree = ast.parse(
        (SCRIPT_DIR / "q32_doubled_period_recurrence_guess.py").read_text(
            encoding="utf-8"
        )
    )
    for statement in tree.body:
        if not isinstance(statement, ast.Assign):
            continue
        if any(
            isinstance(target, ast.Name) and target.id == "RECURRENCE"
            for target in statement.targets
        ):
            return ast.literal_eval(statement.value)
    raise AssertionError("RECURRENCE was not found")


def polynomial_value(coefficients, index):
    return sum(
        ZZ(coefficient) * index^degree
        for degree, coefficient in enumerate(coefficients)
    )


RECURRENCE = stored_candidate()
P = sum(
    sum(
        ZZ(coefficient) * n^degree
        for degree, coefficient in enumerate(row)
    )
    * Sn^shift
    for shift, row in enumerate(RECURRENCE)
)
A4 = 4 * n^5 + 12 * n^4 + 17 * n^3 + 78 * n^2 + 63 * n - 54
B4 = 4 * n^5 + 20 * n^4 + 57 * n^3 + 118 * n^2 + 107 * n + 30
right_factor = A4 / 4 * Sn - B4 / 4
left_factor, factor_remainder = P.quo_rem(right_factor)
assert factor_remainder == 0
assert (P.order(), left_factor.order(), right_factor.order()) == (3, 2, 1)

apery_operator = (
    (n + 1)^3
    - (2 * n + 3)
    * (17 * (n + 1)^2 + 17 * (n + 1) + 5)
    * Sn
    + (n + 2)^3 * Sn^2
)


def sigma(value, amount=1):
    return value(n=n + amount)


def flatten(matrix_value):
    return vector(
        K,
        [
            matrix_value[row, column]
            for row in range(2)
            for column in range(2)
        ],
    )


# Uncouple sigma(G) M_A = M_L G into one scalar rational recurrence.
MA = apery_operator.companion_matrix()
ML = left_factor.companion_matrix()
basis_matrices = []
for index in range(4):
    basis = matrix(K, 2, 2, 0)
    basis[index // 2, index % 2] = 1
    basis_matrices.append(basis)

hom_transition = matrix(
    K,
    4,
    4,
    lambda row, column: flatten(
        ML * basis_matrices[column] * MA.inverse()
    )[row],
)
hom_krylov_rows = [vector(K, [1, 0, 0, 0])]
for _ in range(4):
    hom_krylov_rows.append(
        vector(K, [sigma(entry) for entry in hom_krylov_rows[-1]])
        * hom_transition
    )
hom_krylov = matrix(K, hom_krylov_rows[:4])
assert hom_krylov.det() != 0
hom_relation = hom_krylov.transpose().solve_right(-hom_krylov_rows[4])
hom_operator = OAK(list(hom_relation) + [1]).normalize()
hom_solutions = hom_operator.rational_solutions()
assert len(hom_solutions) == 1

hom_scalar = hom_solutions[0][0]
hom_rhs = vector(K, [sigma(hom_scalar, shift) for shift in range(4)])
gauge_vector = hom_krylov.solve_right(hom_rhs)
gauge = matrix(K, 2, 2, list(gauge_vector))
assert matrix(
    K,
    2,
    2,
    lambda row, column: sigma(gauge[row, column]),
) * MA == ML * gauge
assert gauge.det() != 0


# Generate the two distinguished sequences exactly.
maximum_index = 24
b_values = [ZZ(apery(index)) for index in range(maximum_index + 1)]
j_values = [QQ(45), QQ(225), QQ(3465)]
for index in range(maximum_index - 2):
    recurrence_coefficients = [
        polynomial_value(RECURRENCE[shift], index)
        for shift in range(4)
    ]
    j_values.append(
        -sum(
            recurrence_coefficients[shift] * j_values[index + shift]
            for shift in range(3)
        )
        / recurrence_coefficients[3]
    )

k_values = [
    K(right_factor[1](n=index) * j_values[index + 1]
      + right_factor[0](n=index) * j_values[index])
    for index in range(maximum_index)
]
for index in range(2, 16):
    gauge_at_index = gauge.apply_map(lambda entry: entry(n=index))
    assert gauge_at_index.inverse() * vector(
        K, [k_values[index], k_values[index + 1]]
    ) == -10080 * vector(K, [b_values[index], b_values[index + 1]])


# A compact scalar certificate for the distinguished gauge identity.
D = (n - 1)^2 * n^2 * (n + 2)^3 * (n + 3)^2
P0 = (
    -2136 * n^14
    - 35956 * n^13
    - 275781 * n^12
    - 1295151 * n^11
    - 4204216 * n^10
    - 9914534 * n^9
    - 16999496 * n^8
    - 20566704 * n^7
    - 16905458 * n^6
    - 9365542 * n^5
    - 3782083 * n^4
    - 1149761 * n^3
    - 184974 * n^2
    - 56736 * n
    - 30240
)
P1 = (
    984 * n^14
    + 10292 * n^13
    + 50517 * n^12
    + 159519 * n^11
    + 394364 * n^10
    + 806062 * n^9
    + 1086736 * n^8
    + 510804 * n^7
    - 210002 * n^6
    + 898742 * n^5
    + 2061155 * n^4
    + 572629 * n^3
    - 538602 * n^2
    - 3168 * n
    + 6048
)
assert R(12 * D * (-10080 * gauge[0, 0])) == P0
assert R(12 * D * (-10080 * gauge[0, 1])) == P1


# The extension does not split over Q(n).  Seek V=v0+v1*S satisfying
# right_factor*V = U modulo the Apéry operator, where U is the gauge
# map from an Apéry solution to K.
u0 = K(-10080 * gauge[0, 0])
u1 = K(-10080 * gauge[0, 1])
ar = K(right_factor[1])
br = K(-right_factor[0])
a0, a1, a2 = (K(apery_operator[index]) for index in range(3))
split_transition = matrix(K, 3, 3, 0)
split_transition[1, 0] = -a2 * br / (ar * a0)
split_transition[1, 2] = -a2 * u0 / (ar * a0)
split_transition[0, 0] = a1 / a2 * split_transition[1, 0]
split_transition[0, 1] = br / ar
split_transition[0, 2] = u1 / ar + a1 / a2 * split_transition[1, 2]
split_transition[2, 2] = 1

split_rows = [vector(K, [1, 0, 0])]
for _ in range(3):
    split_rows.append(
        vector(K, [sigma(entry) for entry in split_rows[-1]])
        * split_transition
    )
split_krylov = matrix(K, split_rows[:3])
assert split_krylov.det() != 0
split_relation = split_krylov.transpose().solve_right(-split_rows[3])
split_operator = OAK(list(split_relation) + [1]).normalize()
assert split_operator.rational_solutions() == []


# The unnormalized pullback maps J exactly to b.  Modify it by a left
# multiple of P so that it also annihilates every Apéry solution.
H = -10080 * gauge
H_inverse = H.inverse()
pullback = (
    H_inverse[0, 0] * OAK(right_factor)
    + H_inverse[0, 1] * OAK(Sn) * OAK(right_factor)
)
for index in range(2, 16):
    assert sum(
        K(pullback[shift](n=index)) * j_values[index + shift]
        for shift in range(pullback.order() + 1)
    ) == b_values[index]

P_K = OAK(P)
A_K = OAK(apery_operator)
P_remainder = P_K.quo_rem(A_K)[1]
SP_remainder = (OAK(Sn) * P_K).quo_rem(A_K)[1]
pullback_remainder = pullback.quo_rem(A_K)[1]
crt_matrix = matrix(
    K,
    [
        [P_remainder[0], SP_remainder[0]],
        [P_remainder[1], SP_remainder[1]],
    ],
)
assert crt_matrix.det() != 0
q0, q1 = crt_matrix.solve_right(
    vector(K, [-pullback_remainder[0], -pullback_remainder[1]])
)
origin_operator = pullback + (OAK(q0) + OAK(q1) * OAK(Sn)) * P_K
assert origin_operator.quo_rem(A_K)[1] == 0

# Clear polynomial and rational-number denominators, then remove the
# common integer content shared with the multiplier.
denominator_polynomial = lcm(
    K(origin_operator[shift]).denominator()
    for shift in range(origin_operator.order() + 1)
)
polynomial_coefficients = [
    R(denominator_polynomial * K(origin_operator[shift]))
    for shift in range(origin_operator.order() + 1)
]
constant_denominator = lcm(
    coefficient.denominator()
    for polynomial in polynomial_coefficients
    for coefficient in polynomial.list()
)
integer_coefficients = [
    R(constant_denominator * polynomial)
    for polynomial in polynomial_coefficients
]
multiplier = R(constant_denominator * denominator_polynomial)
integer_content = integer_gcd(
    *[
        int(coefficient)
        for polynomial in integer_coefficients + [multiplier]
        for coefficient in polynomial.list()
    ]
)
integer_coefficients = [
    R(polynomial / integer_content) for polynomial in integer_coefficients
]
multiplier = R(multiplier / integer_content)
integral_operator = sum(
    integer_coefficients[shift] * OAK(Sn)^shift
    for shift in range(len(integer_coefficients))
)
assert integral_operator.quo_rem(A_K)[1] == 0
assert (
    integral_operator - OAK(multiplier) * pullback
).quo_rem(P_K)[1] == 0
assert [coefficient.degree() for coefficient in integer_coefficients] == [
    32, 32, 32, 32, 32
]
assert multiplier.degree() == 30
assert gcd(integer_coefficients) == 1
assert integer_content == 1


def shell_pair(moment, prime, modulus):
    return (
        shell_fast(moment, prime - 2, modulus=modulus)
        + shell_fast(moment, prime, modulus=modulus)
    )


# This congruence holds at targets and non-targets alike.
local_cases = (
    (200, 131),
    (200, 139),
    (200, 181),
    (272, 191),
    (300, 227),
    (321, 193),
)
for global_index, prime in local_cases:
    residue = global_index - prime
    assert residue + origin_operator.order() <= prime - 5
    z_value = sum(
        ZZ(integer_coefficients[shift](n=global_index))
        * shell_pair(global_index + shift, prime, prime)
        for shift in range(origin_operator.order() + 1)
    ) % prime
    old_value = shell_fast(global_index - 1, prime - 1, modulus=prime)
    assert old_value == apery(residue) % prime
    assert z_value == ZZ(multiplier(n=global_index)) * old_value % prime


def newton_carrier(values, lower_node, length):
    return sum(
        (-1)^offset
        * binomial(lower_node + offset, offset)
        * binomial(lower_node + length + 1, length - offset)
        * values[lower_node + offset]
        for offset in range(length + 1)
    )


BLOCKS = (
    (200, 128, 63, (139, 181)),
    (272, 180, 63, (191, 233)),
    (300, 180, 57, (191, 227)),
    (321, 168, 53, (179, 193, 211)),
)


def audit_blocks():
    records = []
    for global_index, block_start, block_length, targets in BLOCKS:
        lower_node = block_start - 1
        stencil_nodes = range(
            lower_node, lower_node + block_length + 1
        )
        required_nodes = sorted(
            {
                node + displacement
                for node in stencil_nodes
                for displacement in (-1, 1)
            }
        )
        moment_shells = {
            shift: shell_batch(global_index + shift, required_nodes)
            for shift in range(origin_operator.order() + 1)
        }
        z_values = {
            node: sum(
                ZZ(integer_coefficients[shift](n=global_index))
                * (
                    moment_shells[shift][node - 1]
                    + moment_shells[shift][node + 1]
                )
                for shift in range(origin_operator.order() + 1)
            )
            for node in stencil_nodes
        }
        old_values = shell_batch(global_index - 1, stencil_nodes)
        old_carrier = newton_carrier(
            old_values, lower_node, block_length
        )
        new_carrier = newton_carrier(
            z_values, lower_node, block_length
        )
        candidate_primes = list(
            prime_range(
                block_start + 1,
                block_start + block_length + 1,
            )
        )
        candidate_product = prod(candidate_primes)
        difference = (
            new_carrier
            - ZZ(multiplier(n=global_index)) * old_carrier
        )
        assert difference % candidate_product == 0
        quotient = difference // candidate_product
        target_product = prod(targets)
        assert gcd(old_carrier, quotient) == 1
        assert gcd(old_carrier, new_carrier) == target_product
        target_residues = tuple(
            (prime, quotient % prime) for prime in targets
        )
        assert all(residue != 0 for _, residue in target_residues)
        records.append(
            (
                global_index,
                target_product,
                target_residues,
                (
                    abs(old_carrier).nbits(),
                    abs(new_carrier).nbits(),
                    abs(quotient).nbits(),
                ),
            )
        )
    return records


block_records = audit_blocks() if "--blocks" in sys.argv else []

print("Q32_DOUBLED_PERIOD_GAUGE_AUDIT=PASS")
print("GAUGE_DETERMINANT_NONZERO", gauge.det() != 0)
print("DISTINGUISHED_PULLBACK_SCALE", -10080)
print("EXTENSION_RATIONAL_SPLIT_SOLUTIONS", 0)
print("ORIGIN_OPERATOR_ORDER", integral_operator.order())
print(
    "ORIGIN_OPERATOR_DEGREES",
    [coefficient.degree() for coefficient in integer_coefficients],
)
print("ORIGIN_MULTIPLIER_DEGREE", multiplier.degree())
print("LOCAL_UNIVERSAL_ALIAS_CHECKS", len(local_cases))
if block_records:
    print("HOSTILE_BLOCK_RECORDS", block_records)
