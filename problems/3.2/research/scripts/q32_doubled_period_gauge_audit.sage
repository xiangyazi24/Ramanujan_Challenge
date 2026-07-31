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


def gcd_many(values):
    out = ZZ(0)
    for value in values:
        out = gcd(out, abs(value))
    return out


def audit_cross_minors():
    records = []
    for global_index, D0, N0, targets in BLOCKS:
        moment = global_index - 1
        target_product = prod(targets)
        maximum_m = min(
            8,
            D0 - moment // 2,
            moment - D0 - N0 + 2,
        )
        lower = D0 - maximum_m + 1
        upper = D0 + N0 + maximum_m - 2
        nodes = list(range(lower, upper + 1))
        required_nodes = list(range(lower - 1, upper + 2))
        old_values = shell_batch(moment, nodes)
        moment_shells = {
            shift: shell_batch(global_index + shift, required_nodes)
            for shift in range(origin_operator.order() + 1)
        }
        new_values = {
            node: sum(
                ZZ(integer_coefficients[shift](n=global_index))
                * (
                    moment_shells[shift][node - 1]
                    + moment_shells[shift][node + 1]
                )
                for shift in range(origin_operator.order() + 1)
            )
            for node in nodes
        }

        history = []
        all_minor_gcd = ZZ(0)
        for multiplicity in range(2, maximum_m + 1):
            length = N0 + multiplicity - 2
            d0 = D0 - multiplicity + 1
            old_carriers = [
                newton_carrier(old_values, d0 + shift, length)
                for shift in range(multiplicity)
            ]
            new_carriers = [
                newton_carrier(new_values, d0 + shift, length)
                for shift in range(multiplicity)
            ]
            primitive_minors = []
            for shift in range(multiplicity - 1):
                prefactor = binomial(D0 + N0 + shift, length)
                determinant = (
                    old_carriers[shift] * new_carriers[shift + 1]
                    - old_carriers[shift + 1] * new_carriers[shift]
                )
                assert determinant % prefactor == 0
                primitive_minor = determinant // prefactor
                assert primitive_minor % target_product == 0
                primitive_minors.append(primitive_minor)
                all_minor_gcd = gcd(all_minor_gcd, abs(primitive_minor))
            minor_gcd = gcd_many(primitive_minors)
            assert minor_gcd % target_product == 0
            nuisance = minor_gcd // target_product
            history.append(
                (
                    multiplicity,
                    abs(nuisance).nbits(),
                    min(abs(value).nbits() for value in primitive_minors),
                    max(abs(value).nbits() for value in primitive_minors),
                )
            )
        coefficient_value_gcd = gcd_many(
            ZZ(coefficient(n=global_index))
            for coefficient in integer_coefficients
        )
        records.append(
            (
                global_index,
                (D0, N0),
                targets,
                maximum_m,
                tuple(history),
                nuisance,
                factor(nuisance),
                (
                    gcd(nuisance, ZZ(multiplier(n=global_index))),
                    gcd(nuisance, coefficient_value_gcd),
                    gcd(nuisance, ZZ(apery(global_index - 1))),
                ),
                all_minor_gcd // target_product,
                factor(all_minor_gcd // target_product),
            )
        )
    return records


cross_minor_records = (
    audit_cross_minors() if "--cross-minors" in sys.argv else []
)


def guess_polynomial_recurrence(
    values,
    modulus,
    maximum_order,
    maximum_degree,
    holdout,
):
    field = GF(modulus)
    for order in range(1, maximum_order + 1):
        equation_count = len(values) - order
        for degree in range(maximum_degree + 1):
            columns = (order + 1) * (degree + 1)
            training_count = equation_count - holdout
            if training_count <= columns:
                continue
            rows = []
            for index in range(training_count):
                powers = [1]
                for _ in range(degree):
                    powers.append(powers[-1] * index % modulus)
                rows.append(
                    [
                        values[index + shift] * powers[power] % modulus
                        for shift in range(order + 1)
                        for power in range(degree + 1)
                    ]
                )
            kernel = matrix(field, rows).right_kernel()
            for candidate in kernel.basis():
                if all(
                    sum(
                        sum(
                            ZZ(candidate[
                                shift * (degree + 1) + power
                            ])
                            * index^power
                            for power in range(degree + 1)
                        )
                        * values[index + shift]
                        for shift in range(order + 1)
                    )
                    % modulus
                    == 0
                    for index in range(training_count, equation_count)
                ):
                    return order, degree
    return None


def audit_ghost_recurrence_guess():
    global_index = 400
    modulus = 1_000_000_007
    nodes = list(range(global_index // 2 + 5, global_index - 1))
    required_nodes = list(range(nodes[0] - 1, nodes[-1] + 2))
    old_values = shell_batch(
        global_index - 1, nodes, modulus=modulus
    )
    moment_shells = {
        shift: shell_batch(
            global_index + shift,
            required_nodes,
            modulus=modulus,
        )
        for shift in range(origin_operator.order() + 1)
    }
    multiplier_value = ZZ(multiplier(n=global_index)) % modulus
    w_values = [
        (
            sum(
                ZZ(integer_coefficients[shift](n=global_index))
                * (
                    moment_shells[shift][node - 1]
                    + moment_shells[shift][node + 1]
                )
                for shift in range(origin_operator.order() + 1)
            )
            - multiplier_value * old_values[node]
        )
        % modulus
        for node in nodes
    ]
    difference_order = 60
    ghost_values = list(w_values)
    for _ in range(difference_order):
        ghost_values = [
            (ghost_values[index + 1] - ghost_values[index]) % modulus
            for index in range(len(ghost_values) - 1)
        ]
    return (
        len(w_values),
        guess_polynomial_recurrence(
            w_values, modulus, 12, 10, 24
        ),
        len(ghost_values),
        guess_polynomial_recurrence(
            ghost_values, modulus, 12, 8, 16
        ),
    )


ghost_guess_record = (
    audit_ghost_recurrence_guess()
    if "--guess-ghost" in sys.argv
    else None
)


def audit_generic_cross_gcds():
    rows = []
    for global_index in (80, 100, 120, 150, 180, 220):
        moment = global_index - 1
        D0 = (13 * global_index) // 20
        N0 = global_index // 5
        d0 = D0 - 2
        length = N0 + 1
        nodes = list(range(d0, d0 + length + 3))
        required_nodes = list(range(nodes[0] - 1, nodes[-1] + 2))
        old_values = shell_batch(moment, nodes)
        moment_shells = {
            shift: shell_batch(global_index + shift, required_nodes)
            for shift in range(origin_operator.order() + 1)
        }
        new_values = {
            node: sum(
                ZZ(integer_coefficients[shift](n=global_index))
                * (
                    moment_shells[shift][node - 1]
                    + moment_shells[shift][node + 1]
                )
                for shift in range(origin_operator.order() + 1)
            )
            for node in nodes
        }
        old_carriers = [
            newton_carrier(old_values, d0 + shift, length)
            for shift in range(3)
        ]
        new_carriers = [
            newton_carrier(new_values, d0 + shift, length)
            for shift in range(3)
        ]
        primitive_minors = []
        for shift in range(2):
            prefactor = binomial(d0 + shift + length + 1, length)
            determinant = (
                old_carriers[shift] * new_carriers[shift + 1]
                - old_carriers[shift + 1] * new_carriers[shift]
            )
            assert determinant % prefactor == 0
            primitive_minors.append(determinant // prefactor)
        common = gcd_many(primitive_minors)
        rows.append(
            (
                global_index,
                (D0, N0),
                abs(common).nbits(),
                factor(common),
                tuple(abs(value).nbits() for value in primitive_minors),
            )
        )
    return rows


generic_cross_records = (
    audit_generic_cross_gcds() if "--scan-cross" in sys.argv else []
)


def full_margin_data(
    global_index, D0, N0, include_state=False, margin=None
):
    """Return the exact shell pair and primitive minors at one margin.

    With ``margin=None`` this uses the maximal admissible margin, as in
    the original full-margin audit.  Passing an integer permits an exact
    comparison of the endpoint carriers attached to the same retained
    core but to different surrounding Newton cells.
    """

    moment = global_index - 1
    available_margin = min(
        D0 - moment // 2,
        moment - D0 - N0 + 2,
    )
    assert available_margin >= 2
    maximum_m = available_margin if margin is None else ZZ(margin)
    assert 2 <= maximum_m <= available_margin
    d0 = D0 - maximum_m + 1
    length = N0 + maximum_m - 2

    # The high differences use W on [d0,d0+length+maximum_m-1].
    nodes = list(range(d0, d0 + length + maximum_m))
    required_nodes = list(range(nodes[0] - 1, nodes[-1] + 2))
    old_values = shell_batch(moment, nodes)
    moment_shells = {
        shift: shell_batch(global_index + shift, required_nodes)
        for shift in range(origin_operator.order() + 1)
    }
    new_values = {
        node: sum(
            ZZ(integer_coefficients[shift](n=global_index))
            * (
                moment_shells[shift][node - 1]
                + moment_shells[shift][node + 1]
            )
            for shift in range(origin_operator.order() + 1)
        )
        for node in nodes
    }
    multiplier_value = ZZ(multiplier(n=global_index))
    w_values = {
        node: new_values[node] - multiplier_value * old_values[node]
        for node in nodes
    }

    old_carriers = [
        newton_carrier(old_values, d0 + shift, length)
        for shift in range(maximum_m)
    ]
    w_carriers = [
        newton_carrier(w_values, d0 + shift, length)
        for shift in range(maximum_m)
    ]
    primitive_minors = []
    for shift in range(maximum_m - 1):
        prefactor = binomial(d0 + shift + length + 1, length)
        determinant = (
            old_carriers[shift] * w_carriers[shift + 1]
            - old_carriers[shift + 1] * w_carriers[shift]
        )
        assert determinant % prefactor == 0
        primitive_minors.append(determinant // prefactor)

    # Compute all Delta^(length+1) W values simultaneously.
    high_differences = [w_values[node] for node in nodes]
    for _ in range(length + 1):
        high_differences = [
            high_differences[index + 1] - high_differences[index]
            for index in range(len(high_differences) - 1)
        ]
    assert len(high_differences) == maximum_m - 1

    old_high_differences = [old_values[node] for node in nodes]
    for _ in range(length + 1):
        old_high_differences = [
            old_high_differences[index + 1]
            - old_high_differences[index]
            for index in range(len(old_high_differences) - 1)
        ]
    assert len(old_high_differences) == maximum_m - 1

    if include_state:
        return (
            maximum_m,
            d0,
            length,
            primitive_minors,
            high_differences,
            old_carriers,
            w_carriers,
            old_high_differences,
        )

    return (
        maximum_m,
        d0,
        length,
        primitive_minors,
        high_differences,
    )


def multi_margin_endpoint_record(global_index, D0, N0):
    """Return one exact multi-margin endpoint-gcd record."""

    moment = global_index - 1
    available_margin = min(
        D0 - moment // 2,
        moment - D0 - N0 + 2,
    )
    assert available_margin >= 2
    candidate_primes = list(
        prime_range(D0 + 1, D0 + N0 + 1)
    )
    targets = tuple(
        prime
        for prime in candidate_primes
        if ZZ(apery(global_index - prime)) % prime == 0
    )
    target_product = prod(targets)

    # Every smaller two-sided window is contained in the maximal one.
    # Computing the distinguished shell arrays once per case makes the
    # dense asymptotic audit practical.
    lower_node = D0 - available_margin + 1
    upper_node = D0 + N0 + available_margin - 2
    nodes = list(range(lower_node, upper_node + 1))
    required_nodes = list(range(lower_node - 1, upper_node + 2))
    old_values = shell_batch(moment, nodes)
    moment_shells = {
        shift: shell_batch(global_index + shift, required_nodes)
        for shift in range(origin_operator.order() + 1)
    }
    new_values = {
        node: sum(
            ZZ(integer_coefficients[shift](n=global_index))
            * (
                moment_shells[shift][node - 1]
                + moment_shells[shift][node + 1]
            )
            for shift in range(origin_operator.order() + 1)
        )
        for node in nodes
    }
    multiplier_value = ZZ(multiplier(n=global_index))
    w_values = {
        node: new_values[node] - multiplier_value * old_values[node]
        for node in nodes
    }

    normalized_endpoints = []
    endpoint_bits = []
    common_weight_bits = []
    target_valuations = []
    for margin in range(2, available_margin + 1):
        d0 = D0 - margin + 1
        length = N0 + margin - 2
        old_left = newton_carrier(old_values, d0, length)
        old_right = newton_carrier(old_values, D0, length)
        w_left = newton_carrier(w_values, d0, length)
        w_right = newton_carrier(w_values, D0, length)
        endpoint = old_left * w_right - w_left * old_right
        weights = [
            binomial(d0 + shift + length + 1, length)
            for shift in range(margin - 1)
        ]
        common_weight = gcd_many(weights)
        assert endpoint % common_weight == 0
        normalized = endpoint // common_weight
        assert normalized % target_product == 0
        normalized_endpoints.append(normalized)
        endpoint_bits.append(abs(normalized).nbits())
        common_weight_bits.append(abs(common_weight).nbits())
        target_valuations.append(
            tuple(
                (
                    valuation(endpoint, prime),
                    valuation(common_weight, prime),
                    valuation(normalized, prime),
                )
                for prime in targets
            )
        )

    running_gcd_bits = []
    running_gcd_values = []
    running_gcd = ZZ(0)
    for endpoint in normalized_endpoints:
        running_gcd = gcd(running_gcd, endpoint)
        running_gcd_values.append(running_gcd)
        running_gcd_bits.append(abs(running_gcd).nbits())
    assert running_gcd % target_product == 0
    residual = running_gcd // target_product
    stabilization_margin = 2 + next(
        index
        for index, value in enumerate(running_gcd_values)
        if value == running_gcd
    )
    reverse_gcd_values = []
    reverse_gcd = ZZ(0)
    for endpoint in reversed(normalized_endpoints):
        reverse_gcd = gcd(reverse_gcd, endpoint)
        reverse_gcd_values.append(reverse_gcd)
    suffix_length = 1 + next(
        index
        for index, value in enumerate(reverse_gcd_values)
        if value == running_gcd
    )
    reverse_gcd_bits = tuple(
        abs(value).nbits() for value in reverse_gcd_values
    )
    two_suffix_residual = reverse_gcd_values[1] // target_product
    assert reverse_gcd_values[1] % target_product == 0
    origin_value_content = gcd_many(
        [ZZ(multiplier(n=global_index))]
        + [
            ZZ(coefficient(n=global_index))
            for coefficient in integer_coefficients
        ]
    )
    assert reverse_gcd_values[1] % origin_value_content == 0
    two_suffix_primitive = (
        reverse_gcd_values[1] // origin_value_content
    )
    return (
        global_index,
        (D0, N0),
        available_margin,
        targets,
        tuple(endpoint_bits),
        tuple(common_weight_bits),
        tuple(running_gcd_bits),
        abs(running_gcd).nbits(),
        abs(residual).nbits(),
        factor(residual) if abs(residual).nbits() <= 256 else None,
        stabilization_margin,
        suffix_length,
        reverse_gcd_bits,
        abs(reverse_gcd_values[1]).nbits(),
        abs(two_suffix_residual).nbits(),
        (
            factor(two_suffix_residual)
            if abs(two_suffix_residual).nbits() <= 256
            else None
        ),
        factor(origin_value_content),
        abs(two_suffix_primitive).nbits(),
        (
            factor(two_suffix_primitive)
            if abs(two_suffix_primitive).nbits() <= 256
            else None
        ),
        tuple(target_valuations),
    )


def audit_multi_margin_endpoint_gcd(dense=False):
    """Intersect normalized endpoint carriers over every valid margin."""

    cases = [
        (80, 52, 16),
        (120, 78, 24),
        (160, 104, 32),
        (200, 130, 40),
        (200, 128, 63),
        (272, 180, 63),
        (300, 180, 57),
        (321, 168, 53),
    ]
    if dense:
        cases.extend(
            (
                global_index,
                (13 * global_index) // 20,
                global_index // 5,
            )
            for global_index in range(60, 501, 20)
        )
    records = []
    for global_index, D0, N0 in dict.fromkeys(cases):
        records.append(
            multi_margin_endpoint_record(
                global_index, D0, N0
            )
        )
    return tuple(records)


def audit_full_margin_curvature_geometry():
    """Measure the exact boundary and curvature terms from (68.99)."""

    cases = (
        (80, 52, 16),
        (120, 78, 24),
        (160, 104, 32),
        (200, 130, 40),
        (200, 128, 63),
        (272, 180, 63),
        (300, 180, 57),
        (321, 168, 53),
    )
    records = []
    for global_index, D0, N0 in cases:
        (
            maximum_m,
            d0,
            length,
            primitive_minors,
            w_high_differences,
            old_carriers,
            w_carriers,
            old_high_differences,
        ) = full_margin_data(
            global_index, D0, N0, include_state=True
        )
        packet_vectors = [
            (
                (-1) ** (length + 1)
                * old_high_differences[shift],
                (-1) ** (length + 1)
                * w_high_differences[shift],
            )
            for shift in range(maximum_m - 1)
        ]
        carrier_vectors = list(zip(old_carriers, w_carriers))
        weights = [
            binomial(d0 + shift + length + 1, length)
            for shift in range(maximum_m - 1)
        ]
        for shift, edge in enumerate(primitive_minors):
            assert edge == (
                packet_vectors[shift][0]
                * carrier_vectors[shift][1]
                - packet_vectors[shift][1]
                * carrier_vectors[shift][0]
            )

        endpoint = (
            carrier_vectors[0][0] * carrier_vectors[-1][1]
            - carrier_vectors[0][1] * carrier_vectors[-1][0]
        )
        visible = sum(
            weights[shift] * primitive_minors[shift]
            for shift in range(maximum_m - 1)
        )
        curvature = sum(
            weights[left]
            * weights[right]
            * (
                packet_vectors[right][0]
                * packet_vectors[left][1]
                - packet_vectors[right][1]
                * packet_vectors[left][0]
            )
            for right in range(maximum_m - 1)
            for left in range(right)
        )
        assert endpoint == visible + curvature
        common_weight = gcd_many(weights)
        assert endpoint % common_weight == 0

        candidate_primes = list(
            prime_range(D0 + 1, D0 + N0 + 1)
        )
        targets = tuple(
            prime
            for prime in candidate_primes
            if ZZ(apery(global_index - prime)) % prime == 0
        )
        target_product = prod(targets)
        assert endpoint % (target_product**2) == 0
        records.append(
            (
                global_index,
                (D0, N0),
                maximum_m,
                length,
                targets,
                (
                    abs(endpoint).nbits(),
                    abs(visible).nbits(),
                    abs(curvature).nbits(),
                ),
                (
                    abs(common_weight).nbits(),
                    abs(endpoint // common_weight).nbits(),
                ),
                (
                    min(abs(value).nbits()
                        for value in primitive_minors),
                    max(abs(value).nbits()
                        for value in primitive_minors),
                ),
                abs(endpoint // (target_product**2)).nbits(),
            )
        )
    return tuple(records)


def audit_full_margin_candidate_separation():
    """Test exact separation only on prime nodes in the retained core.

    For a non-target candidate prime q, (68.60) says that q divides all
    primitive cross-minors if and only if every displayed high
    difference vanishes modulo q.  Unlike the height of the whole gcd,
    this candidate-restricted condition ignores near-n and half-scale
    Kummer nuisances outside the core.
    """

    central_indices = (
        range(60, 201)
        if "--candidate-dense" in sys.argv
        else range(60, 401, 20)
    )
    central_cases = [
        (
            global_index,
            (13 * global_index) // 20,
            global_index // 5,
        )
        for global_index in central_indices
    ]
    hostile_cases = [
        (global_index, D0, N0)
        for global_index, D0, N0, _ in BLOCKS
        if global_index <= 400
    ]
    cases = list(dict.fromkeys(central_cases + hostile_cases))

    records = []
    total_candidates = 0
    total_targets = 0
    total_exceptions = 0
    total_first_failures = 0
    total_last_failures = 0
    total_endpoint_pair_failures = 0
    total_active_boundary_failures = 0
    for global_index, D0, N0 in cases:
        (
            maximum_m,
            d0,
            length,
            primitive_minors,
            high_differences,
        ) = full_margin_data(global_index, D0, N0)
        candidate_primes = list(prime_range(D0 + 1, D0 + N0 + 1))
        candidate_product = prod(candidate_primes)
        targets = tuple(
            prime
            for prime in candidate_primes
            if ZZ(apery(global_index - prime)) % prime == 0
        )
        exceptions = tuple(
            prime
            for prime in candidate_primes
            if prime not in targets
            and all(value % prime == 0 for value in high_differences)
        )
        first_failures = tuple(
            prime
            for prime in candidate_primes
            if prime not in targets
            and high_differences[0] % prime == 0
        )
        last_failures = tuple(
            prime
            for prime in candidate_primes
            if prime not in targets
            and high_differences[-1] % prime == 0
        )
        endpoint_pair_failures = tuple(
            prime
            for prime in candidate_primes
            if prime not in targets
            and high_differences[0] % prime == 0
            and high_differences[-1] % prime == 0
        )
        left_margin = D0 - (global_index - 1) // 2
        right_margin = (global_index - 1) - D0 - N0 + 2
        active_indices = []
        if left_margin <= right_margin:
            active_indices.append(0)
        if right_margin <= left_margin:
            active_indices.append(-1)
        active_boundary_failures = tuple(
            prime
            for prime in candidate_primes
            if prime not in targets
            and all(
                high_differences[index] % prime == 0
                for index in active_indices
            )
        )
        common_candidate_part = gcd(
            candidate_product, gcd_many(primitive_minors)
        )
        all_minor_gcd = gcd_many(primitive_minors)
        endpoint_minor_gcd = gcd(
            abs(primitive_minors[0]),
            abs(primitive_minors[-1]),
        )
        origin_value_content = gcd_many(
            [ZZ(multiplier(n=global_index))]
            + [
                ZZ(coefficient(n=global_index))
                for coefficient in integer_coefficients
            ]
        )
        assert all_minor_gcd % origin_value_content == 0
        assert endpoint_minor_gcd % origin_value_content == 0
        boundary_shift_product = prod(
            global_index + shift for shift in range(-6, 7)
        )
        all_residual = all_minor_gcd // (
            prod(targets) * origin_value_content
        )
        endpoint_residual = endpoint_minor_gcd // (
            prod(targets) * origin_value_content
        )
        all_residual_radical = prod(
            prime for prime, _ in factor(all_residual)
        )
        endpoint_residual_radical = prod(
            prime for prime, _ in factor(endpoint_residual)
        )
        assert boundary_shift_product % all_residual_radical == 0
        assert boundary_shift_product % endpoint_residual_radical == 0
        expected = prod(targets + exceptions)
        assert common_candidate_part == expected

        # Independently check the integral congruence (68.60), including
        # its sign, for every candidate and every translated minor.
        for prime in candidate_primes:
            residue = ZZ(apery(global_index - prime)) % prime
            for minor, difference in zip(
                primitive_minors, high_differences
            ):
                assert (
                    minor
                    - (-1)^length * residue * difference
                ) % prime == 0

        total_candidates += len(candidate_primes)
        total_targets += len(targets)
        total_exceptions += len(exceptions)
        total_first_failures += len(first_failures)
        total_last_failures += len(last_failures)
        total_endpoint_pair_failures += len(endpoint_pair_failures)
        total_active_boundary_failures += len(active_boundary_failures)
        records.append(
            (
                global_index,
                (D0, N0),
                maximum_m,
                length,
                len(candidate_primes),
                targets,
                exceptions,
                common_candidate_part,
                first_failures,
                last_failures,
                endpoint_pair_failures,
                tuple(active_indices),
                active_boundary_failures,
                abs(all_minor_gcd).nbits(),
                abs(endpoint_minor_gcd).nbits(),
                factor(all_minor_gcd // prod(targets)),
                factor(endpoint_minor_gcd // prod(targets)),
                (
                    min(abs(value).nbits() for value in high_differences),
                    max(abs(value).nbits() for value in high_differences),
                ),
                factor(origin_value_content),
                factor(all_residual),
                factor(endpoint_residual),
            )
        )

    return (
        total_candidates,
        total_targets,
        total_exceptions,
        total_first_failures,
        total_last_failures,
        total_endpoint_pair_failures,
        total_active_boundary_failures,
        tuple(records),
    )


full_margin_candidate_records = (
    audit_full_margin_candidate_separation()
    if (
        "--candidate-scan" in sys.argv
        or "--candidate-dense" in sys.argv
    )
    else None
)

full_margin_curvature_records = (
    audit_full_margin_curvature_geometry()
    if "--curvature-scan" in sys.argv
    else None
)

multi_margin_endpoint_records = (
    audit_multi_margin_endpoint_gcd(
        dense="--multi-margin-dense" in sys.argv
    )
    if (
        "--multi-margin-scan" in sys.argv
        or "--multi-margin-dense" in sys.argv
    )
    else None
)

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
if "--origin-factors" in sys.argv:
    print("ORIGIN_MULTIPLIER_FACTOR", factor(multiplier))
    for shift, coefficient in enumerate(integer_coefficients):
        print(
            "ORIGIN_COEFFICIENT_FACTOR",
            shift,
            factor(coefficient),
        )
print("LOCAL_UNIVERSAL_ALIAS_CHECKS", len(local_cases))
if block_records:
    print("HOSTILE_BLOCK_RECORDS", block_records)
if cross_minor_records:
    print("CROSS_MINOR_RECORDS", cross_minor_records)
if ghost_guess_record is not None:
    print("GHOST_RECURRENCE_GUESS", ghost_guess_record)
if generic_cross_records:
    print("GENERIC_CROSS_GCDS", generic_cross_records)
if full_margin_candidate_records is not None:
    print(
        "FULL_MARGIN_CANDIDATE_SEPARATION",
        full_margin_candidate_records,
    )
if full_margin_curvature_records is not None:
    print(
        "FULL_MARGIN_CURVATURE_GEOMETRY",
        full_margin_curvature_records,
    )
if multi_margin_endpoint_records is not None:
    if "--multi-margin-dense" in sys.argv:
        checkpoint_margins = (2, 3, 4, 5, 6, 8, 12, 16, 24, 32, 48, 64)
        displayed_multi_margin_records = tuple(
            (
                record[0],
                record[1],
                record[2],
                record[3],
                tuple(
                    (margin, record[6][margin - 2])
                    for margin in checkpoint_margins
                    if margin <= record[2]
                ),
                record[7],
                record[8],
                record[9],
                record[10],
                record[11],
                tuple(
                    (length, record[12][length - 1])
                    for length in (1, 2, 3, 4, 6, 8, 12, 16)
                    if length <= len(record[12])
                ),
                (record[13], record[14], record[15]),
                (record[16], record[17], record[18]),
            )
            for record in multi_margin_endpoint_records
        )
    else:
        displayed_multi_margin_records = multi_margin_endpoint_records
    print(
        "MULTI_MARGIN_ENDPOINT_GCD",
        displayed_multi_margin_records,
    )
