#!/usr/bin/env sage
"""Audit the fixed origin/shift-one ghost minor at the first divided digit.

For a candidate prime ``p`` in a Newton core, every origin-cancelled or
common-annihilator ghost carrier is divisible by ``p``.  Dividing the two
endpoint carriers by ``p`` and reducing modulo ``p`` is equivalent, up to
the common unit ``Q_I / p``, to dividing by the complete core primorial.

This script computes those first divided digits directly modulo ``p^2``.
It audits the fixed minor

    det(R_{-1,m}, R_{1,m}),

where ``R_{-1,m}`` is the origin-cancelled column and ``R_{1,m}`` is the
shift-one common-annihilator ghost.  It also verifies the closed formula
for the divided Newton functional in ``divided_newton_formula``.

Examples:

    sage q32_fixed_minor_jet_audit.sage --central 300
    sage q32_fixed_minor_jet_audit.sage --random 40 --seed 320
    sage q32_fixed_minor_jet_audit.sage \
        --central-all 200 --all-ghosts --allow-failures

All shell arithmetic is performed modulo ``p^2``; no large exact shell
integers are constructed.
"""

from pathlib import Path
import argparse
import random
import sys

from ore_algebra import *


HERE = Path(__file__).resolve().parent
load(str(HERE / "q32_doubled_period_gauge_audit.sage"))


def primitive_integral_coefficients(operator):
    coefficients = [K(operator[i]) for i in range(operator.order() + 1)]
    polynomial_denominator = lcm(value.denominator() for value in coefficients)
    polynomials = [R(polynomial_denominator * value) for value in coefficients]
    rational_denominator = lcm(
        coefficient.denominator()
        for polynomial in polynomials
        for coefficient in polynomial.list()
    )
    polynomials = [R(rational_denominator * value) for value in polynomials]
    content = gcd(
        [
            ZZ(coefficient)
            for polynomial in polynomials
            for coefficient in polynomial.list()
        ]
    )
    answer = tuple(R(value / content) for value in polynomials)
    assert gcd(answer) == 1
    return answer


common_operator = P.lclm(apery_operator)
ghost_coefficients = primitive_integral_coefficients(common_operator)
integral_common_operator = sum(
    ghost_coefficients[shift] * OAK(Sn)^shift
    for shift in range(len(ghost_coefficients))
)
twisted_apery_operator = sum(
    K(apery_operator[shift])
    / K(multiplier(n=n + shift))
    * OAK(Sn)^shift
    for shift in range(apery_operator.order() + 1)
)
ore_quotient, ore_remainder = (
    (twisted_apery_operator * integral_operator).quo_rem(
        integral_common_operator
    )
)
assert ore_remainder == 0
assert ore_quotient.order() == 1


def margin_data(index, D, N):
    moment = index - 1
    maximum_margin = min(
        D - moment // 2,
        moment - D - N + 2,
    )
    assert maximum_margin >= 2
    return maximum_margin, (maximum_margin - 1, maximum_margin)


def newton_mod(values, lower_node, length, modulus):
    return ZZ(
        sum(
            (-1)^offset
            * binomial(lower_node + offset, offset)
            * binomial(
                lower_node + length + 1,
                length - offset,
            )
            * values[lower_node + offset]
            for offset in range(length + 1)
        )
        % modulus
    )


def divided_newton_formula(values, lower_node, length, prime):
    """Return ``G_{d,L}(F)/p mod p`` from the first divided node digit.

    This applies when ``L<p``, ``p-1=d+j`` is in the stencil, and
    ``F_{p-1}=0 mod p``.  Put ``s=L-j``.  The off-node coefficients are
    the once-divided Newton weights:

      i<j: (-1)^(j-i-1) C(j,i) s!(j-i-1)!/(s+j-i)!,
      i=j+t: (-1)^t C(s,t) (t-1)!j!/(j+t)!.
    """

    j = prime - 1 - lower_node
    s = length - j
    assert 0 <= j <= length < prime
    assert values[prime - 1] % prime == 0
    answer = (values[prime - 1] // prime) % prime

    for i in range(j):
        numerator = (
            (-1)^(j - i - 1)
            * binomial(j, i)
            * factorial(s)
            * factorial(j - i - 1)
        )
        denominator = factorial(s + j - i)
        answer += (
            numerator
            * inverse_mod(denominator, prime)
            * values[lower_node + i]
        )

    for t in range(1, s + 1):
        numerator = (
            (-1)^t
            * binomial(s, t)
            * factorial(t - 1)
            * factorial(j)
        )
        denominator = factorial(j + t)
        answer += (
            numerator
            * inverse_mod(denominator, prime)
            * values[prime - 1 + t]
        )

    return ZZ(answer % prime)


def rational_mod(value, prime):
    value = QQ(value)
    numerator = ZZ(value.numerator()) % prime
    denominator = ZZ(value.denominator()) % prime
    if denominator == 0:
        return None
    return ZZ(numerator * inverse_mod(denominator, prime) % prime)


def apery_mod(index, prime):
    return ZZ(apery(index) % prime)


def build_node_functions(index, D, N, modulus):
    maximum_margin, margins = margin_data(index, D, N)

    terminal_d = D - maximum_margin + 1
    terminal_length = N + maximum_margin - 2
    lower_node = terminal_d
    upper_node = D + terminal_length
    nodes = tuple(range(lower_node, upper_node + 1))
    neighbours = tuple(range(lower_node - 1, upper_node + 2))

    y_values = shell_batch(index - 1, nodes, modulus=modulus)
    maximum_ghost_shift = 4 if "--all-ghosts" in sys.argv else 1
    maximum_time_index = max(index + 6, index + maximum_ghost_shift + 5)
    time_shells = {
        time_index: shell_batch(
            time_index,
            neighbours,
            modulus=modulus,
        )
        for time_index in range(index, maximum_time_index + 1)
    }

    def shell_pair(time_index, node):
        return (
            time_shells[time_index][node - 1]
            + time_shells[time_index][node + 1]
        ) % modulus

    shifted_y = {
        0: y_values,
        1: {node: time_shells[index][node] for node in nodes},
        2: {node: time_shells[index + 1][node] for node in nodes},
    }
    origin_shifts = {
        time_shift: {
            node: ZZ(
                (
                    sum(
                        ZZ(
                            integer_coefficients[operator_shift](
                                n=index + time_shift
                            )
                        )
                        * shell_pair(
                            index + time_shift + operator_shift,
                            node,
                        )
                        for operator_shift in range(
                            len(integer_coefficients)
                        )
                    )
                    - ZZ(multiplier(n=index + time_shift))
                    * shifted_y[time_shift][node]
                )
                % modulus
            )
            for node in nodes
        }
        for time_shift in range(3)
    }
    ghost_shifts = {
        time_shift: {
            node: ZZ(
                sum(
                    ZZ(
                        ghost_coefficients[operator_shift](
                            n=index + time_shift
                        )
                    )
                    * shell_pair(
                        index + time_shift + operator_shift,
                        node,
                    )
                    for operator_shift in range(
                        len(ghost_coefficients)
                    )
                )
                % modulus
            )
            for node in nodes
        }
        for time_shift in range(maximum_ghost_shift + 1)
    }
    apery_defect = {
        node: ZZ(
            sum(
                ZZ(apery_operator[time_shift](n=index))
                * shifted_y[time_shift][node]
                for time_shift in range(3)
            )
            % modulus
        )
        for node in nodes
    }
    return (
        margins,
        nodes,
        origin_shifts[0],
        ghost_shifts[1],
        {
            "origin_shifts": origin_shifts,
            "ghost_zero": ghost_shifts[0],
            "ghost_shifts": ghost_shifts,
            "apery_defect": apery_defect,
        },
    )


def audit_candidate(index, D, N, prime, node_functions=None):
    assert D < prime <= D + N
    assert index - prime + 9 <= prime - 5
    modulus = prime^2
    if node_functions is None:
        (
            margins,
            nodes,
            origin,
            ghost_one,
            auxiliary,
        ) = build_node_functions(
            index,
            D,
            N,
            modulus,
        )
    else:
        margins, nodes, origin, ghost_one, auxiliary = node_functions

    assert origin[prime - 1] % prime == 0
    assert ghost_one[prime - 1] % prime == 0

    residue = index - prime
    records = []
    ore_elimination_records = []
    all_ghost_records = []
    for margin in margins:
        d = D - margin + 1
        length = N + margin - 2
        assert length < prime
        columns = []
        for values in (origin, ghost_one):
            left = newton_mod(values, d, length, modulus)
            right = newton_mod(values, D, length, modulus)
            assert left % prime == 0
            assert right % prime == 0
            divided = (left // prime % prime, right // prime % prime)
            formula = (
                divided_newton_formula(values, d, length, prime),
                divided_newton_formula(values, D, length, prime),
            )
            assert divided == formula
            columns.append(divided)

        determinant = ZZ(
            (
                columns[0][0] * columns[1][1]
                - columns[0][1] * columns[1][0]
            )
            % prime
        )
        records.append((margin, determinant, tuple(columns)))

        if "--all-ghosts" in sys.argv:
            all_columns = [columns[0]]
            all_columns.extend(
                (
                    divided_newton_formula(
                        auxiliary["ghost_shifts"][shift],
                        d,
                        length,
                        prime,
                    ),
                    divided_newton_formula(
                        auxiliary["ghost_shifts"][shift],
                        D,
                        length,
                        prime,
                    ),
                )
                for shift in range(5)
            )
            all_minors = tuple(
                ZZ(
                    (
                        all_columns[left][0]
                        * all_columns[right][1]
                        - all_columns[left][1]
                        * all_columns[right][0]
                    )
                    % prime
                )
                for left in range(6)
                for right in range(left + 1, 6)
            )
            all_ghost_records.append(
                (
                    margin,
                    any(all_minors),
                    tuple(
                        index
                        for index, value in enumerate(all_minors)
                        if value
                    ),
                )
            )

        alpha = tuple(
            rational_mod(
                twisted_apery_operator[shift](n=residue),
                prime,
            )
            for shift in range(3)
        )
        q_values = tuple(
            rational_mod(ore_quotient[shift](n=residue), prime)
            for shift in range(2)
        )
        if None not in alpha + q_values:
            def divided_vector(values):
                return vector(
                    GF(prime),
                    [
                        divided_newton_formula(
                            values,
                            lower,
                            length,
                            prime,
                        )
                        for lower in (d, D)
                    ],
                )

            origin_vectors = tuple(
                divided_vector(auxiliary["origin_shifts"][shift])
                for shift in range(3)
            )
            ghost_zero_vector = divided_vector(
                auxiliary["ghost_zero"]
            )
            ghost_one_vector = divided_vector(ghost_one)
            defect_vector = divided_vector(
                auxiliary["apery_defect"]
            )
            left_side = sum(
                (
                    GF(prime)(alpha[shift])
                    * origin_vectors[shift]
                    for shift in range(3)
                ),
                vector(GF(prime), [0, 0]),
            ) + defect_vector
            right_side = (
                GF(prime)(q_values[0]) * ghost_zero_vector
                + GF(prime)(q_values[1]) * ghost_one_vector
            )
            assert left_side == right_side

            def wedge(left, right):
                return ZZ(
                    left[0] * right[1] - left[1] * right[0]
                )

            wedge_residuals = (
                ZZ(
                    GF(prime)(alpha[1])
                    * wedge(origin_vectors[0], origin_vectors[1])
                ),
                ZZ(
                    GF(prime)(alpha[2])
                    * wedge(origin_vectors[0], origin_vectors[2])
                ),
                wedge(origin_vectors[0], defect_vector),
                ZZ(
                    -GF(prime)(q_values[0])
                    * wedge(origin_vectors[0], ghost_zero_vector)
                ),
            )
            assert (
                sum(wedge_residuals)
                - GF(prime)(q_values[1]) * determinant
            ) % prime == 0
            ore_elimination_records.append(
                (
                    margin,
                    alpha,
                    q_values,
                    tuple(value % prime for value in wedge_residuals),
                )
            )

    q0 = rational_mod(ore_quotient[0](n=residue), prime)
    q1 = rational_mod(ore_quotient[1](n=residue), prime)
    record = {
        "n": index,
        "core": (D, N),
        "p": prime,
        "r": residue,
        "target": apery_mod(residue, prime) == 0,
        "b_next": apery_mod(residue + 1, prime),
        "ore_q": (q0, q1),
        "records": tuple(records),
        "ore_elimination": tuple(ore_elimination_records),
        "all_ghost_rank_two": tuple(all_ghost_records),
    }
    if "--profiles" in sys.argv:
        record["support_mod_p"] = {
            "origin": tuple(
                node - (prime - 1)
                for node in nodes
                if origin[node] % prime
            ),
            "ghost_one": tuple(
                node - (prime - 1)
                for node in nodes
                if ghost_one[node] % prime
            ),
        }
    return record


def audit_block(index, D, N):
    candidates = tuple(
        prime
        for prime in prime_range(D + 1, D + N + 1)
        if index - prime + 9 <= prime - 5
    )
    if not candidates:
        return ()
    modulus = prod(prime^2 for prime in candidates)
    node_functions = build_node_functions(index, D, N, modulus)
    return tuple(
        audit_candidate(index, D, N, prime, node_functions)
        for prime in candidates
    )


def central_cases(limit, step=20):
    return tuple(
        (index, (13 * index) // 20, index // 5)
        for index in range(60, limit + 1, step)
    )


def random_cases(count, seed, maximum_index):
    generator = random.Random(seed)
    cases = []
    while len(cases) < count:
        index = generator.randrange(80, maximum_index + 1)
        N = generator.randrange(max(8, index // 12), max(9, index // 4))
        lower_D = (index - 1) // 2 + 3
        upper_D = index - N - 2
        if lower_D > upper_D:
            continue
        D = generator.randrange(lower_D, upper_D + 1)
        maximum_margin, _ = margin_data(index, D, N)
        if maximum_margin < 2:
            continue
        candidates = tuple(
            prime
            for prime in prime_range(D + 1, D + N + 1)
            if index - prime + 9 <= prime - 5
        )
        if not candidates:
            continue
        prime = generator.choice(candidates)
        cases.append((index, D, N, prime))
    return tuple(cases)


def random_blocks(count, seed, maximum_index):
    generator = random.Random(seed)
    blocks = []
    while len(blocks) < count:
        index = generator.randrange(80, maximum_index + 1)
        N = generator.randrange(max(8, index // 12), max(9, index // 4))
        lower_D = (index - 1) // 2 + 3
        upper_D = index - N - 2
        if lower_D > upper_D:
            continue
        D = generator.randrange(lower_D, upper_D + 1)
        maximum_margin, _ = margin_data(index, D, N)
        if maximum_margin < 2:
            continue
        if not any(
            index - prime + 9 <= prime - 5
            for prime in prime_range(D + 1, D + N + 1)
        ):
            continue
        blocks.append((index, D, N))
    return tuple(blocks)


def balanced_blocks(count, seed, maximum_index):
    """Random cores whose left and right first-cell margins differ by <= 1."""

    generator = random.Random(seed)
    blocks = []
    while len(blocks) < count:
        index = generator.randrange(80, maximum_index + 1)
        N = generator.randrange(max(8, index // 12), max(9, index // 4))
        midpoint = (index - 1) // 2
        # Solve D-midpoint = index+1-D-N up to the nearest integer.
        D = (index + 1 + midpoint - N) // 2
        if abs((D - midpoint) - (index + 1 - D - N)) > 1:
            continue
        if margin_data(index, D, N)[0] < 2:
            continue
        if not any(
            index - prime + 9 <= prime - 5
            for prime in prime_range(D + 1, D + N + 1)
        ):
            continue
        blocks.append((index, D, N))
    return tuple(blocks)


parser = argparse.ArgumentParser(add_help=False)
parser.add_argument("--central", type=int, default=0)
parser.add_argument("--central-all", type=int, default=0)
parser.add_argument("--random", type=int, default=0)
parser.add_argument("--random-blocks", type=int, default=0)
parser.add_argument("--balanced-blocks", type=int, default=0)
parser.add_argument("--seed", type=int, default=320)
parser.add_argument("--max-n", type=int, default=800)
parser.add_argument("--case", nargs=4, type=int, action="append")
arguments, _ = parser.parse_known_args()

records = []
work = []
if arguments.central:
    for index, D, N in central_cases(arguments.central):
        records.extend(audit_block(index, D, N))
if arguments.central_all:
    for index, D, N in central_cases(arguments.central_all, step=1):
        records.extend(audit_block(index, D, N))
if arguments.random:
    work.extend(
        random_cases(arguments.random, arguments.seed, arguments.max_n)
    )
if arguments.random_blocks:
    for block in random_blocks(
        arguments.random_blocks,
        arguments.seed,
        arguments.max_n,
    ):
        records.extend(audit_block(*block))
if arguments.balanced_blocks:
    for block in balanced_blocks(
        arguments.balanced_blocks,
        arguments.seed,
        arguments.max_n,
    ):
        records.extend(audit_block(*block))
if arguments.case:
    work.extend(tuple(case) for case in arguments.case)
if not records and not work:
    work = [
        (200, 128, 63, 139),
        (200, 128, 63, 181),
        (321, 168, 53, 179),
        (321, 168, 53, 193),
        (321, 168, 53, 211),
    ]

records.extend(audit_candidate(*case) for case in work)
failures = []
rank_failures = []
for record in records:
    if any(determinant == 0 for _, determinant, _ in record["records"]):
        failures.append(record)
        print("FIXED_MINOR_FAILURE", record)
    else:
        print("FIXED_MINOR_PASS", record)
    rank_failures.extend(
        (record["n"], record["core"], record["p"], margin)
        for margin, rank_two, _ in record["all_ghost_rank_two"]
        if not rank_two
    )

print(
    "Q32_FIXED_MINOR_JET_AUDIT",
    {
        "candidate_incidences": len(records),
        "margin_incidences": 2 * len(records),
        "failures": len(failures),
        "six_column_rank_failures": len(rank_failures),
    },
)
assert not failures or "--allow-failures" in sys.argv
assert not rank_failures
