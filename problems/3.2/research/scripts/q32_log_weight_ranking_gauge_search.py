#!/usr/bin/env python3
"""Search the tangent-gauge affine space for one-sided transition support.

The zero-block equations leave many free tangent parameters.  This
script builds the complete linear transition map for a finite Laurent
gauge basis and first asks the necessary self-edge question:

    can all self shifts with dz <= 0 (or dz >= 0) be cancelled?

If so, a z-directed ranking remains possible; if neither side can be
cancelled, the script returns exact left-kernel obstructions for this
gauge basis.
"""

from pathlib import Path
from math import gcd
from itertools import combinations

import numpy as np
import sympy as sp
from sympy.polys.matrices import DomainMatrix


HELPER = Path(__file__).with_name("q32_log_weight_zero_block_gauge.py")
prefix = HELPER.read_text().split("original =")[0]
namespace = {"__file__": str(HELPER)}
exec(prefix, namespace)

x, y, z = namespace["x"], namespace["y"], namespace["z"]
a, b, c = namespace["a"], namespace["b"], namespace["c"]
fields = namespace["fields"]
pair_fields = namespace["pair_fields"]
theta = namespace["theta"]
canonical_division = namespace["canonical_division"]
laurent_ledger = namespace["laurent_ledger"]


def scale_field(scalar, field):
    return tuple(sp.expand(scalar * entry) for entry in field)


def coefficient_vector(coefficient):
    polynomial = sp.Poly(sp.expand(coefficient), a, b, c)
    return (
        polynomial.coeff_monomial(a),
        polynomial.coeff_monomial(b),
        polynomial.coeff_monomial(c),
        polynomial.coeff_monomial(1),
    )


def transition_map(field):
    beta = sum(
        theta(field[index], (x, y, z)[index])
        for index in range(3)
    )
    coefficient = sp.expand(
        a * field[0] + b * field[1] + c * field[2] + beta
    )
    quotients, _ = canonical_division(coefficient)
    out = {}
    for target, quotient in enumerate(quotients):
        for shift, value in laurent_ledger(quotient).items():
            vector = coefficient_vector(value)
            for coordinate, entry in enumerate(vector):
                if entry != 0:
                    out[(target, shift, coordinate)] = entry
    return out


def polynomial_monomials(degree):
    return [
        x**i * y**j * z**k
        for i in range(degree + 1)
        for j in range(degree + 1 - i)
        for k in range(degree + 1 - i - j)
    ]


monomials = polynomial_monomials(3) + [
    z**power for power in range(-10, 0)
]
decorated = [
    (name, monomial, scale_field(monomial, field))
    for name, field in zip(("yx", "yz", "xz"), pair_fields)
    for monomial in monomials
]
print("GAUGE_BASIS", len(decorated))

basis_maps = []
for index, (_, _, field) in enumerate(decorated):
    basis_maps.append(transition_map(field))
    if (index + 1) % 30 == 0:
        print(" BUILT", index + 1)
original_maps = [transition_map(field) for field in fields]

universe = sorted(
    set().union(
        *(set(mapping) for mapping in basis_maps),
        *(set(mapping) for mapping in original_maps),
    )
)
row_index = {key: index for index, key in enumerate(universe)}
matrix_all = sp.zeros(len(universe), len(decorated))
for column, mapping in enumerate(basis_maps):
    for key, value in mapping.items():
        matrix_all[row_index[key], column] = value
original_vectors = []
for mapping in original_maps:
    vector = sp.zeros(len(universe), 1)
    for key, value in mapping.items():
        vector[row_index[key], 0] = value
    original_vectors.append(vector)
effect_rank = DomainMatrix.from_Matrix(matrix_all).rank()
print("UNIVERSE_ROWS", len(universe), "EFFECT_RANK", effect_rank)


def exact_solve(source, selected_keys):
    indices = [row_index[key] for key in selected_keys]
    matrix = matrix_all.extract(indices, range(matrix_all.cols))
    rhs = -original_vectors[source].extract(indices, [0])
    augmented = matrix.row_join(rhs)
    reduced_domain, pivots = DomainMatrix.from_Matrix(augmented).rref()
    rank = sum(pivot < matrix.cols for pivot in pivots)
    if matrix.cols in pivots:
        return None, rank, None
    reduced = reduced_domain.to_Matrix()
    chosen = [sp.Rational(0)] * matrix.cols
    for row, pivot in enumerate(pivots):
        chosen[pivot] = reduced[row, matrix.cols]
    chosen = tuple(chosen)
    assert matrix * sp.Matrix(chosen) == rhs
    return chosen, rank, None


for source, source_name in enumerate(("y", "x", "z")):
    for direction_name, bad_predicate in (
        ("KEEP_DZ_POSITIVE", lambda dz: dz <= 0),
        ("KEEP_DZ_NEGATIVE", lambda dz: dz >= 0),
    ):
        selected = [
            key
            for key in universe
            if key[0] == source and bad_predicate(key[1][2])
        ]
        solution, rank, obstruction = exact_solve(source, selected)
        print(
            "SELF_TEST", source_name, direction_name,
            "ROWS", len(selected), "RANK", rank,
            "EXISTS", solution is not None,
        )
        if solution is not None:
            nonzero = [
                (decorated[index][0], decorated[index][1], coefficient)
                for index, coefficient in enumerate(solution)
                if coefficient != 0
            ]
            print(" NONZERO", nonzero)
        elif obstruction is not None:
            vector, pairing = obstruction
            support = [
                (selected[index], value)
                for index, value in enumerate(vector)
                if value != 0
            ]
            print(" OBSTRUCTION_PAIRING", pairing)
            print(" OBSTRUCTION_SUPPORT", support[:12])


# Compress the unavoidable y self-chain to one negative and one positive
# pure-z shift while retaining the complete zero-block cancellation.
zero_block_keys = [
    key for key in universe if key[1] == (0, 0, 0)
]
pure_z_y_keys = [
    key
    for key in universe
    if key[0] == 0
    and key[1][0] == 0
    and key[1][1] == 0
    and key[1][2] != 0
]
negative_shifts = sorted({key[1][2] for key in pure_z_y_keys if key[1][2] < 0})
positive_shifts = sorted({key[1][2] for key in pure_z_y_keys if key[1][2] > 0})
print("Y_PURE_Z_RANGE", negative_shifts, positive_shifts)
pair_solution = None
for negative in negative_shifts:
    for positive in positive_shifts:
        selected = zero_block_keys + [
            key
            for key in pure_z_y_keys
            if key[1][2] not in (negative, positive)
        ]
        solution, rank, _ = exact_solve(0, selected)
        if solution is None:
            continue
        pair_solution = (negative, positive, solution, rank)
        break
    if pair_solution is not None:
        break

print("Y_TWO_SIDED_PAIR_EXISTS", pair_solution is not None)
if pair_solution is not None:
    negative, positive, solution, rank = pair_solution
    print(" Y_PAIR", negative, positive, "RANK", rank)
    nonzero = [
        (decorated[index][0], decorated[index][1], coefficient)
        for index, coefficient in enumerate(solution)
        if coefficient != 0
    ]
    print(" Y_PAIR_GAUGE", nonzero)
    resulting = original_vectors[0] + matrix_all * sp.Matrix(solution)
    surviving_pure = [
        (key, sp.factor(resulting[row_index[key], 0]))
        for key in pure_z_y_keys
        if resulting[row_index[key], 0] != 0
    ]
    print(" Y_PAIR_SURVIVING", surviving_pure)


def float_consistent(source, selected_keys):
    indices = [row_index[key] for key in selected_keys]
    left = np.array(
        matrix_all.extract(indices, range(matrix_all.cols)).tolist(),
        dtype=float,
    )
    right = np.array(
        (-original_vectors[source].extract(indices, [0])).tolist(),
        dtype=float,
    ).reshape(-1)
    solution, _, _, _ = np.linalg.lstsq(left, right, rcond=1e-11)
    residual = np.linalg.norm(left @ solution - right, ord=np.inf)
    scale = max(1.0, np.linalg.norm(right, ord=np.inf))
    return residual <= 1e-8 * scale


if pair_solution is None:
    all_pure_shifts = negative_shifts + positive_shifts
    compressed_solution = None
    for allowed_count in range(3, 7):
        tested = 0
        float_survivors = 0
        for allowed in combinations(all_pure_shifts, allowed_count):
            if not any(shift < 0 for shift in allowed):
                continue
            if not any(shift > 0 for shift in allowed):
                continue
            selected = zero_block_keys + [
                key
                for key in pure_z_y_keys
                if key[1][2] not in allowed
            ]
            tested += 1
            if not float_consistent(0, selected):
                continue
            float_survivors += 1
            solution, rank, _ = exact_solve(0, selected)
            if solution is not None:
                compressed_solution = (allowed, solution, rank)
                break
        print(
            "Y_ALLOWED_SEARCH", allowed_count,
            "TESTED", tested,
            "FLOAT_SURVIVORS", float_survivors,
            "EXACT_EXISTS", compressed_solution is not None,
        )
        if compressed_solution is not None:
            break

    if compressed_solution is not None:
        allowed, solution, rank = compressed_solution
        print(" Y_MIN_ALLOWED", allowed, "RANK", rank)
        nonzero = [
            (decorated[index][0], decorated[index][1], coefficient)
            for index, coefficient in enumerate(solution)
            if coefficient != 0
        ]
        print(" Y_COMPRESSED_GAUGE", nonzero)
        resulting = original_vectors[0] + matrix_all * sp.Matrix(solution)
        surviving_pure = [
            (key, sp.factor(resulting[row_index[key], 0]))
            for key in pure_z_y_keys
            if resulting[row_index[key], 0] != 0
        ]
        print(" Y_COMPRESSED_SURVIVING", surviving_pure)


def compress_source_pure_z(source, full_zero_row):
    source_name = ("y", "x", "z")[source]
    if full_zero_row:
        base_keys = zero_block_keys
    else:
        base_keys = [
            key
            for key in zero_block_keys
            if key[0] == source
        ]
    pure_keys = [
        key
        for key in universe
        if key[0] == source
        and key[1][0] == 0
        and key[1][1] == 0
        and key[1][2] != 0
    ]
    shifts = sorted({key[1][2] for key in pure_keys})
    print(source_name.upper() + "_PURE_Z_RANGE", shifts)
    answer = None
    for allowed_count in range(2, 7):
        tested = 0
        float_survivors = 0
        for allowed in combinations(shifts, allowed_count):
            if not any(shift < 0 for shift in allowed):
                continue
            if not any(shift > 0 for shift in allowed):
                continue
            selected = base_keys + [
                key for key in pure_keys if key[1][2] not in allowed
            ]
            tested += 1
            if not float_consistent(source, selected):
                continue
            float_survivors += 1
            solution, rank, _ = exact_solve(source, selected)
            if solution is not None:
                answer = (allowed, solution, rank, pure_keys)
                break
        print(
            source_name.upper() + "_ALLOWED_SEARCH",
            allowed_count,
            "TESTED", tested,
            "FLOAT_SURVIVORS", float_survivors,
            "EXACT_EXISTS", answer is not None,
        )
        if answer is not None:
            break
    if answer is None:
        return
    allowed, solution, rank, pure_keys = answer
    print(source_name.upper() + "_MIN_ALLOWED", allowed, "RANK", rank)
    nonzero = [
        (decorated[index][0], decorated[index][1], coefficient)
        for index, coefficient in enumerate(solution)
        if coefficient != 0
    ]
    print(source_name.upper() + "_COMPRESSED_GAUGE", nonzero)
    resulting = original_vectors[source] + matrix_all * sp.Matrix(solution)
    surviving = [
        (key, sp.factor(resulting[row_index[key], 0]))
        for key in pure_keys
        if resulting[row_index[key], 0] != 0
    ]
    print(source_name.upper() + "_COMPRESSED_SURVIVING", surviving)


compress_source_pure_z(1, full_zero_row=True)
compress_source_pure_z(2, full_zero_row=False)

print("Q32_LOG_WEIGHT_RANKING_GAUGE_SEARCH=PASS")
