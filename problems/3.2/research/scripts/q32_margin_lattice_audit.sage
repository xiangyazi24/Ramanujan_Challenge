#!/usr/bin/env sage
"""Compare fixed-length and diagonal-margin Newton exterior lattices.

This is a universal coefficient audit.  It does not use the Apéry shell
values.  Every determinant is expanded in the Plücker coordinates

    X_{i,j} = Y_i Z_j - Y_j Z_i.

The script determines whether the two target-preserving carrier families
from Sections 68.6 and 68.11 are rationally or integrally equivalent.
"""

from itertools import combinations
from math import gcd as integer_gcd


def gcd_many(values):
    out = 0
    for value in values:
        out = integer_gcd(out, int(value))
    return ZZ(abs(out))


def carrier_row(lower, length, nodes):
    """Coefficient row of G_{lower,length} in the global node basis."""

    positions = {node: index for index, node in enumerate(nodes)}
    row = vector(ZZ, [0] * len(nodes))
    for offset in range(length + 1):
        node = lower + offset
        row[positions[node]] = (
            (-1) ** offset
            * binomial(lower + offset, offset)
            * binomial(lower + length + 1, length - offset)
        )
    return row


def wedge_row(left, right):
    """Plücker coefficient row of det(left*X,right*X)."""

    return vector(
        ZZ,
        [
            left[i] * right[j] - left[j] * right[i]
            for i, j in combinations(range(len(left)), 2)
        ],
    )


def primitive_divide(row, divisor):
    assert divisor > 0
    assert all(entry % divisor == 0 for entry in row)
    return vector(ZZ, [entry // divisor for entry in row])


def family_matrices(D, N, maximum_margin):
    lower = D - maximum_margin + 1
    upper = D + N + maximum_margin - 2
    nodes = tuple(range(lower, upper + 1))

    # Section 68.6: one maximal length and every adjacent translation.
    maximal_length = N + maximum_margin - 2
    fixed_rows = []
    for shift in range(maximum_margin - 1):
        d = lower + shift
        left = carrier_row(d, maximal_length, nodes)
        right = carrier_row(d + 1, maximal_length, nodes)
        prefactor = binomial(d + maximal_length + 1, maximal_length)
        fixed_rows.append(
            primitive_divide(wedge_row(left, right), prefactor)
        )

    # Section 68.11: every diagonal margin from 2 to maximum_margin.
    diagonal_rows = []
    for margin in range(2, maximum_margin + 1):
        d = D - margin + 1
        length = N + margin - 2
        left = carrier_row(d, length, nodes)
        right = carrier_row(D, length, nodes)
        common_weight = gcd_many(
            binomial(d + shift + length + 1, length)
            for shift in range(margin - 1)
        )
        diagonal_rows.append(
            primitive_divide(wedge_row(left, right), common_weight)
        )

    return matrix(ZZ, fixed_rows), matrix(ZZ, diagonal_rows)


def change_of_rows(source, target):
    """Return T over QQ with target=T*source, or None."""

    source_q = matrix(QQ, source)
    target_q = matrix(QQ, target)
    if source_q.rank() != source_q.nrows():
        return None
    pivots = source_q.pivots()
    square = source_q.matrix_from_columns(pivots)
    assert square.det() != 0
    transform = target_q.matrix_from_columns(pivots) * square.inverse()
    if transform * source_q != target_q:
        return None
    return transform


def denominator_lcm(matrix_value):
    return lcm(entry.denominator() for entry in matrix_value.list())


def audit_case(D, N, maximum_margin):
    fixed, diagonal = family_matrices(D, N, maximum_margin)
    fixed_rank = fixed.rank()
    diagonal_rank = diagonal.rank()
    combined_rank = block_matrix([[fixed], [diagonal]]).rank()
    intersection_rank = fixed_rank + diagonal_rank - combined_rank
    fixed_to_diagonal = change_of_rows(fixed, diagonal)
    diagonal_to_fixed = change_of_rows(diagonal, fixed)

    record = {
        "parameters": (D, N, maximum_margin),
        "shape": fixed.dimensions(),
        "ranks": (fixed_rank, diagonal_rank),
        "combined_rank": combined_rank,
        "intersection_rank": intersection_rank,
        "same_rational_space": (
            fixed_to_diagonal is not None
            and diagonal_to_fixed is not None
        ),
    }
    nodes = tuple(
        range(
            D - maximum_margin + 1,
            D + N + maximum_margin - 1,
        )
    )
    pairs = tuple(combinations(range(len(nodes)), 2))
    target_star_records = []
    for prime in prime_range(D + 1, D + N + 1):
        maximal_length = N + maximum_margin - 2
        if maximal_length >= prime:
            continue
        marked_index = nodes.index(prime - 1)
        star_columns = tuple(
            column
            for column, pair in enumerate(pairs)
            if marked_index in pair
        )
        nonstar_columns = tuple(
            column
            for column, pair in enumerate(pairs)
            if marked_index not in pair
        )
        combined_mod = block_matrix([[fixed], [diagonal]]).change_ring(
            GF(prime)
        )
        assert combined_mod.matrix_from_columns(nonstar_columns).is_zero()
        target_star_records.append(
            (
                prime,
                len(star_columns),
                combined_mod.matrix_from_columns(star_columns).rank(),
            )
        )
    record["target_star_records"] = tuple(target_star_records)
    if fixed_to_diagonal is not None:
        record["fixed_to_diagonal_denominator"] = denominator_lcm(
            fixed_to_diagonal
        )
        record["fixed_to_diagonal_det"] = factor(
            fixed_to_diagonal.det()
        )
    if diagonal_to_fixed is not None:
        record["diagonal_to_fixed_denominator"] = denominator_lcm(
            diagonal_to_fixed
        )
        record["diagonal_to_fixed_det"] = factor(
            diagonal_to_fixed.det()
        )
    return record


CASES = (
    (8, 3, 3),
    (10, 4, 4),
    (13, 5, 5),
    (17, 6, 6),
    (19, 7, 7),
)

records = tuple(audit_case(*case) for case in CASES)
print("MARGIN_LATTICE_RECORDS", records)
print("Q32_MARGIN_LATTICE_AUDIT=PASS")
