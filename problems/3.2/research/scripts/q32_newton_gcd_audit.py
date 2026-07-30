#!/usr/bin/env python3
"""Exact audits for the adjacent-Newton gcd route in Sections 53--55.

The default run verifies the universal interpolation identities and the
hostile n=321 collapse, including the local Smith form of a stencil
rectangle.  Pass ``--extended`` to reproduce the five-row minimal-stencil
residual table, and ``--blocks`` to reproduce the optimal two-carrier
block table; those runs are substantially slower.
"""

from argparse import ArgumentParser
from itertools import combinations
from math import comb, gcd, prod

from q32_cartier_packet_audit import newton_weight, shell_fast


def carrier_from_values(values, d, length):
    return sum(
        newton_weight(d, length, i) * values[d + i]
        for i in range(length + 1)
    )


def coefficient_row(d, length, lower, upper):
    """Return one Newton-stencil coefficient row on a fixed node interval."""

    row = [0] * (upper - lower + 1)
    for i in range(length + 1):
        row[d + i - lower] = newton_weight(d, length, i)
    return row


def rank_mod(rows, prime):
    """Return the row rank of an integer matrix modulo ``prime``."""

    matrix = [[entry % prime for entry in row] for row in rows]
    rank = 0
    columns = len(matrix[0]) if matrix else 0
    for column in range(columns):
        pivot = next(
            (
                index
                for index in range(rank, len(matrix))
                if matrix[index][column]
            ),
            None,
        )
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        inverse = pow(matrix[rank][column], -1, prime)
        matrix[rank] = [
            entry * inverse % prime for entry in matrix[rank]
        ]
        for index in range(len(matrix)):
            if index == rank or not matrix[index][column]:
                continue
            multiple = matrix[index][column]
            matrix[index] = [
                (
                    matrix[index][other]
                    - multiple * matrix[rank][other]
                )
                % prime
                for other in range(columns)
            ]
        rank += 1
    return rank


def determinant(matrix):
    """Return the determinant of a small square integer matrix."""

    if len(matrix) == 1:
        return matrix[0][0]
    return sum(
        (-1) ** column
        * matrix[0][column]
        * determinant(
            [
                row[:column] + row[column + 1 :]
                for row in matrix[1:]
            ]
        )
        for column in range(len(matrix))
    )


def determinantal_divisor(matrix, size):
    """Return the gcd of all ``size``-minors of ``matrix``."""

    divisor = 0
    for row_indices in combinations(range(len(matrix)), size):
        for column_indices in combinations(range(len(matrix[0])), size):
            minor = [
                [matrix[row][column] for column in column_indices]
                for row in row_indices
            ]
            divisor = gcd(divisor, abs(determinant(minor)))
    return divisor


def forward_difference(values, order):
    current = list(values)
    for _ in range(order):
        current = [
            current[i + 1] - current[i]
            for i in range(len(current) - 1)
        ]
    assert len(current) == len(values) - order
    return current


def audit_rectangle_smith(D, R, A, B, prime):
    """Audit the target-local Smith profile of a stencil rectangle.

    The independent boundary consists of the bottom edge
    ``G_{D+t,A}``, followed by the right edge ``G_{D+R,L}``.
    If ``prime - 1`` is in every stencil, all boundary rows coincide
    modulo ``prime``.  Dividing each successive edge difference by its
    single forced prime factor produces a full-rank matrix modulo the
    prime.  Hence the local Smith form is
    ``diag(1, prime, ..., prime)``.
    """

    assert R <= A <= B < prime
    node = prime - 1
    assert D + R <= node <= D + A
    lower = D
    upper = D + R + B

    boundary = [
        coefficient_row(D + shift, A, lower, upper)
        for shift in range(R + 1)
    ]
    boundary.extend(
        coefficient_row(D + R, length, lower, upper)
        for length in range(A + 1, B + 1)
    )
    boundary_rank = R + B - A + 1
    assert len(boundary) == boundary_rank

    coordinate = [0] * (upper - lower + 1)
    coordinate[node - lower] = 1
    assert all(
        [entry % prime for entry in row] == coordinate
        for row in boundary
    )
    assert rank_mod(boundary, prime) == 1

    divided = [boundary[0]]
    for shift in range(R):
        difference = [
            boundary[shift][column] - boundary[shift + 1][column]
            for column in range(len(coordinate))
        ]
        assert all(entry % prime == 0 for entry in difference)
        divided.append([entry // prime for entry in difference])

    right_edge_start = R
    for offset in range(B - A):
        lower_row = boundary[right_edge_start + offset]
        upper_row = boundary[right_edge_start + offset + 1]
        difference = [
            upper_row[column] - lower_row[column]
            for column in range(len(coordinate))
        ]
        assert all(entry % prime == 0 for entry in difference)
        divided.append([entry // prime for entry in difference])

    assert rank_mod(divided, prime) == boundary_rank


def audit_rectangle_identities():
    sequence = {
        d: (d + 2) ** 7 + 3 ** d + 11 * d
        for d in range(1, 30)
    }
    for d in range(2, 10):
        for length in range(0, 7):
            current = carrier_from_values(sequence, d, length)
            horizontal = (
                current
                - carrier_from_values(sequence, d + 1, length)
            )
            delta = forward_difference(
                [
                    sequence[node]
                    for node in range(d, d + length + 2)
                ],
                length + 1,
            )[0]
            assert horizontal == (
                (-1) ** (length + 1)
                * comb(d + length + 1, length)
                * delta
            )

            vertical = (
                carrier_from_values(sequence, d, length + 1)
                - current
            )
            assert vertical == (
                (-1) ** (length + 1)
                * comb(d + length + 1, length + 1)
                * delta
            )
            assert (d + 1) * horizontal == (length + 1) * vertical

    # The common-node interval is [166, 211].  The check includes both
    # target and non-target primes: the Smith statement is formal and
    # independent of the shell values.
    for prime in (167, 173, 179, 181, 191, 193, 197, 199, 211):
        audit_rectangle_smith(161, 5, 50, 50, prime)

    # In the coordinates X=G_{d,L-1}, U=Delta^L Y_d, and
    # V=Delta^L Y_{d+1}, the four-value square has Smith invariants
    # 1, gcd(A,C), B.
    for d in range(2, 13):
        for length in range(1, 8):
            sign = (-1) ** length
            A = comb(d + length, length)
            C = comb(d + length, length - 1)
            B = A + C
            content = gcd(A, C)
            matrix = [
                [1, 0, 0],
                [1, -sign * C, 0],
                [1, sign * A, 0],
                [1, -sign * C, sign * B],
            ]
            assert determinantal_divisor(matrix, 1) == 1
            assert determinantal_divisor(matrix, 2) == content
            assert determinantal_divisor(matrix, 3) == content * B


def audit_universal_identities():
    # Deliberately non-polynomial data: none of these checks relies on a
    # finite-degree interpolation accident.
    sequence = {
        d: (d + 2) ** 7 + 3 ** d + 11 * d
        for d in range(1, 30)
    }

    for d in range(2, 10):
        for length in range(1, 8):
            g_dl = carrier_from_values(sequence, d, length)
            g_d_lm1 = carrier_from_values(sequence, d, length - 1)
            g_dp1_lm1 = carrier_from_values(sequence, d + 1, length - 1)
            g_dp1_l = carrier_from_values(sequence, d + 1, length)

            delta_l = forward_difference(
                [sequence[j] for j in range(d, d + length + 1)],
                length,
            )[0]
            delta_lp1 = forward_difference(
                [sequence[j] for j in range(d, d + length + 2)],
                length + 1,
            )[0]

            assert g_dl - g_d_lm1 == (
                (-1) ** length * comb(d + length, length) * delta_l
            )
            assert g_dl - g_dp1_lm1 == (
                (-1) ** length
                * comb(d + length + 1, length)
                * delta_l
            )
            assert g_dl - g_dp1_l == (
                (-1) ** (length + 1)
                * comb(d + length + 1, length)
                * delta_lp1
            )
            assert length * g_dl == (
                (d + length + 1) * g_d_lm1
                - (d + 1) * g_dp1_lm1
            )


def minimal_stencil_row(n, targets):
    """Return the exact two-window residual data for one target row."""

    moment = n - 1
    left_target_node = min(targets) - 1
    right_target_node = max(targets) - 1
    length = right_target_node - left_target_node + 1
    d = left_target_node - 1

    values = {
        node: shell_fast(moment, node)
        for node in range(d, d + length + 2)
    }
    first = carrier_from_values(values, d, length)
    second = carrier_from_values(values, d + 1, length)
    delta = forward_difference(
        [values[node] for node in range(d, d + length + 2)],
        length + 1,
    )[0]
    binomial_content = comb(d + length + 1, length)

    assert first - second == (
        (-1) ** (length + 1) * binomial_content * delta
    )

    target_product = prod(targets)
    carrier_gcd = gcd(abs(first), abs(second))
    assert carrier_gcd % target_product == 0
    return (
        length,
        carrier_gcd // target_product,
        gcd(abs(first), abs(delta)),
    )


def audit_hostile_collapse():
    moment, length = 320, 50
    values = {
        node: shell_fast(moment, node) for node in range(161, 217)
    }
    carriers = [
        carrier_from_values(values, d, length)
        for d in range(161, 167)
    ]
    target_product = 179 * 193 * 211

    assert gcd(abs(carriers[0]), abs(carriers[1])) == target_product
    running = 0
    for value in carriers:
        running = gcd(running, abs(value))
    assert running == target_product
    assert {len(str(abs(value))) for value in carriers} <= {528, 529}
    return running


def audit_extended_table():
    rows = (
        (200, (139, 181), (43, 11, 11)),
        (272, (191, 233), (43, 5, 1)),
        (300, (191, 227), (37, 15, 15)),
        (321, (179, 193, 211), (33, 111, 111)),
        (755, (593, 733), (141, 20075, 5)),
    )
    for n, targets, expected in rows:
        actual = minimal_stencil_row(n, targets)
        assert actual == expected, (n, actual, expected)
        print("ROW", n, actual)


def audit_optimal_block_table():
    """Reproduce the exact two-carrier block gcds in Section 56."""

    rows = (
        (200, 128, 63, (139, 181), 47),
        (272, 180, 63, (191, 233), 1),
        (300, 180, 57, (191, 227), 1),
        (321, 168, 53, (179, 193, 211), 43),
        (755, 582, 161, (593, 733), 275),
    )
    for n, D, length, targets, expected_nuisance in rows:
        moment = n - 1
        values = {
            node: shell_fast(moment, node)
            for node in range(D - 1, D + length + 1)
        }
        left = carrier_from_values(values, D - 1, length)
        right = carrier_from_values(values, D, length)
        block_gcd = gcd(abs(left), abs(right))
        target_product = prod(targets)
        assert block_gcd == expected_nuisance * target_product
        print(
            "BLOCK",
            n,
            (D, length),
            "digits",
            (len(str(abs(left))), len(str(abs(right)))),
            "nuisance",
            expected_nuisance,
        )


def main():
    parser = ArgumentParser()
    parser.add_argument("--extended", action="store_true")
    parser.add_argument("--blocks", action="store_true")
    args = parser.parse_args()

    audit_universal_identities()
    audit_rectangle_identities()
    hostile_gcd = audit_hostile_collapse()
    if args.extended:
        audit_extended_table()
    if args.blocks:
        audit_optimal_block_table()
    print("PASS: universal Newton identities; hostile gcd", hostile_gcd)


if __name__ == "__main__":
    main()
