#!/usr/bin/env python3
"""Audit the fixed-length Newton-carrier curvature obstruction.

For a fixed carrier length ``L`` set

    A_d(F) = G_{d,L}(F),
    J_d(F) = sum_q (-1)^q binom(L+1,q) F[d+q],
    P_d    = binom(d+L+1,L).

The beta-Pade difference law is

    A_d(F) - A_{d+1}(F) = P_d J_d(F).

For two integer arrays Y,Z let a_d=(A_d(Y),A_d(Z)) and
j_d=(J_d(Y),J_d(Z)).  This script independently checks the exact
endpoint-curvature identity

    det(a_0,a_r)
      = sum_s P_s E_s
        + sum_{t<s} P_t P_s det(j_s,j_t),

where E_s=det(a_s,a_{s+1})/P_s.  It also realizes an explicit
integer-array counterexample showing that divisibility of every
adjacent normalized minor need not imply endpoint divisibility.

Thus an endpoint collapse cannot follow from the universal
Newton/Pascal kernel alone; any useful collapse must use additional
Apéry-shell structure.
"""

from __future__ import annotations

from fractions import Fraction
from math import comb, gcd
from random import Random


def determinant(left: tuple[int, int], right: tuple[int, int]) -> int:
    return left[0] * right[1] - left[1] * right[0]


def carrier(values: list[int], d: int, length: int) -> int:
    return sum(
        (-1) ** i
        * comb(d + i, i)
        * comb(d + length + 1, length - i)
        * values[d + i]
        for i in range(length + 1)
    )


def high_difference(values: list[int], d: int, length: int) -> int:
    return sum(
        (-1) ** q * comb(length + 1, q) * values[d + q]
        for q in range(length + 2)
    )


def realize_packets(
    d0: int, length: int, initial_carrier: int, packets: list[int]
) -> list[int]:
    """Realize arbitrary integral ``A_0,J_0,...,J_{r-1}`` data."""

    upper = d0 + len(packets) + length
    values = [0] * (upper + 1)
    for shift, target in enumerate(packets):
        d = d0 + shift
        known = sum(
            (-1) ** q * comb(length + 1, q) * values[d + q]
            for q in range(length + 1)
        )
        values[d + length + 1] = (
            target - known
        ) * (-1) ** (length + 1)
        assert high_difference(values, d, length) == target

    # The Newton row has coefficient sum one, while every high
    # difference kills constants.
    correction = initial_carrier - carrier(values, d0, length)
    for index in range(d0, upper + 1):
        values[index] += correction
    assert carrier(values, d0, length) == initial_carrier
    assert [
        high_difference(values, d0 + shift, length)
        for shift in range(len(packets))
    ] == packets
    return values


def audit_random_arrays() -> int:
    random = Random(20260730)
    checks = 0
    for length in range(0, 8):
        for d0 in range(1, 8):
            for r in range(2, 9):
                upper = d0 + r + length
                for _ in range(4):
                    y = [random.randrange(-10**5, 10**5) for _ in range(upper + 1)]
                    z = [random.randrange(-10**5, 10**5) for _ in range(upper + 1)]
                    a = [
                        (
                            carrier(y, d0 + s, length),
                            carrier(z, d0 + s, length),
                        )
                        for s in range(r + 1)
                    ]
                    j = [
                        (
                            high_difference(y, d0 + s, length),
                            high_difference(z, d0 + s, length),
                        )
                        for s in range(r)
                    ]
                    weights = [
                        comb(d0 + s + length + 1, length)
                        for s in range(r)
                    ]
                    edges: list[int] = []
                    for s in range(r):
                        difference = (
                            a[s][0] - a[s + 1][0],
                            a[s][1] - a[s + 1][1],
                        )
                        assert difference == (
                            weights[s] * j[s][0],
                            weights[s] * j[s][1],
                        )
                        raw = determinant(a[s], a[s + 1])
                        assert raw % weights[s] == 0
                        edge = raw // weights[s]
                        assert edge == determinant(j[s], a[s])
                        assert edge == determinant(j[s], a[s + 1])
                        edges.append(edge)

                    curvature = sum(
                        weights[t]
                        * weights[s]
                        * determinant(j[s], j[t])
                        for s in range(r)
                        for t in range(s)
                    )
                    assert determinant(a[0], a[r]) == sum(
                        weights[s] * edges[s] for s in range(r)
                    ) + curvature
                    checks += 1
    return checks


def audit_realization_and_counterexamples() -> int:
    checks = 0
    for d0, length, r, prime in [
        (1, 0, 2, 101),
        (2, 5, 3, 1009),
        (7, 4, 8, 10007),
    ]:
        weights = [
            comb(d0 + s + length + 1, length) for s in range(r)
        ]
        assert weights[0] % prime and weights[1] % prime
        left_inverse = pow(weights[0], -1, prime)
        right_inverse = pow(weights[1], -1, prime)
        packets = (
            [(left_inverse, 0), (0, -right_inverse)]
            + [(0, 0)] * (r - 2)
        )
        a0 = (
            prime + weights[0] * left_inverse,
            prime,
        )

        y = realize_packets(
            d0, length, a0[0], [packet[0] for packet in packets]
        )
        z = realize_packets(
            d0, length, a0[1], [packet[1] for packet in packets]
        )
        path = [
            (
                carrier(y, d0 + s, length),
                carrier(z, d0 + s, length),
            )
            for s in range(r + 1)
        ]
        edges = [
            determinant(path[s], path[s + 1]) // weights[s]
            for s in range(r)
        ]
        assert edges[:2] == [
            prime * left_inverse,
            prime * right_inverse,
        ]
        assert edges[2:] == [0] * (r - 2)
        common = 0
        for edge in edges:
            common = gcd(common, abs(edge))
        assert common % prime == 0
        assert determinant(path[0], path[-1]) % prime == 1
        checks += 1
    return checks


def rational_rank(matrix: list[list[int]]) -> int:
    work = [[Fraction(entry) for entry in row] for row in matrix]
    row = 0
    for column in range(len(work[0]) if work else 0):
        pivot = next(
            (candidate for candidate in range(row, len(work))
             if work[candidate][column]),
            None,
        )
        if pivot is None:
            continue
        work[row], work[pivot] = work[pivot], work[row]
        scale = work[row][column]
        work[row] = [entry / scale for entry in work[row]]
        for candidate in range(len(work)):
            if candidate == row:
                continue
            factor = work[candidate][column]
            if factor:
                work[candidate] = [
                    left - factor * right
                    for left, right in zip(work[candidate], work[row])
                ]
        row += 1
        if row == len(work):
            break
    return row


def audit_curvature_rank() -> int:
    for size in range(1, 25):
        skew = [
            [
                0 if row == column else (1 if row > column else -1)
                for column in range(size)
            ]
            for row in range(size)
        ]
        expected = size if size % 2 == 0 else size - 1
        assert rational_rank(skew) == expected
    return 24


def main() -> None:
    random_checks = audit_random_arrays()
    counterexamples = audit_realization_and_counterexamples()
    rank_checks = audit_curvature_rank()
    print(
        "CROSS_CURVATURE_AUDIT_OK",
        {
            "random_array_checks": random_checks,
            "integral_counterexamples": counterexamples,
            "skew_rank_checks": rank_checks,
        },
    )


if __name__ == "__main__":
    main()
