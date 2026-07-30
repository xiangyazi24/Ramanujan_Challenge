#!/usr/bin/env python3
"""Exact audits for the first-cell shell decomposition in Section 57.

The default run checks the decomposition, endpoint identity, and
finite-difference splitting on small, selected large, and deterministic
random inputs.  Pass ``--rows`` to reproduce the five expensive
primitive-residual computations in (57.13).
"""

from argparse import ArgumentParser
from math import comb, gcd, prod
from random import Random

from q32_cartier_packet_audit import apery, shell_fast
from q32_newton_gcd_audit import carrier_from_values, forward_difference


def C(n, k):
    return comb(n, k) if 0 <= k <= n else 0


def first_cell_parts(moment, node):
    """Return b_M, the long core, two boundary sums, and their total."""

    assert moment // 2 < node <= moment
    residue = moment - node
    core = 0
    low = 0
    for t in range(moment + 1):
        outer = C(moment, t)
        upper = 2 * moment - t
        central = C(upper, moment)
        global_term = C(upper, residue)
        core += outer**2 * (
            2 * central * global_term + global_term**2
        )
        if t <= residue:
            low_term = C(upper, residue - t)
            x_term = C(moment, residue - t)
            low += outer * (
                outer
                * (
                    2 * (central + global_term) * low_term
                    + low_term**2
                )
                + x_term
                * (central + global_term + low_term) ** 2
            )

    high = sum(
        C(moment, k)
        * C(moment, residue - k)
        * (C(moment + k, k) + C(moment + k, residue)) ** 2
        for k in range(residue + 1)
    )
    origin = apery(moment)
    return origin, core, low, high, origin + core + low + high


def long_core_coefficient(moment, residue):
    """Return [z^r] Kcal_M(z) from the two polynomial families."""

    first = 2 * sum(
        C(moment, k) ** 2
        * C(moment + k, k)
        * C(moment + k, residue)
        for k in range(moment + 1)
    )
    second = sum(
        C(moment, k) ** 2 * C(moment + k, residue) ** 2
        for k in range(moment + 1)
    )
    return first + second


def zeta2_apery(moment):
    return sum(
        C(moment, k) ** 2 * C(moment + k, k)
        for k in range(moment + 1)
    )


def endpoint_value(moment):
    central = C(2 * moment, moment)
    return (
        apery(moment)
        + 2 * zeta2_apery(moment)
        + central**2
        + 7 * central
        + 11
    )


def audit_formula(moment, node):
    parts = first_cell_parts(moment, node)
    assert shell_fast(moment, node) == parts[-1]
    assert parts[1] == long_core_coefficient(moment, moment - node)


def audit_difference(moment, node, order):
    rows = [
        first_cell_parts(moment, node + shift)
        for shift in range(order + 1)
    ]
    actual = forward_difference(
        [shell_fast(moment, node + shift) for shift in range(order + 1)],
        order,
    )[0]
    split = sum(
        forward_difference([row[index] for row in rows], order)[0]
        for index in (1, 2, 3)
    )
    assert actual == split


def audit_default():
    for moment in range(1, 35):
        assert shell_fast(moment, moment) == endpoint_value(moment)
        for node in range(moment // 2 + 1, moment + 1):
            audit_formula(moment, node)
            for order in range(1, min(5, moment - node) + 1):
                audit_difference(moment, node, order)

    for moment in (199, 271, 299, 320, 754):
        nodes = sorted(
            {
                moment // 2 + 1,
                moment // 2 + 2,
                3 * moment // 4,
                moment - 7,
                moment - 1,
                moment,
            }
        )
        for node in nodes:
            if moment // 2 < node <= moment:
                audit_formula(moment, node)
        assert shell_fast(moment, moment) == endpoint_value(moment)
        print("FORMULA", moment, nodes)

    random = Random(5711)
    for _ in range(40):
        moment = random.randrange(20, 180)
        node = random.randrange(moment // 2 + 1, moment + 1)
        audit_formula(moment, node)
        if node < moment:
            order = random.randrange(1, min(7, moment - node + 1))
            audit_difference(moment, node, order)


def audit_rows():
    rows = (
        (200, (139, 181), 5),
        (272, (191, 233), 385),
        (300, (191, 227), 1),
        (321, (179, 193, 211), 1),
        (755, (593, 733), 85),
    )
    for n, targets, expected_nuisance in rows:
        moment = n - 1
        node = moment // 2 + 1
        length = max(targets) - 1 - node
        values = {
            index: shell_fast(moment, index)
            for index in range(node, node + length + 2)
        }
        left = carrier_from_values(values, node, length)
        right = carrier_from_values(values, node + 1, length)
        high_difference = forward_difference(
            [
                values[index]
                for index in range(node, node + length + 2)
            ],
            length + 1,
        )[0]
        target_product = prod(targets)
        assert gcd(abs(left), abs(high_difference)) == 1
        assert (
            gcd(abs(left), abs(right))
            == expected_nuisance * target_product
        )
        print(
            "ROW",
            n,
            "length",
            length,
            "nuisance",
            expected_nuisance,
        )


def main():
    parser = ArgumentParser()
    parser.add_argument("--rows", action="store_true")
    args = parser.parse_args()

    audit_default()
    if args.rows:
        audit_rows()
    print("PASS: first-cell decomposition and primitive residuals")


if __name__ == "__main__":
    main()
