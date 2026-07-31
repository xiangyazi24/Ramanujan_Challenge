#!/usr/bin/env python3
"""Exact directional-cofactor audit for the terminal Y/W packet.

The theorem checked here is

    rad(g_raw(n)) | rad(D_dir(n)),

where ``g_raw`` is the gcd of the 28 raw 2-by-2 minors and ``D_dir``
is the gcd of 560 directional cofactors of the actual common-
convolution observability matrix.  The additional divisibility by
the short product (n+2)...(n+7) is deliberately labelled as a
conjectural regression assertion.

The origin operator is read from the exact JSON dump produced by
``/tmp/p32_dump_origin_operator.sage``.  This avoids loading
``ore_algebra`` during a dense coefficient scan.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from q32_cartier_packet_audit import shell_batch


def load_origin_coefficients(path: Path):
    prefix = "ORIGIN_OPERATOR_JSON "
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(prefix):
            payload = json.loads(line[len(prefix) :])
            return payload["coefficients"], payload["multiplier"]
    raise RuntimeError(f"no {prefix.strip()} row in {path}")


def polynomial_value(coefficients, value):
    output = 0
    for coefficient in reversed(coefficients):
        output = output * value + coefficient
    return output


def gcd_many(values):
    output = 0
    for value in values:
        output = math.gcd(output, abs(value))
    return output


def radical_divides(left, right):
    """Test rad(left) | rad(right) without factoring either integer."""

    remaining = abs(left)
    right = abs(right)
    while remaining > 1:
        common = math.gcd(remaining, right)
        if common == 1:
            return False
        remaining //= common
    return True


def bareiss_determinant(rows):
    """Fraction-free determinant over the integers."""

    matrix = [list(row) for row in rows]
    size = len(matrix)
    sign = 1
    previous = 1
    for pivot_index in range(size - 1):
        if matrix[pivot_index][pivot_index] == 0:
            swap = next(
                (
                    row
                    for row in range(pivot_index + 1, size)
                    if matrix[row][pivot_index] != 0
                ),
                None,
            )
            if swap is None:
                return 0
            matrix[pivot_index], matrix[swap] = (
                matrix[swap],
                matrix[pivot_index],
            )
            sign = -sign
        pivot = matrix[pivot_index][pivot_index]
        for row in range(pivot_index + 1, size):
            for column in range(pivot_index + 1, size):
                matrix[row][column] = (
                    matrix[row][column] * pivot
                    - matrix[row][pivot_index]
                    * matrix[pivot_index][column]
                ) // previous
        previous = pivot
    return sign * matrix[-1][-1]


def terminal_record(index, origin_coefficients, width=7):
    moment = index - 1
    nodes = list(range(moment - width, moment + 1))
    required = list(range(moment - width - 1, moment + 2))
    assert required[0] >= 1

    y_values = shell_batch(moment, nodes)
    u_values = [
        polynomial_value(row, index) for row in origin_coefficients
    ]
    shifted = [
        shell_batch(index + shift, required)
        for shift in range(len(origin_coefficients))
    ]
    shell_columns = [
        {
            node: shifted[shift][node - 1]
            + shifted[shift][node + 1]
            for node in nodes
        }
        for shift in range(len(origin_coefficients))
    ]
    z_values = {
        node: sum(
            u_values[shift] * shell_columns[shift][node]
            for shift in range(len(origin_coefficients))
        )
        for node in nodes
    }

    ordered_nodes = list(reversed(nodes))
    raw_vectors = [
        (y_values[node], z_values[node]) for node in ordered_nodes
    ]
    raw_minors = [
        left[0] * right[1] - left[1] * right[0]
        for left_index, left in enumerate(raw_vectors)
        for right in raw_vectors[left_index + 1 :]
    ]
    raw_content = gcd_many(raw_minors)

    observability = [
        [
            *(
                shell_columns[shift][node]
                for shift in range(len(origin_coefficients))
            ),
            -y_values[node],
        ]
        for node in ordered_nodes
    ]

    directional_certificates = []
    for selection in itertools.combinations(range(width + 1), 5):
        selected = [observability[row] for row in selection]
        cofactors = [
            (-1) ** column
            * bareiss_determinant(
                [
                    row[:column] + row[column + 1 :]
                    for row in selected
                ]
            )
            for column in range(6)
        ]
        directional_certificates.extend(
            u_values[left] * cofactors[right]
            - u_values[right] * cofactors[left]
            for left in range(5)
            for right in range(left + 1, 5)
        )
    directional_content = gcd_many(directional_certificates)
    assert radical_divides(raw_content, directional_content)

    origin_content = gcd_many(u_values)
    boundary_product = math.prod(index + offset for offset in range(2, 8))
    conjectural_bound = origin_content * boundary_product
    return {
        "raw_content": raw_content,
        "directional_content": directional_content,
        "origin_content": origin_content,
        "boundary_product": boundary_product,
        "conjectural_raw_divisibility": (
            conjectural_bound % raw_content == 0
        ),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--start", type=int, default=10)
    parser.add_argument("--stop", type=int, default=40)
    parser.add_argument(
        "--operator-dump",
        type=Path,
        default=Path("/tmp/p32_origin_operator_dump.out"),
    )
    parser.add_argument(
        "--assert-short-bound",
        action="store_true",
        help=(
            "assert the still-conjectural raw divisibility "
            "by origin content times (n+2)...(n+7)"
        ),
    )
    parser.add_argument(
        "--origin-resultant",
        action="store_true",
        help="compute the fixed resultant Res(u_0,u_3) with SymPy",
    )
    args = parser.parse_args()
    origin_coefficients, _ = load_origin_coefficients(
        args.operator_dump
    )
    origin_resultant = None
    if args.origin_resultant:
        import sympy

        variable = sympy.symbols("n")
        polynomials = [
            sympy.Poly.from_list(
                list(reversed(row)),
                gens=variable,
                domain=sympy.ZZ,
            )
            for row in origin_coefficients
        ]
        assert sympy.gcd(polynomials[0], polynomials[3]).degree() == 0
        origin_resultant = abs(
            int(sympy.resultant(polynomials[0], polynomials[3]))
        )
        assert origin_resultant
        print(
            "ORIGIN_RESULTANT",
            "pair",
            (0, 3),
            "bits",
            origin_resultant.bit_length(),
            "digits",
            len(str(origin_resultant)),
        )

    failures = []
    for index in range(args.start, args.stop + 1):
        record = terminal_record(index, origin_coefficients)
        if not record["conjectural_raw_divisibility"]:
            failures.append(index)
        if origin_resultant is not None:
            assert (
                origin_resultant % record["origin_content"] == 0
            )
        print(
            "DIRECTIONAL_COFACTOR",
            "n",
            index,
            "raw_bits",
            record["raw_content"].bit_length(),
            "directional_bits",
            record["directional_content"].bit_length(),
            "origin_bits",
            record["origin_content"].bit_length(),
            "short_bound",
            record["conjectural_raw_divisibility"],
        )
    if args.assert_short_bound:
        assert not failures, failures
    print("THEOREM_RAD_RAW_DIVIDES_RAD_DIRECTIONAL=PASS")
    print("CONJECTURAL_SHORT_BOUND_FAILURES", failures)


if __name__ == "__main__":
    main()
