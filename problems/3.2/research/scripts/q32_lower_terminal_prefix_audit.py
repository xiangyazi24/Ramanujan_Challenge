#!/usr/bin/env python3
"""Exact lower-packet prefix audit after top-half Freshman reduction.

For ``n = p + r`` and terminal offsets ``0 <= t <= width``, the
Freshman decomposition reduces the actual terminal Y/W vectors modulo
``p`` to ``X_t(r)`` below (away from the bounded borrow strip).  This
script forms the same signed-binomial transform and cumulative Newton
vectors as the characteristic-zero terminal packet, but with

    (-1)^L binom(n,L) = (-1)^L binom(r,L)  (mod p).

It records the gcd of the reduced exterior-prefix determinants.  This
is diagnostic only.  In particular, this bounded endpoint packet is
not the central selector packet: genuine target primes need not divide
its prefix determinants.  A theorem would still need both an exact
symbolic Bezout bound in the residue variable and a target-preserving
transport from the central packet.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from q32_cartier_packet_audit import LAMBDA, coefficient


def load_origin_operator(path: Path):
    prefix = "ORIGIN_OPERATOR_JSON "
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(prefix):
            payload = json.loads(line[len(prefix) :])
            return payload["coefficients"], payload["multiplier"]
    raise RuntimeError(f"no {prefix.strip()} row in {path}")


def polynomial_value(coefficients, value):
    output = 0
    for entry in reversed(coefficients):
        output = output * value + entry
    return output


def gcd_many(values):
    output = 0
    for value in values:
        output = math.gcd(output, abs(value))
    return output


def weighted_correlation(moment, exponent):
    return sum(
        weight
        * coefficient(
            moment,
            exponent * point[0],
            exponent * point[1],
            exponent * point[2],
        )
        for point, weight in LAMBDA.items()
    )


def lower_raw_vector(
    residue, terminal_offset, origin_coefficients, multiplier
):
    exponent = residue - 1 - terminal_offset
    y_value = weighted_correlation(residue - 1, exponent)
    z_value = 0
    for shift, polynomial in enumerate(origin_coefficients):
        scalar = polynomial_value(polynomial, residue)
        z_value += scalar * (
            weighted_correlation(residue + shift, exponent - 1)
            + weighted_correlation(residue + shift, exponent + 1)
        )
    w_value = (
        z_value
        - polynomial_value(multiplier, residue) * y_value
    )
    return (y_value, w_value)


def determinant(left, right):
    return left[0] * right[1] - left[1] * right[0]


def add_scaled(left, scalar, right):
    return (
        left[0] + scalar * right[0],
        left[1] + scalar * right[1],
    )


def record(residue, origin_coefficients, multiplier, width=7):
    raw = [
        lower_raw_vector(
            residue, offset, origin_coefficients, multiplier
        )
        for offset in range(width + 1)
    ]
    raw_content = gcd_many(
        determinant(raw[left], raw[right])
        for left in range(width + 1)
        for right in range(left + 1, width + 1)
    )

    # B_L = sum_{t=0}^L (-1)^t binom(L,t) X_t.
    differences = []
    for order in range(width + 1):
        differences.append(
            (
                sum(
                    (-1) ** offset
                    * math.comb(order, offset)
                    * raw[offset][0]
                    for offset in range(order + 1)
                ),
                sum(
                    (-1) ** offset
                    * math.comb(order, offset)
                    * raw[offset][1]
                    for offset in range(order + 1)
                ),
            )
        )

    cumulative = [raw[0]]
    exterior = []
    for order in range(1, width + 1):
        exterior.append(
            determinant(differences[order], cumulative[-1])
        )
        scalar = (-1) ** order * math.comb(residue, order)
        cumulative.append(
            add_scaled(cumulative[-1], scalar, differences[order])
        )

    prefix_content = gcd_many(exterior)
    separator_contents = [
        gcd_many(vector) for vector in cumulative[:-1]
    ]
    return {
        "raw_content": raw_content,
        "prefix_content": prefix_content,
        "ratio": prefix_content // math.gcd(
            prefix_content, raw_content
        ),
        "separator_contents": separator_contents,
    }


def factor_small(value):
    factors = []
    remaining = abs(value)
    prime = 2
    while prime * prime <= remaining:
        if remaining % prime == 0:
            exponent = 0
            while remaining % prime == 0:
                remaining //= prime
                exponent += 1
            factors.append((prime, exponent))
        prime = 3 if prime == 2 else prime + 2
    if remaining > 1:
        factors.append((remaining, 1))
    return factors


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--start", type=int, default=12)
    parser.add_argument("--stop", type=int, default=80)
    parser.add_argument("--width", type=int, default=7)
    parser.add_argument(
        "--operator-dump",
        type=Path,
        default=Path("/tmp/p32_origin_operator_dump.out"),
    )
    args = parser.parse_args()
    coefficients, multiplier = load_origin_operator(
        args.operator_dump
    )
    for residue in range(args.start, args.stop + 1):
        entry = record(
            residue, coefficients, multiplier, args.width
        )
        print(
            "LOWER_PREFIX",
            "r",
            residue,
            "raw",
            factor_small(entry["raw_content"]),
            "prefix",
            factor_small(entry["prefix_content"]),
            "ratio",
            factor_small(entry["ratio"]),
            "h",
            tuple(
                factor_small(value)
                for value in entry["separator_contents"]
            ),
        )


if __name__ == "__main__":
    main()
