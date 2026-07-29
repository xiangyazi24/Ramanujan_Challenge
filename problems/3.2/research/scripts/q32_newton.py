#!/usr/bin/env python3
"""Exact height probe for the folded Newton certificates in Q451.

For J=floor((n-1)/3), form the degree-J Newton interpolant

    F_J(X) = sum_{k=0}^J (Delta^k A)_0 binom(X,k)

of the Apéry numbers A_0,...,A_J, and evaluate it at the two arguments
corresponding to the direct and reflected q=1 branches.
"""

from __future__ import annotations

import math


SAMPLES = (30, 60, 120, 240, 480, 720, 1000)


def apery_numbers(limit: int) -> list[int]:
    values = [1]
    if limit == 0:
        return values
    values.append(5)
    for j in range(1, limit):
        numerator = (
            (34 * j**3 + 51 * j**2 + 27 * j + 5) * values[j]
            - j**3 * values[j - 1]
        )
        denominator = (j + 1) ** 3
        assert numerator % denominator == 0
        values.append(numerator // denominator)
    return values


def forward_differences(values: list[int]) -> list[int]:
    row = values[:]
    coefficients = []
    while row:
        coefficients.append(row[0])
        row = [b - a for a, b in zip(row, row[1:])]
    return coefficients


def evaluate_newton(coefficients: list[int], x: int) -> int:
    total = 0
    binomial = 1
    for k, coefficient in enumerate(coefficients):
        if k:
            binomial = binomial * (x - k + 1) // k
        total += coefficient * binomial
    return total


def log_abs(value: int) -> float:
    value = abs(value)
    assert value
    shift = max(0, value.bit_length() - 53)
    return math.log(value >> shift) + shift * math.log(2)


def main() -> None:
    max_j = max((n - 1) // 3 for n in SAMPLES)
    apery = apery_numbers(max_j)
    print("n J log|F_J(n)|/n log|F_J(-n-1)|/n sum")
    for n in SAMPLES:
        j = (n - 1) // 3
        coefficients = forward_differences(apery[: j + 1])
        positive = log_abs(evaluate_newton(coefficients, n)) / n
        negative = log_abs(evaluate_newton(coefficients, -n - 1)) / n
        print(f"{n} {j} {positive:.9f} {negative:.9f} {positive + negative:.9f}")


if __name__ == "__main__":
    main()
