#!/usr/bin/env python3
"""Modular coboundary scan for the 4F3 trajectory found from the author code."""

from fractions import Fraction as F
import sys

import sympy as sp
from ramanujantools import Position
from ramanujantools.cmf import pFq

from p25_meijer_gauge import (
    PRIME,
    challenge,
    inv,
    inverse_transpose,
    matrix_rank,
    rat,
)


n = sp.Symbol("n")
x = sp.symbols("x:5")
y = sp.symbols("y:4")
base = Position(
    {
        x[0]: 4,
        x[1]: 4,
        x[2]: -sp.Rational(1, 2),
        x[3]: 0,
        y[0]: sp.Rational(7, 2),
        y[1]: sp.Rational(7, 2),
        y[2]: 4,
    }
)


def source_expressions(middle_sign):
    trajectory = Position(
        {
            x[0]: 0,
            x[1]: 0,
            x[2]: 0,
            x[3]: 0,
            y[0]: -2,
            y[1]: middle_sign,
            y[2]: 2,
        }
    )
    return pFq(4, 3, 1).trajectory_matrix(trajectory, base, n)


def mod_rational(value):
    value = sp.Rational(value)
    return int(value.p) % PRIME * inv(int(value.q)) % PRIME


def compile_entry(value):
    numerator, denominator = sp.cancel(value).as_numer_denom()
    return (
        [mod_rational(c) for c in sp.Poly(numerator, n).all_coeffs()],
        [mod_rational(c) for c in sp.Poly(denominator, n).all_coeffs()],
    )


def compile_matrix(matrix):
    return [[compile_entry(matrix[i, j]) for j in range(3)] for i in range(3)]


def evaluate_poly(coefficients, value):
    answer = 0
    for coefficient in coefficients:
        answer = (answer * value + coefficient) % PRIME
    return answer


def evaluate_matrix(compiled, value):
    residue = rat(value)
    answer = []
    for row in compiled:
        answer_row = []
        for numerator, denominator in row:
            den = evaluate_poly(denominator, residue)
            if not den:
                raise ZeroDivisionError
            answer_row.append(evaluate_poly(numerator, residue) * inv(den) % PRIME)
        answer.append(answer_row)
    return answer


def transpose(matrix):
    return [[matrix[j][i] for j in range(3)] for i in range(3)]


def transform(matrix, variant):
    if variant == 0:
        return matrix
    if variant == 1:
        return transpose(matrix)
    inverse_t = inverse_transpose(matrix)
    return transpose(inverse_t) if variant == 2 else inverse_t


def build_system(degree, samples, left_function, right_function):
    width = degree + 1
    columns = 9 * width
    rows = []
    for value in samples:
        left = left_function(value)
        right = right_function(value)
        powers = [pow(value % PRIME, r, PRIME) for r in range(width)]
        next_powers = [pow((value + 1) % PRIME, r, PRIME) for r in range(width)]
        for i in range(3):
            for j in range(3):
                row = [0] * columns
                for k in range(3):
                    for r in range(width):
                        # left(n) U(n+1) = U(n) right(n)
                        index = (k * 3 + j) * width + r
                        row[index] = (row[index] + left[i][k] * next_powers[r]) % PRIME
                        index = (i * 3 + k) * width + r
                        row[index] = (row[index] - powers[r] * right[k][j]) % PRIME
                rows.append(row)
    return columns - matrix_rank(rows, columns)


def main():
    maximum_degree = int(sys.argv[1]) if len(sys.argv) > 1 else 16
    signs = (int(sys.argv[2]),) if len(sys.argv) > 2 else (1, -1)
    source_variants = (int(sys.argv[3]),) if len(sys.argv) > 3 else range(4)
    target_variants = (int(sys.argv[4]),) if len(sys.argv) > 4 else range(4)
    shifts = (int(sys.argv[5]),) if len(sys.argv) > 5 else range(-2, 13)
    for middle_sign in signs:
        print("BUILD", middle_sign, flush=True)
        compiled = compile_matrix(source_expressions(middle_sign))
        for source_variant in source_variants:
            for target_variant in target_variants:
                print("VARIANTS", middle_sign, source_variant, target_variant, flush=True)
                for shift_integer in shifts:
                    shift = F(shift_integer)

                    def source(value, shift=shift):
                        return transform(evaluate_matrix(compiled, F(value) + shift), source_variant)

                    def target(value):
                        return transform(challenge(value), target_variant)

                    degree = maximum_degree
                    start = max(20, 5 - shift_integer)
                    samples = range(start, start + degree + 8)
                    try:
                        nullity = build_system(degree, samples, target, source)
                    except ZeroDivisionError:
                        continue
                    print("degree", shift, degree, nullity, flush=True)
                    if nullity:
                        print("HIT", shift, degree, nullity, flush=True)
                        return


if __name__ == "__main__":
    main()
