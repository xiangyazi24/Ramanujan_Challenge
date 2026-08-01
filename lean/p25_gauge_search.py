#!/usr/bin/env sage -python
"""Temporary search for a rational 3F2/challenge gauge in Problem 2.5."""

import sys

import sympy as sp
from sage.all import GF, Matrix as SageMatrix

def challenge_matrix(n):
    entries = [
        (-2 * n - 5)
        * (n + 3) ** 2
        * (136 * n**4 + 1424 * n**3 + 5548 * n**2 + 9551 * n + 6141),
        384 * n**6
        + 6384 * n**5
        + 44168 * n**4
        + 162698 * n**3
        + 336377 * n**2
        + 369933 * n
        + 169011,
        -(480 * n**4 + 4980 * n**3 + 19210 * n**2 + 32690 * n + 20730),
        (n + 2) ** 2
        * (n + 3) ** 2
        * (4 * n + 10)
        * (48 * n**3 + 386 * n**2 + 1017 * n + 879),
        (n + 2) ** 2
        * (
            -272 * n**5
            - 3848 * n**4
            - 21732 * n**3
            - 61184 * n**2
            - 85761 * n
            - 47808
        ),
        (n + 2) ** 2 * (320 * n**3 + 2540 * n**2 + 6610 * n + 5640),
        (-4 * n - 10)
        * (n + 2) ** 2
        * (n + 3) ** 2
        * (32 * n**4 + 302 * n**3 + 1037 * n**2 + 1530 * n + 813),
        (n + 2) ** 2
        * (
            192 * n**6
            + 2984 * n**5
            + 19116 * n**4
            + 64452 * n**3
            + 120256 * n**2
            + 117279 * n
            + 46476
        ),
        (n + 2) ** 2
        * (
            -16 * n**5
            - 408 * n**4
            - 2912 * n**3
            - 8884 * n**2
            - 12254 * n
            - 6240
        ),
    ]
    delta = -2 * (n + 2) ** 2 * (n + 3) ** 2 * (2 * n + 5) * (2 * n + 7) ** 2
    return sp.Matrix(3, 3, entries) / delta


def pfq_matrix():
    n = sp.symbols("n")
    d0 = 131072 * (n + 1) ** 4 * (2 * n + 1) ** 3 * (4 * n + 1)
    d1 = 131072 * (n + 1) ** 4 * (2 * n + 1) ** 4 * (4 * n + 1)
    entries = [
        (4 * n + 5)
        * (5382144*n**7 + 21127168*n**6 + 34369024*n**5 + 29840256*n**4
           + 14790336*n**3 + 4118312*n**2 + 578120*n + 29163) / d0,
        (4 * n + 3) * (4 * n + 5)
        * (2220032*n**7 + 8093696*n**6 + 12238336*n**5 + 9900160*n**4
           + 4591616*n**3 + 1204376*n**2 + 160816*n + 7843) / (2*d0),
        (4 * n + 3) * (4 * n + 5)
        * (3702784*n**8 + 16572416*n**7 + 32495616*n**6 + 36204544*n**5
           + 24775552*n**4 + 10490976*n**3 + 2624504*n**2 + 342324*n + 16529) / (4*d0),
        (4 * n + 5)
        * (7585792*n**7 + 29718528*n**6 + 48417280*n**5 + 42366720*n**4
           + 21407968*n**3 + 6209784*n**2 + 949140*n + 58047) / d1,
        (4 * n + 3) * (4 * n + 5)
        * (3162112*n**7 + 11520000*n**6 + 17475072*n**5 + 14278656*n**4
           + 6768992*n**3 + 1854536*n**2 + 270012*n + 15895) / (2*d1),
        (4 * n + 3) * (4 * n + 5)
        * (5177344*n**8 + 23011328*n**7 + 44912640*n**6 + 50018816*n**5
           + 34462336*n**4 + 14873728*n**3 + 3876824*n**2 + 550120*n + 31949) / (4*d1),
        (4 * n + 5)
        * (555008*n**6 + 2081280*n**5 + 3192576*n**4 + 2561440*n**3
           + 1132792*n**2 + 261740*n + 24705) / (d1/2),
        (4 * n + 3) * (4 * n + 5)
        * (231424*n**6 + 804352*n**5 + 1143808*n**4 + 852128*n**3
           + 351080*n**2 + 75972*n + 6761) / d1,
        (4 * n + 3) * (4 * n + 5)
        * (385024*n**7 + 1650688*n**6 + 3076608*n**5 + 3221888*n**4
           + 2034624*n**3 + 768632*n**2 + 159512*n + 13939) / (2*d1),
    ]
    return n, sp.Matrix(3, 3, entries) / 4


def compile_matrix(matrix, symbol):
    result = []
    for value in matrix:
        numerator, denominator = sp.cancel(value).as_numer_denom()
        numerator_coefficients = [
            int(coefficient) for coefficient in sp.Poly(numerator, symbol).all_coeffs()
        ]
        denominator_coefficients = [
            int(coefficient) for coefficient in sp.Poly(denominator, symbol).all_coeffs()
        ]
        result.append((numerator_coefficients, denominator_coefficients))
    return result


def evaluate_coefficients(coefficients, value, field):
    result = field(0)
    for coefficient in coefficients:
        result = result * field(value) + field(coefficient)
    return result


def evaluate_compiled(compiled, value, field):
    entries = []
    for numerator, denominator in compiled:
        entries.append(
            evaluate_coefficients(numerator, value, field)
            / evaluate_coefficients(denominator, value, field)
        )
    return [entries[0:3], entries[3:6], entries[6:9]]


def evaluation_system(left_values, right_values, degree, samples, field):
    columns = 9 * (degree + 1)
    rows = []
    for value in samples:
        left_field = left_values[value]
        right_field = right_values[value]
        for i in range(3):
            for j in range(3):
                row = [field(0)] * columns
                for degree_index in range(degree + 1):
                    next_power = field(value + 1) ** degree_index
                    current_power = field(value) ** degree_index
                    for a in range(3):
                        # left * U(n+1)
                        variable = degree_index * 9 + a * 3 + j
                        row[variable] += left_field[i][a] * next_power
                    for b in range(3):
                        # -U(n) * right
                        variable = degree_index * 9 + i * 3 + b
                        row[variable] -= current_power * right_field[b][j]
                rows.append(row)
    return SageMatrix(field, rows)


def main():
    max_degree = int(sys.argv[1]) if len(sys.argv) > 1 else 20
    n, pfq = pfq_matrix()
    challenge = challenge_matrix(n)
    field = GF(2147483647)
    variants = {
        "B_to_A": (pfq, challenge),
        "A_to_B": (challenge, pfq),
        "Bt_to_At": (pfq.T, challenge.T),
        "At_to_Bt": (challenge.T, pfq.T),
    }
    for name, (left, right) in variants.items():
        print(name, flush=True)
        left_compiled = compile_matrix(left, n)
        right_compiled = compile_matrix(right, n)
        left_values = {
            value: evaluate_compiled(left_compiled, value, field)
            for value in range(1, max_degree + 16)
        }
        right_values = {
            value: evaluate_compiled(right_compiled, value, field)
            for value in range(1, max_degree + 16)
        }
        for degree in range(max_degree + 1):
            sample_count = degree + 14
            system = evaluation_system(
                left_values, right_values, degree, range(1, sample_count + 1), field
            )
            nullity = system.ncols() - system.rank()
            print(degree, nullity, flush=True)
            if nullity:
                break


if __name__ == "__main__":
    main()
