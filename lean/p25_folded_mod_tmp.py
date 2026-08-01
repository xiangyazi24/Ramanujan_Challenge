#!/usr/bin/env python3
"""Modular low-degree scan for a folded challenge/Wilson quotient."""

import sys

import numpy as np
import sympy as s


PRIME = 1_000_003
r = s.symbols("r")


def challenge_at(x):
    entries = [
        (-2*x-5)*(x+3)**2*(136*x**4+1424*x**3+5548*x**2+9551*x+6141),
        384*x**6+6384*x**5+44168*x**4+162698*x**3+336377*x**2+369933*x+169011,
        -(480*x**4+4980*x**3+19210*x**2+32690*x+20730),
        (x+2)**2*(x+3)**2*(4*x+10)*(48*x**3+386*x**2+1017*x+879),
        (x+2)**2*(-272*x**5-3848*x**4-21732*x**3-61184*x**2-85761*x-47808),
        (x+2)**2*(320*x**3+2540*x**2+6610*x+5640),
        (-4*x-10)*(x+2)**2*(x+3)**2*(32*x**4+302*x**3+1037*x**2+1530*x+813),
        (x+2)**2*(192*x**6+2984*x**5+19116*x**4+64452*x**3+120256*x**2+117279*x+46476),
        (x+2)**2*(-16*x**5-408*x**4-2912*x**3-8884*x**2-12254*x-6240),
    ]
    delta = -2*(x+2)**2*(x+3)**2*(2*x+5)*(2*x+7)**2
    return s.Matrix(3, 3, entries) / delta


def wilson_at(x):
    d0 = 8*(x+3)**2*(4*x+11)**2
    d1 = d0*(4*x+13)**2
    return s.Matrix([
        [
            (4*x+9)**2*(40*x**2+228*x+325)/d0,
            (1536*x**4+16512*x**3+66496*x**2+118896*x+79641)/d0,
        ],
        [
            (4*x+9)**2*(1536*x**4+18432*x**3+82880*x**2+165504*x+123841)/d1,
            (59392*x**6+1012736*x**5+7184384*x**4+27140352*x**3+
             57583336*x**2+65059404*x+30580677)/d1,
        ],
    ])


def compile_matrix(value):
    answer = []
    for entry in value:
        numerator, denominator = s.cancel(entry).as_numer_denom()
        answer.append((
            [int(c) % PRIME for c in s.Poly(numerator, r).all_coeffs()],
            [int(c) % PRIME for c in s.Poly(denominator, r).all_coeffs()],
        ))
    return answer


def horner(coefficients, point):
    answer = 0
    for coefficient in coefficients:
        answer = (answer*point + coefficient) % PRIME
    return answer


def evaluate(compiled, rows, columns, point):
    values = []
    for numerator, denominator in compiled:
        values.append(horner(numerator, point)*pow(horner(denominator, point), PRIME-2, PRIME) % PRIME)
    return np.array(values, dtype=np.int64).reshape((rows, columns))


def rank_mod(matrix):
    matrix = matrix.copy() % PRIME
    row = 0
    for column in range(matrix.shape[1]):
        nonzero = np.flatnonzero(matrix[row:, column])
        if not len(nonzero):
            continue
        pivot = row + int(nonzero[0])
        matrix[[row, pivot]] = matrix[[pivot, row]]
        matrix[row] = matrix[row]*pow(int(matrix[row, column]), PRIME-2, PRIME) % PRIME
        for start in range(0, matrix.shape[0], 200):
            indices = np.arange(start, min(start+200, matrix.shape[0]))
            indices = indices[indices != row]
            factors = matrix[indices, column]
            active = indices[factors != 0]
            if len(active):
                matrix[active] = (matrix[active] -
                    matrix[active, column, None]*matrix[row]) % PRIME
        row += 1
        if row == matrix.shape[0]:
            break
    return row


def scan(phase, q_expression, maximum_degree):
    left = challenge_at(2*r+phase)*challenge_at(2*r+phase+1)
    right = wilson_at(r+phase).T
    cleft = compile_matrix(left)
    cright = compile_matrix(right)
    qnum = [int(c) % PRIME for c in s.Poly(s.expand(q_expression), r).all_coeffs()]
    samples = list(range(1, maximum_degree+20))
    left_values = {x: evaluate(cleft, 3, 3, x) for x in samples}
    right_values = {x: evaluate(cright, 2, 2, x) for x in samples}
    for degree in range(maximum_degree+1):
        columns = 6*(degree+1)
        points = samples[:(columns+5)//6+2]
        rows = []
        for point in points:
            old_powers = [pow(point, k, PRIME) for k in range(degree+1)]
            new_powers = [pow(point+1, k, PRIME) for k in range(degree+1)]
            qold = horner(qnum, point)
            qnew = horner(qnum, point+1)
            left_value = left_values[point]
            right_value = right_values[point]
            for i in range(3):
                for j in range(2):
                    row = np.zeros(columns, dtype=np.int64)
                    for power in range(degree+1):
                        base = 6*power
                        for a in range(3):
                            row[base+2*a+j] += qold*int(left_value[i, a])*new_powers[power]
                        for b in range(2):
                            row[base+2*i+b] -= qnew*int(right_value[b, j])*old_powers[power]
                    rows.append(row % PRIME)
        data = np.array(rows, dtype=np.int64)
        nullity = columns-rank_mod(data)
        print("phase", phase, "degree", degree, "nullity", nullity, flush=True)
        if nullity:
            return degree
    return None


maximum = int(sys.argv[1]) if len(sys.argv) > 1 else 24
q_candidates = [
    1,
    (r+1)*(r+2)*(2*r+1)*(2*r+3)*(4*r+3)*(4*r+5)*(4*r+7)*(4*r+9),
    (r+1)**2*(r+2)**2*(2*r+1)**2*(2*r+3)**2*(4*r+3)**2*(4*r+5)**2*(4*r+7)**2*(4*r+9)**2,
]
for phase in (0, 1):
    for q_value in q_candidates:
        print("Q", s.factor(q_value), flush=True)
        if scan(phase, q_value, maximum) is not None:
            raise SystemExit
