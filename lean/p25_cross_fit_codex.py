#!/usr/bin/env python3
"""Exact low-degree fit test for the Wilson/challenge cross-products."""

from fractions import Fraction as F
from math import comb

import sympy as sp


def gb(x, m):
    answer = F(1)
    for i in range(m):
        answer *= x - i
        answer /= i + 1
    return answer


def wilson_pair(index):
    z = 4 * index
    a = F(z - 1, 2)
    denominator = sum(
        F(comb(index, j)) * gb(a, j) * gb(a + j, j)
        for j in range(index + 1)
    )
    correction = F(z, 4) * sum(
        F(comb(index, j))
        * sum(
            gb(a + j, j - k)
            * gb(a - k, j - k)
            * F((-1) ** (k - 1), k * k * comb(j, k) ** 2)
            for k in range(1, j + 1)
        )
        for j in range(1, index + 1)
    )
    partial = sum(F((-1) ** k, (2 * k + 1) ** 2) for k in range(2 * index))
    return denominator, denominator * partial + correction / F(2 * z)


def challenge_matrix(n):
    return [
        [
            -(2*n+5)*(n+3)**2*(136*n**4+1424*n**3+5548*n**2+9551*n+6141),
            384*n**6+6384*n**5+44168*n**4+162698*n**3+336377*n**2+369933*n+169011,
            -(480*n**4+4980*n**3+19210*n**2+32690*n+20730),
        ],
        [
            (n+2)**2*(n+3)**2*(4*n+10)*(48*n**3+386*n**2+1017*n+879),
            (n+2)**2*(-272*n**5-3848*n**4-21732*n**3-61184*n**2-85761*n-47808),
            (n+2)**2*(320*n**3+2540*n**2+6610*n+5640),
        ],
        [
            -(4*n+10)*(n+2)**2*(n+3)**2*(32*n**4+302*n**3+1037*n**2+1530*n+813),
            (n+2)**2*(192*n**6+2984*n**5+19116*n**4+64452*n**3+120256*n**2+117279*n+46476),
            (n+2)**2*(-16*n**5-408*n**4-2912*n**3-8884*n**2-12254*n-6240),
        ],
    ]


count = 16
wilson = [wilson_pair(i + 2) for i in range(count)]
p = [F(30921), F(-32972), F(8240)]
q = [F(33750), F(-36000), F(9000)]
signs = [1, -1, 1]
gauge = F(1)
differences = [[], [], []]
for n in range(count):
    u, v = wilson[n]
    for j in range(3):
        differences[j].append(F(signs[j]) * (p[j] * u - q[j] * v) / gauge)
    matrix = challenge_matrix(n)
    p = [sum(p[i] * matrix[i][j] for i in range(3)) for j in range(3)]
    q = [sum(q[i] * matrix[i][j] for i in range(3)) for j in range(3)]
    gauge *= F(-2 * (n+2)**2 * (n+3)**2 * (2*n+5) * (2*n+7)**2)

x = sp.symbols("n")
for column in range(3):
    ratios = [
        differences[column][i + 1] / differences[column][i]
        for i in range(count - 1)
    ]
    print("column", column, "signs", {value > 0 for value in differences[column]},
          "first", differences[column][0], flush=True)
    found = False
    for total in range(11):
        for numerator_degree in range(total + 1):
            denominator_degree = total - numerator_degree
            unknowns = total + 2
            rows = []
            for i, value in enumerate(ratios[: unknowns + 2]):
                rational = sp.Rational(value.numerator, value.denominator)
                rows.append(
                    [sp.Integer(i) ** k for k in range(numerator_degree + 1)]
                    + [-rational * sp.Integer(i) ** k
                       for k in range(denominator_degree + 1)]
                )
            kernel = sp.Matrix(rows).nullspace()
            if len(kernel) != 1:
                continue
            vector = kernel[0]
            numerator = sum(
                vector[k] * x**k for k in range(numerator_degree + 1)
            )
            denominator = sum(
                vector[numerator_degree + 1 + k] * x**k
                for k in range(denominator_degree + 1)
            )
            if denominator != 0 and all(
                sp.cancel(
                    numerator.subs(x, i) / denominator.subs(x, i)
                    - sp.Rational(value.numerator, value.denominator)
                ) == 0
                for i, value in enumerate(ratios)
            ):
                print("fit", numerator_degree, denominator_degree,
                      sp.factor(numerator / denominator), flush=True)
                found = True
                break
        if found:
            break
    if not found:
        print("no fit of total degree <= 10", flush=True)


def fit_ratio_sequence(name, values):
    ratios = [values[i + 1] / values[i] for i in range(len(values) - 1)]
    print(name, "signs", {value > 0 for value in values}, "first", values[0],
          flush=True)
    for total in range(11):
        for numerator_degree in range(total + 1):
            denominator_degree = total - numerator_degree
            unknowns = total + 2
            rows = []
            for i, value in enumerate(ratios[: unknowns + 2]):
                rational = sp.Rational(value.numerator, value.denominator)
                rows.append(
                    [sp.Integer(i) ** k for k in range(numerator_degree + 1)]
                    + [-rational * sp.Integer(i) ** k
                       for k in range(denominator_degree + 1)]
                )
            kernel = sp.Matrix(rows).nullspace()
            if len(kernel) != 1:
                continue
            vector = kernel[0]
            numerator = sum(
                vector[k] * x**k for k in range(numerator_degree + 1)
            )
            denominator = sum(
                vector[numerator_degree + 1 + k] * x**k
                for k in range(denominator_degree + 1)
            )
            if denominator != 0 and all(
                sp.cancel(
                    numerator.subs(x, i) / denominator.subs(x, i)
                    - sp.Rational(value.numerator, value.denominator)
                ) == 0
                for i, value in enumerate(ratios)
            ):
                print("  ratio fit", numerator_degree, denominator_degree,
                      sp.factor(numerator / denominator), flush=True)
                return
    print("  no ratio fit of total degree <= 10", flush=True)


def fit_rational_values(name, values, maximum_total=12):
    print(name, "direct rational fit", flush=True)
    for total in range(maximum_total + 1):
        for numerator_degree in range(total + 1):
            denominator_degree = total - numerator_degree
            rows = []
            for i, value in enumerate(values):
                rational = sp.Rational(value.numerator, value.denominator)
                rows.append(
                    [sp.Integer(i) ** k for k in range(numerator_degree + 1)]
                    + [-rational * sp.Integer(i) ** k
                       for k in range(denominator_degree + 1)]
                )
            kernel = sp.Matrix(rows).nullspace()
            if len(kernel) != 1:
                continue
            vector = kernel[0]
            numerator = sum(vector[k] * x**k for k in range(numerator_degree + 1))
            denominator = sum(
                vector[numerator_degree + 1 + k] * x**k
                for k in range(denominator_degree + 1)
            )
            if denominator != 0 and all(
                sp.cancel(
                    numerator.subs(x, i) / denominator.subs(x, i)
                    - sp.Rational(value.numerator, value.denominator)
                ) == 0
                for i, value in enumerate(values)
            ):
                print("  fit", numerator_degree, denominator_degree,
                      sp.factor(numerator / denominator), flush=True)
                return
    print("  no direct fit", flush=True)


p = [F(30921), F(32972), F(8240)]
q = [F(33750), F(36000), F(9000)]
absolute_gauge = F(1)
projection_zero = []
projection_one = []
for n in range(count):
    u, v = wilson[n]
    qn = [value / absolute_gauge for value in q]
    pn = [value / absolute_gauge for value in p]
    determinant = qn[0] * pn[1] - qn[1] * pn[0]
    projection_zero.append((u * pn[1] - v * qn[1]) / determinant)
    projection_one.append((v * qn[0] - u * pn[0]) / determinant)
    positive = challenge_matrix(n)
    positive = [
        [F(signs[i] * signs[j] * -positive[i][j]) for j in range(3)]
        for i in range(3)
    ]
    p = [sum(p[i] * positive[i][j] for i in range(3)) for j in range(3)]
    q = [sum(q[i] * positive[i][j] for i in range(3)) for j in range(3)]
    absolute_gauge *= F(2 * (n+2)**2 * (n+3)**2 * (2*n+5) * (2*n+7)**2)

fit_ratio_sequence("projection coefficient 0", projection_zero)
fit_ratio_sequence("projection coefficient 1", projection_one)


def wilson_pair_at(index, z):
    a = F(z - 1, 2)
    denominator = sum(
        F(comb(index, j)) * gb(a, j) * gb(a + j, j)
        for j in range(index + 1)
    )
    correction = F(z, 4) * sum(
        F(comb(index, j))
        * sum(
            gb(a + j, j - k)
            * gb(a - k, j - k)
            * F((-1) ** (k - 1), k * k * comb(j, k) ** 2)
            for k in range(1, j + 1)
        )
        for j in range(1, index + 1)
    )
    half_z = z // 2
    partial = sum(F((-1) ** k, (2 * k + 1) ** 2) for k in range(half_z))
    return denominator, (
        denominator * partial + F((-1) ** half_z) * correction / F(2 * z)
    )


wilson_intervals = [
    [wilson_pair_at(m, 4 * m), wilson_pair_at(m, 4 * m + 2)]
    for m in range(2, count + 3)
]
transition_entries = [[], [], [], []]
for n in range(count):
    old = wilson_intervals[n]
    new = wilson_intervals[n + 1]
    old_det = old[0][0] * old[1][1] - old[0][1] * old[1][0]
    inverse = [
        [old[1][1] / old_det, -old[0][1] / old_det],
        [-old[1][0] / old_det, old[0][0] / old_det],
    ]
    transition = [
        [sum(new[i][k] * inverse[k][j] for k in range(2)) for j in range(2)]
        for i in range(2)
    ]
    for i in range(2):
        for j in range(2):
            transition_entries[2 * i + j].append(transition[i][j])

for i, values in enumerate(transition_entries):
    fit_ratio_sequence(f"Wilson interval transition entry {i // 2},{i % 2}", values)
    fit_rational_values(f"Wilson interval transition entry {i // 2},{i % 2}", values)


def delannoy_term(n, k):
    if k > n:
        return F(0)
    return F(2**k * comb(2 * k, k) * comb(n, k) * comb(n + k, k))


def decompose(values):
    coefficients = []
    for n, value in enumerate(values):
        remainder = value - sum(
            coefficients[k] * delannoy_term(n, k) for k in range(n)
        )
        coefficients.append(remainder / delannoy_term(n, n))
    return coefficients


raw_p = [F(30921), F(-32972), F(8240)]
raw_q = [F(33750), F(-36000), F(9000)]
normalizing_gauge = F(1)
first_p = []
first_q = []
for n in range(count):
    first_p.append(raw_p[0] / normalizing_gauge)
    first_q.append(raw_q[0] / normalizing_gauge)
    raw_matrix = challenge_matrix(n)
    raw_p = [sum(raw_p[i] * raw_matrix[i][j] for i in range(3)) for j in range(3)]
    raw_q = [sum(raw_q[i] * raw_matrix[i][j] for i in range(3)) for j in range(3)]
    normalizing_gauge *= F(-2 * (n+2)**2 * (n+3)**2 * (2*n+5) * (2*n+7)**2)

f_coefficients = decompose(first_q)
g_coefficients = decompose(first_p)
for k in range(8):
    target_ratio = g_coefficients[k] / f_coefficients[k]
    comparisons = []
    for offset in range(5):
        m = k + offset + 1
        u, v = wilson_pair_at(m, 4 * m)
        comparisons.append(target_ratio - v / u)
    print("Delannoy/Wilson", k, comparisons, flush=True)

lower_cross = []
upper_cross = []
for k in range(len(f_coefficients)):
    m = k + 1
    lower_u, lower_v = wilson_pair_at(m, 4 * m)
    upper_u, upper_v = wilson_pair_at(m, 4 * m + 2)
    lower_cross.append(g_coefficients[k] * lower_u - f_coefficients[k] * lower_v)
    upper_cross.append(g_coefficients[k] * upper_u - f_coefficients[k] * upper_v)
fit_ratio_sequence("Delannoy/lower Wilson cross", lower_cross)
fit_ratio_sequence("Delannoy/upper Wilson cross", upper_cross)
pair_cross = [
    g_coefficients[k + 1] * f_coefficients[k]
    - g_coefficients[k] * f_coefficients[k + 1]
    for k in range(len(f_coefficients) - 1)
]
fit_ratio_sequence("consecutive Delannoy ratio cross", pair_cross)
print("lower cross signs", [value > 0 for value in lower_cross], flush=True)
print("upper cross signs", [value > 0 for value in upper_cross], flush=True)

# Six-dimensional Wilson/challenge cross state, with challenge columns 1,2
# sign-flipped so every observed coordinate is positive.
raw_p = [F(30921), F(-32972), F(8240)]
raw_q = [F(33750), F(-36000), F(9000)]
for n in range(12):
    rows = wilson_intervals[n]
    state = []
    for u, v in rows:
        state.append([
            F(signs[j] * (1 if j == 0 else -1))
            * (raw_p[j] * u - raw_q[j] * v)
            for j in range(3)
        ])
    base = state[0][0]
    print("cross state", n,
          [[float(value / base) for value in row] for row in state], flush=True)
    raw_matrix = challenge_matrix(n)
    raw_p = [sum(raw_p[i] * raw_matrix[i][j] for i in range(3)) for j in range(3)]
    raw_q = [sum(raw_q[i] * raw_matrix[i][j] for i in range(3)) for j in range(3)]
