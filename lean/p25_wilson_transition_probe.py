#!/usr/bin/env python3
"""Recover the rational transition of the two Wilson--Pade bounds."""

from fractions import Fraction as F
from math import comb

import sympy as sp


def gb(x, m):
    answer = F(1)
    for i in range(m):
        answer *= x - i
        answer /= i + 1
    return answer


def pair(m, z):
    a = F(z - 1, 2)
    u = sum(F(comb(m, j)) * gb(a, j) * gb(a + j, j)
            for j in range(m + 1))
    correction = F(z, 4) * sum(
        F(comb(m, j)) * sum(
            gb(a + j, j - k) * gb(a - k, j - k)
            * F((-1) ** (k - 1), k * k * comb(j, k) ** 2)
            for k in range(1, j + 1))
        for j in range(1, m + 1))
    half = z // 2
    partial = sum(F((-1) ** k, (2 * k + 1) ** 2) for k in range(half))
    v = u * partial + F((-1) ** half) * correction / F(2 * z)
    return u, v


def bounds(m):
    return [list(pair(m, 4 * m)), list(pair(m, 4 * m + 2))]


def transition(m):
    old = bounds(m)
    new = bounds(m + 1)
    det = old[0][0] * old[1][1] - old[0][1] * old[1][0]
    inverse = [
        [old[1][1] / det, -old[0][1] / det],
        [-old[1][0] / det, old[0][0] / det],
    ]
    return [[sum(new[i][k] * inverse[k][j] for k in range(2))
             for j in range(2)] for i in range(2)]


x = sp.symbols("m")
count = 48
values = [[[] for _ in range(2)] for _ in range(2)]
for m in range(1, count + 1):
    matrix = transition(m)
    for i in range(2):
        for j in range(2):
            value = matrix[i][j]
            values[i][j].append(sp.Rational(value.numerator, value.denominator))

for i in range(2):
    for j in range(2):
        data = [(sp.Integer(m + 1), value)
                for m, value in enumerate(values[i][j])]
        print("entry", i, j, flush=True)
        found = False
        for numerator_degree in range(0, 32):
            sample_count = min(len(data), 2 * numerator_degree + 12)
            if sample_count <= numerator_degree:
                continue
            candidate = sp.factor(sp.rational_interpolate(
                data[:sample_count], numerator_degree, X=x))
            if all(sp.cancel(candidate.subs(x, point) - value) == 0
                   for point, value in data):
                numerator, denominator = sp.cancel(candidate).as_numer_denom()
                print("degrees", sp.degree(numerator, x), sp.degree(denominator, x))
                print(candidate)
                found = True
                break
        if not found:
            print("not found")
