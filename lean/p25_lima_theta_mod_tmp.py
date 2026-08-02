#!/usr/bin/env python3.12
"""Modular recurrence search for standard-Catalan bracket coefficients."""

import sys

sys.path.insert(0, "/private/tmp/p25-py312-deps")
from flint import nmod_mat


P = 1_000_000_007
COUNT = 900


def inv(value):
    return pow(value % P, P-2, P)


def matrix(n):
    raw = [
        (-2*n-5)*(n+3)**2*(136*n**4+1424*n**3+5548*n**2+9551*n+6141),
        384*n**6+6384*n**5+44168*n**4+162698*n**3+336377*n**2+369933*n+169011,
        -(480*n**4+4980*n**3+19210*n**2+32690*n+20730),
        (n+2)**2*(n+3)**2*(4*n+10)*(48*n**3+386*n**2+1017*n+879),
        (n+2)**2*(-272*n**5-3848*n**4-21732*n**3-61184*n**2-85761*n-47808),
        (n+2)**2*(320*n**3+2540*n**2+6610*n+5640),
        (-4*n-10)*(n+2)**2*(n+3)**2*(32*n**4+302*n**3+1037*n**2+1530*n+813),
        (n+2)**2*(192*n**6+2984*n**5+19116*n**4+64452*n**3+120256*n**2+117279*n+46476),
        (n+2)**2*(-16*n**5-408*n**4-2912*n**3-8884*n**2-12254*n-6240),
    ]
    d = -2*(n+2)**2*(n+3)**2*(2*n+5)*(2*n+7)**2
    di = inv(d)
    return [[raw[3*i+j] % P*di % P for j in range(3)] for i in range(3)]


def row_step(row, value):
    return [sum(row[i]*value[i][j] for i in range(3)) % P for j in range(3)]


def choose_table(count):
    table = [[0]*(count+1) for _ in range(count+1)]
    table[0][0] = 1
    for n in range(1, count+1):
        table[n][0] = table[n][n] = 1
        for k in range(1, n):
            table[n][k] = (table[n-1][k-1]+table[n-1][k]) % P
    return table


choose = choose_table(2*COUNT+2)


def bsummand(n, k):
    return (pow(2, k, P)*choose[2*k][k]*choose[n][k]*choose[n+k][k]) % P


q = [33750, -36000, 9000]
p = [30921, -32972, 8240]
qvalues = []
pvalues = []
for n in range(COUNT+1):
    qvalues.append(q[0] % P)
    pvalues.append(p[0] % P)
    value = matrix(n)
    q = row_step(q, value)
    p = row_step(p, value)


def decompose(values):
    result = []
    for n in range(len(values)):
        residual = values[n]-sum(result[k]*bsummand(n, k) for k in range(n))
        result.append(residual % P*inv(bsummand(n, n)) % P)
    return result


f = decompose(qvalues)
g = decompose(pvalues)


def term(k):
    return ((-1 if k % 2 else 1)*inv((2*k+1)**2)) % P


partial = 0
A = []
B = []
for k in range(COUNT+1):
    value = (g[k]-partial*f[k]) % P*inv(term(k)) % P
    A.append(value)
    B.append((f[k]-value) % P)
    partial = (partial+term(k)) % P
print("data", len(A), flush=True)


def verify(sequence, order, degree, vector, start):
    for k in range(start, len(sequence)-order):
        powers = [1]
        for _ in range(degree):
            powers.append(powers[-1]*k % P)
        total = sum(vector[j*(degree+1)+d]*powers[d]*sequence[k+j]
                    for j in range(order+1) for d in range(degree+1)) % P
        if total:
            return False, k
    return True, None


for label, sequence in (("A", A), ("B", B)):
    print("search", label, flush=True)
    found = False
    grid = [(order, degree) for order in range(6, 19)
            for degree in (10, 20, 30, 40, 50, 60, 70, 80)]
    for order, degree in grid:
        columns = (order+1)*(degree+1)
        if columns+20+order >= len(sequence):
            continue
        rows = []
        for k in range(columns+12):
            powers = [1]
            for _ in range(degree):
                powers.append(powers[-1]*k % P)
            rows.append([sequence[k+j]*powers[d] % P
                         for j in range(order+1) for d in range(degree+1)])
        kernel, nullity = nmod_mat(rows, P).nullspace()
        if not nullity:
            continue
        vector = [int(kernel[i, 0]) for i in range(columns)]
        good, failure = verify(sequence, order, degree, vector, len(rows))
        print("candidate", order, degree, nullity, good, failure, flush=True)
        if good:
            found = True
            break
    if not found:
        print("none", flush=True)
