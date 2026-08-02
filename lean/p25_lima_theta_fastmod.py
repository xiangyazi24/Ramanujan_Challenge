#!/usr/bin/env python3.12
"""Fast modular recurrence search for the Lima coefficient split."""

import sys

sys.path.insert(0, "/private/tmp/p25-py312-deps")
from flint import nmod_mat


P = 1_000_000_007
COUNT = 500


def inv(value):
    return pow(value % P, P - 2, P)


def transition(n):
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
    delta = -2*(n+2)**2*(n+3)**2*(2*n+5)*(2*n+7)**2
    scale = inv(delta)
    return [[raw[3*i+j] % P * scale % P for j in range(3)] for i in range(3)]


def row_step(row, value):
    return [sum(row[i]*value[i][j] for i in range(3)) % P for j in range(3)]


choose = [[0]*(2*COUNT+2) for _ in range(2*COUNT+2)]
choose[0][0] = 1
for n in range(1, len(choose)):
    choose[n][0] = choose[n][n] = 1
    for k in range(1, n):
        choose[n][k] = (choose[n-1][k-1] + choose[n-1][k]) % P


def basis(n, k):
    return pow(2, k, P)*choose[2*k][k]*choose[n][k]*choose[n+k][k] % P


q = [33750, -36000, 9000]
p = [30921, -32972, 8240]
qvalues = []
pvalues = []
for n in range(COUNT + 1):
    qvalues.append(q[0] % P)
    pvalues.append(p[0] % P)
    value = transition(n)
    q = row_step(q, value)
    p = row_step(p, value)


def decompose(values):
    answer = []
    for n, value in enumerate(values):
        residual = value - sum(answer[k]*basis(n, k) for k in range(n))
        answer.append(residual % P * inv(basis(n, n)) % P)
    return answer


f = decompose(qvalues)
g = decompose(pvalues)


def term(k):
    return ((-1 if k % 2 else 1)*(3*k+2)*pow(8, k, P)
            * inv(2*(2*k+1)**3*choose[2*k][k]**3)) % P


partial = 0
A = []
B = []
for k in range(COUNT + 1):
    value = (g[k] - partial*f[k]) % P * inv(term(k)) % P
    A.append(value)
    B.append((f[k] - value) % P)
    partial = (partial + term(k)) % P
print("data ready", flush=True)


def recurrence_matrix(sequence, order, degree, extra=8):
    columns = (order + 1)*(degree + 1)
    rows = []
    for k in range(columns + extra):
        powers = [1]
        for _ in range(degree):
            powers.append(powers[-1]*k % P)
        rows.append([sequence[k+j]*powers[d] % P
                     for j in range(order+1) for d in range(degree+1)])
    return nmod_mat(rows, P)


def verify(sequence, order, degree, vector, start):
    for k in range(start, len(sequence)-order):
        powers = [1]
        for _ in range(degree):
            powers.append(powers[-1]*k % P)
        total = sum(vector[j*(degree+1)+d]*powers[d]*sequence[k+j]
                    for j in range(order+1) for d in range(degree+1)) % P
        if total:
            return False
    return True


for label, sequence in (("A", A), ("B", B)):
    print("search", label, flush=True)
    hit = None
    for order in range(1, 13):
        for degree in range(0, 81, 10):
            columns = (order+1)*(degree+1)
            if columns+20+order >= len(sequence):
                continue
            matrix = recurrence_matrix(sequence, order, degree)
            nullity = matrix.ncols() - matrix.rank()
            print(order, degree, nullity, flush=True)
            if nullity:
                low = max(0, degree-9)
                for exact_degree in range(low, degree+1):
                    exact = recurrence_matrix(sequence, order, exact_degree)
                    kernel, exact_nullity = exact.nullspace()
                    if not exact_nullity:
                        continue
                    vector = [int(kernel[i, 0]) for i in range(exact.ncols())]
                    if verify(sequence, order, exact_degree, vector, exact.nrows()):
                        hit = (order, exact_degree, vector)
                        break
                break
        if hit:
            break
    print("HIT", label, None if hit is None else hit[:2], flush=True)
    if hit:
        order, degree, vector = hit
        for shift in range(order+1):
            coefficients = vector[shift*(degree+1):(shift+1)*(degree+1)]
            while coefficients and coefficients[-1] == 0:
                coefficients.pop()
            print("coefficient", shift, coefficients, flush=True)
