#!/usr/bin/env python3
"""Temporary modular scan for recurrences of Wilson/challenge cross-products."""

import sys
import numpy as np


PRIME = 1_000_000_007
COUNT = int(sys.argv[1]) if len(sys.argv) > 1 else 520


def inv(value):
    return pow(value % PRIME, PRIME - 2, PRIME)


def wilson_coefficients(m):
    a4 = 12265 + 29296*m + 26176*m**2 + 10368*m**3 + 1536*m**4
    c4 = 313 + 1904*m + 4288*m**2 + 4224*m**3 + 1536*m**4
    b10 = (
        111992515 + 1144683736*m + 5147619352*m**2
        + 13412393984*m**3 + 22433518592*m**4
        + 25185342464*m**5 + 19235018752*m**6
        + 9876373504*m**7 + 3265527808*m**8
        + 628359168*m**9 + 53477376*m**10
    )
    a = 4*(m+1)**2*(4*m+1)**2*(4*m+3)**2*a4
    c = 4*(m+2)**2*(4*m+5)**2*(4*m+7)**2*c4
    return a % PRIME, b10 % PRIME, c % PRIME


def challenge_matrix(n):
    entries = [
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
    return [[entries[3*i+j] % PRIME * scale % PRIME for j in range(3)]
            for i in range(3)]


def generate():
    u = [1, 19 * inv(4) % PRIME]
    v = [0, 313 * inv(72) % PRIME]
    for m in range(COUNT):
        a, b, c = wilson_coefficients(m)
        u.append((b*u[-1] - a*u[-2]) * inv(c) % PRIME)
        v.append((b*v[-1] - a*v[-2]) * inv(c) % PRIME)
    p = [30921 % PRIME, -32972 % PRIME, 8240 % PRIME]
    q = [33750 % PRIME, -36000 % PRIME, 9000 % PRIME]
    sequences = [[] for _ in range(3)]
    for n in range(COUNT):
        for j in range(3):
            sequences[j].append((p[j]*u[n+2] - q[j]*v[n+2]) % PRIME)
        matrix = challenge_matrix(n)
        p = [sum(p[i]*matrix[i][j] for i in range(3)) % PRIME for j in range(3)]
        q = [sum(q[i]*matrix[i][j] for i in range(3)) % PRIME for j in range(3)]
    return sequences


def nullity(rows, columns):
    rows = np.array(rows, dtype=np.int64)
    rank = 0
    for column in range(columns):
        candidates = np.flatnonzero(rows[rank:, column])
        if not len(candidates):
            continue
        pivot = rank + int(candidates[0])
        rows[[rank, pivot]] = rows[[pivot, rank]]
        scale = inv(int(rows[rank, column]))
        rows[rank, column:] = rows[rank, column:] * scale % PRIME
        factors = rows[rank + 1:, column].copy()
        active = np.flatnonzero(factors)
        if len(active):
            indices = rank + 1 + active
            rows[indices, column:] = (
                rows[indices, column:]
                - factors[active, None] * rows[rank, column:]
            ) % PRIME
        rank += 1
        if rank == len(rows):
            break
    return columns - rank


def scan(sequence, order, degree):
    columns = (order + 1) * (degree + 1)
    equation_count = min(len(sequence) - order, columns + 8)
    rows = []
    for n in range(equation_count):
        powers = [1]
        for _ in range(degree):
            powers.append(powers[-1] * n % PRIME)
        rows.append([sequence[n+shift] * powers[power] % PRIME
                     for shift in range(order + 1)
                     for power in range(degree + 1)])
    return nullity(rows, columns)


sequences = generate()
print("generated", COUNT, flush=True)
for order in range(3, 9):
    for degree in range(0, 71, 2):
        columns = (order + 1) * (degree + 1)
        if columns + order + 8 > COUNT:
            break
        dimension = scan(sequences[0], order, degree)
        if dimension:
            print("FOUND", order, degree, dimension, flush=True)
            raise SystemExit
    print("order", order, "none", flush=True)
print("none", flush=True)
