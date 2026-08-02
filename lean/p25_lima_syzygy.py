#!/usr/bin/env python3.12
"""Modular search for polynomial syzygies of the Lima endpoint-error rows."""

import sys

sys.path.insert(0, "/private/tmp/p25-py312-deps")
from flint import nmod_mat


PRIME = 1_000_000_007
TERMS = 1100


def div(a, b):
    return a % PRIME * pow(b % PRIME, PRIME - 2, PRIME) % PRIME


def challenge(n):
    raw = [
        (2*n+5)*(n+3)**2*(136*n**4+1424*n**3+5548*n**2+9551*n+6141),
        384*n**6+6384*n**5+44168*n**4+162698*n**3+336377*n**2+369933*n+169011,
        480*n**4+4980*n**3+19210*n**2+32690*n+20730,
        (n+2)**2*(n+3)**2*(4*n+10)*(48*n**3+386*n**2+1017*n+879),
        (n+2)**2*(272*n**5+3848*n**4+21732*n**3+61184*n**2+85761*n+47808),
        (n+2)**2*(320*n**3+2540*n**2+6610*n+5640),
        (4*n+10)*(n+2)**2*(n+3)**2*(32*n**4+302*n**3+1037*n**2+1530*n+813),
        (n+2)**2*(192*n**6+2984*n**5+19116*n**4+64452*n**3+120256*n**2+117279*n+46476),
        (n+2)**2*(16*n**5+408*n**4+2912*n**3+8884*n**2+12254*n+6240),
    ]
    gauge = 2*(n+2)**2*(n+3)**2*(2*n+5)*(2*n+7)**2
    inverse = pow(gauge % PRIME, PRIME - 2, PRIME)
    return [[raw[3*i+j] % PRIME * inverse % PRIME for j in range(3)]
            for i in range(3)]


def ratio(m):
    return -div((m+1)**3*(3*m+5), (2*m+3)**3*(3*m+2)) % PRIME


def row_mul(row, matrix):
    return [sum(row[i]*matrix[i][j] for i in range(3)) % PRIME
            for j in range(3)]


term = 1
partial = 0
for m in range(5):
    if m < 4:
        partial = (partial + term) % PRIME
    if m < 4:
        term = term * ratio(m) % PRIME
c4 = term

p = [30921, 32972, 8240]
q = [33750, 36000, 9000]
lower = [(partial*q[j] - p[j]) % PRIME for j in range(3)]
tail = [c4*q[j] % PRIME for j in range(3)]
upper = [(-lower[j] - tail[j]) % PRIME for j in range(3)]
histories = {"lower": [], "upper": []}

for n in range(TERMS):
    histories["lower"].append(lower[:])
    histories["upper"].append(upper[:])
    matrix = challenge(n)
    m = 2*n + 4
    odd_ratio = ratio(m)
    two_ratio = odd_ratio * ratio(m + 1) % PRIME
    lower_force = (1 + odd_ratio) % PRIME
    upper_force = (-odd_ratio - two_ratio) % PRIME
    lower = row_mul([(lower[i] + lower_force*tail[i]) % PRIME
                     for i in range(3)], matrix)
    upper = row_mul([(upper[i] + upper_force*tail[i]) % PRIME
                     for i in range(3)], matrix)
    tail = [(two_ratio*value) % PRIME for value in row_mul(tail, matrix)]


def verify(history, degree, vector, start):
    for n in range(start, len(history)):
        powers = [1]
        for _ in range(degree):
            powers.append(powers[-1] * n % PRIME)
        total = sum(history[n][j] * vector[j*(degree+1)+d] * powers[d]
                    for j in range(3) for d in range(degree+1)) % PRIME
        if total:
            return False, n
    return True, None


for name, history in histories.items():
    print("search", name, flush=True)
    for degree in list(range(0, 21)) + [24, 28, 32, 40, 48, 56, 64, 80,
                                        96, 112, 128, 160, 192, 224, 256,
                                        288, 320]:
        columns = 3*(degree+1)
        rows = []
        for n in range(columns + 12):
            powers = [1]
            for _ in range(degree):
                powers.append(powers[-1] * n % PRIME)
            rows.append([history[n][j] * powers[d] % PRIME
                         for j in range(3) for d in range(degree+1)])
        kernel, nullity = nmod_mat(rows, PRIME).nullspace()
        print("degree", degree, "nullity", nullity, flush=True)
        if not nullity:
            continue
        for column in range(nullity):
            vector = [int(kernel[i, column]) for i in range(columns)]
            good, failure = verify(history, degree, vector, len(rows))
            if good:
                print("FOUND", name, degree, flush=True)
                for j in range(3):
                    coefficients = vector[j*(degree+1):(j+1)*(degree+1)]
                    while coefficients and not coefficients[-1]:
                        coefficients.pop()
                    print(j, coefficients, flush=True)
                raise SystemExit
            print("candidate failed", failure, flush=True)
