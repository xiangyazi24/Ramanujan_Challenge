#!/usr/bin/env python3.12
"""Modular search for g_k - C_k f_k as a rational Ore associate of f."""

import sys

sys.path.insert(0, "/private/tmp/p25-py312-deps")
from flint import nmod_mat


PRIME = 1_000_000_007
COUNT = 1400


def inverse(value):
    return pow(value % PRIME, PRIME - 2, PRIME)


def matrix_entries(n):
    return [
        [(-2*n-5)*(n+3)**2*(136*n**4+1424*n**3+5548*n**2+9551*n+6141),
         384*n**6+6384*n**5+44168*n**4+162698*n**3+336377*n**2+369933*n+169011,
         -(480*n**4+4980*n**3+19210*n**2+32690*n+20730)],
        [(n+2)**2*(n+3)**2*(4*n+10)*(48*n**3+386*n**2+1017*n+879),
         (n+2)**2*(-272*n**5-3848*n**4-21732*n**3-61184*n**2-85761*n-47808),
         (n+2)**2*(320*n**3+2540*n**2+6610*n+5640)],
        [(-4*n-10)*(n+2)**2*(n+3)**2*(32*n**4+302*n**3+1037*n**2+1530*n+813),
         (n+2)**2*(192*n**6+2984*n**5+19116*n**4+64452*n**3+120256*n**2+117279*n+46476),
         (n+2)**2*(-16*n**5-408*n**4-2912*n**3-8884*n**2-12254*n-6240)],
    ]


def gauge(n):
    return -2*(n+2)**2*(n+3)**2*(2*n+5)*(2*n+7)**2


def row_step(row, matrix, divisor):
    scale = inverse(divisor)
    return [sum(row[i]*matrix[i][j] for i in range(3)) % PRIME * scale % PRIME
            for j in range(3)]


qrow = [33750, -36000, 9000]
prow = [30921, -32972, 8240]
qvalues = []
pvalues = []
for n in range(COUNT):
    qvalues.append(qrow[0] % PRIME)
    pvalues.append(prow[0] % PRIME)
    transition = matrix_entries(n)
    qrow = row_step(qrow, transition, gauge(n))
    prow = row_step(prow, transition, gauge(n))


fvalues = []
gvalues = []
basis_row = [1]
for n in range(COUNT):
    if n:
        basis_row = [basis_row[k]*(n+k) % PRIME*inverse(n-k) % PRIME
                     for k in range(n)]
        diagonal = pow(2, n, PRIME)
        central = 1
        for j in range(1, n+1):
            central = central*(n+j) % PRIME*inverse(j) % PRIME
        basis_row.append(diagonal*central*central % PRIME)
    fq = qvalues[n]
    gp = pvalues[n]
    for k in range(n):
        fq = (fq-fvalues[k]*basis_row[k]) % PRIME
        gp = (gp-gvalues[k]*basis_row[k]) % PRIME
    diagonal_inverse = inverse(basis_row[n])
    fvalues.append(fq*diagonal_inverse % PRIME)
    gvalues.append(gp*diagonal_inverse % PRIME)
    if n and n % 200 == 0:
        print("data", n, flush=True)


partial = 0
target = []
for k in range(COUNT):
    target.append((gvalues[k]-partial*fvalues[k]) % PRIME)
    term = (-1 if k % 2 else 1)*inverse((2*k+1)**2)
    partial = (partial+term) % PRIME


def verify(order, degree, vector, start):
    blocks = order+2
    for k in range(start, COUNT-order):
        values = [target[k]]+[fvalues[k+shift] for shift in range(order+1)]
        total = 0
        powers = [1]
        for _ in range(degree):
            powers.append(powers[-1]*k % PRIME)
        for block in range(blocks):
            for power in range(degree+1):
                total += (vector[block*(degree+1)+power]
                          * powers[power]*values[block])
        if total % PRIME:
            return False, k
    return True, None


for order in range(0, 8):
    blocks = order+2
    print("order", order, flush=True)
    maximum_degree = min(180, (COUNT-order-40)//blocks-1)
    for degree in range(0, maximum_degree+1):
        unknowns = blocks*(degree+1)
        if degree > 20 and degree % 5:
            continue
        rows = []
        for k in range(unknowns+12):
            values = [target[k]]+[fvalues[k+shift] for shift in range(order+1)]
            power_values = [1]
            for _ in range(degree):
                power_values.append(power_values[-1]*k % PRIME)
            rows.append([values[block]*power_values[power] % PRIME
                         for block in range(blocks)
                         for power in range(degree+1)])
        kernel, nullity = nmod_mat(rows, PRIME).nullspace()
        if not nullity:
            continue
        print("candidate", order, degree, "nullity", nullity, flush=True)
        for column in range(nullity):
            vector = [int(kernel[row, column]) for row in range(unknowns)]
            # The first polynomial must be nonzero, otherwise this is only a
            # recurrence among shifts of f.
            if not any(vector[:degree+1]):
                continue
            good, failure = verify(order, degree, vector, unknowns+12)
            print("verified", good, failure, flush=True)
            if good:
                print("FOUND", order, degree, flush=True)
                for block in range(blocks):
                    coefficients = vector[block*(degree+1):(block+1)*(degree+1)]
                    while coefficients and not coefficients[-1]:
                        coefficients.pop()
                    print(block, coefficients, flush=True)
                raise SystemExit

print("none", flush=True)
