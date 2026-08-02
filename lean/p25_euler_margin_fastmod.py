#!/usr/bin/env python3.12
"""Guess polynomial recurrences for the rational Euler-squeeze margins."""

import sys

sys.path.insert(0, "/private/tmp/p25-py312-deps")
from flint import nmod_mat


P = 1_000_000_007
COUNT = 900


def inv(value):
    return pow(value % P, P - 2, P)


def positive_transition(n):
    entries = [
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
    denominator = 2*(n+2)**2*(n+3)**2*(2*n+5)*(2*n+7)**2
    scale = inv(denominator)
    return [[entries[3*i+j] % P * scale % P for j in range(3)]
            for i in range(3)]


def row_step(row, matrix):
    return [sum(row[i]*matrix[i][j] for i in range(3)) % P
            for j in range(3)]


def euler_step(a, odd, partial, k):
    partial = (partial + a*odd) % P
    a = a*(k+1) % P*inv(2*k+3) % P
    odd = (odd + inv(2*k+3)) % P
    return a, odd, partial


p = [30921, 32972, 8240]
q = [33750, 36000, 9000]
a, odd, partial = inv(2), 1, 0
k = 0
for _ in range(12):
    a, odd, partial = euler_step(a, odd, partial, k)
    k += 1

lower = []
upper = []
all_errors = [[], [], []]
for n in range(COUNT):
    for j in range(3):
        all_errors[j].append((q[j]*partial-p[j]) % P)
    lower.append((q[2]*partial-p[2]) % P)
    upper.append((p[0]-q[0]*(partial+2*a*odd)) % P)
    value = positive_transition(n)
    p = row_step(p, value)
    q = row_step(q, value)
    for _ in range(6):
        a, odd, partial = euler_step(a, odd, partial, k)
        k += 1

print("data ready", flush=True)


def recurrence_matrix(sequence, order, degree, rows):
    data = []
    for k in range(rows):
        powers = [1]
        for _ in range(degree):
            powers.append(powers[-1]*k % P)
        data.append([sequence[k+j]*powers[d] % P
                     for j in range(order+1) for d in range(degree+1)])
    return nmod_mat(data, P)


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


def search(label, sequence):
    print("search", label, flush=True)
    for order in range(1, 19):
        for degree in range(0, 121, 10):
            columns = (order+1)*(degree+1)
            rows = columns+12
            if rows+order >= len(sequence):
                continue
            matrix = recurrence_matrix(sequence, order, degree, rows)
            nullity = matrix.ncols()-matrix.rank()
            print(order, degree, nullity, flush=True)
            if not nullity:
                continue
            for exact_degree in range(max(0, degree-9), degree+1):
                exact_columns = (order+1)*(exact_degree+1)
                exact_rows = exact_columns+12
                exact = recurrence_matrix(sequence, order, exact_degree,
                                          exact_rows)
                kernel, exact_nullity = exact.nullspace()
                if not exact_nullity:
                    continue
                for column in range(exact_nullity):
                    vector = [int(kernel[i, column])
                              for i in range(exact.ncols())]
                    if verify(sequence, order, exact_degree, vector,
                              exact.nrows()):
                        print("HIT", label, order, exact_degree, flush=True)
                        for shift in range(order+1):
                            coefficients = vector[
                                shift*(exact_degree+1):(shift+1)*(exact_degree+1)]
                            while coefficients and coefficients[-1] == 0:
                                coefficients.pop()
                            print("coefficient", shift, coefficients, flush=True)
                        return order, exact_degree, vector
            raise RuntimeError("coarse nullspace did not refine")
    print("NONE", label, flush=True)
    return None


search("lower", lower)
search("upper", upper)
