#!/usr/bin/env python3.12
"""Fast modular recurrence search for the two Lima bracket margins."""

import sys

sys.path.insert(0, "/private/tmp/p25-py312-deps")
from flint import nmod_mat


prime = 1000000007
terms = 1000


def divide(a, b):
    return a % prime * pow(b % prime, prime-2, prime) % prime


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
    inverse = pow(gauge % prime, prime-2, prime)
    return [[raw[3*i+j] % prime * inverse % prime for j in range(3)] for i in range(3)]


def ratio(m):
    return -divide((m+1)**3*(3*m+5), (2*m+3)**3*(3*m+2)) % prime


def row_mul(row, C):
    return [sum(row[i]*C[i][j] for i in range(3)) % prime for j in range(3)]


current_term = 1
partial = 0
for m in range(5):
    if m < 4:
        partial = (partial+current_term) % prime
    if m < 4:
        current_term = current_term*ratio(m) % prime
c4 = current_term

p = [30921, 32972, 8240]
q = [33750, 36000, 9000]
lower = [(q[j]*partial-p[j]) % prime for j in range(3)]
tail = [q[j]*c4 % prime for j in range(3)]
upper = [(-lower[j]-tail[j]) % prime for j in range(3)]
lower_values, upper_values = [], []

for n in range(terms):
    lower_values.append(lower[2])
    upper_values.append(upper[0])
    C = challenge(n)
    m = 2*n+4
    a = ratio(m)
    s = a*ratio(m+1) % prime
    r = (1+a) % prime
    t = (-a-s) % prime
    lower = row_mul([(lower[i]+r*tail[i]) % prime for i in range(3)], C)
    upper = row_mul([(upper[i]+t*tail[i]) % prime for i in range(3)], C)
    tail = [(s*value) % prime for value in row_mul(tail, C)]


def verify(values, order, degree, vector, start):
    for n in range(start, len(values)-order):
        power = [1]
        for _ in range(degree):
            power.append(power[-1]*n % prime)
        total = 0
        for shift in range(order+1):
            polynomial = sum(vector[shift*(degree+1)+d]*power[d]
                             for d in range(degree+1)) % prime
            total = (total+polynomial*values[n+shift]) % prime
        if total:
            return False, n
    return True, None


def search(values, name):
    print("search", name, flush=True)
    for order in [4, 5, 6]:
        for degree in [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120]:
            unknowns = (order+1)*(degree+1)
            if unknowns+12+order >= len(values):
                continue
            rows = []
            for n in range(unknowns+8):
                powers = [1]
                for _ in range(degree):
                    powers.append(powers[-1]*n % prime)
                rows.append([values[n+shift]*powers[d] % prime
                             for shift in range(order+1) for d in range(degree+1)])
            kernel, nullity = nmod_mat(rows, prime).nullspace()
            print("order", order, "degree", degree, "nullity", nullity, flush=True)
            if nullity:
                vector = [int(kernel[i, 0]) for i in range(unknowns)]
                good, failure = verify(values, order, degree, vector, unknowns+8)
                print("verified", good, failure, flush=True)
                if good:
                    return order, degree, vector
    return None


for values, name in [(lower_values, "lower"), (upper_values, "upper")]:
    answer = search(values, name)
    if answer:
        order, degree, vector = answer
        print("FOUND", name, order, degree, flush=True)
        for shift in range(order+1):
            coefficients = vector[shift*(degree+1):(shift+1)*(degree+1)]
            while coefficients and not coefficients[-1]:
                coefficients.pop()
            print(shift, len(coefficients)-1, coefficients, flush=True)
