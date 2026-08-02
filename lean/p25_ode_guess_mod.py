#!/usr/bin/env python3.12
"""Guess low-order polynomial ODEs for the normalized P2.5 generating series."""

import sys

sys.path.insert(0, "/private/tmp/p25-py312-deps")
from flint import nmod_mat


PRIME = 1_000_000_007
TERMS = 500


def inv(value):
    return pow(value % PRIME, PRIME - 2, PRIME)


def transition(n):
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


def series(seed):
    row = [value % PRIME for value in seed]
    values = []
    for n in range(TERMS):
        values.append(row[0])
        matrix = transition(n)
        row = [sum(row[i] * matrix[i][j] for i in range(3)) % PRIME
               for j in range(3)]
    return values


def falling(value, count):
    answer = 1
    for shift in range(count):
        answer = answer * (value - shift) % PRIME
    return answer


def search(value_sets, name):
    if value_sets and isinstance(value_sets[0], int):
        value_sets = [value_sets]
    print("search", name, flush=True)
    for order in range(1, 14):
        for degree in range(0, 41):
            unknowns = (order + 1) * (degree + 1)
            if unknowns + order + degree + 12 >= TERMS:
                continue
            rows = []
            # Coefficient of z^m in sum p_i(z) D^i F.
            rows_per_series = (unknowns + 8 + len(value_sets) - 1) // len(value_sets)
            for values in value_sets:
                for m in range(rows_per_series):
                    row = []
                    for derivative in range(order + 1):
                        for power in range(degree + 1):
                            source = m - power + derivative
                            if source < derivative or source >= len(values):
                                row.append(0)
                            else:
                                row.append(values[source]
                                           * falling(source, derivative) % PRIME)
                    rows.append(row)
            kernel, nullity = nmod_mat(rows, PRIME).nullspace()
            if not nullity:
                continue
            for column in range(nullity):
                vector = [int(kernel[index, column]) for index in range(unknowns)]
                good = True
                for values in value_sets:
                    for m in range(rows_per_series, TERMS - order):
                        total = 0
                        cursor = 0
                        for derivative in range(order + 1):
                            for power in range(degree + 1):
                                source = m - power + derivative
                                if derivative <= source < len(values):
                                    total += (vector[cursor] * values[source]
                                              * falling(source, derivative))
                                cursor += 1
                        if total % PRIME:
                            good = False
                            break
                    if not good:
                        break
                if good:
                    print("FOUND", order, degree, "nullity", nullity, flush=True)
                    try:
                        import sympy as sp
                        z = sp.symbols("z")
                        for derivative in range(order + 1):
                            coefficients = vector[
                                derivative*(degree+1):(derivative+1)*(degree+1)]
                            polynomial = sum(coefficient*z**power
                                             for power, coefficient in enumerate(coefficients))
                            print("D", derivative,
                                  sp.factor(polynomial, modulus=PRIME), flush=True)
                    except Exception as error:
                        print("factor display failed", error, flush=True)
                    return order, degree, vector
    print("none", flush=True)
    return None


q_values = series([33750, -36000, 9000])
p_values = series([30921, -32972, 8240])
search(q_values, "q")
search(p_values, "p")
search([q_values, p_values], "common")
