#!/usr/bin/env python3
"""Probe the Delannoy coefficient ratios against the Lima--Guillera series."""

from fractions import Fraction as F
from math import comb

import mpmath as mp
import sympy as sp


def matrix(n):
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


def delta(n):
    return -2*(n+2)**2*(n+3)**2*(2*n+5)*(2*n+7)**2


def basis(n, k):
    return F(2**k*comb(2*k, k)*comb(n, k)*comb(n+k, k))


def coefficients(count):
    qrow = [F(33750), F(-36000), F(9000)]
    prow = [F(30921), F(-32972), F(8240)]
    qvalues, pvalues = [], []
    for n in range(count):
        qvalues.append(qrow[0])
        pvalues.append(prow[0])
        raw, d = matrix(n), F(delta(n))
        qrow = [sum(qrow[i]*raw[i][j]/d for i in range(3)) for j in range(3)]
        prow = [sum(prow[i]*raw[i][j]/d for i in range(3)) for j in range(3)]

    def invert(values):
        answer = []
        for n, value in enumerate(values):
            answer.append((value-sum(answer[k]*basis(n, k) for k in range(n)))
                          / basis(n, n))
        return answer

    return invert(qvalues), invert(pvalues)


def lima_term(m):
    return F((-1)**m * (3*m+2) * 8**m,
             2 * (2*m+1)**3 * comb(2*m, m)**3)


count = 36
f, g = coefficients(count)
partials = [F(0)]
for m in range(4*count+20):
    partials.append(partials[-1] + lima_term(m))

mp.mp.dps = 100
G = mp.catalan

print("epsilon / Lima tails")
for k in list(range(8)) + [12, 20, 30, 35]:
    r = mp.mpf(g[k].numerator)/g[k].denominator / (mp.mpf(f[k].numerator)/f[k].denominator)
    eps = r-G
    items = []
    for slope, shift in [(1, 0), (1, 1), (1, 2), (1, 3), (2, 1), (3, 3)]:
        m = slope*k+shift
        a = partials[m]
        tail = G-mp.mpf(a.numerator)/a.denominator
        items.append((slope, shift, mp.nstr(eps/tail, 12)))
    print(k, mp.nstr(eps, 12), items)

print("\ncorrection quotient fits for g-f*A_(k+s)")
x = sp.symbols("k")


def mod_rank(rows, prime):
    rows = [[entry % prime for entry in row] for row in rows]
    rank = 0
    columns = len(rows[0])
    for column in range(columns):
        pivot = next((i for i in range(rank, len(rows)) if rows[i][column]), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        inverse = pow(rows[rank][column], -1, prime)
        rows[rank] = [(entry*inverse) % prime for entry in rows[rank]]
        for i in range(len(rows)):
            if i != rank and rows[i][column]:
                scale = rows[i][column]
                rows[i] = [(rows[i][j]-scale*rows[rank][j]) % prime
                           for j in range(columns)]
        rank += 1
    return rank


def possible_degrees(values, maximum=12):
    prime = 2305843009213693951
    points = []
    for k, value in enumerate(values[:30]):
        residue = value.numerator % prime * pow(value.denominator % prime, -1, prime) % prime
        points.append((k, residue))
    answers = []
    for pdegree in range(maximum+1):
        for qdegree in range(maximum+1):
            unknowns = pdegree+qdegree+2
            if len(points) < unknowns+2:
                continue
            rows = []
            for k, value in points:
                rows.append([pow(k, j, prime) for j in range(pdegree+1)] +
                            [(-value*pow(k, j, prime)) % prime for j in range(qdegree+1)])
            if mod_rank(rows, prime) < unknowns:
                answers.append((pdegree, qdegree))
    return answers


for shift in range(0, 10):
    d = [g[k]-f[k]*partials[k+shift] for k in range(count)]
    ratios = [d[k+1]/d[k] for k in range(count-1) if d[k] != 0]
    print("shift", shift, "possible degrees", possible_degrees(ratios))

print("\naffine combinations of adjacent Lima partial sums")
# Solve r_k = u_k A_(k+s) + (1-u_k) A_(k+s+2), then test u_k.
for shift in range(0, 8):
    weights = []
    for k in range(count):
        r = g[k]/f[k]
        a, b = partials[k+shift], partials[k+shift+2]
        weights.append((r-b)/(a-b))
    print("shift", shift, "weight possible degrees", possible_degrees(weights),
          "range", float(min(weights)), float(max(weights)))
