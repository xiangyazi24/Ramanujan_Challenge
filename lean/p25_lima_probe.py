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
for shift in range(0, 10):
    d = [g[k]-f[k]*partials[k+shift] for k in range(count)]
    ratios = [d[k+1]/d[k] for k in range(count-1) if d[k] != 0]
    points = [(sp.Integer(k), sp.Rational(ratios[k].numerator, ratios[k].denominator))
              for k in range(min(28, len(ratios)))]
    found = None
    for numerator_degree in range(0, 14):
        for sample_count in [14, 18, 22]:
            if sample_count > len(points) or sample_count <= numerator_degree:
                continue
            candidate = sp.cancel(sp.rational_interpolate(
                points[:sample_count], numerator_degree, X=x))
            if all(sp.cancel(candidate.subs(x, point)-value) == 0
                   for point, value in points):
                found = candidate
                break
        if found is not None:
            break
    print("shift", shift, "fit", found)

print("\naffine combinations of adjacent Lima partial sums")
# Solve r_k = u_k A_(k+s) + (1-u_k) A_(k+s+2), then test u_k.
for shift in range(0, 8):
    weights = []
    for k in range(count):
        r = g[k]/f[k]
        a, b = partials[k+shift], partials[k+shift+2]
        weights.append((r-b)/(a-b))
    points = [(sp.Integer(k), sp.Rational(v.numerator, v.denominator))
              for k, v in enumerate(weights[:28])]
    found = None
    for numerator_degree in range(0, 14):
        candidate = sp.cancel(sp.rational_interpolate(points[:20], numerator_degree, X=x))
        if all(sp.cancel(candidate.subs(x, point)-value) == 0 for point, value in points):
            found = candidate
            break
    print("shift", shift, "weight fit", found,
          "range", float(min(weights)), float(max(weights)))
