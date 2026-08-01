#!/usr/bin/env python3
"""Compare Delannoy coefficient ratios with exact Pade approximants to Catalan."""

from fractions import Fraction as F
from functools import lru_cache
from math import comb


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


def delannoy_coefficients(count=12):
    rows = [[F(1), F(0), F(0)], [F(0), F(1), F(0)], [F(0), F(0), F(1)]]
    q0 = [F(33750), F(-36000), F(9000)]
    p0 = [F(30921), F(-32972), F(8240)]
    qvalues, pvalues = [], []
    for n in range(count):
        qvalues.append(sum(q0[i]*rows[i][0] for i in range(3)))
        pvalues.append(sum(p0[i]*rows[i][0] for i in range(3)))
        raw, d = matrix(n), F(delta(n))
        transition = [[F(raw[i][j], d) for j in range(3)] for i in range(3)]
        rows = [[sum(row[i]*transition[i][j] for i in range(3)) for j in range(3)]
                for row in rows]
    def invert(values):
        answer = []
        for n, value in enumerate(values):
            answer.append((value-sum(answer[k]*basis(n, k) for k in range(n)))
                          / basis(n, n))
        return answer
    return invert(qvalues), invert(pvalues)


@lru_cache(maxsize=None)
def pade_value(L, M, z=F(1)):
    """Value at z of the [L/M] Pade approximant to sum (-z)^j/(2j+1)^2."""
    a = [F((-1) ** j, (2 * j + 1) ** 2) for j in range(L + M + 1)]
    # Gaussian elimination for q_1,...,q_M.
    aug = []
    for n in range(L + 1, L + M + 1):
        aug.append([a[n-j] for j in range(1, M + 1)] + [-a[n]])
    for col in range(M):
        pivot = next(row for row in range(col, M) if aug[row][col])
        aug[col], aug[pivot] = aug[pivot], aug[col]
        scale = aug[col][col]
        aug[col] = [v / scale for v in aug[col]]
        for row in range(M):
            if row != col and aug[row][col]:
                scale = aug[row][col]
                aug[row] = [aug[row][j] - scale * aug[col][j]
                            for j in range(M + 1)]
    q = [F(1)] + [aug[j][-1] for j in range(M)]
    p = [sum(q[j] * a[n-j] for j in range(min(n, M) + 1))
         for n in range(L + 1)]
    pz = sum(v * z**j for j, v in enumerate(p))
    qz = sum(v * z**j for j, v in enumerate(q))
    return pz / qz


f, g = delannoy_coefficients()
targets = [g[k] / f[k] for k in range(12)]
for k, target in enumerate(targets[:6]):
    matches = []
    for L in range(13):
        for M in range(1, 13):
            try:
                if pade_value(L, M) == target:
                    matches.append((L, M))
            except (StopIteration, ZeroDivisionError):
                pass
    print("k", k, "target", target, "matches", matches)

print("nearest standard Pade indices")
for k, target in enumerate(targets[:8]):
    candidates = []
    for L in range(15):
        for M in range(1, 15):
            try:
                v = pade_value(L, M)
                candidates.append((abs(float(v - target)), L, M, v))
            except (StopIteration, ZeroDivisionError):
                pass
    print(k, sorted(candidates)[:4])
