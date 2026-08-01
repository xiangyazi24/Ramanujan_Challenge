#!/usr/bin/env python3
"""Exact Padé-table probe for F(z)=sum (-z)^k/(2k+1)^2."""

from fractions import Fraction as F


def solve(a, b):
    n = len(b)
    for col in range(n):
        pivot = next(row for row in range(col, n) if a[row][col])
        a[col], a[pivot] = a[pivot], a[col]
        b[col], b[pivot] = b[pivot], b[col]
        scale = a[col][col]
        a[col] = [x / scale for x in a[col]]
        b[col] /= scale
        for row in range(n):
            if row == col:
                continue
            scale = a[row][col]
            if scale:
                a[row] = [x - scale * y for x, y in zip(a[row], a[col])]
                b[row] -= scale * b[col]
    return b


def pade(l, m):
    c = [F((-1) ** k, (2 * k + 1) ** 2) for k in range(l + m + 1)]
    # sum_{j=1}^m q_j c_{k-j} = -c_k, k=l+1,...,l+m
    qtail = solve(
        [[c[k - j] for j in range(1, m + 1)] for k in range(l + 1, l + m + 1)],
        [-c[k] for k in range(l + 1, l + m + 1)],
    ) if m else []
    q = [F(1)] + qtail
    p = [sum(q[j] * c[k-j] for j in range(min(k, m) + 1)) for k in range(l + 1)]
    return sum(p), sum(q), p, q


targets = [(30921,33750), (32972,36000), (8240,9000)]
for l in range(13):
    for m in range(1, 13):
        try:
            p, q, _, _ = pade(l,m)
        except (StopIteration, ZeroDivisionError):
            continue
        ratio=p/q
        hits=[]
        for x,y in targets:
            if ratio == F(x,y): hits.append((x,y))
        if hits:
            print('HIT',l,m,hits)
        if abs(float(ratio)-0.915965594177219) < 0.0003 and l<=m+2 and m<=l+2:
            print(l,m,ratio,float(ratio)-0.915965594177219)
