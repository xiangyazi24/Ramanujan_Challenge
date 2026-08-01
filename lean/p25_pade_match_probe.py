#!/usr/bin/env python3
"""Compare Delannoy coefficient ratios with exact Pade approximants to Catalan."""

from contextlib import redirect_stdout
from fractions import Fraction as F
from io import StringIO

with redirect_stdout(StringIO()):
    import p25_delta2_probe as data


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


targets = [data.g[k] / data.f[k] for k in range(12)]
for k, target in enumerate(targets[:6]):
    matches = []
    for L in range(25):
        for M in range(1, 25):
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
