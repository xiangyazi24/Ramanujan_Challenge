#!/usr/bin/env python3
"""Empirical check of (HM)_k: sum_{m<X^2} (K_X(m))_k vs X^2 * lambda_X^k,
using the ACTUAL Apery zero sets Z_p for p in (X, 2X]."""
import numpy as np

def primes_upto(n):
    s = np.ones(n + 1, dtype=bool); s[:2] = False
    for i in range(2, int(n**0.5) + 1):
        if s[i]: s[i*i::i] = False
    return np.nonzero(s)[0]

def apery_zeros(p):
    """r in [0,p) with b_r == 0 mod p. O(p) with one pass (batch inverse of cubes)."""
    # iterate cleared form: y_{m+1} = P(m) y_m - m^6 y_{m-1} with y_m = (m!)^3 b_m
    # then b_m = 0 iff y_m = 0 (m! unit for m < p)
    y0, y1 = 1, 5  # b_0, b_1 -> cleared: Y_0 = 1, Y_1 = 1^3*5 = 5
    zeros = []
    if y0 % p == 0: zeros.append(0)
    if y1 % p == 0: zeros.append(1)
    Ym1, Y = 1 % p, 5 % p
    for m in range(1, p - 1):
        P = (34*m*m*m + 51*m*m + 27*m + 5) % p
        Yn = (P * Y - pow(m, 6, p) * Ym1) % p
        Ym1, Y = Y, Yn
        if Y == 0: zeros.append(m + 1)
    return zeros

def check(X):
    ps = primes_upto(2 * X)
    ps = ps[(ps > X) & (ps >= 7)]
    M = X * X
    K = np.zeros(M, dtype=np.int8)
    lam = 0.0
    zc = {}
    for p in ps:
        Z = apery_zeros(int(p))
        zc[int(p)] = len(Z)
        lam += len(Z) / p
        for r in Z:
            K[r::p] += 1
    out = [X, len(ps), lam, int(K.max())]
    for k in (2, 3, 4):
        Kf = K.astype(np.float64)
        fk = np.ones(M)
        for j in range(k):
            fk *= (Kf - j)
        fk[K < k] = 0.0
        Sk = fk.sum()
        pred = M * lam**k
        out.append(Sk / pred if pred > 0 else float('nan'))
    return out

print(f"{'X':>6} {'#p':>5} {'lambda':>8} {'maxK':>4} {'R2':>8} {'R3':>8} {'R4':>8}")
for X in (128, 256, 512, 1024, 2048):
    r = check(X)
    print(f"{r[0]:>6} {r[1]:>5} {r[2]:>8.4f} {r[3]:>4} {r[4]:>8.3f} {r[5]:>8.3f} {r[6]:>8.3f}")
