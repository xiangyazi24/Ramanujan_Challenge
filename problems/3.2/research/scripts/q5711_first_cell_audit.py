#!/usr/bin/env python3
"""Exact standard-library audit for Q5711.

Checks the first quotient-cell decomposition of C_M(d), endpoint collapse at
d=M, finite-difference splitting, Newton residual gcds on the requested rows,
and deterministic hostile random cases.  No third-party packages are used.
"""
from __future__ import annotations

from functools import lru_cache
from math import comb, gcd, isqrt
from random import Random


def C(n: int, k: int) -> int:
    return comb(n, k) if 0 <= k <= n else 0


@lru_cache(maxsize=None)
def apery(n: int) -> int:
    if n == 0:
        return 1
    a0, a1 = 1, 5
    if n == 1:
        return a1
    for m in range(1, n):
        num = (34 * m**3 + 51 * m**2 + 27 * m + 5) * a1 - m**3 * a0
        den = (m + 1) ** 3
        assert num % den == 0
        a0, a1 = a1, num // den
    return a1


def apery_mod(r: int, p: int) -> int:
    if r == 0:
        return 1
    b0, b1 = 1, 5 % p
    if r == 1:
        return b1
    for m in range(1, r):
        num = ((34 * m**3 + 51 * m**2 + 27 * m + 5) * b1 - m**3 * b0) % p
        den = pow(m + 1, 3, p)
        b0, b1 = b1, num * pow(den, -1, p) % p
    return b1


def primes_upto(n: int) -> list[int]:
    s = bytearray(b"\x01") * (n + 1)
    if n >= 0:
        s[0] = 0
    if n >= 1:
        s[1] = 0
    for p in range(2, isqrt(n) + 1):
        if s[p]:
            s[p * p : n + 1 : p] = b"\x00" * (((n - p * p) // p) + 1)
    return [p for p in range(2, n + 1) if s[p]]


@lru_cache(maxsize=None)
def shell(M: int, d: int) -> int:
    """Exact C_M(d), specialized to a=1 when d>M/2."""
    if d > M:
        return apery(M)
    a = M // d
    total = 0
    for t in range(M + 1):
        X = sum(C(M, t + d * u) for u in range(-a, a + 1))
        N = 2 * M - t
        Z = sum(C(N, M + d * v) for v in range(-a, a + 1))
        total += C(M, t) * X * Z * Z
    return total


def first_cell_parts(M: int, d: int) -> tuple[int, int, int, int, int]:
    """Return b_M, full core, low boundary, high boundary, and total.

    Put r=M-d.  The identity is
      C_M(M-r) = b_M + core_M(r) + low_M(r) + high_M(r).
    """
    assert M // 2 < d <= M
    r = M - d
    core = 0
    low = 0
    for t in range(M + 1):
        A = C(M, t)
        N = 2 * M - t
        B = C(N, M)
        P = C(N, r)
        core += A * A * (2 * B * P + P * P)
        if t <= r:
            Q = C(N, r - t)
            U = C(M, r - t)
            low += A * (
                A * (2 * (B + P) * Q + Q * Q)
                + U * (B + P + Q) ** 2
            )
    high = 0
    for k in range(r + 1):
        A = C(M, k)
        U = C(M, r - k)
        B = C(M + k, k)
        P = C(M + k, r)
        high += A * U * (B + P) ** 2
    bM = apery(M)
    return bM, core, low, high, bM + core + low + high


def zeta2_apery(M: int) -> int:
    return sum(C(M, k) ** 2 * C(M + k, k) for k in range(M + 1))


def endpoint_formula(M: int) -> int:
    B = C(2 * M, M)
    return apery(M) + 2 * zeta2_apery(M) + B * B + 7 * B + 11


def delta_values(values: list[int]) -> int:
    row = values[:]
    while len(row) > 1:
        row = [row[i + 1] - row[i] for i in range(len(row) - 1)]
    return row[0]


def delta(M: int, d: int, order: int) -> int:
    return delta_values([shell(M, d + i) for i in range(order + 1)])


def weight(d: int, L: int, i: int) -> int:
    return (-1) ** i * C(d + i, i) * C(d + L + 1, L - i)


def G(M: int, d: int, L: int) -> int:
    return sum(weight(d, L, i) * shell(M, d + i) for i in range(L + 1))


def factor_trial(n: int, limit: int = 2_000_000) -> tuple[dict[int, int], int]:
    n = abs(n)
    fac: dict[int, int] = {}
    for p in primes_upto(min(limit, isqrt(n) if n else 0)):
        if p * p > n:
            break
        if n % p == 0:
            e = 0
            while n % p == 0:
                n //= p
                e += 1
            fac[p] = e
    return fac, n


def fmt_factor(n: int) -> str:
    fac, rem = factor_trial(n)
    out = [str(p) if e == 1 else f"{p}^{e}" for p, e in fac.items()]
    if rem != 1:
        out.append(str(rem) if rem < 10**18 else f"C{len(str(rem))}")
    return " * ".join(out) if out else str(abs(n))


def targets(n: int) -> list[tuple[int, int, int]]:
    out = []
    for q in primes_upto(n):
        if q <= isqrt(n):
            continue
        a, r = divmod(n, q)
        if 1 <= r <= q - 2 and apery_mod(r, q) == 0:
            out.append((q, a, r))
    return out


def adjacent_pairs(n: int) -> list[tuple[int, int, int]]:
    out = []
    ts = targets(n)
    for a in sorted({a for _, a, _ in ts}):
        qs = sorted(q for q, aa, _ in ts if aa == a)
        out.extend((q, ell, a) for q, ell in zip(qs, qs[1:]))
    return out


def audit_formula(M: int, d: int) -> None:
    direct = shell(M, d)
    parts = first_cell_parts(M, d)
    assert parts[-1] == direct, (M, d)


def audit_split_difference(M: int, d: int, order: int) -> None:
    assert M // 2 < d and d + order <= M
    direct = delta(M, d, order)
    rows = [first_cell_parts(M, d + i) for i in range(order + 1)]
    # b_M is constant and disappears; the other three pieces split exactly.
    split = sum(delta_values([row[j] for row in rows]) for j in (1, 2, 3))
    assert direct == split


def audit_pair(n: int, q: int, ell: int, a: int) -> None:
    M = n - a
    d = q - 1
    L = ell - q
    assert a == 1, (n, q, ell, a)
    gd = G(M, d, L)
    high = delta(M, d, L + 1)
    g1 = G(M, d + 1, L)
    B = C(d + L + 1, L)
    assert gd - g1 == (-1) ** (L + 1) * B * high

    node_targets = []
    for i in range(L + 1):
        p = d + i + 1
        if p in primes_upto(d + L + 1) and p > L and apery_mod(n - p, p) == 0:
            node_targets.append(p)
    target_product = 1
    for p in node_targets:
        target_product *= p
    R = gcd(abs(gd), abs(high))
    assert R % target_product == 0

    width = min(5, M - (d + L))
    boundary = [G(M, d + s, L) for s in range(width + 1)]
    boundary += [delta(M, d + s, L + 1) for s in range(width + 1)]
    RG = 0
    for x in boundary:
        RG = gcd(RG, abs(x))

    print(
        "PAIR",
        f"n={n}", f"M={M}", f"q={q}", f"ell={ell}", f"d={d}", f"L={L}",
        f"targets={node_targets}",
    )
    print(
        "  digits(G,high,R,rectangle)=",
        (len(str(abs(gd))), len(str(abs(high))), len(str(R)), len(str(RG))),
    )
    print("  R=", R, "factor=", fmt_factor(R), "R/targets=", R // target_product)
    print("  rectangle=", RG, "factor=", fmt_factor(RG), "rectangle/targets=", RG // target_product)


def main() -> None:
    # Full small audit.
    for M in range(1, 35):
        assert shell(M, M) == endpoint_formula(M)
        for d in range(M // 2 + 1, M + 1):
            audit_formula(M, d)
            max_order = min(5, M - d)
            for order in range(1, max_order + 1):
                audit_split_difference(M, d, order)

    # Requested large M values: hostile points throughout the first cell.
    for M in (199, 271, 299, 320, 754):
        ds = sorted({M // 2 + 1, M // 2 + 2, (3 * M) // 4, M - 7, M - 1, M})
        for d in ds:
            if M // 2 < d <= M:
                audit_formula(M, d)
        assert shell(M, M) == endpoint_formula(M)
        for d in ds:
            if M // 2 < d and d + 4 <= M:
                audit_split_difference(M, d, 4)
        print("FORMULA_PASS", M, "sample_d=", ds)

    # Deterministic hostile random checks.
    rng = Random(5711)
    random_rows = []
    for _ in range(40):
        M = rng.randrange(20, 180)
        d = rng.randrange(M // 2 + 1, M + 1)
        audit_formula(M, d)
        order = rng.randrange(1, min(7, M - d + 1)) if d < M else 0
        if order:
            audit_split_difference(M, d, order)
        random_rows.append((M, d, order))
    print("RANDOM_PASS", random_rows)

    # Requested rows M=n-1.
    for n in (200, 272, 300, 321, 755):
        ts = targets(n)
        ps = adjacent_pairs(n)
        print("ROW", n, "targets=", ts, "pairs=", ps)
        for q, ell, a in ps:
            audit_pair(n, q, ell, a)

    print("ALL_Q5711_CHECKS_PASS")


if __name__ == "__main__":
    main()
