#!/usr/bin/env python3
"""Dependency-free exact audit for Q5729.

For 20 <= M <= LIMIT and H=floor(M^(1/3)), write
    M = 2*P + eps, eps in {0,1},
    s_j = eps + 2*j, p_j = P-j, 0 <= j < H.
The moving diagonal is p_j prime and p_j | b_{s_j}.

The script verifies, for every candidate node in the range:
  * Apéry--Lucas synchronization b_M == 73*b_s (mod p);
  * reflection b_{p-1-s} == b_s (mod p);
  * divisibility of the parity Newton carrier Q_H(P);
  * positivity/lower bound Delta^k b_eps >= b_eps*72^k;
  * the exact two-node resultant for every pair of hits;
  * the elementary size obstruction b_s >= p at every hit.

Only Python's standard library and exact integer arithmetic are used.
"""

from __future__ import annotations

from math import comb, log

LIMIT = 5000


def icbrt(n: int) -> int:
    lo, hi = 0, 1
    while hi**3 <= n:
        hi *= 2
    while lo + 1 < hi:
        mid = (lo + hi) // 2
        if mid**3 <= n:
            lo = mid
        else:
            hi = mid
    return lo


def primes_upto(n: int) -> list[bool]:
    isprime = bytearray(b"\x01") * (n + 1)
    if n >= 0:
        isprime[0] = 0
    if n >= 1:
        isprime[1] = 0
    p = 2
    while p * p <= n:
        if isprime[p]:
            isprime[p * p : n + 1 : p] = b"\x00" * (((n - p * p) // p) + 1)
        p += 1
    return [bool(x) for x in isprime]


def apery_numbers(n: int) -> list[int]:
    if n == 0:
        return [1]
    b = [1, 5]
    for k in range(1, n):
        num = (34 * k**3 + 51 * k**2 + 27 * k + 5) * b[k] - k**3 * b[k - 1]
        den = (k + 1) ** 3
        q, r = divmod(num, den)
        assert r == 0
        b.append(q)
    return b[: n + 1]


def forward_differences(values: list[int]) -> list[int]:
    out = []
    level = values[:]
    while level:
        out.append(level[0])
        level = [level[i + 1] - level[i] for i in range(len(level) - 1)]
    return out


def parity_newton_carrier(P: int, values: list[int]) -> tuple[int, list[int]]:
    diffs = forward_differences(values)
    carrier = sum(comb(P, k) * diffs[k] for k in range(len(diffs)))
    return carrier, diffs


def main() -> None:
    b = apery_numbers(LIMIT)
    isprime = primes_upto(LIMIT // 2 + 10)

    node_checks = 0
    reflection_checks = 0
    lucas_checks = 0
    carrier_checks = 0
    pair_checks = 0
    positivity_checks = 0
    total_hits = 0
    hit_rows = 0
    max_hits = 0
    max_weight_over_H = 0.0
    max_carrier_log_over_HlogM = 0.0
    rows_with_hits: list[tuple[int, list[tuple[int, int, int]]]] = []

    for M in range(20, LIMIT + 1):
        H = icbrt(M)
        eps = M & 1
        P = (M - eps) // 2
        values = [b[eps + 2 * j] for j in range(H)]
        carrier, diffs = parity_newton_carrier(P, values)

        for k, dk in enumerate(diffs):
            assert dk >= b[eps] * 72**k, (M, k, dk)
            positivity_checks += 1

        lower = comb(P, H - 1) * b[eps] * 72 ** (H - 1)
        assert carrier >= lower
        if H > 1 and M > 1:
            ratio = log(carrier) / (H * log(M))
            max_carrier_log_over_HlogM = max(max_carrier_log_over_HlogM, ratio)

        hits: list[tuple[int, int, int]] = []
        for j in range(H):
            s = eps + 2 * j
            p = P - j
            if p <= s or p >= len(isprime) or not isprime[p]:
                continue
            node_checks += 1

            assert b[M] % p == 73 * (b[s] % p) % p, (M, s, p)
            lucas_checks += 1

            t = p - 1 - s
            assert 0 <= t < p
            assert b[t] % p == b[s] % p, (M, s, p, t)
            reflection_checks += 1

            hit = b[s] % p == 0
            if p != 73:
                assert hit == (b[M] % p == 0), (M, s, p)
            if hit:
                assert b[s] >= p
                assert carrier % p == 0, (M, s, p)
                carrier_checks += 1
                hits.append((j, s, p))

        for x in range(len(hits)):
            j, s, p = hits[x]
            for y in range(x + 1, len(hits)):
                k, t, q = hits[y]
                E = (P - k) * values[j] - (P - j) * values[k]
                assert E % (p * q) == 0, (M, hits[x], hits[y])
                pair_checks += 1

        if hits:
            hit_rows += 1
            total_hits += len(hits)
            max_hits = max(max_hits, len(hits))
            weight = sum(log(p) for _, _, p in hits)
            max_weight_over_H = max(max_weight_over_H, weight / H)
            rows_with_hits.append((M, hits))

    print("Q5729_DIAGONAL_AUDIT=PASS")
    print("LIMIT", LIMIT)
    print("NODE_CHECKS", node_checks)
    print("LUCAS_CHECKS", lucas_checks)
    print("REFLECTION_CHECKS", reflection_checks)
    print("POSITIVITY_CHECKS", positivity_checks)
    print("CARRIER_HIT_CHECKS", carrier_checks)
    print("PAIR_RESULTANT_CHECKS", pair_checks)
    print("HIT_ROWS", hit_rows, "TOTAL_HITS", total_hits, "MAX_HITS_IN_ONE_ROW", max_hits)
    print("MAX_WEIGHT_OVER_H", max_weight_over_H)
    print("MAX_LOG_CARRIER_OVER_H_LOG_M", max_carrier_log_over_HlogM)
    print("ROWS_WITH_HITS_BEGIN")
    for row in rows_with_hits:
        print(row)
    print("ROWS_WITH_HITS_END")


if __name__ == "__main__":
    main()
