#!/usr/bin/env python3
"""Exact standard-library audit for Q5727's local target-alias no-go.

Checks:
* Newton evaluation at -1 is the matched-shell projector modulo every prime
  node in a short stencil;
* the exact two-endpoint coefficient module;
* the cross-weighted attempt that is coefficient-forced at one endpoint and
  alias-forced at the other;
* rank <= 1 modulo each preserved endpoint and the resulting minor factors.
"""
from itertools import combinations
from math import comb, gcd


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    d = 2
    while d * d <= n:
        if n % d == 0:
            return False
        d += 1
    return True


def valuation(n: int, p: int) -> int:
    n = abs(n)
    if n == 0:
        return 10**9
    e = 0
    while n % p == 0:
        n //= p
        e += 1
    return e


def newton_row(a: int, b: int) -> list[int]:
    """Coefficients of interpolation at -1 on nodes a,...,b."""
    m = b - a
    return [(-1) ** i * comb(a + i, i) * comb(b + 1, m - i)
            for i in range(m + 1)]


def weighted_row(a: int, b: int, center: int) -> list[int]:
    return [(a + i - center) * c for i, c in enumerate(newton_row(a, b))]


def rank_mod(A: list[list[int]], p: int) -> int:
    if not A:
        return 0
    B = [[x % p for x in row] for row in A]
    m, n = len(B), len(B[0])
    r = 0
    for c in range(n):
        pivot = next((i for i in range(r, m) if B[i][c]), None)
        if pivot is None:
            continue
        B[r], B[pivot] = B[pivot], B[r]
        z = pow(B[r][c], -1, p)
        B[r] = [(z * x) % p for x in B[r]]
        for i in range(m):
            if i != r and B[i][c]:
                z = B[i][c]
                B[i] = [(x - z * y) % p for x, y in zip(B[i], B[r])]
        r += 1
    return r


def det2(row1: list[int], row2: list[int], c1: int, c2: int) -> int:
    return row1[c1] * row2[c2] - row1[c2] * row2[c1]


def audit_interval(q: int, ell: int, left: int, right: int) -> None:
    assert is_prime(q) and is_prime(ell) and q < ell
    assert left <= q - 1 < ell - 1 <= right
    assert right - left < q
    # In the P3.2 short-margin range this also excludes second multiples.
    assert right + 1 < 2 * q

    row = newton_row(left, right)
    iq, ie = q - 1 - left, ell - 1 - left

    # Exact matched-node projector modulo each endpoint.
    for p, ip in ((q, iq), (ell, ie)):
        assert row[ip] % p == 1
        assert all(c % p == 0 for j, c in enumerate(row) if j != ip)

    # Exact pair-preserving coordinate module:
    # ell at q-coordinate, q at ell-coordinate, q*ell elsewhere.
    assert row[iq] % ell == 0 and row[iq] % q != 0
    assert row[ie] % q == 0 and row[ie] % ell != 0
    assert all(c % (q * ell) == 0
               for j, c in enumerate(row) if j not in (iq, ie))

    # The minimal CRT cross attempt.
    rq = weighted_row(left, right, q - 1)
    re = weighted_row(left, right, ell - 1)
    L = ell - q
    assert all(c % q == 0 for c in rq)       # coefficient-forced at q
    assert re[iq] % q == (-L) % q            # alias at q
    assert rq[ie] % ell == L % ell           # alias at ell
    assert all(c % ell == 0 for c in re)     # coefficient-forced at ell

    A = [rq, re]
    assert rank_mod(A, q) == 1
    assert rank_mod(A, ell) == 1
    D = det2(rq, re, iq, ie)
    assert valuation(D, q) == 1
    assert valuation(D, ell) == 1

    # Any two pair-preserving Newton rows have every 2x2 minor divisible q*ell.
    rows = []
    for s in range(0, min(5, q - (right - left))):
        a = left - s
        b = right + s
        if b - a >= q or a < 0 or b + 1 >= 2 * q:
            continue
        rows.append(newton_row(a, b))
    # Embed rows into their common coordinate interval.
    lo = min(left - i for i in range(len(rows))) if rows else left
    hi = max(right + i for i in range(len(rows))) if rows else right
    embedded = []
    for i, r in enumerate(rows):
        a = left - i
        v = [0] * (hi - lo + 1)
        for j, c in enumerate(r):
            v[a + j - lo] = c
        embedded.append(v)
    for x, y in combinations(embedded, 2):
        assert rank_mod([x, y], q) <= 1
        assert rank_mod([x, y], ell) <= 1
        for c1, c2 in combinations(range(len(x)), 2):
            assert det2(x, y, c1, c2) % (q * ell) == 0


def main() -> None:
    cases = [
        (179, 193, 174, 198),
        (193, 211, 188, 216),
        (101, 113, 97, 117),
        (137, 149, 132, 154),
    ]
    for case in cases:
        audit_interval(*case)
    print("Q5727 local coefficient-module audit: PASS", len(cases))


if __name__ == "__main__":
    main()
