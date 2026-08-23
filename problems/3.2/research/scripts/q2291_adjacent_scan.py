#!/usr/bin/env python3
"""Exact dependency-free audit for Q2291.

Checks the adjacent identity, the factorial-cube gcd bound, and searches for
counterexamples to the claimed propagation from K(r,r-1) to all later
nonmate/noncentral spectral nodes.
"""
from __future__ import annotations

from math import comb, factorial, gcd, isqrt


def primes_upto(n: int) -> list[int]:
    sieve = bytearray(b"\x01") * (n + 1)
    sieve[:2] = b"\x00\x00"
    for q in range(2, isqrt(n) + 1):
        if sieve[q]:
            sieve[q*q:n+1:q] = b"\x00" * (((n - q*q)//q) + 1)
    return [q for q in range(2, n + 1) if sieve[q]]


def D(n: int, k: int) -> int:
    if k < 0 or k > n:
        return 0
    return comb(n, k) * comb(n + k, k)


def apery(n: int) -> int:
    return sum(D(n, k) ** 2 for k in range(n + 1))


def K(r: int, s: int) -> int:
    return sum(D(r, k) * D(s, k) for k in range(min(r, s) + 1))


def verify_adjacent(limit: int = 100) -> None:
    for r in range(1, limit + 1):
        br = apery(r)
        bm = apery(r - 1)
        kr = K(r, r - 1)
        assert 6 * kr == br + bm
        mr = r**3 * kr
        assert factorial(r) ** 3 % gcd(mr, br) == 0
    print(f"adjacent identity and gcd(M_r,b_r)|(r!)^3: r=1..{limit} VERIFIED")


def rows_mod_p(p: int) -> list[list[int]]:
    rows: list[list[int]] = []
    for n in range(p):
        row = [1]
        cur = 1
        for k in range(n):
            cur = (
                cur
                * (n - k)
                * (n + k + 1)
                * pow(k + 1, -2, p)
            ) % p
            row.append(cur)
        rows.append(row)
    return rows


def dot_rows(a: list[int], b: list[int], p: int) -> int:
    return sum(x * y for x, y in zip(a, b)) % p


def search_propagation(pmax: int = 1000) -> None:
    first_any = None
    first_s_noncentral = None
    first_both_noncentral = None
    first_strict_lower_half = None
    first_marked = None
    for p in primes_upto(pmax):
        if p < 5:
            continue
        rows = rows_mod_p(p)
        vals = [dot_rows(row, row, p) for row in rows]
        central = (p - 1) // 2
        for r in range(1, p - 1):
            if not (p > r + 1 and vals[r] == 0):
                continue
            kadj = dot_rows(rows[r], rows[r - 1], p)
            assert kadj != 0
            mate = p - 1 - r
            for s in range(r + 1, p - 1):
                if s == mate:
                    continue
                if dot_rows(rows[r], rows[s], p) != 0:
                    continue
                row = {
                    "p": p,
                    "r": r,
                    "s": s,
                    "mate": mate,
                    "central": central,
                    "b_r_mod_p": vals[r],
                    "b_s_mod_p": vals[s],
                    "K_adj_mod_p": kadj,
                    "K_rs_mod_p": 0,
                }
                if first_any is None:
                    first_any = row
                if s != central and first_s_noncentral is None:
                    first_s_noncentral = row
                if r != central and s != central and first_both_noncentral is None:
                    first_both_noncentral = row
                if r < s < central and first_strict_lower_half is None:
                    first_strict_lower_half = row
                if vals[s] == 0 and first_marked is None:
                    first_marked = row
    print("first later nonmate zero (central allowed):", first_any)
    print("first later nonmate zero (s noncentral):", first_s_noncentral)
    print("first later nonmate zero (r,s both noncentral):", first_both_noncentral)
    print("first later nonmate zero (strict lower half r<s<M):", first_strict_lower_half)
    print("first later marked zero b_s=0:", first_marked)
    if first_strict_lower_half is not None:
        p = first_strict_lower_half["p"]
        r = first_strict_lower_half["r"]
        s = first_strict_lower_half["s"]
        br = apery(r)
        bs = apery(s)
        ka = K(r, r - 1)
        ks = K(r, s)
        assert br % p == 0 and bs % p == first_strict_lower_half["b_s_mod_p"]
        assert ka % p == first_strict_lower_half["K_adj_mod_p"] and ks % p == 0
        print("exact strict-lower-half witness:", {
            "p": p,
            "r": r,
            "s": s,
            "b_r": br,
            "b_r_over_p": br // p,
            "b_s": bs,
            "K_adj": ka,
            "K_rs": ks,
            "K_rs_over_p": ks // p,
        })


if __name__ == "__main__":
    verify_adjacent(100)
    search_propagation(1000)
