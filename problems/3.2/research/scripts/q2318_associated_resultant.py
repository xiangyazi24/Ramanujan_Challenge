#!/usr/bin/env python3
"""Exact finite-field audit for Q2318.

The script uses the monic Racah recurrence

    P_(n+1) = (T+a_n) P_n - beta_n P_(n-1),
    a_n=(n^2+n+1)/2,
    beta_n=n^6/(4(4n^2-1)).

For h=s-r, the transfer/associated block used in Q2306 is
S_(h-1)^(r), not the literal S_h^(r):

    S_0^(r)=1,
    S_1^(r)=T+a_(r+1),
    S_(j+1)^(r)=(T+a_(r+j+1))S_j^(r)-beta_(r+j+1)S_(j-1)^(r).

All arithmetic is exact in F_p and uses only Python's standard library.
"""
from __future__ import annotations

import json
from itertools import combinations
from math import isqrt


# ---------- prime and polynomial utilities ----------

def primes_upto(n: int) -> list[int]:
    sieve = bytearray(b"\x01") * (n + 1)
    sieve[:2] = b"\x00\x00"
    for q in range(2, isqrt(n) + 1):
        if sieve[q]:
            sieve[q * q:n + 1:q] = b"\x00" * (((n - q * q) // q) + 1)
    return [q for q in range(2, n + 1) if sieve[q]]


def trim(a: list[int], p: int) -> list[int]:
    a = [x % p for x in a] or [0]
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return a


def add(a: list[int], b: list[int], p: int) -> list[int]:
    return trim([
        (a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)
        for i in range(max(len(a), len(b)))
    ], p)


def sub(a: list[int], b: list[int], p: int) -> list[int]:
    return trim([
        (a[i] if i < len(a) else 0) - (b[i] if i < len(b) else 0)
        for i in range(max(len(a), len(b)))
    ], p)


def scale(a: list[int], c: int, p: int) -> list[int]:
    return trim([c * x for x in a], p)


def mul(a: list[int], b: list[int], p: int) -> list[int]:
    a, b = trim(a, p), trim(b, p)
    if a == [0] or b == [0]:
        return [0]
    z = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            z[i + j] = (z[i + j] + x * y) % p
    return trim(z, p)


def divmod_poly(a: list[int], b: list[int], p: int) -> tuple[list[int], list[int]]:
    a, b = trim(a, p), trim(b, p)
    if b == [0]:
        raise ZeroDivisionError
    q = [0] * max(1, len(a) - len(b) + 1)
    ib = pow(b[-1], -1, p)
    while a != [0] and len(a) >= len(b):
        d = len(a) - len(b)
        c = a[-1] * ib % p
        q[d] = c
        for j, y in enumerate(b):
            a[d + j] = (a[d + j] - c * y) % p
        a = trim(a, p)
    return trim(q, p), a


def monic(a: list[int], p: int) -> list[int]:
    a = trim(a, p)
    if a == [0]:
        return [0]
    return scale(a, pow(a[-1], -1, p), p)


def gcd_poly(a: list[int], b: list[int], p: int) -> list[int]:
    a, b = trim(a, p), trim(b, p)
    while b != [0]:
        a, b = b, divmod_poly(a, b, p)[1]
    return monic(a, p)


def eval_poly(a: list[int], x: int, p: int) -> int:
    z = 0
    for c in reversed(a):
        z = (z * x + c) % p
    return z


def resultant(f: list[int], g: list[int], p: int) -> int:
    """Euclidean resultant, with Res(f,g) convention."""
    f, g = trim(f, p), trim(g, p)
    m, n = len(f) - 1, len(g) - 1
    if n == 0:
        return pow(g[0], m, p)
    if m < n:
        return ((-1 if (m * n) & 1 else 1) * resultant(g, f, p)) % p
    _, h = divmod_poly(f, g, p)
    if h == [0]:
        return 0
    k = len(h) - 1
    return (
        (-1 if (m * n) & 1 else 1)
        * pow(g[-1], m - k, p)
        * resultant(g, h, p)
    ) % p


# ---------- Apéry values and Racah polynomials ----------

def apery_half(p: int) -> list[int]:
    m = (p - 1) // 2
    if m == 0:
        return [1]
    vals = [1, 5 % p]
    for n in range(1, m):
        middle = (34 * n**3 + 51 * n**2 + 27 * n + 5) % p
        num = (middle * vals[n] - pow(n, 3, p) * vals[n - 1]) % p
        vals.append(num * pow(pow(n + 1, 3, p), -1, p) % p)
    return vals[:m + 1]


def aa(n: int, p: int) -> int:
    return (n * n + n + 1) * pow(2, -1, p) % p


def beta(n: int, p: int) -> int:
    den = 4 * (4 * n * n - 1) % p
    return pow(n, 6, p) * pow(den, -1, p) % p


def build_P(max_n: int, p: int) -> list[list[int]]:
    ps = [[1]]
    if max_n == 0:
        return ps
    ps.append([aa(0, p), 1])
    for n in range(1, max_n):
        nxt = sub(mul([aa(n, p), 1], ps[n], p), scale(ps[n - 1], beta(n, p), p), p)
        assert len(nxt) - 1 == n + 1 and nxt[-1] == 1
        ps.append(nxt)
    return ps


def S(base_r: int, j: int, p: int) -> list[int]:
    """S_j^(base_r), degree j."""
    if j < 0:
        return [0]
    if j == 0:
        return [1]
    prev = [1]
    cur = [aa(base_r + 1, p), 1]
    for k in range(1, j):
        n = base_r + k + 1
        prev, cur = cur, sub(mul([aa(n, p), 1], cur, p), scale(prev, beta(n, p), p), p)
    return cur


def transfer_block(r: int, s: int, p: int) -> list[int]:
    return S(r, s - r - 1, p)


def literal_C_h(r: int, s: int, p: int) -> list[int]:
    """Literal C_h under C_0=1,C_1=linear; degree h, not Q2306's block."""
    return S(r, s - r, p)


def row_D(n: int, p: int) -> list[int]:
    row = [1]
    cur = 1
    for k in range(n):
        cur = cur * (n - k) * (n + k + 1) % p
        cur = cur * pow((k + 1) ** 2, -1, p) % p
        row.append(cur)
    return row


def K_mod(r: int, s: int, p: int) -> int:
    a, b = row_D(r, p), row_D(s, p)
    return sum(x * y for x, y in zip(a, b)) % p


def lam(n: int, p: int) -> int:
    return n * (n + 1) % p


def beta_product(lo: int, hi: int, p: int) -> int:
    z = 1
    for n in range(lo, hi + 1):
        z = z * beta(n, p) % p
    return z


def verify_transfer(ps: list[list[int]], r: int, s: int, p: int) -> None:
    h = s - r
    c = S(r, h - 1, p)
    shifted = S(r + 1, h - 2, p)
    rhs = sub(mul(c, ps[r + 1], p), scale(mul(shifted, ps[r], p), beta(r + 1, p), p), p)
    assert rhs == ps[s]

    if h >= 3:
        lhs = sub(
            mul(S(r, h - 2, p), S(r + 1, h - 2, p), p),
            mul(S(r, h - 1, p), S(r + 1, h - 3, p), p),
            p,
        )
        assert lhs == [beta_product(r + 2, s - 1, p)]


# ---------- audits ----------

def marked_pair_rows(pmax: int = 1000) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for p in primes_upto(pmax):
        if p < 7:
            continue
        vals = apery_half(p)
        m = (p - 1) // 2
        zeros = [n for n in range(1, m) if vals[n] == 0]
        if len(zeros) < 2:
            continue
        ps = build_P(max(zeros), p)
        for r, s in combinations(zeros, 2):
            verify_transfer(ps, r, s, p)
            c = transfer_block(r, s, p)
            c_lit = literal_C_h(r, s, p)
            g = gcd_poly(ps[r], c, p)
            block_res = resultant(ps[r], c, p)
            adj_res = resultant(ps[r], ps[r + 1], p)
            full_res = resultant(ps[r], ps[s], p)
            assert full_res == adj_res * block_res % p
            kval = K_mod(r, s, p)
            lr, ls = lam(r, p), lam(s, p)
            assert vals[r] == vals[s] == 0
            # Monic P-values differ from K by unit leading scalars, so only zero/nonzero is compared.
            assert (eval_poly(ps[r], ls, p) == 0) == (kval == 0)
            assert (eval_poly(c, ls, p) == 0) == (kval == 0)
            assert (eval_poly(c, lr, p) == 0) == (kval == 0)
            rows.append({
                "p": p,
                "r": r,
                "s": s,
                "h": s - r,
                "lambda_r": lr,
                "lambda_s": ls,
                "K_mod_p": kval,
                "associated_degree": len(c) - 1,
                "associated_resultant": block_res,
                "literal_C_h_degree": len(c_lit) - 1,
                "literal_C_h_resultant": resultant(ps[r], c_lit, p),
                "adjacent_resultant": adj_res,
                "full_monic_resultant": full_res,
                "gcd_degree": len(g) - 1,
                "C_at_lambda_r": eval_poly(c, lr, p),
                "C_at_lambda_s": eval_poly(c, ls, p),
            })
    return rows


def first_actual_block_collision(pmax: int = 251) -> dict[str, object] | None:
    """Find a stable Racah triple with Res(P_r,S_(h-1)^r)=0."""
    for p in primes_upto(pmax):
        if p < 7:
            continue
        m = (p - 1) // 2
        ps = build_P(m - 1, p)
        for s in range(3, m):
            for r in range(1, s - 1):
                c = transfer_block(r, s, p)
                g = gcd_poly(ps[r], c, p)
                if len(g) > 1:
                    roots = [x for x in range(p) if eval_poly(g, x, p) == 0]
                    verify_transfer(ps, r, s, p)
                    return {
                        "p": p,
                        "r": r,
                        "s": s,
                        "h": s - r,
                        "resultant": 0,
                        "gcd": g,
                        "base_field_common_roots": roots,
                        "all_beta_units": all(beta(n, p) != 0 for n in range(1, s)),
                    }
    return None


def first_bs_only_collision(pmax: int = 1000) -> dict[str, object] | None:
    """Find b_s=0, b_r!=0 but associated resultant zero, if present."""
    for p in primes_upto(pmax):
        if p < 7:
            continue
        vals = apery_half(p)
        m = (p - 1) // 2
        zero_s = [s for s in range(3, m) if vals[s] == 0]
        if not zero_s:
            continue
        ps = build_P(max(zero_s), p)
        for s in zero_s:
            for r in range(1, s - 1):
                if vals[r] == 0:
                    continue
                c = transfer_block(r, s, p)
                g = gcd_poly(ps[r], c, p)
                if len(g) > 1:
                    roots = [x for x in range(p) if eval_poly(g, x, p) == 0]
                    return {
                        "p": p,
                        "r": r,
                        "s": s,
                        "b_r": vals[r],
                        "b_s": vals[s],
                        "resultant": 0,
                        "gcd": g,
                        "base_field_common_roots": roots,
                        "K_mod_p": K_mod(r, s, p),
                    }
    return None


def main() -> None:
    marked = marked_pair_rows(1000)
    collision = first_actual_block_collision(251)
    bs_only = first_bs_only_collision(1000)

    fixture = next(row for row in marked if (row["p"], row["r"], row["s"]) == (181, 19, 47))
    assert fixture["associated_resultant"] == 19
    assert fixture["associated_degree"] == 27

    summary = {
        "marked_scan_bound": 1000,
        "marked_pair_count": len(marked),
        "all_marked_associated_resultants_nonzero": all(row["associated_resultant"] != 0 for row in marked),
        "all_marked_K_nonzero": all(row["K_mod_p"] != 0 for row in marked),
        "marked_pairs": marked,
        "first_stable_associated_collision": collision,
        "first_b_s_zero_only_collision": bs_only,
        "indexing_fixture": fixture,
    }
    print("ANSWER Q2318 7cd8ebec")
    print(json.dumps(summary, sort_keys=True, separators=(",", ":")))


if __name__ == "__main__":
    main()
