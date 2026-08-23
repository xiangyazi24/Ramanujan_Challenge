#!/usr/bin/env python3
from __future__ import annotations

from itertools import combinations, product
from math import comb, isqrt
import json


def primes_upto(n: int) -> list[int]:
    z = bytearray(b"\x01") * (n + 1)
    z[:2] = b"\x00\x00"
    for q in range(2, isqrt(n) + 1):
        if z[q]:
            z[q*q:n+1:q] = b"\x00" * (((n - q*q)//q) + 1)
    return [q for q in range(2, n + 1) if z[q]]


def apery_P(n: int) -> int:
    return 34*n**3 + 51*n**2 + 27*n + 5


def apery_mod(p: int, last: int) -> list[int]:
    if last == 0:
        return [1]
    b = [1, 5 % p]
    for n in range(1, last):
        den = pow(n + 1, 3, p)
        assert den
        b.append((apery_P(n)*b[n] - n**3*b[n-1]) * pow(den, -1, p) % p)
    return b


def trim(a: list[int]) -> list[int]:
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return a


def add(a: list[int], b: list[int], p: int) -> list[int]:
    n = max(len(a), len(b))
    c = [0] * n
    for i in range(n):
        c[i] = ((a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)) % p
    return trim(c)


def sub(a: list[int], b: list[int], p: int) -> list[int]:
    n = max(len(a), len(b))
    c = [0] * n
    for i in range(n):
        c[i] = ((a[i] if i < len(a) else 0) - (b[i] if i < len(b) else 0)) % p
    return trim(c)


def scale(a: list[int], c: int, p: int) -> list[int]:
    return trim([(c*x) % p for x in a])


def mul(a: list[int], b: list[int], p: int) -> list[int]:
    c = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            c[i+j] = (c[i+j] + x*y) % p
    return trim(c)


def divmod_poly(a: list[int], b: list[int], p: int) -> tuple[list[int], list[int]]:
    a = trim(a[:])
    b = trim(b[:])
    if b == [0]:
        raise ZeroDivisionError
    if len(a) < len(b):
        return [0], a
    q = [0] * (len(a) - len(b) + 1)
    ib = pow(b[-1], -1, p)
    while len(a) >= len(b) and a != [0]:
        sh = len(a) - len(b)
        c = a[-1] * ib % p
        q[sh] = c
        for j, y in enumerate(b):
            a[sh+j] = (a[sh+j] - c*y) % p
        trim(a)
    return trim(q), trim(a)


def gcd_poly(a: list[int], b: list[int], p: int) -> list[int]:
    a, b = trim(a[:]), trim(b[:])
    while b != [0]:
        _, r = divmod_poly(a, b, p)
        a, b = b, r
    ia = pow(a[-1], -1, p)
    return scale(a, ia, p)


def resultant(f: list[int], g: list[int], p: int) -> int:
    f, g = trim(f[:]), trim(g[:])
    if f == [0] or g == [0]:
        return 0
    m, n = len(f) - 1, len(g) - 1
    if n == 0:
        return pow(g[0], m, p)
    if m == 0:
        return pow(f[0], n, p)
    if m < n:
        return ((-1 if (m*n) & 1 else 1) * resultant(g, f, p)) % p
    _, r = divmod_poly(f, g, p)
    if r == [0]:
        return 0
    k = len(r) - 1
    return ((-1 if (m*n) & 1 else 1)
            * pow(g[-1], m-k, p)
            * resultant(g, r, p)) % p


def evaluate(f: list[int], x: int, p: int) -> int:
    y = 0
    for c in reversed(f):
        y = (y*x + c) % p
    return y


def alpha(n: int, p: int) -> int:
    return (n*n + n + 1) * pow(2, -1, p) % p


def beta(n: int, p: int) -> int:
    den = 4 * (4*n*n - 1) % p
    assert den, (p, n)
    return pow(n, 6, p) * pow(den, -1, p) % p


def monic_racah_rows(p: int, M: int) -> list[list[int]]:
    P = [[1], [pow(2, -1, p), 1]]
    for n in range(1, M):
        nxt = sub(mul([alpha(n, p), 1], P[n], p), scale(P[n-1], beta(n, p), p), p)
        assert len(nxt) == n + 2 and nxt[-1] == 1
        P.append(nxt)
    return P


def associated_block(p: int, r: int, h: int) -> list[int]:
    # C_h = S_{h-1}^{(r)} in Q2306 notation.
    assert h >= 1
    if h == 1:
        return [1]
    prev = [1]  # S_0
    cur = [alpha(r + 1, p), 1]  # S_1
    if h == 2:
        return cur
    for j in range(1, h - 1):
        n = r + j + 1
        nxt = sub(mul([alpha(n, p), 1], cur, p), scale(prev, beta(n, p), p), p)
        prev, cur = cur, nxt
    assert len(cur) == h
    return cur


def D(n: int, k: int) -> int:
    return comb(n, k) * comb(n + k, k)


def K(r: int, s: int, p: int) -> int:
    return sum(D(r, k)*D(s, k) for k in range(min(r, s)+1)) % p


def bal(x: int, p: int) -> int:
    x %= p
    return x if x <= p//2 else x-p


def inv(x: int, p: int) -> int:
    return pow(x % p, -1, p)


def canonical_pair(pair: tuple[int, int], p: int) -> tuple[int, int]:
    a, b = pair
    a = min(a, p-1-a)
    b = min(b, p-1-b)
    return tuple(sorted((a, b)))


def main() -> None:
    rows = []
    full_nonmate_count = 0
    relevant_full_sets = {}

    for p in primes_upto(500):
        if p < 5:
            continue
        M = (p-1)//2
        bfull = apery_mod(p, p-2)
        fullzeros = [n for n in range(1, p-1) if bfull[n] == 0]
        fullpairs = [(r, s) for r, s in combinations(fullzeros, 2) if r+s+1 != p]
        full_nonmate_count += len(fullpairs)

        lowerzeros = [n for n in fullzeros if n <= M]
        if len(lowerzeros) < 2:
            continue
        Pr = monic_racah_rows(p, M)
        cpairs = list(combinations(lowerzeros, 2))
        orbit_counts = {}
        for pair in fullpairs:
            c = canonical_pair(pair, p)
            orbit_counts[c] = orbit_counts.get(c, 0) + 1
        relevant_full_sets[p] = fullzeros

        for r, s in cpairs:
            h = s-r
            C = associated_block(p, r, h)
            R = resultant(Pr[r], C, p)
            g = gcd_poly(Pr[r], C, p)
            kval = K(r, s, p)
            v = (r+s+1) % p
            delta = h*v % p
            assert bfull[r] == bfull[s] == 0
            assert evaluate(Pr[r], s*(s+1) % p, p) != 0
            assert len(g) == 1 and R != 0
            # Q2306 calibration.
            if (p, r, s) == (181, 19, 47):
                assert R == 19
            # Verify the block factorization modulo P_r.
            rem = divmod_poly(sub(Pr[s], mul(Pr[r+1], C, p), p), Pr[r], p)[1]
            assert rem == [0]
            adj = resultant(Pr[r], Pr[r+1], p)
            whole = resultant(Pr[r], Pr[s], p)
            assert whole == adj*R % p
            rows.append({
                "p": p, "r": r, "s": s, "h": h, "v": v,
                "delta": delta, "rs": r*s % p, "K": kval,
                "R": R, "R_bal": bal(R, p),
                "R_over_r": R*inv(r, p) % p,
                "R_over_s": R*inv(s, p) % p,
                "R_over_h": R*inv(h, p) % p,
                "R_over_v": R*inv(v, p) % p,
                "R_over_delta": R*inv(delta, p) % p,
                "R_over_K": R*inv(kval, p) % p,
                "orbit_count": orbit_counts.get((r, s), 0),
                "degC": len(C)-1,
            })

    # Test a useful bank of candidate formulae.
    def candidate_values(row: dict) -> dict[str, int]:
        p = row["p"]
        r, s, h, v, d, k = (row[x] for x in ("r", "s", "h", "v", "delta", "K"))
        vals = {
            "1": 1, "-1": -1, "r": r, "-r": -r,
            "s": s, "-s": -s, "h": h, "-h": -h,
            "v=r+s+1": v, "-v": -v, "r*s": r*s,
            "delta=h*v": d, "K": k, "-K": -k,
            "K^2": k*k, "r/h": r*inv(h,p),
            "h/r": h*inv(r,p), "r/v": r*inv(v,p),
            "v/r": v*inv(r,p), "r*s/(h*v)": r*s*inv(d,p),
        }
        return {name: value % p for name, value in vals.items()}

    exact_matches = []
    for name in candidate_values(rows[0]):
        if all(candidate_values(row)[name] == row["R"] for row in rows):
            exact_matches.append(name)

    # Search monomials +/- r^a s^b h^c v^d K^e, exponents -2..2.
    monomial_matches = []
    bases = ("r", "s", "h", "v", "K")
    for exps in product(range(-2, 3), repeat=len(bases)):
        if all(e == 0 for e in exps):
            continue
        for sign in (1, -1):
            ok = True
            for row in rows:
                p = row["p"]
                val = sign % p
                for base, e in zip(bases, exps):
                    x = row[base] % p
                    val = val * (pow(x, e, p) if e >= 0 else pow(inv(x,p), -e, p)) % p
                if val != row["R"]:
                    ok = False
                    break
            if ok:
                monomial_matches.append({"sign": sign, "exponents": dict(zip(bases, exps))})

    print("Q2328_JSON_BEGIN")
    print(json.dumps({
        "canonical_pair_count": len(rows),
        "full_nonmate_unordered_pair_count": full_nonmate_count,
        "relevant_full_zero_sets": relevant_full_sets,
        "rows": rows,
        "simple_exact_matches": exact_matches,
        "monomial_matches_exponents_-2_to_2": monomial_matches,
    }, sort_keys=True))
    print("Q2328_JSON_END")


if __name__ == "__main__":
    main()
