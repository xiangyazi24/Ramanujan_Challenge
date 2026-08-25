#!/usr/bin/env python3
from __future__ import annotations

import math
import struct
from collections import defaultdict
from pathlib import Path

NMAX = 300_000
PAIR = struct.Struct("<II")
DATA = Path("problems/3.2/data_zp_pairs.bin")


def sieve(n: int) -> bytearray:
    s = bytearray(b"\x01") * (n + 1)
    s[:2] = b"\x00\x00"
    for p in range(2, math.isqrt(n) + 1):
        if s[p]:
            s[p*p:n+1:p] = b"\x00" * (((n - p*p) // p) + 1)
    return s


def P(m: int) -> int:
    return 34*m**3 + 51*m**2 + 27*m + 5


def apery_companion_mod(q: int, mmax: int):
    assert mmax + 1 < q
    b = [0] * (mmax + 2)
    a = [0] * (mmax + 2)
    b[0], a[0] = 1, 0
    if mmax + 1 >= 1:
        b[1], a[1] = 5 % q, 6 % q
    for m in range(1, mmax + 1):
        inv = pow(m + 1, -3, q)
        b[m+1] = (P(m) * b[m] - m**3 * b[m-1]) * inv % q
        a[m+1] = (P(m) * a[m] - m**3 * a[m-1]) * inv % q
    for m in range(0, mmax + 1):
        lhs = (a[m+1] * b[m] - a[m] * b[m+1]) % q
        rhs = 6 * pow(m + 1, -3, q) % q
        assert lhs == rhs, (q, m, lhs, rhs)
    return a, b


def forward_continuant_mod(s: int, t: int, q: int) -> int:
    """Delta_t times solution u_t with u_s=0,u_(s+1)=1."""
    assert s < t < q
    if t == s + 1:
        return 1
    km1, k = 0, 1
    for m in range(s + 1, t):
        km1, k = k, (P(m) * k - pow(m, 6, q) * km1) % q
    return k


def read_pairs(path: Path):
    raw = path.read_bytes()
    assert len(raw) % PAIR.size == 0
    out = []
    prev = None
    for off in range(0, len(raw), PAIR.size):
        pair = PAIR.unpack_from(raw, off)
        if prev is not None:
            assert prev < pair
        prev = pair
        out.append(pair)
    return out


def rho(n: int) -> int:
    for r in range(1, 5):
        if (r - (1 - n)) % 4 == 0:
            return r
    raise AssertionError


def main():
    prime = sieve(2_000_000)
    pairs = read_pairs(DATA)
    target = set()
    for p, r in pairs:
        assert prime[p] and 0 <= r < p
        if p <= NMAX:
            j = min(r, p - 1 - r)
            target.add((p, j))

    events = defaultdict(dict)
    for p, j in sorted(target):
        nd = p + j
        if nd <= NMAX and 3*j < nd - 1:
            assert p == nd - j
            events[(nd, j)]["D"] = p
        nr = 2*p - 1 - j
        if nr <= NMAX and 3*j < nr - 1:
            assert p == (nr + 1 + j) // 2
            events[(nr, j)]["R"] = p

    repeated = []
    by_n = defaultdict(list)
    for (n, j), e in sorted(events.items()):
        if set(e) != {"D", "R"}:
            continue
        pd, pr = e["D"], e["R"]
        if pd == pr:
            continue
        assert 1 <= j and 3*j < n - 1
        assert pd == n - j
        assert 2*pr == n + 1 + j
        assert prime[pd] and prime[pr]
        assert n + j == 1 or (n + j) % 4 == 1
        assert pd > (2*n + 1) / 3
        assert n / 2 < pr < (2*n + 1) / 3
        M = pd * pr
        assert 2*M == n*(n+1) - j*(j+1)
        row = {"n": n, "j": j, "pD": pd, "pR": pr, "M": M}
        repeated.append(row)
        by_n[n].append(row)

    # Exact reuse, spacing, and quadratic-node Vandermonde audits.
    all_assigned = []
    spacing_checks = 0
    for n, rows in by_n.items():
        rows.sort(key=lambda z: z["j"])
        rr = rho(n)
        for idx, row in enumerate(rows):
            assert row["j"] % 4 == rr % 4
            if idx:
                assert row["j"] - rows[idx-1]["j"] >= 4
            all_assigned += [row["pD"], row["pR"]]
        for x in range(len(rows)):
            for y in range(x + 1, len(rows)):
                a, b = rows[x], rows[y]
                assert math.gcd(a["M"], b["M"]) == 1
                diff = b["j"]*(b["j"]+1) - a["j"]*(a["j"]+1)
                assert diff == 2*(a["M"] - b["M"])
                assert math.gcd(a["M"], diff) == 1
                assert math.gcd(b["M"], diff) == 1
                spacing_checks += 1
    assert len(all_assigned) == len(set(all_assigned))

    # Actual recurrence, companion, cross-hit, and short-continuant audit.
    modular_cache = {}
    cross_hits = []
    ordered_tests = 0
    continuant_tests = 0
    for n, rows in sorted(by_n.items()):
        js = [r["j"] for r in rows]
        mmax = max(js) + 1
        for row in rows:
            j = row["j"]
            for qkind in ("pD", "pR"):
                q = row[qkind]
                if (q, mmax) not in modular_cache:
                    modular_cache[(q, mmax)] = apery_companion_mod(q, mmax)
                a, b = modular_cache[(q, mmax)]
                assert b[j] == 0
                assert a[j] != 0 and b[j+1] != 0
                for k in js:
                    det = (a[j]*b[k] - b[j]*a[k]) % q
                    hit = b[k] == 0
                    assert (det == 0) == hit
                    assert hit == ((q, k) in target)
                    ordered_tests += 1
                    if k != j and hit:
                        cross_hits.append((n, j, q, k, qkind))
                    if j < k:
                        K = forward_continuant_mod(j, k, q)
                        assert (K == 0) == hit
                        continuant_tests += 1

    # Sieve local-root count audit for every repeated n and odd ell <= 1000.
    root_count_tests = 0
    small_prime = sieve(1000)
    for n in by_n:
        r0 = rho(n)
        cD = n - r0
        cR2 = n + 1 + r0  # 2*L_R = cR2 + 4t
        for ell in range(3, 1001, 2):
            if not small_prime[ell]:
                continue
            roots = set()
            for t in range(ell):
                if (cD - 4*t) % ell == 0 or (cR2 + 4*t) % ell == 0:
                    roots.add(t)
            expected = 1 if (2*n + 1) % ell == 0 else 2
            assert len(roots) == expected, (n, ell, len(roots), expected)
            root_count_tests += 1

    print("Q3448 repeated-node audit")
    print(f"binary_records={len(pairs)}")
    print(f"folded_target_pairs_p_le_{NMAX}={len(target)}")
    print(f"repeated_nodes={len(repeated)}")
    print(f"n_with_repeated_nodes={len(by_n)}")
    print(f"max_repeated_nodes_at_one_n={max((len(v) for v in by_n.values()), default=0)}")
    print(f"assigned_prime_reuse=0")
    print(f"quadratic_spacing_pairs_checked={spacing_checks}")
    print(f"ordered_cross_hit_tests={ordered_tests}")
    print(f"forward_continuant_tests={continuant_tests}")
    print(f"sieve_root_count_tests={root_count_tests}")
    print(f"cross_hits={cross_hits}")
    print("rows:")
    for n in sorted(by_n):
        print(n, [(r["j"], r["pD"], r["pR"]) for r in by_n[n]])


if __name__ == "__main__":
    main()
