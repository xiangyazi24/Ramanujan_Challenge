#!/usr/bin/env python3
"""Q5714: exact pair-preserving Newton-margin audit.

Only the Python standard library is used.  The script:
* verifies the full left/right Pascal calculus for H_{s,t}=G_{d-s,L+s+t};
* verifies rational rank S+T+1 and q/ell-local primitive rank;
* computes actual-shell rectangle gcds for n=200,272,300,321,755;
* scans target-blind square-margin rules on all short prime pairs in those rows.
"""
from __future__ import annotations

from functools import lru_cache
from math import comb, gcd, isqrt, log
from random import Random


def C(n: int, k: int) -> int:
    return comb(n, k) if 0 <= k <= n else 0


def primes_upto(n: int) -> list[int]:
    if n < 2:
        return []
    s = bytearray(b"\x01") * (n + 1)
    s[0:2] = b"\x00\x00"
    for p in range(2, isqrt(n) + 1):
        if s[p]:
            s[p * p : n + 1 : p] = b"\x00" * (((n - p * p) // p) + 1)
    return [p for p in range(2, n + 1) if s[p]]


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


@lru_cache(None)
def shell(M: int, d: int) -> int:
    """C_M(d), by the exact one-fold cyclic-binomial formula."""
    assert d >= 1
    a = M // d
    out = 0
    for t in range(M + 1):
        X = sum(C(M, M - t + d * u) for u in range(-a, a + 1))
        Z = sum(C(2 * M - t, M - t + d * v) for v in range(-a, a + 1))
        out += C(M, t) * X * Z * Z
    return out


def delta(M: int, d: int, k: int) -> int:
    return sum((-1) ** (k - i) * C(k, i) * shell(M, d + i) for i in range(k + 1))


def weight(d: int, L: int, i: int) -> int:
    return (-1) ** i * C(d + i, i) * C(d + L + 1, L - i)


def G(M: int, d: int, L: int) -> int:
    return sum(weight(d, L, i) * shell(M, d + i) for i in range(L + 1))


def Hmargin(M: int, d: int, L: int, s: int, t: int) -> int:
    return G(M, d - s, L + s + t)


def valuation(n: int, p: int) -> int:
    n = abs(n)
    if n == 0:
        return 10**9
    e = 0
    while n % p == 0:
        n //= p
        e += 1
    return e


def prod(xs) -> int:
    z = 1
    for x in xs:
        z *= x
    return z


# Deterministic primality for 64-bit integers, used only to certify displayed factors.
def is_prime64(n: int) -> bool:
    if n < 2:
        return False
    small = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)
    for p in small:
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in (2, 325, 9375, 28178, 450775, 9780504, 1795265022):
        if a % n == 0:
            continue
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def pollard_rho(n: int) -> int | None:
    if n % 2 == 0:
        return 2
    if n % 3 == 0:
        return 3
    rng = Random(n ^ 0x9E3779B97F4A7C15)
    for _ in range(32):
        c = rng.randrange(1, n - 1)
        x = rng.randrange(2, n - 1)
        y = x
        g = 1
        for _ in range(250000):
            x = (x * x + c) % n
            y = (y * y + c) % n
            y = (y * y + c) % n
            g = gcd(abs(x - y), n)
            if g == 1:
                continue
            if g != n:
                return g
            break
    return None


def factor_rigorous(n: int) -> tuple[list[tuple[int, int]], list[int]]:
    """Return certified <2^64 prime factors and exact unresolved cofactors."""
    n = abs(n)
    if n in (0, 1):
        return [], [] if n == 1 else [0]
    fs: list[int] = []
    for p in primes_upto(100000):
        while n % p == 0:
            fs.append(p)
            n //= p
        if p * p > n:
            break
    stack = [] if n == 1 else [n]
    unresolved: list[int] = []
    while stack:
        m = stack.pop()
        if m == 1:
            continue
        if m < (1 << 64) and is_prime64(m):
            fs.append(m)
            continue
        f = pollard_rho(m)
        if f is None:
            unresolved.append(m)
        else:
            stack.extend((f, m // f))
    counts: dict[int, int] = {}
    for p in fs:
        counts[p] = counts.get(p, 0) + 1
    return sorted(counts.items()), sorted(unresolved)


def fmt_factor(n: int) -> str:
    if n in (0, 1):
        return str(n)
    fs, rem = factor_rigorous(n)
    parts = [str(p) if e == 1 else f"{p}^{e}" for p, e in fs]
    parts += [f"[{r}; exact unresolved cofactor, {len(str(r))} digits]" for r in rem]
    return " * ".join(parts) if parts else str(n)


def top_targets(n: int) -> list[int]:
    b = apery_numbers(n // 2)
    return [q for q in primes_upto(n - 1) if 2 * q > n and b[n - q] % q == 0]


def available_margins(M: int, d: int, L: int, cap: int = 10) -> tuple[int, int]:
    # Stay in the a=1 quotient cell and retain length < q.
    q = d + 1
    s_full = d - (M // 2 + 1)
    t_full = M - (d + L)
    s = min(cap, s_full)
    t = min(cap, t_full)
    while s >= 0 and L + s + t >= q:
        if s >= t:
            s -= 1
        else:
            t -= 1
    assert s >= 0 and t >= 0
    return s, t


def coeff_vector(d: int, L: int, S: int, T: int, s: int, t: int) -> list[int]:
    """Coefficients on Delta^r Y_{d-S}, 0<=r<=L+S+T."""
    a0 = d - S
    a = d - s
    h = a - a0
    m = L + s + t
    N = L + S + T
    out = []
    for r in range(N + 1):
        z = 0
        for k in range(max(0, r - m), min(h, r) + 1):
            z += (-1) ** (r - k) * C(a + r - k, r - k) * C(h, k)
        out.append(z)
    return out


def rank_mod(A: list[list[int]], p: int) -> int:
    if not A:
        return 0
    a = [[x % p for x in row] for row in A]
    m, n = len(a), len(a[0])
    r = 0
    for c in range(n):
        pivot = next((i for i in range(r, m) if a[i][c]), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        inv = pow(a[r][c], -1, p)
        a[r] = [(x * inv) % p for x in a[r]]
        for i in range(m):
            if i != r and a[i][c]:
                z = a[i][c]
                a[i] = [(x - z * y) % p for x, y in zip(a[i], a[r])]
        r += 1
        if r == m:
            break
    return r


def rectangle_algebra_audit(M: int, d: int, L: int, S: int, T: int) -> None:
    q, ell = d + 1, d + L + 1
    assert L + S + T < q < ell
    core = [p for p in primes_upto(ell) if q <= p <= ell]
    Pi = prod(core)
    table = [[Hmargin(M, d, L, s, t) for t in range(T + 1)] for s in range(S + 1)]

    for s in range(S + 1):
        for t in range(T + 1):
            m = L + s + t
            if s:
                lhs = table[s][t] - table[s - 1][t]
                rhs = (-1) ** m * C(ell + t, m) * delta(M, d - s, m)
                assert lhs == rhs
            if t:
                lhs = table[s][t] - table[s][t - 1]
                rhs = (-1) ** m * C(ell + t - 1, m) * delta(M, d - s, m)
                assert lhs == rhs
            if s and t:
                c = gcd(q - s, m)
                assert (m // c) * table[s][t] + ((q - s) // c) * table[s - 1][t] == ((ell + t) // c) * table[s][t - 1]

    base = coeff_vector(d, L, S, T, 0, 0)
    P: list[list[int]] = []
    for s in range(S + 1):
        for t in range(T + 1):
            if s == 0 and t == 0:
                continue
            row = coeff_vector(d, L, S, T, s, t)
            diff = [row[r] - base[r] for r in range(L + 1, L + S + T + 1)]
            assert all(x % Pi == 0 for x in diff)
            P.append(diff)
    r = S + T
    assert rank_mod(P, 1000000007) == r
    assert rank_mod([[x // q for x in row] for row in P], q) == r
    assert rank_mod([[x // ell for x in row] for row in P], ell) == r
    for p in core:
        assert rank_mod([[x // p for x in row] for row in P], p) == r

    # Independent boundary diagonal after row differencing.
    As = [C(ell, L + s) for s in range(1, S + 1)]
    Bs = [C(ell + t - 1, L + t) for t in range(1, T + 1)]
    for z in As + Bs:
        for p in core:
            assert valuation(z, p) == 1
        assert gcd(z // Pi, Pi) == 1


def family_table(M: int, d: int, L: int, S: int, T: int) -> list[list[int]]:
    return [[Hmargin(M, d, L, s, t) for t in range(T + 1)] for s in range(S + 1)]


def gcd_grid(table: list[list[int]]) -> list[list[int]]:
    S, T = len(table) - 1, len(table[0]) - 1
    out = [[0] * (T + 1) for _ in range(S + 1)]
    for s in range(S + 1):
        for t in range(T + 1):
            g = abs(table[s][t])
            if s:
                g = gcd(g, out[s - 1][t])
            if t:
                g = gcd(g, out[s][t - 1])
            out[s][t] = g
    return out


def minimal_unit_points(D: list[list[int]], target_product: int) -> list[tuple[int, int]]:
    pts = []
    S, T = len(D) - 1, len(D[0]) - 1
    for s in range(S + 1):
        for t in range(T + 1):
            if D[s][t] != target_product:
                continue
            if any(D[i][j] == target_product for i in range(s + 1) for j in range(t + 1) if (i, j) != (s, t)):
                continue
            pts.append((s, t))
    return pts


def analyze_interval(n: int, q: int, ell: int, preserved_targets: list[int], label: str) -> None:
    M, d, L = n - 1, q - 1, ell - q
    S, T = available_margins(M, d, L, 10)
    assert L + S + T < q
    rectangle_algebra_audit(M, d, L, S, T)
    H = family_table(M, d, L, S, T)
    D = gcd_grid(H)
    P = prod(preserved_targets)
    assert all(x % P == 0 for row in H for x in row)
    print("=" * 96)
    print(f"INTERVAL {label} n={n} M={M} q={q} ell={ell} L={L} Smax={S} Tmax={T}")
    print(f"preserved_targets={preserved_targets} target_product={P}")
    print(f"H00_digits={len(str(abs(H[0][0])))} D00_digits={len(str(D[0][0]))}")
    print(f"minimal_unit_points={minimal_unit_points(D, P)}")
    first = None
    for total in range(S + T + 1):
        for s in range(S + 1):
            t = total - s
            if 0 <= t <= T and D[s][t] == P:
                first = (s, t)
                break
        if first is not None:
            break
    print(f"first_by_total_margin={first}")
    for k in range(min(S, T) + 1):
        N = D[k][k] // P
        print(f"square k={k}: gcd={D[k][k]} nuisance={N} nuisance_factor={fmt_factor(N)}")
    Nfull = D[S][T] // P
    print(f"full_rectangle_gcd={D[S][T]}")
    print(f"persistent_nuisance={Nfull}")
    print(f"persistent_nuisance_factor={fmt_factor(Nfull)}")
    # Print only strict gcd drops, a compact exact history.
    last = None
    changes = []
    for total in range(S + T + 1):
        for s in range(S + 1):
            t = total - s
            if 0 <= t <= T:
                z = D[s][t] // P
                if z != last:
                    changes.append((s, t, z))
                    last = z
    print(f"total-order_nuisance_changes={changes}")


def candidate_ledger(n: int, A: int = 4, kmax: int = 3) -> None:
    M = n - 1
    ps = [p for p in primes_upto(n - 1) if 2 * p > n]
    gapmax = max(2, int(A * log(n)))
    pairs = [(q, e) for i, q in enumerate(ps) for e in ps[i + 1 :] if e - q <= gapmax]
    print("-" * 96)
    print(f"TARGET_BLIND_LEDGER n={n} gapmax={gapmax} candidate_pairs={len(pairs)} H=n^(1/3)={n**(1/3):.9f}")
    for k in range(kmax + 1):
        total_log = 0.0
        nontrivial = 0
        used = 0
        max_digits = 1
        for q, ell in pairs:
            d, L = q - 1, ell - q
            Sfull, Tfull = available_margins(M, d, L, k)
            if Sfull < k or Tfull < k:
                continue
            table = family_table(M, d, L, k, k)
            D = gcd_grid(table)[k][k]
            total_log += log(max(2, D))
            nontrivial += D > 1
            max_digits = max(max_digits, len(str(D)))
            used += 1
        print(f"rule square({k},{k}): used={used} nontrivial={nontrivial} sum_log={total_log:.12f} sum_log/H={total_log/(n**(1/3)):.12f} max_digits={max_digits}")


def main() -> None:
    rows = (200, 272, 300, 321, 755)
    for n in rows:
        ts = top_targets(n)
        print(f"ROW n={n} top_half_targets={ts}")
        for q, ell in zip(ts, ts[1:]):
            analyze_interval(n, q, ell, [q, ell], f"adjacent-{q}-{ell}")
        if n == 321:
            assert ts == [179, 193, 211]
            analyze_interval(n, 179, 211, ts, "three-target-cover-179-211")

    for n in rows:
        candidate_ledger(n)
    print("Q5714_AUDIT=PASS")


if __name__ == "__main__":
    main()
