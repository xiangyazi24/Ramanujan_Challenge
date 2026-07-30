#!/usr/bin/env python3
"""Exact standard-library audit for Q5719.

The scan uses the Q5715 block generator, but only retains blocks satisfying
Q5715's own extended quotient-cell assertion.  It enumerates every outside
common ghost rho>N through n<=limit.  Under rho>s+a it evaluates the correct
three-variable Cartier-Hadamard scalar.  It evaluates the unreduced shell in
the non-separated cases and in deterministic cross-checks, and compares every
zero scalar with both Newton carriers modulo rho.
"""
from __future__ import annotations

from argparse import ArgumentParser
from functools import lru_cache
from math import comb, isqrt


def primes_upto(n: int) -> list[int]:
    if n < 2:
        return []
    sieve = bytearray(b"\x01") * (n + 1)
    sieve[:2] = b"\x00\x00"
    for p in range(2, isqrt(n) + 1):
        if sieve[p]:
            sieve[p*p:n+1:p] = b"\x00" * (((n-p*p)//p)+1)
    return [p for p in range(2, n + 1) if sieve[p]]


PRIMES = primes_upto(2000)


def prime_divisors(n: int) -> list[int]:
    out: list[int] = []
    x = n
    for p in PRIMES:
        if p*p > x:
            break
        if x % p == 0:
            out.append(p)
            while x % p == 0:
                x //= p
        if x == 1:
            break
    if x > 1:
        out.append(x)
    return out


def cell_bounds(n: int, a: int) -> tuple[int, int]:
    M = n - a
    lo = M // (a + 1) + 1
    hi = M // a
    lo = max(lo, isqrt(n))
    return lo, hi


def blocks_for_cell(n: int, a: int) -> list[tuple[int, int]]:
    """Literal Q5715 block generator."""
    lo, hi = cell_bounds(n, a)
    length = hi - lo + 1
    if length < 2:
        return []
    h = max(2, round(n ** (1 / 3)))
    scales = sorted({2, 3, 4, 5, 8, h, min(length, 2*h), length})
    out: list[tuple[int, int]] = []
    for N in scales:
        if N > length:
            continue
        step = max(1, N // 2)
        starts = list(range(lo, hi - N + 2, step))
        if starts[-1] != hi - N + 1:
            starts.append(hi - N + 1)
        out.extend((D, N) for D in starts)
    return sorted(set(out))


def in_scaled_P(v: tuple[int, int, int], a: int) -> bool:
    x, y, z = v
    return (
        -a <= x <= a and -a <= y <= a and -a <= z <= a
        and x - y <= a and x - z <= a
    )


def polytope_points(a: int):
    for x in range(-a, a + 1):
        for y in range(-a, a + 1):
            for z in range(-a, a + 1):
                if x - y <= a and x - z <= a:
                    yield x, y, z


@lru_cache(maxsize=None)
def factorial_tables(p: int) -> tuple[tuple[int, ...], tuple[int, ...]]:
    fact = [1] * p
    for i in range(1, p):
        fact[i] = fact[i-1] * i % p
    invfact = [1] * p
    invfact[p-1] = pow(fact[p-1], p-2, p)
    for i in range(p-1, 0, -1):
        invfact[i-1] = invfact[i] * i % p
    return tuple(fact), tuple(invfact)


def small_binom_mod(n: int, k: int, p: int) -> int:
    if k < 0 or k > n:
        return 0
    fact, invfact = factorial_tables(p)
    return fact[n] * invfact[k] % p * invfact[n-k] % p


def binom_mod(n: int, k: int, p: int) -> int:
    if k < 0 or k > n:
        return 0
    out = 1
    while n or k:
        ni, ki = n % p, k % p
        if ki > ni:
            return 0
        out = out * small_binom_mod(ni, ki, p) % p
        n //= p
        k //= p
    return out


@lru_cache(maxsize=None)
def coeff_mod(n: int, u: int, v: int, w: int, p: int) -> int:
    out = 0
    for t in range(n + 1):
        out += (
            binom_mod(n, t, p)
            * binom_mod(n, t-u, p)
            * binom_mod(2*n-t, n-v, p)
            * binom_mod(2*n-t, n-w, p)
        )
        out %= p
    return out


@lru_cache(maxsize=None)
def shell_mod(M: int, d: int, p: int) -> int:
    a = M // d
    rowM = [binom_mod(M, j, p) for j in range(M + 1)]
    out = 0
    for t in range(M + 1):
        xp = 0
        yp = 0
        for u in range(-a, a + 1):
            j = M - t + d*u
            if 0 <= j <= M:
                xp += rowM[j]
            if 0 <= j <= 2*M-t:
                yp += binom_mod(2*M-t, j, p)
        out = (out + rowM[t] * (xp % p) * (yp % p) ** 2) % p
    return out


def weight_mod(d: int, L: int, j: int, p: int) -> int:
    return ((-1 if j & 1 else 1)
            * binom_mod(d+j, j, p)
            * binom_mod(d+L+1, L-j, p)) % p


def carrier_mod(M: int, d: int, L: int, p: int) -> int:
    return sum(weight_mod(d, L, j, p) * shell_mod(M, d+j, p)
               for j in range(L+1)) % p


def ghost_hadamard_mod(M: int, rho: int, k: int) -> int:
    alpha, s = divmod(M, rho)
    d = k*rho - 1
    a = M // d
    out = 0
    for kap in polytope_points(a):
        kk = tuple(k*x for x in kap)
        neg = tuple(-x for x in kap)
        if not in_scaled_P(kk, alpha) or not in_scaled_P(neg, s):
            continue
        out += coeff_mod(alpha, *kk, rho) * coeff_mod(s, *neg, rho)
        out %= rho
    return out


def enumerate_geometry(limit: int):
    raw: list[dict] = []
    core_blocks = 0
    valid_blocks = 0
    invalid_blocks = 0
    for n in range(20, limit + 1):
        for a in range(1, isqrt(n) + 2):
            M = n - a
            for D, N in blocks_for_cell(n, a):
                if not all(M // d == a for d in range(D, D+N)):
                    continue
                core_blocks += 1
                extended_ok = all(d > 0 and M // d == a
                                  for d in range(D-1, D+N+1))
                if not extended_ok:
                    invalid_blocks += 1
                    continue
                valid_blocks += 1
                for m in range(D+1, D+N+1):
                    for rho in prime_divisors(m):
                        if rho <= N:
                            continue
                        k = m // rho
                        if k < 2:  # k=1 is a candidate prime, not outside ghost
                            continue
                        alpha, s = divmod(M, rho)
                        d = m - 1
                        assert M // d == a
                        raw.append({
                            "n": n, "a": a, "M": M, "D": D, "N": N,
                            "rho": rho, "k": k, "d": d,
                            "alpha": alpha, "s": s,
                            "sep": rho > s + a,
                        })
    return raw, core_blocks, valid_blocks, invalid_blocks


def ranges(rows: list[dict]) -> dict:
    keys = ["n", "a", "M", "D", "N", "rho", "k", "d", "alpha", "s"]
    if not rows:
        return {}
    return {key: [min(r[key] for r in rows), max(r[key] for r in rows)]
            for key in keys}


def main() -> None:
    ap = ArgumentParser()
    ap.add_argument("--limit", type=int, default=2000)
    ap.add_argument("--evaluate", action="store_true")
    args = ap.parse_args()

    raw, core_blocks, valid_blocks, invalid_blocks = enumerate_geometry(args.limit)
    unique_map: dict[tuple, dict] = {}
    occurrence_map: dict[tuple, list[dict]] = {}
    for r in raw:
        key = (r["n"], r["a"], r["rho"], r["k"], r["d"])
        unique_map.setdefault(key, r)
        occurrence_map.setdefault(key, []).append(r)
    unique = list(unique_map.values())
    sep = [r for r in unique if r["sep"]]
    nonsep = [r for r in unique if not r["sep"]]

    print("LITERAL_Q5715_FIRST_ASSERT_FAILURE n=20 a=1 D=10 N=2 M=19 because 19//9=2")
    print("GEOMETRY", {
        "core_blocks": core_blocks,
        "extended_valid_blocks": valid_blocks,
        "extended_invalid_blocks": invalid_blocks,
        "raw_outside_common_ghost_occurrences": len(raw),
        "unique_outside_common_ghosts": len(unique),
        "separated_unique": len(sep),
        "nonseparated_unique": len(nonsep),
    })
    print("RANGES_ALL", ranges(unique))
    print("RANGES_SEPARATED", ranges(sep))
    print("RANGES_NONSEPARATED", ranges(nonsep))
    print("MAX_K_EXAMPLE", max(unique, key=lambda r: r["k"]) if unique else None)
    print("FIRST_A1_K_GT_2", next((r for r in unique if r["a"] == 1 and r["k"] > 2), None))
    print("FIRST_NONSEPARATED", nonsep[0] if nonsep else None)

    if not args.evaluate:
        return

    zeros: list[dict] = []
    shell_crosschecks = 0
    hadamard_checks = 0
    carrier_checks = 0
    for idx, original in enumerate(unique, 1):
        r = dict(original)
        if r["sep"]:
            scalar = ghost_hadamard_mod(r["M"], r["rho"], r["k"])
            hadamard_checks += 1
            # Deterministic direct-shell checks, plus every actual zero.
            direct = idx <= 200 or idx % 1000 == 0 or scalar == 0
            if direct:
                sh = shell_mod(r["M"], r["d"], r["rho"])
                assert sh == scalar, ("H mismatch", r, sh, scalar)
                shell_crosschecks += 1
        else:
            scalar = shell_mod(r["M"], r["d"], r["rho"])
            shell_crosschecks += 1
        r["scalar"] = scalar
        if scalar == 0:
            key = (r["n"], r["a"], r["rho"], r["k"], r["d"])
            for q in occurrence_map[key]:
                A = carrier_mod(q["M"], q["D"]-1, q["N"], q["rho"])
                B = carrier_mod(q["M"], q["D"], q["N"], q["rho"])
                assert A == scalar and B == scalar, ("carrier mismatch", q, A, B)
                carrier_checks += 1
            r["block_occurrences"] = len(occurrence_map[key])
            zeros.append(r)
        if idx % 1000 == 0:
            print("PROGRESS", idx, "of", len(unique), "zeros", len(zeros), flush=True)

    sepzeros = [r for r in zeros if r["sep"]]
    nonsepzeros = [r for r in zeros if not r["sep"]]
    print("EVALUATION", {
        "hadamard_evaluations": hadamard_checks,
        "direct_shell_crosschecks": shell_crosschecks,
        "carrier_occurrence_checks": carrier_checks,
        "unique_scalar_zeros": len(zeros),
        "separated_zeros": len(sepzeros),
        "nonseparated_zeros": len(nonsepzeros),
    })
    print("ZERO_RANGES", ranges(zeros))
    print("SEPARATED_ZERO_RANGES", ranges(sepzeros))
    print("NONSEPARATED_ZERO_RANGES", ranges(nonsepzeros))
    print("ZERO_ROWS_BEGIN")
    for r in zeros:
        print(r)
    print("ZERO_ROWS_END")


if __name__ == "__main__":
    main()
