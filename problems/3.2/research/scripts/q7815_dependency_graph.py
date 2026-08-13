#!/usr/bin/env python3
"""Exact X=512 audit for the Apéry zero-event overlap graphs.

No floating-point arithmetic is used for graph construction.  We distinguish:
  * unrestricted CRT graph: every pair of nonempty zero sets is adjacent;
  * truncated co-occurrence graph on a specified integer interval;
  * pairwise-dependence graph on the uniform measure of that interval.

The latter two are often confused with a Janson dependency graph.  Pairwise
independence is only a necessary condition for a nonedge in a genuine
standard dependency graph, not a sufficient one.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
from math import isqrt

X = 512


def primes_upto(n: int) -> list[int]:
    sieve = bytearray(b"\x01") * (n + 1)
    sieve[0:2] = b"\x00\x00"
    for p in range(2, isqrt(n) + 1):
        if sieve[p]:
            sieve[p * p : n + 1 : p] = b"\x00" * (((n - p * p) // p) + 1)
    return [p for p in range(2, n + 1) if sieve[p]]


def P(n: int) -> int:
    return 34 * n**3 + 51 * n**2 + 27 * n + 5


def apery_row_mod_p(p: int) -> list[int]:
    b = [0] * p
    b[0] = 1
    b[1] = 5 % p
    for n in range(1, p - 1):
        den = pow(n + 1, 3, p)
        assert den != 0
        num = (P(n) * b[n] - pow(n, 3, p) * b[n - 1]) % p
        b[n + 1] = num * pow(den, -1, p) % p
    return b


def count_residue(lo: int, hi: int, p: int, r: int) -> int:
    """Count n in [lo,hi) with n == r mod p."""
    first = lo + ((r - lo) % p)
    if first >= hi:
        return 0
    return 1 + (hi - 1 - first) // p


def crt_lift(r: int, p: int, s: int, q: int) -> int:
    """The unique x in [0,pq) with x=r mod p and x=s mod q."""
    return r + p * (((s - r) * pow(p, -1, q)) % q)


def count_pair(lo: int, hi: int, p: int, zp: tuple[int, ...], q: int, zq: tuple[int, ...]) -> int:
    mod = p * q
    total = 0
    for r in zp:
        for s in zq:
            x = crt_lift(r, p, s, q)
            k = (lo - x) // mod
            y = x + k * mod
            if y < lo:
                y += mod
            if y < hi:
                total += 1 + (hi - 1 - y) // mod
    return total


def summarize_interval(name: str, lo: int, hi: int, primes: list[int], zeros: dict[int, tuple[int, ...]]) -> None:
    N = hi - lo
    active = [p for p in primes if zeros[p]]
    counts = {
        p: sum(count_residue(lo, hi, p, r) for r in zeros[p])
        for p in primes
    }

    co_adj = {p: set() for p in primes}
    dep_adj = {p: set() for p in primes}
    overlap_hist = Counter()
    exact_independent_pairs: list[tuple[int, int, int, int, int]] = []
    active_nonoverlap: list[tuple[int, int]] = []

    for i, p in enumerate(primes):
        for q in primes[i + 1 :]:
            J = count_pair(lo, hi, p, zeros[p], q, zeros[q])
            overlap_hist[J] += 1
            if J > 0:
                co_adj[p].add(q)
                co_adj[q].add(p)
            elif zeros[p] and zeros[q]:
                active_nonoverlap.append((p, q))

            # Pairwise dependence under uniform n in [lo,hi).
            if J * N != counts[p] * counts[q]:
                dep_adj[p].add(q)
                dep_adj[q].add(p)
            else:
                exact_independent_pairs.append((p, q, J, counts[p], counts[q]))

    def graph_stats(adj: dict[int, set[int]], vertices: list[int]) -> tuple[int, int, Fraction, int, int]:
        degrees = [len(adj[p] & set(vertices)) for p in vertices]
        edges = sum(degrees) // 2
        possible = len(vertices) * (len(vertices) - 1) // 2
        density = Fraction(edges, possible) if possible else Fraction(0, 1)
        return min(degrees, default=0), max(degrees, default=0), density, edges, sum(degrees)

    co_active = graph_stats(co_adj, active)
    dep_active = graph_stats(dep_adj, active)
    co_all = graph_stats(co_adj, primes)
    dep_all = graph_stats(dep_adj, primes)

    # Direct incidence array for factorial moments and maximum K.
    kvals = [0] * N
    for p in primes:
        for r in zeros[p]:
            first = lo + ((r - lo) % p)
            for n in range(first, hi, p):
                kvals[n - lo] += 1
    F1 = sum(kvals)
    F2 = sum(k * (k - 1) for k in kvals)
    F3 = sum(k * (k - 1) * (k - 2) for k in kvals)
    max_k = max(kvals, default=0)
    max_rows = [lo + i for i, k in enumerate(kvals) if k == max_k]

    print(f"INTERVAL {name} lo={lo} hi={hi} N={N}")
    print(f"  prime_count={len(primes)} active_count={len(active)} unrestricted_CRT_max_degree={max(0,len(active)-1)}")
    print(f"  active_primes={active}")
    print(f"  zero_sizes={[(p,len(zeros[p])) for p in active]}")
    print(
        "  cooccurrence_active "
        f"min_degree={co_active[0]} max_degree={co_active[1]} "
        f"edges={co_active[3]} density={float(co_active[2]):.12f} "
        f"density_exact={co_active[2]}"
    )
    print(
        "  pairwise_dependence_active "
        f"min_degree={dep_active[0]} max_degree={dep_active[1]} "
        f"edges={dep_active[3]} density={float(dep_active[2]):.12f} "
        f"density_exact={dep_active[2]}"
    )
    print(
        "  cooccurrence_all "
        f"min_degree={co_all[0]} max_degree={co_all[1]} edges={co_all[3]} density={float(co_all[2]):.12f}"
    )
    print(
        "  pairwise_dependence_all "
        f"min_degree={dep_all[0]} max_degree={dep_all[1]} edges={dep_all[3]} density={float(dep_all[2]):.12f}"
    )
    print(f"  active_nonoverlap_count={len(active_nonoverlap)} active_nonoverlap_pairs={active_nonoverlap}")
    print(f"  exact_pairwise_independence_count={len(exact_independent_pairs)}")
    print(f"  exact_pairwise_independence_pairs={exact_independent_pairs}")
    print(f"  overlap_histogram={sorted(overlap_hist.items())}")
    print(f"  cooccurrence_degrees_active={[(p,len(co_adj[p] & set(active))) for p in active]}")
    print(f"  pairwise_dependence_degrees_active={[(p,len(dep_adj[p] & set(active))) for p in active]}")
    print(f"  F1={F1} F2={F2} F3={F3} maxK={max_k} max_rows={max_rows[:30]} max_row_count={len(max_rows)}")
    print()


def main() -> None:
    primes = [p for p in primes_upto(2 * X) if X < p <= 2 * X]
    zeros: dict[int, tuple[int, ...]] = {}
    for p in primes:
        row = apery_row_mod_p(p)
        zp = tuple(i for i, x in enumerate(row) if x == 0)
        # Exact regression checks from the proved structure.
        assert 0 not in zp and p - 1 not in zp
        assert all((p - 1 - r) in zp for r in zp)
        assert all((r + 1) not in zp for r in zp if r + 1 < p)
        zeros[p] = zp

    print(f"Q7815 X={X}")
    print(f"window_primes={primes}")
    print(f"prime_count={len(primes)}")
    print(f"active_count={sum(bool(zeros[p]) for p in primes)}")
    print(f"total_zero_count={sum(len(zeros[p]) for p in primes)}")
    print()

    M = X * X
    summarize_interval("HM3_prefix", 0, M, primes, zeros)
    summarize_interval("exact_dyadic_shell", M + 1, 2 * M + 1, primes, zeros)
    summarize_interval("half_open_dyadic_shell", M, 2 * M, primes, zeros)
    summarize_interval("double_prefix", 0, 2 * M, primes, zeros)
    print("Q7815_DEPENDENCY_GRAPH_SUCCESS")


if __name__ == "__main__":
    main()
