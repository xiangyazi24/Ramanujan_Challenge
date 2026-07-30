#!/usr/bin/env python3
"""Exact audit for Q5708: adjacent Newton carriers and their residual gcds.

The script uses only exact integer arithmetic.  It audits the five requested
rows and a wider n<=1000 scan.  The shell formula is the one-fold formula in
Q32_SEPARATION_ANALYSIS.md, Sections 48--51.
"""

from __future__ import annotations

from collections import defaultdict
from functools import lru_cache
from math import comb, gcd, isqrt, log

try:
    from sympy import factorint, isprime
except Exception:
    factorint = None
    isprime = None


def C(n: int, k: int) -> int:
    return comb(n, k) if 0 <= k <= n else 0


def primes_upto(N: int) -> list[int]:
    sieve = bytearray(b"\x01") * (N + 1)
    if N >= 0:
        sieve[0] = 0
    if N >= 1:
        sieve[1] = 0
    for p in range(2, isqrt(N) + 1):
        if sieve[p]:
            sieve[p * p : N + 1 : p] = b"\x00" * (((N - p * p) // p) + 1)
    return [p for p in range(2, N + 1) if sieve[p]]


PRIMES = primes_upto(200000)


def apery_mod(r: int, p: int) -> int:
    return sum((C(r, k) * C(r + k, k)) ** 2 for k in range(r + 1)) % p


@lru_cache(maxsize=None)
def apery(n: int) -> int:
    return sum((C(n, k) * C(n + k, k)) ** 2 for k in range(n + 1))


@lru_cache(maxsize=None)
def coefficient(M: int, u: int, v: int, w: int) -> int:
    out = 0
    for t in range(M + 1):
        out += C(M, t) * C(M, t - u) * C(2 * M - t, M - v) * C(2 * M - t, M - w)
    return out


@lru_cache(maxsize=None)
def shell(M: int, d: int) -> int:
    """C_M(d), with the quotient a=floor(M/d) determined exactly."""
    a = M // d
    out = 0
    for t in range(M + 1):
        xp = sum(C(M, M - t + d * u) for u in range(-a, a + 1))
        yp = sum(C(2 * M - t, M - t + d * v) for v in range(-a, a + 1))
        out += C(M, t) * xp * yp * yp
    return out


def delta(M: int, d: int, order: int) -> int:
    return sum((-1) ** (order - i) * C(order, i) * shell(M, d + i) for i in range(order + 1))


def weight(d: int, L: int, i: int) -> int:
    return (-1) ** i * C(d + i, i) * C(d + L + 1, L - i)


def G(M: int, d: int, L: int) -> int:
    return sum(weight(d, L, i) * shell(M, d + i) for i in range(L + 1))


def valuation(x: int, p: int) -> int:
    x = abs(x)
    e = 0
    while x and x % p == 0:
        x //= p
        e += 1
    return e


def short_factor(n: int, trial_limit: int = 200000) -> tuple[dict[int, int], int, str]:
    """Exact trial factorization plus primality status of the remaining cofactor."""
    n = abs(n)
    fac: dict[int, int] = {}
    if n in (0, 1):
        return fac, n, "unit" if n == 1 else "zero"
    for p in PRIMES:
        if p > trial_limit or p * p > n:
            break
        if n % p == 0:
            e = 0
            while n % p == 0:
                n //= p
                e += 1
            fac[p] = e
    if n == 1:
        return fac, 1, "complete"
    if isprime is not None and isprime(n):
        fac[n] = fac.get(n, 0) + 1
        return fac, 1, "complete"
    # Ask factorint only for moderate cofactors; with limit this remains bounded.
    if factorint is not None and n.bit_length() <= 512:
        try:
            ff = factorint(n, limit=trial_limit, use_ecm=True)
            rem = 1
            for p, e in ff.items():
                p = int(p)
                if isprime is not None and isprime(p):
                    fac[p] = fac.get(p, 0) + int(e)
                else:
                    rem *= p ** int(e)
            if rem == 1:
                return dict(sorted(fac.items())), 1, "complete"
            return dict(sorted(fac.items())), rem, "composite cofactor"
        except Exception as exc:
            return dict(sorted(fac.items())), n, f"unfactored ({type(exc).__name__})"
    return dict(sorted(fac.items())), n, "composite cofactor"


def fmt_factor(n: int) -> str:
    fac, rem, status = short_factor(n)
    pieces = [f"{p}^{e}" if e != 1 else str(p) for p, e in fac.items()]
    if rem not in (0, 1):
        pieces.append(f"C{len(str(rem))}={rem}")
    if not pieces:
        pieces = [str(abs(n))]
    return " * ".join(pieces) + f" [{status}]"


def targets(n: int) -> list[tuple[int, int, int]]:
    """Large-prime targets (q,a,r), q>sqrt(n), grouped by fixed quotient."""
    out = []
    for q in primes_upto(n):
        if q <= isqrt(n):
            continue
        a, r = divmod(n, q)
        if 1 <= r <= q - 2 and apery_mod(r, q) == 0:
            out.append((q, a, r))
    return out


def adjacent_same_cell(n: int) -> list[tuple[int, int, int]]:
    by_a: dict[int, list[int]] = defaultdict(list)
    for q, a, _ in targets(n):
        by_a[a].append(q)
    out = []
    for a, qs in by_a.items():
        qs.sort()
        for q, ell in zip(qs, qs[1:]):
            out.append((q, ell, a))
    return sorted(out)


def primitive_relation(d: int, L: int) -> tuple[int, int, int, int]:
    t = gcd(d + 1, L)
    u = (d + 1) // t
    v = L // t
    return t, u, v, u + v


def packet_content(M: int, d: int) -> tuple[int, int, int]:
    """gcd of origin+nonorigin coefficients and nonorigin coefficients."""
    a = M // d
    all_g = 0
    radial_g = 0
    count = 0
    for x in range(-a, a + 1):
        for y in range(-a, a + 1):
            for z in range(-a, a + 1):
                c = coefficient(M, d * x, d * y, d * z)
                if c == 0:
                    continue
                count += 1
                all_g = gcd(all_g, c)
                if (x, y, z) != (0, 0, 0):
                    radial_g = gcd(radial_g, c)
    return all_g, radial_g, count


def mobius(n: int) -> int:
    if n == 1:
        return 1
    mu = 1
    x = n
    for p in primes_upto(isqrt(n) + 1):
        if p * p > x:
            break
        if x % p == 0:
            x //= p
            mu = -mu
            if x % p == 0:
                return 0
            while x % p == 0:
                x //= p
    if x > 1:
        mu = -mu
    return mu


def primitive_radial_packet(M: int, d: int) -> int:
    """P_M(d)=sum_{gcd(|nu|)=d} c_M(nu), by finite Mobius inversion."""
    bM = apery(M)
    return sum(mobius(k) * (shell(M, k * d) - bM) for k in range(1, M // d + 1))


def analyze_pair(n: int, q: int, ell: int, a: int, packet_audit: bool = False) -> dict[str, object]:
    M = n - a
    d = q - 1
    L = ell - q
    assert M // d == a and M // (d + L + 1) == a
    vals = [shell(M, d + i) for i in range(L + 2)]
    gdL = G(M, d, L)
    g1L = G(M, d + 1, L)
    gdLm = G(M, d, L - 1) if L else vals[0]
    g1Lm = G(M, d + 1, L - 1) if L else vals[1]
    high = delta(M, d, L + 1)
    B = C(d + L + 1, L)
    assert gdL - g1L == (-1) ** (L + 1) * B * high
    R = gcd(gdL, high)
    GG = gcd(gdL, g1L)
    G4 = gcd(gcd(abs(gdL), abs(g1L)), gcd(abs(gdLm), abs(g1Lm)))
    t, u, v, uv = primitive_relation(d, L)
    assert v * gdL + u * g1Lm == uv * gdLm

    # Smith content of the 4x3 Pascal evaluation matrix.
    A = C(d + L, L)
    Cc = C(d + L, L - 1)
    Ap = C(d + L + 1, L)
    smith3 = Ap * gcd(A, Cc)
    assert gcd(A, Cc) * u == A
    assert gcd(A, Cc) * v == Cc

    radial = [primitive_radial_packet(M, d + i) for i in range(L + 2)]
    # Exact Mobius reconstruction at each shell.
    bM = apery(M)
    for i in range(L + 2):
        dd = d + i
        assert shell(M, dd) == bM + sum(primitive_radial_packet(M, k * dd) for k in range(1, M // dd + 1))

    result: dict[str, object] = {
        "n": n, "q": q, "ell": ell, "a": a, "M": M, "d": d, "L": L,
        "R": R, "GG": GG, "G4": G4, "high": high, "B": B,
        "R_factor": fmt_factor(R), "GG_factor": fmt_factor(GG), "G4_factor": fmt_factor(G4),
        "target_v_R": (valuation(R, q), valuation(R, ell)),
        "target_v_GG": (valuation(GG, q), valuation(GG, ell)),
        "target_v_G4": (valuation(G4, q), valuation(G4, ell)),
        "primitive_relation": (u, v, uv),
        "unit_pivot": min(u, v) == 1,
        "smith3": smith3,
        "radial_gcd": gcd(*[abs(x) for x in radial]) if radial else 0,
        "radial_endpoint_coefficients": (1, (-1) ** (L + 1)),
        "digits": (len(str(abs(gdL))), len(str(abs(high))), len(str(abs(R))), len(str(abs(GG)))),
    }
    if packet_audit:
        result["packet_content_d"] = packet_content(M, d)
        result["packet_content_ell"] = packet_content(M, d + L)
    return result


def print_result(r: dict[str, object]) -> None:
    print(f"PAIR n={r['n']} q={r['q']} ell={r['ell']} a={r['a']} M={r['M']} d={r['d']} L={r['L']}")
    print(f"  digits(G,high,R,gcdGG)={r['digits']}")
    print(f"  R={r['R_factor']}")
    print(f"  gcd(G_dL,G_d+1,L)={r['GG_factor']}")
    print(f"  gcd4={r['G4_factor']}")
    print(f"  valuations R(q,ell)={r['target_v_R']} GG(q,ell)={r['target_v_GG']} G4(q,ell)={r['target_v_G4']}")
    print(f"  primitive Pascal relation (u,v,u+v)={r['primitive_relation']} unit_pivot={r['unit_pivot']}")
    print(f"  Smith maximal-minor content={r['smith3']}")
    print(f"  gcd primitive radial packet values={r['radial_gcd']}")
    if 'packet_content_d' in r:
        print(f"  coefficient packet contents at d,ell-node={r['packet_content_d']}, {r['packet_content_ell']}")


def main() -> None:
    requested = [200, 272, 300, 321, 755]
    print("=== REQUESTED ROW TARGETS ===")
    selected_results = []
    for n in requested:
        ts = targets(n)
        pairs = adjacent_same_cell(n)
        print(f"n={n} targets={ts}")
        print(f"n={n} adjacent_same_cell={pairs}")
        for q, ell, a in pairs:
            r = analyze_pair(n, q, ell, a, packet_audit=(n in (200, 321)))
            selected_results.append(r)
            print_result(r)

    print("=== WIDER TARGET SCAN n<=1000, gap<=ceil(10 log n) ===")
    candidate_pairs = []
    for n in range(20, 1001):
        gap_bound = max(2, int(10 * log(n) + 0.999999))
        for q, ell, a in adjacent_same_cell(n):
            if ell - q <= gap_bound:
                candidate_pairs.append((n, q, ell, a))
    print(f"candidate_pair_count={len(candidate_pairs)}")
    print("candidate_pairs=", candidate_pairs)

    # Exact residual audit for every feasible pair with L<=40; this is the
    # range that completes comfortably while still extending well beyond the
    # five named rows.  Larger candidates are listed above for reproduction.
    audited = []
    for tup in candidate_pairs:
        n, q, ell, a = tup
        if ell - q <= 40:
            rr = analyze_pair(n, q, ell, a, packet_audit=False)
            audited.append(rr)
    print(f"audited_pair_count={len(audited)}")
    counts = defaultdict(int)
    for rr in audited:
        key = (
            rr["target_v_R"][0] > 0,
            rr["target_v_R"][1] > 0,
            rr["target_v_GG"][0] > 0,
            rr["target_v_GG"][1] > 0,
            rr["target_v_G4"][0] > 0,
            rr["target_v_G4"][1] > 0,
            rr["unit_pivot"],
        )
        counts[key] += 1
    print("classification_counts=")
    for k in sorted(counts):
        print(" ", k, counts[k])

    print("=== AUDITED PAIRS ===")
    for rr in audited:
        print_result(rr)


if __name__ == "__main__":
    main()
