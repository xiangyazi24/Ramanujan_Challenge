#!/usr/bin/env python3
"""Q8345 exact verifier: far physical one-label obstruction.

Standard library only; exact integer arithmetic only.

This verifier is for finite identities and falsification/regression.  No finite
scan is used as an asymptotic proof.  The accompanying English proofs are in
chatgpt_q8345_far_physical_one_label_obstruction.md.
"""

from __future__ import annotations

from fractions import Fraction
from math import gcd


def P(n: int) -> int:
    return 34 * n**3 + 51 * n**2 + 27 * n + 5


def tau(n: int) -> int:
    return n * (n + 1)


def apery_row_mod(p: int) -> list[int]:
    """b_0,...,b_{p-1} modulo p from the exact Apéry recurrence."""
    assert p >= 3
    b = [0] * p
    b[0] = 1
    b[1] = 5 % p
    for n in range(1, p - 1):
        numerator = (P(n) * b[n] - n**3 * b[n - 1]) % p
        denominator = (n + 1) ** 3 % p
        b[n + 1] = numerator * pow(denominator, -1, p) % p
    return b


def zero_set(p: int) -> list[int]:
    row = apery_row_mod(p)
    return [r for r, value in enumerate(row) if value == 0]


def fold(p: int, r: int) -> int:
    assert 0 <= r < p
    return min(r, p - 1 - r)


def signed_lift(m: int, p: int) -> tuple[int, int, int, int, str]:
    """Return (r,h,s,c,side), with m=c*p+s and tau(s)=tau(h)."""
    r = m % p
    if 2 * r <= p - 1:
        h = r
        s = h
        side = "direct"
    else:
        h = p - 1 - r
        s = -h - 1
        side = "reflected"
    assert -p <= s <= p - 1
    assert (m - s) % p == 0
    c = (m - s) // p
    assert m == c * p + s
    assert tau(s) == tau(h)
    return r, h, s, c, side


def kummer_order(p: int, h: int) -> int:
    return (p - 1) // gcd(p - 1, h)


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2
    return True


def shell_prime_divisors(n: int, X: int) -> list[int]:
    return [p for p in range(X + 1, 2 * X + 1) if is_prime(p) and n % p == 0]


def euler_phi(n: int) -> int:
    out = n
    d = 2
    x = n
    while d * d <= x:
        if x % d == 0:
            while x % d == 0:
                x //= d
            out -= out // d
        d += 1
    if x > 1:
        out -= out // x
    return out


def character_order_count(p: int, D: int) -> int:
    """Number of exponents r mod p-1 whose character order is <=D."""
    n = p - 1
    count = 0
    for r in range(n):
        order = n // gcd(n, r)
        if order <= D:
            count += 1
    by_divisors = sum(euler_phi(d) for d in range(1, D + 1) if n % d == 0)
    assert count == by_divisors
    return count


def graph_edges(signed_positions: dict[int, int], g: int) -> list[tuple[int, int]]:
    labels = sorted(signed_positions)
    out = []
    for i, q in enumerate(labels):
        for ell in labels[i + 1 :]:
            if abs(signed_positions[q] - signed_positions[ell]) == g:
                out.append((q, ell))
    return out


def graph_is_path_forest(signed_positions: dict[int, int], g: int) -> bool:
    edges = graph_edges(signed_positions, g)
    deg = {q: 0 for q in signed_positions}
    for a, b in edges:
        deg[a] += 1
        deg[b] += 1
    if any(d > 2 for d in deg.values()):
        return False
    # Every edge changes the integer coordinate by +-g.  A cycle would have a
    # smallest coordinate vertex with two larger neighbours at the same point.
    # The explicit finite check below additionally verifies forest edge counts.
    vertices = {q for e in edges for q in e}
    if not vertices:
        return True
    seen = set()
    for root in vertices:
        if root in seen:
            continue
        stack = [root]
        comp_v = set()
        comp_e_twice = 0
        while stack:
            u = stack.pop()
            if u in comp_v:
                continue
            comp_v.add(u)
            seen.add(u)
            for a, b in edges:
                if a == u or b == u:
                    comp_e_twice += 1
                    v = b if a == u else a
                    if v not in comp_v:
                        stack.append(v)
        comp_e = comp_e_twice // 2
        if comp_e != len(comp_v) - 1:
            return False
    return True


def verify_actual_apery_edge() -> None:
    X = 128
    m = 321
    data = {p: signed_lift(m, p) for p in (179, 193, 211)}

    assert data[179] == (142, 36, -37, 2, "reflected")
    assert data[193] == (128, 64, -65, 2, "reflected")
    assert data[211] == (110, 100, -101, 2, "reflected")

    # Exact Apéry zero sets in the two endpoints used for the counterexample.
    z179 = zero_set(179)
    z193 = zero_set(193)
    z211 = zero_set(211)
    assert z179 == [36, 142]
    assert 64 in z193 and 128 in z193
    assert 100 in z211 and 110 in z211

    row179 = apery_row_mod(179)
    assert row179[36] == 0
    assert row179[64] != 0  # Cross-characteristic edge is not a q=179 gap-28 return.

    r179, h179, s179, c179, _ = data[179]
    r193, h193, s193, c193, _ = data[193]
    assert abs(s179 - s193) == 28
    assert abs(h179 - h193) == 28
    assert c179 * 179 + 28 == c193 * 193 == 386
    assert shell_prime_divisors(c179 * 179 + 28, X) == [193]
    assert shell_prime_divisors(c179 * 179 - 28, X) == []

    # Triangular carrier is one-label: tau(m)-tau(h)=c*p*(m+s+1).
    for p, (_, h, s, c, _) in data.items():
        assert tau(m) - tau(h) == c * p * (m + s + 1)

    # For every eta<1/15, D at X=128 is <=2: exponent <1/2,
    # sqrt(128)<12 and log(128)>4 imply D<3.  We verify the retained orders.
    assert kummer_order(179, 36) == 89
    assert kummer_order(193, 64) == 3
    assert kummer_order(211, 100) == 21
    assert min(kummer_order(p, data[p][1]) for p in data) > 2

    # Exact D<=2 isolation for the two endpoints of the g=28 edge.
    for p in (179, 193):
        folded = sorted({fold(p, r) for r in zero_set(p)})
        h = data[p][1]
        assert all(abs(h - t) > 2 for t in folded if t != h)

    # The three signed lifts make the path-forest geometry literal.
    signed = {p: data[p][2] for p in data}
    for g in (28, 36, 64):
        assert graph_is_path_forest(signed, g)
        edges = graph_edges(signed, g)
        assert len(edges) == 1
        incident = {v for e in edges for v in e}
        assert len(incident) == 2


def verify_partner_uniqueness() -> None:
    # Exact size argument tested on a broad finite grid: N<=X^2+X-1 cannot
    # contain two distinct prime divisors both >X.
    for X in range(8, 80):
        shell = [p for p in range(X + 1, 2 * X + 1) if is_prime(p)]
        limit = X * X + X - 1
        for i, p in enumerate(shell):
            for q in shell[i + 1 :]:
                assert p * q > limit


def verify_graph_combinatorics() -> None:
    # Exhaust all subsets of a small integer interval and all gaps.
    positions = list(range(-5, 6))
    for mask in range(1 << len(positions)):
        chosen = [positions[i] for i in range(len(positions)) if (mask >> i) & 1]
        labels = {i: s for i, s in enumerate(chosen)}
        for g in range(1, 7):
            assert graph_is_path_forest(labels, g)
            edges = graph_edges(labels, g)
            E = len(edges)
            incident = {v for e in edges for v in e}
            V = len(incident)
            if V:
                assert V / 2 <= E < V
            # Exact rectangle/chain forcing inequality used in the report.
            for T in (4, 5, 6, 7, 8):
                if V > T:
                    pair_edges = E * (E - 1) // 2
                    assert V - T <= Fraction(16 * pair_edges, T)


def verify_character_order_count() -> None:
    for p in (179, 193, 211, 257):
        for D in range(1, 15):
            count = character_order_count(p, D)
            assert count <= D * (D + 1) // 2


def verify_exponent_ledger() -> None:
    # One exact sample eta demonstrates the algebraic exponent identities.
    eta = Fraction(1, 60)
    aD = Fraction(7, 15) + eta / 2
    aG = 2 * aD
    aL = Fraction(14, 15) + eta
    aT = 2 * aL - 1
    alpha_star = aT / aG
    assert aG == Fraction(19, 20)
    assert aT == Fraction(9, 10)
    assert alpha_star == Fraction(18, 19) < 1
    # Quadratic fixed-gap height misses the threshold by exactly one X-power
    # at the exponent level: 2*aG-aT = 1.
    assert 2 * aG - aT == 1
    # Linear height still misses for eta<1/15.
    assert aG - aT == Fraction(1, 20) > 0


def main() -> None:
    verify_actual_apery_edge()
    verify_partner_uniqueness()
    verify_graph_combinatorics()
    verify_character_order_count()
    verify_exponent_ledger()
    print("actual_edge=(X,m,q,ell,g)=(128,321,179,193,28)")
    print("Z_179=[36,142], so 179 does not divide b_64")
    print("orders=(89,3,21), all > D because D<=2 at X=128")
    print("Q8345_FAR_PHYSICAL_OBSTRUCTION=PASS")


if __name__ == "__main__":
    main()
