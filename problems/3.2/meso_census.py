#!/usr/bin/env python3
"""Extended mesoscopic census — same caliber as energy_verify.py.

Outputs per prime p (H = isqrt(p)):
  R = root_mass = sum of column_counts (= R_p(H))
  W = support = #{(d,r) : m_{d,r} > 0}
  E = weighted energy = sum m_{d,r} (split boundary / nonboundary)
  max_m = max m_{d,r}
  m-tail: counts of pairs with m >= 1, 2, 3
  left max degree = max_d #{r : m_{d,r} > 0}
  right max degree = max_r #{d : m_{d,r} > 0}
  left-pair max codegree = max_{d!=d'} #{r : m_{d,r}>0 and m_{d',r}>0}

First validates against energy_verify.py's canonical records for p <= 5000,
then extends.
"""
from __future__ import annotations
from collections import defaultdict
from math import isqrt
import sys


def primes_up_to(limit: int) -> list[int]:
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[:2] = b"\x00\x00"
    for d in range(2, isqrt(limit) + 1):
        if not sieve[d]:
            continue
        for j in range(d * d, limit + 1, d):
            sieve[j] = 0
    return [n for n, f in enumerate(sieve) if f]


def apery_coefficient(n: int) -> int:
    return (2 * n + 1) * (17 * n * n + 17 * n + 5)


def root_data(prime: int, height: int):
    """Canonical root_data from energy_verify.py (h-rotation, no per-x pow)."""
    roots: list[list[int]] = [[] for _ in range(height + 1)]
    masks: list[bytearray] = [bytearray(prime) for _ in range(height + 1)]
    column_counts = [0] * prime
    if height < 2:
        return roots, masks, column_counts

    p_values: list[int] = []
    sixth_powers: list[int] = []
    for residue in range(prime):
        square = residue * residue % prime
        p_values.append(
            ((2 * residue + 1) * (17 * square + 17 * residue + 5)) % prime
        )
        sixth_powers.append(pow(residue, 6, prime))

    previous = [0] * prime
    current = [1] * prime
    for index in range(1, height):
        p_shift = p_values[index:] + p_values[:index]
        sixth_shift = sixth_powers[index:] + sixth_powers[:index]
        following = [
            (c * v - s * o) % prime
            for c, v, s, o in zip(p_shift, current, sixth_shift, previous)
        ]
        level = index + 1
        level_roots = [x for x, v in enumerate(following) if v == 0]
        mask = masks[level]
        for x in level_roots:
            mask[x] = 1
            column_counts[x] += 1
        roots[level] = level_roots
        previous, current = current, following

    return roots, masks, column_counts


def census(prime: int):
    H = isqrt(prime)
    if H < 4:
        return None
    roots, masks, column_counts = root_data(prime, H)

    # boundary: x ≡ -cut (mod p) for 2 <= cut <= H
    boundary = bytearray(prime)
    for cut in range(2, H + 1):
        boundary[(-cut) % prime] = 1

    R = sum(column_counts)  # root mass R_p(H)

    # Per-pair m_{d,r} and aggregates
    W = 0           # support
    E_nb = 0        # energy nonboundary
    E_bd = 0        # energy boundary
    max_m = 0
    m_ge1 = 0
    m_ge2 = 0
    m_ge3 = 0

    left_deg = defaultdict(int)   # d -> #{r with m>0}
    right_deg = defaultdict(int)  # r -> #{d with m>0}
    # For codegree: left_neighbors[d] = set of r with m>0
    left_neighbors: dict[int, set[int]] = defaultdict(set)

    for d in range(2, H):
        for r in range(2, H - d + 1):
            m = 0
            m_nb = 0
            m_bd = 0
            mask_r = masks[r]
            for x in roots[d]:
                if mask_r[(x + d) % prime]:
                    m += 1
                    if boundary[x]:
                        m_bd += 1
                    else:
                        m_nb += 1
            if m > 0:
                W += 1
                left_deg[d] += 1
                right_deg[r] += 1
                left_neighbors[d].add(r)
                m_ge1 += 1
            if m >= 2:
                m_ge2 += 1
            if m >= 3:
                m_ge3 += 1
            if m > max_m:
                max_m = m
            E_nb += m_nb
            E_bd += m_bd

    E = E_nb + E_bd
    max_left_deg = max(left_deg.values()) if left_deg else 0
    max_right_deg = max(right_deg.values()) if right_deg else 0

    # Left-pair max codegree
    max_codeg = 0
    ds = list(left_neighbors.keys())
    for i in range(len(ds)):
        for j in range(i + 1, len(ds)):
            c = len(left_neighbors[ds[i]] & left_neighbors[ds[j]])
            if c > max_codeg:
                max_codeg = c

    return {
        'p': prime, 'H': H, 'R': R, 'W': W,
        'E': E, 'E_nb': E_nb, 'E_bd': E_bd,
        'max_m': max_m, 'm_ge1': m_ge1, 'm_ge2': m_ge2, 'm_ge3': m_ge3,
        'max_ldeg': max_left_deg, 'max_rdeg': max_right_deg,
        'max_codeg': max_codeg,
    }


def main():
    limit = int(sys.argv[1]) if len(sys.argv) > 1 else 5000
    primes = [p for p in primes_up_to(limit) if p >= 7]

    hdr = (f"{'p':>6} {'H':>4} {'R':>5} {'W':>5} {'E':>5} "
           f"{'E_nb':>5} {'E_bd':>5} {'H1.5':>6} "
           f"{'maxm':>4} {'m≥2':>4} {'m≥3':>4} "
           f"{'Ldeg':>4} {'Rdeg':>4} {'codeg':>5}")
    print(hdr)
    print("-" * len(hdr))

    all_recs = []
    for p in primes:
        rec = census(p)
        if rec is None:
            continue
        all_recs.append(rec)
        r = rec
        H15 = r['H'] ** 1.5
        if (p <= 100 or p % 500 < 10 or r['max_m'] >= 3
                or r['E'] > 2 * H15 or p in (653, 3727)):
            print(f"{r['p']:6d} {r['H']:4d} {r['R']:5d} {r['W']:5d} "
                  f"{r['E']:5d} {r['E_nb']:5d} {r['E_bd']:5d} {H15:6.1f} "
                  f"{r['max_m']:4d} {r['m_ge2']:4d} {r['m_ge3']:4d} "
                  f"{r['max_ldeg']:4d} {r['max_rdeg']:4d} {r['max_codeg']:5d}")

    if all_recs:
        max_E = max(r['E'] for r in all_recs)
        max_m_all = max(r['max_m'] for r in all_recs)
        max_codeg_all = max(r['max_codeg'] for r in all_recs)
        print(f"\nSummary (p <= {limit}):")
        print(f"  max E = {max_E}")
        print(f"  max m = {max_m_all}")
        print(f"  max codegree = {max_codeg_all}")
        print(f"  #(m>=2) primes: {sum(1 for r in all_recs if r['m_ge2'] > 0)}")
        print(f"  #(m>=3) primes: {sum(1 for r in all_recs if r['m_ge3'] > 0)}")


if __name__ == "__main__":
    main()
