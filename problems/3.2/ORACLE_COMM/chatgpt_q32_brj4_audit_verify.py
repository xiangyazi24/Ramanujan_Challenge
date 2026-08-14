#!/usr/bin/env python3
"""Exact stdlib verifier for the Q8369 BRJ4 hostile audit.

This script uses only the Python standard library.  It has four deliberately
separate roles:

THEOREM SELF-TESTS
    Check the algebraic BRJ/Vandermonde rewrites and the fixed-g direct/
    reflected geometry identities on finite integer examples.  These examples
    are *formal algebra tests only* and are not presented as Apéry data or as
    counterexamples.

FINITE EVIDENCE
    Recompute the integer Apéry numbers from their three-term recurrence,
    enumerate all defining pairs p | b_h for 5 <= p <= --max-prime, group them
    by the direct-chart common value m=p+h, and verify BRJ by exact integer
    arithmetic.  With the default bound 500 the script pins the census quoted
    in chatgpt_q32_brj4_audit.md.

CONDITIONAL TARGET
    Print the exact [WFQA16] residue alignment that would convert BRJ into a
    target-selective scalar divisible by R.

NO-GO SCOPE
    Exhaustively check over several small finite fields that when the chart
    coefficient c is a unit, z is freely and uniquely solvable from
    V*B-c*z-Gamma=0.  The p=5/c=5 exceptional degeneration is tested
    separately.

No third-party package, network access, randomness, or synthetic Apéry values
are used.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from math import gcd, prod
from typing import Dict, Iterable, List, Sequence, Tuple


Pair = Tuple[int, int, int]  # (p, h, m=p+h)
ShortPair = Tuple[int, int]  # (p, h)
EvidenceRow = Tuple[int, int, int, int, int, int, int, int]


EXPECTED_MULTI_500: Dict[int, Tuple[ShortPair, ...]] = {
    200: ((139, 61), (181, 19)),
    272: ((191, 81), (233, 39)),
    300: ((191, 109), (227, 73)),
    321: ((179, 142), (193, 128), (211, 110)),
}

EXPECTED_ROWS_500 = {
    (200, 139, 61, 17, 29, 13, 42, 0),
    (200, 181, 19, 170, 149, 79, 139, 0),
    (321, 179, 142, 33, 97, 158, 90, 0),
    (321, 193, 128, 187, 178, 43, 134, 0),
    (321, 211, 110, 118, 44, 17, 154, 0),
    (272, 191, 81, 47, 134, 158, 42, 0),
    (272, 233, 39, 195, 13, 133, 191, 0),
    (300, 191, 109, 169, 26, 33, 36, 0),
    (300, 227, 73, 186, 4, 94, 191, 0),
}


def apery_numbers(limit: int) -> List[int]:
    """Return [b_0,...,b_limit] by exact division in the Apéry recurrence."""
    if limit < 0:
        raise ValueError("limit must be nonnegative")
    if limit == 0:
        return [1]

    b = [1, 5]
    for n in range(1, limit):
        coefficient = 34 * n**3 + 51 * n**2 + 27 * n + 5
        numerator = coefficient * b[n] - n**3 * b[n - 1]
        denominator = (n + 1) ** 3
        assert numerator % denominator == 0, (n, numerator, denominator)
        b.append(numerator // denominator)
    return b


def primes_up_to(n: int) -> List[int]:
    if n < 2:
        return []
    sieve = [True] * (n + 1)
    sieve[0] = sieve[1] = False
    p = 2
    while p * p <= n:
        if sieve[p]:
            for q in range(p * p, n + 1, p):
                sieve[q] = False
        p += 1
    return [p for p in range(2, n + 1) if sieve[p]]


def defining_pairs(max_prime: int, b: Sequence[int]) -> List[Pair]:
    """Actual direct-chart pairs (p,h,m) with prime p and p | b_h."""
    out: List[Pair] = []
    for p in primes_up_to(max_prime):
        if p < 5:
            continue
        for h in range(p):
            if b[h] % p != 0:
                continue
            m = p + h
            # Direct Apéry-Lucas coefficient b_1=5; this is checked, not assumed
            # in the quotient arithmetic below.
            assert (b[m] - 5 * b[h]) % p == 0, (p, h, m)
            out.append((p, h, m))
    return out


def group_by_m(pairs: Iterable[Pair]) -> Dict[int, List[Pair]]:
    groups: Dict[int, List[Pair]] = defaultdict(list)
    for pair in pairs:
        groups[pair[2]].append(pair)
    for group in groups.values():
        group.sort()
    return dict(groups)


def audit_direct_group(m: int, group: Sequence[Pair], b: Sequence[int]) -> List[EvidenceRow]:
    """Check exact direct-chart BRJ and return residues for one common-m group."""
    primes = [p for p, _, _ in group]
    assert len(primes) == len(set(primes)), (m, group)

    R = prod(primes)
    assert b[m] % R == 0, (m, R)
    B = b[m] // R

    rows: List[EvidenceRow] = []
    for p, h, mm in group:
        assert mm == m
        assert b[h] % p == 0
        z = b[h] // p

        gamma_numerator = b[m] - 5 * b[h]
        assert gamma_numerator % p == 0
        gamma = (gamma_numerator // p) % p

        V = prod(q - p for q, _, _ in group if q != p)
        # Since m=p+h=q+h_j on the direct chart, q-p=h-h_j exactly.
        V_from_nodes = prod(h - hh for q, hh, _ in group if q != p)
        assert V == V_from_nodes, (m, p, h, V, V_from_nodes)

        # Distinct prime factors make V a unit mod p (also true for singleton V=1).
        assert gcd(V, p) == 1

        # The key claim: R/p is merely rewritten as the oriented Vandermonde.
        assert (R // p - V) % p == 0, (m, p, R // p, V)

        brj = (B * V - 5 * z - gamma) % p
        assert brj == 0, (m, p, h, brj)

        rows.append((m, p, h, B % p, z % p, gamma, V % p, brj))
    return rows


def finite_evidence(max_prime: int) -> Tuple[List[Pair], Dict[int, Tuple[ShortPair, ...]], List[EvidenceRow]]:
    # m=p+h <= 2*max_prime-1, so this bound is safely sufficient.
    b = apery_numbers(2 * max_prime)
    pairs = defining_pairs(max_prime, b)
    groups = group_by_m(pairs)

    all_rows: List[EvidenceRow] = []
    for m in sorted(groups):
        rows = audit_direct_group(m, groups[m], b)
        if len(groups[m]) >= 2:
            all_rows.extend(rows)

    multi = {
        m: tuple((p, h) for p, h, _ in groups[m])
        for m in sorted(groups)
        if len(groups[m]) >= 2
    }

    if max_prime == 500:
        assert len(pairs) == 95, len(pairs)
        assert multi == EXPECTED_MULTI_500, (multi, EXPECTED_MULTI_500)
        assert set(all_rows) == EXPECTED_ROWS_500, set(all_rows) ^ EXPECTED_ROWS_500

    return pairs, multi, all_rows


def verify_unit_linear_freeness() -> int:
    """Finite-field self-test of z=c^{-1}(VB-Gamma) when c is a unit."""
    cases = 0
    for p in (7, 11, 13, 17):
        c = 5 % p
        c_inv = pow(c, -1, p)
        for V in range(p):
            for B in range(p):
                for gamma in range(p):
                    z = ((V * B - gamma) * c_inv) % p
                    assert (V * B - c * z - gamma) % p == 0
                    cases += 1
    return cases


def verify_p5_exception() -> None:
    """At p=5 the coefficient of z vanishes, so z is not eliminated."""
    p = 5
    c = 5 % p
    assert c == 0
    V = 1
    B = 0
    gamma = 0
    residues = {(V * B - c * z - gamma) % p for z in range(p)}
    assert residues == {0}


def vandermonde(primes: Sequence[int], index: int) -> int:
    p = primes[index]
    return prod(q - p for j, q in enumerate(primes) if j != index)


def verify_fixed_g_geometry() -> None:
    """Formal integer tests of the direct/reflected oriented geometry identities."""
    # Formal algebra tests only: these are not asserted to be Apéry packets.
    g = 2
    m = 100

    direct_h = (10, 20, 30)
    direct_p = tuple((m - h) // g for h in direct_h)
    assert all(m - h == g * p for p, h in zip(direct_p, direct_h))
    n = len(direct_h)
    for i, h in enumerate(direct_h):
        V = vandermonde(direct_p, i)
        node_v = prod(h - hh for j, hh in enumerate(direct_h) if j != i)
        assert g ** (n - 1) * V == node_v

    reflected_h = (1, 4, 7)
    # m=(g+1)p-1-h, so h == -m-1 (mod g+1).
    assert all((m + h + 1) % (g + 1) == 0 for h in reflected_h)
    reflected_p = tuple((m + h + 1) // (g + 1) for h in reflected_h)
    for i, h in enumerate(reflected_h):
        V = vandermonde(reflected_p, i)
        node_v = prod(hh - h for j, hh in enumerate(reflected_h) if j != i)
        assert (g + 1) ** (n - 1) * V == node_v


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--max-prime",
        type=int,
        default=500,
        help="largest prime in the actual Apéry finite census (default: 500)",
    )
    parser.add_argument(
        "--show-groups",
        action="store_true",
        help="print all common-m groups of size at least two",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.max_prime < 5:
        raise SystemExit("--max-prime must be at least 5")

    verify_fixed_g_geometry()
    free_cases = verify_unit_linear_freeness()
    verify_p5_exception()
    pairs, multi, rows = finite_evidence(args.max_prime)

    print("THEOREM")
    print("  PASS: R/p_i == prod_{j!=i}(p_j-p_i) (mod p_i) is checked inside every actual direct group.")
    print("  PASS: fixed-g direct and reflected oriented Vandermonde geometry identities pass exact integer tests.")
    print("  NOTE: formal fixed-g geometry tests are algebra self-tests, not synthetic Apéry counterexamples.")

    print("FINITE EVIDENCE")
    print(f"  prime bound: {args.max_prime}")
    print(f"  actual defining pairs p|b_h: {len(pairs)}")
    print(f"  common-m groups of size >=2: {len(multi)}")
    print(f"  largest common-m group: {max((len(v) for v in multi.values()), default=1)}")
    print(f"  BRJ rows checked in multi-node groups: {len(rows)}")
    if args.max_prime == 500:
        print("  PASS: exact 95-pair / 4-multigroup / 9-row pinned census matches the audit report.")
    if args.show_groups:
        for m, group in multi.items():
            print(f"  m={m}: {group}")
        for row in sorted(rows):
            print("  row", row)

    print("CONDITIONAL TARGET")
    print("  [WFQA16]: c_g*z_i + Gamma_i == V_i*Q_T (mod p_i) for all 16 packet nodes,")
    print("  with Q_T chosen from packet/face data independently of CRT representatives.")
    print("  BRJ then gives R | (B-Q_T); the additional height 0 < |B-Q_T| < R is contradictory.")
    print("  Zero target Q_T=0 gives R|B and hence R^2|b_m.")

    print("NO-GO SCOPE")
    print(f"  PASS: {free_cases} finite-field cases verify unique local z solvability when c=5 is a unit.")
    print("  PASS: p=5 is separately detected as the coefficient-degenerate exception.")
    print("  Scope: this does not cover any new p-adic relation that already constrains z_i or B,")
    print("  nor a pre-reduction integral/height identity coupling different primes.")


if __name__ == "__main__":
    main()
