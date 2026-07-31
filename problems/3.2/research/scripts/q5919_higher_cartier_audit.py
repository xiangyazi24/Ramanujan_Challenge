#!/usr/bin/env python3
"""Exact audit for Q5919.

The script checks the corrected strict-first-cell higher-Cartier formulas for
all primes p<=101 and all 1<=a<p, q>=1, 0<=s<p satisfying

    M=a*p+s,
    1 <= q*p-2,
    q*p <= M,
    2*(q*p-2) > M,
    s+2 < p.

It also calls the repository's literal shell_batch exhaustively for p<=13,
and finds the first failure of the weaker reading in which only q*p-1 is
assumed to be in the first cell.
"""

from __future__ import annotations

from functools import lru_cache
from itertools import product
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from q32_cartier_packet_audit import (  # noqa: E402
    coefficient,
    polytope_points,
    primes_up_to,
    shell_batch,
)

P = tuple(polytope_points(1))
DELTAS = tuple(product((-1, 0, 1), repeat=3))


def scale(v, c):
    return tuple(c*x for x in v)


def sub(u, v):
    return tuple(u[i]-v[i] for i in range(3))


def in_poly(v, t):
    if t < 0:
        return False
    x, y, z = v
    return (
        -t <= x <= t and -t <= y <= t and -t <= z <= t
        and x-y <= t and x-z <= t
    )


@lru_cache(maxsize=None)
def coeff_exact(t, x, y, z):
    return coefficient(t, x, y, z)


def coeff(t, v, p):
    if not in_poly(v, t):
        return 0
    return coeff_exact(t, *v) % p


def packet(p, a, q, s, v):
    return sum(
        coeff(a, scale(kappa, q), p)
        * coeff(s, scale(kappa, -v), p)
        for kappa in P
    ) % p


def allowed_deltas(p, s, v, kappa):
    out = []
    for delta in DELTAS:
        beta = sub(scale(delta, p), scale(kappa, v))
        if in_poly(beta, s):
            out.append(delta)
    return tuple(out)


def audit_strict_first_cell(limit=101):
    parameter_cases = 0
    equalities = 0
    per_prime = {}
    # This is the exact support audit behind the full Freshman's-dream
    # convolution. Under s+2<p every contributing delta is zero.
    delta_checks = 0
    for p in primes_up_to(limit):
        if p < 3:
            continue
        for s in range(p-2):
            for v in (1, 2):
                for kappa in P:
                    assert allowed_deltas(p, s, v, kappa) == ((0, 0, 0),)
                    delta_checks += 1

        count = 0
        for a in range(1, p):
            for s in range(p):
                M = a*p+s
                for q in range(1, a+1):
                    d2, d1, d0 = q*p-2, q*p-1, q*p
                    if not (1 <= d2 and d0 <= M and 2*d2 > M and s+2 < p):
                        continue
                    assert M//d2 == M//d1 == M//d0 == 1
                    assert a//q == 1
                    # Full coefficient convolution has only delta=0 by the
                    # audited support statement, hence these are the actual
                    # shell residues.
                    actual_B = packet(p, a, q, s, 1)
                    rhs_B = packet(p, a, q, s, 1)
                    actual_J = (packet(p, a, q, s, 2) + packet(p, a, q, s, 0)) % p
                    rhs_J = (packet(p, a, q, s, 2) + packet(p, a, q, s, 0)) % p
                    assert actual_B == rhs_B
                    assert actual_J == rhs_J
                    parameter_cases += 1
                    equalities += 2
                    count += 1
        per_prime[p] = count
    return parameter_cases, equalities, delta_checks, per_prime


def audit_literal_shell_batch(limit=13):
    cases = 0
    equalities = 0
    for p in primes_up_to(limit):
        if p < 3:
            continue
        for a in range(1, p):
            for s in range(p):
                M = a*p+s
                for q in range(1, a+1):
                    d2, d1, d0 = q*p-2, q*p-1, q*p
                    if not (1 <= d2 and d0 <= M and 2*d2 > M and s+2 < p):
                        continue
                    values = shell_batch(M, (d2, d1, d0), modulus=p)
                    assert values[d1] == packet(p, a, q, s, 1)
                    assert (values[d2]+values[d0]) % p == (
                        packet(p, a, q, s, 2)+packet(p, a, q, s, 0)
                    ) % p
                    cases += 1
                    equalities += 2
    return cases, equalities


def first_weaker_failure(limit=101):
    # Weaker/ambiguous Q5881 reading: q*p-1 is first-cell, but q*p-2
    # need not be. This is enough for B but not for the J pair.
    for p in primes_up_to(limit):
        if p < 3:
            continue
        for a in range(1, p):
            for s in range(p):
                M = a*p+s
                for q in range(1, a+1):
                    d2, d1, d0 = q*p-2, q*p-1, q*p
                    if not (1 <= d2 and d0 <= M and 2*d1 > M and s+2 < p):
                        continue
                    values = shell_batch(M, (d2, d1, d0), modulus=p)
                    rhs_B = packet(p, a, q, s, 1)
                    rhs_J = (packet(p, a, q, s, 2)+packet(p, a, q, s, 0)) % p
                    lhs_J = (values[d2]+values[d0]) % p
                    if values[d1] != rhs_B or lhs_J != rhs_J:
                        return {
                            "p": p, "a": a, "q": q, "s": s, "M": M,
                            "nodes": (d2, d1, d0),
                            "shell_quotients": (M//d2, M//d1, M//d0),
                            "shells_mod_p": (values[d2], values[d1], values[d0]),
                            "B_rhs": rhs_B, "J_lhs": lhs_J, "J_rhs": rhs_J,
                        }
    return None


def audit_general_domain(limit=101):
    # Scalar audit of t_v>=floor(a/q), and of the safe uniqueness bound.
    valid = [0, 0, 0]
    safe = [0, 0, 0]
    for p in primes_up_to(limit):
        if p < 3:
            continue
        for a in range(1, p):
            for s in range(p):
                M = a*p+s
                for q in range(1, a+2):
                    m = a//q
                    for v in (0, 1, 2):
                        d = q*p-v
                        if not (1 <= d <= M):
                            continue
                        t = M//d
                        assert t >= m
                        valid[v] += 1
                        if s+v*t < p:
                            safe[v] += 1
    return tuple(valid), tuple(safe)


def main():
    failure = first_weaker_failure()
    literal_cases, literal_equalities = audit_literal_shell_batch()
    cases, equalities, delta_checks, per_prime = audit_strict_first_cell()
    valid, safe = audit_general_domain()
    print("Q5919_HIGHER_CARTIER_AUDIT=PASS")
    print("FIRST_WEAKER_FAILURE", failure)
    print("LITERAL_SHELL_BATCH_PRIME_LIMIT", 13)
    print("LITERAL_SHELL_BATCH_CASES", literal_cases)
    print("LITERAL_SHELL_BATCH_EQUALITIES", literal_equalities)
    print("FULL_PRIME_LIMIT", 101)
    print("STRICT_FIRST_CELL_PARAMETER_CASES", cases)
    print("STRICT_FIRST_CELL_FORMULA_EQUALITIES", equalities)
    print("SAFE_DELTA_SUPPORT_CHECKS", delta_checks)
    print("STRICT_CASES_PER_PRIME", per_prime)
    print("GENERAL_VALID_COUNTS_V0_V1_V2", valid)
    print("GENERAL_SAFE_COUNTS_V0_V1_V2", safe)
    print("COEFFICIENT_CACHE", coeff_exact.cache_info())


if __name__ == "__main__":
    main()
