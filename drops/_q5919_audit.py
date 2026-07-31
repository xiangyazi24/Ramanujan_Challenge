#!/usr/bin/env python3
"""Temporary exact audit for Q5919.

The literal repository shell_batch is used on an exhaustive low-prime
subrange and on the first-failure search.  The full p<=101 audit uses an
exact Freshman's-dream convolution of the same coefficient definition;
this avoids shell_batch's quadratic-in-M binomial-row construction.
"""

from __future__ import annotations

from functools import lru_cache
from itertools import product
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "problems" / "3.2" / "research" / "scripts"
sys.path.insert(0, str(SCRIPTS))

from q32_cartier_packet_audit import (  # noqa: E402
    apery,
    coefficient,
    polytope_points,
    primes_up_to,
    shell_batch,
)

P1 = tuple(polytope_points(1))
DELTAS = tuple(product((-1, 0, 1), repeat=3))


def scale(v, c):
    return tuple(c * x for x in v)


def add(u, v):
    return tuple(u[i] + v[i] for i in range(3))


def sub(u, v):
    return tuple(u[i] - v[i] for i in range(3))


def in_dilate(v, t):
    if t < 0:
        return False
    x, y, z = v
    return (
        -t <= x <= t
        and -t <= y <= t
        and -t <= z <= t
        and x - y <= t
        and x - z <= t
    )


@lru_cache(maxsize=None)
def coeff_exact(t, x, y, z):
    return coefficient(t, x, y, z)


def coeff(t, v, p):
    if not in_dilate(v, t):
        return 0
    return coeff_exact(t, *v) % p


def rhs_general(p, a, s, q, v):
    m = a // q
    total = 0
    for kappa in polytope_points(m):
        total += coeff(a, scale(kappa, q), p) * coeff(s, scale(kappa, -v), p)
    return total % p


def shell_by_full_cartier_convolution(p, a, s, q, v):
    """Exact C_(ap+s)(qp-v) mod p, independent of uniqueness.

    Since s<p and kappa lies in tP, any contributing
    delta=q*kappa-mu satisfies |delta_i|<=1 under the safe cases audited
    here, so DELTAS is complete.  We assert this scalar bound separately.
    """
    M = a * p + s
    d = q * p - v
    assert 1 <= d <= M
    t = M // d
    total = 0
    for kappa in polytope_points(t):
        target = scale(kappa, d)
        # p*mu+beta=target, delta=q*kappa-mu.
        for delta in DELTAS:
            mu = sub(scale(kappa, q), delta)
            beta = sub(target, scale(mu, p))
            if in_dilate(mu, a) and in_dilate(beta, s):
                total += coeff(a, mu, p) * coeff(s, beta, p)
    return total % p


def strict_triple_cases(limit=101):
    checks = 0
    parameter_cases = 0
    per_prime = {}
    for p in primes_up_to(limit):
        if p < 3:
            continue
        count = 0
        for a in range(1, p):
            for s in range(p):
                M = a * p + s
                for q in range(1, a + 1):
                    d2, d1, d0 = q * p - 2, q * p - 1, q * p
                    if not (1 <= d2 and d0 <= M and 2 * d2 > M and s + 2 < p):
                        continue
                    parameter_cases += 1
                    count += 1
                    assert M // d2 == M // d1 == M // d0 == 1
                    assert a // q == 1
                    b_actual = shell_by_full_cartier_convolution(p, a, s, q, 1)
                    b_rhs = rhs_general(p, a, s, q, 1)
                    assert b_actual == b_rhs, ("B", p, a, q, s, b_actual, b_rhs)
                    j_actual = (
                        shell_by_full_cartier_convolution(p, a, s, q, 2)
                        + shell_by_full_cartier_convolution(p, a, s, q, 0)
                    ) % p
                    j_rhs = (rhs_general(p, a, s, q, 2) + rhs_general(p, a, s, q, 0)) % p
                    assert j_actual == j_rhs, ("J", p, a, q, s, j_actual, j_rhs)
                    checks += 2
        per_prime[p] = count
    return parameter_cases, checks, per_prime


def literal_shell_batch_checks(limit=13):
    checks = 0
    cases = 0
    for p in primes_up_to(limit):
        if p < 3:
            continue
        for a in range(1, p):
            for s in range(p):
                M = a * p + s
                for q in range(1, a + 1):
                    d2, d1, d0 = q * p - 2, q * p - 1, q * p
                    if not (1 <= d2 and d0 <= M and 2 * d2 > M and s + 2 < p):
                        continue
                    values = shell_batch(M, (d2, d1, d0), modulus=p)
                    b_rhs = rhs_general(p, a, s, q, 1)
                    j_rhs = (rhs_general(p, a, s, q, 2) + rhs_general(p, a, s, q, 0)) % p
                    assert values[d1] == b_rhs
                    assert (values[d2] + values[d0]) % p == j_rhs
                    cases += 1
                    checks += 2
    return cases, checks


def first_central_only_failure(limit=101):
    """Search the weaker reading used in Q5881: only qp-1 is first-cell."""
    for p in primes_up_to(limit):
        if p < 3:
            continue
        for a in range(1, p):
            for s in range(p):
                M = a * p + s
                for q in range(1, a + 1):
                    d2, d1, d0 = q * p - 2, q * p - 1, q * p
                    if not (1 <= d2 and d0 <= M and 2 * d1 > M and s + 2 < p):
                        continue
                    values = shell_batch(M, (d2, d1, d0), modulus=p)
                    b_rhs = rhs_general(p, a, s, q, 1)
                    j_rhs = (rhs_general(p, a, s, q, 2) + rhs_general(p, a, s, q, 0)) % p
                    if values[d1] != b_rhs or (values[d2] + values[d0]) % p != j_rhs:
                        return {
                            "p": p,
                            "a": a,
                            "q": q,
                            "s": s,
                            "M": M,
                            "nodes": (d2, d1, d0),
                            "quotients": (M // d2, M // d1, M // d0),
                            "shells_mod_p": (values[d2], values[d1], values[d0]),
                            "B_rhs": b_rhs,
                            "J_lhs": (values[d2] + values[d0]) % p,
                            "J_rhs": j_rhs,
                        }
    return None


def general_safe_domain_checks(limit=101):
    """Exhaust scalar domain and uniqueness inequalities for v=0,1,2."""
    safe = [0, 0, 0]
    domain = [0, 0, 0]
    for p in primes_up_to(limit):
        if p < 3:
            continue
        for a in range(1, p):
            for s in range(p):
                M = a * p + s
                for q in range(1, a + 2):
                    m = a // q
                    for v in (0, 1, 2):
                        d = q * p - v
                        if not (1 <= d <= M):
                            continue
                        t = M // d
                        assert t >= m
                        domain[v] += 1
                        if s + v * t < p:
                            safe[v] += 1
                            # Exhaust the only possible small deltas for every shell ray.
                            for kappa in polytope_points(t):
                                for delta in DELTAS:
                                    mu = sub(scale(kappa, q), delta)
                                    beta = sub(scale(kappa, d), scale(mu, p))
                                    if in_dilate(mu, a) and in_dilate(beta, s):
                                        assert delta == (0, 0, 0), (p, a, q, s, v, t, kappa, delta)
                                        assert in_dilate(kappa, m)
    return tuple(domain), tuple(safe)


def main():
    failure = first_central_only_failure()
    literal_cases, literal_checks = literal_shell_batch_checks()
    parameter_cases, formula_checks, per_prime = strict_triple_cases()
    domain_counts, safe_counts = general_safe_domain_checks()
    print("Q5919_HIGHER_CARTIER_AUDIT=PASS")
    print("FIRST_CENTRAL_ONLY_FAILURE", failure)
    print("LITERAL_SHELL_BATCH_LIMIT", 13)
    print("LITERAL_SHELL_BATCH_CASES", literal_cases)
    print("LITERAL_SHELL_BATCH_EQUALITIES", literal_checks)
    print("FULL_PRIME_LIMIT", 101)
    print("STRICT_FIRST_CELL_PARAMETER_CASES", parameter_cases)
    print("STRICT_FIRST_CELL_FORMULA_EQUALITIES", formula_checks)
    print("STRICT_CASES_PER_PRIME", per_prime)
    print("GENERAL_DOMAIN_COUNTS_V0_V1_V2", domain_counts)
    print("GENERAL_SAFE_COUNTS_V0_V1_V2", safe_counts)
    print("COEFFICIENT_CACHE", coeff_exact.cache_info())


if __name__ == "__main__":
    main()
