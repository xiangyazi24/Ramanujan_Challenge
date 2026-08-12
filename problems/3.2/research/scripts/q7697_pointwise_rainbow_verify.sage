#!/usr/bin/env sage
"""Exact verifier for q7697_pointwise_rainbow_correlation.sage.

The verifier checks more than the local zero sets:

  * three banked exact Apéry zero sets;
  * every pair intersection J_{p,q}(M) against a fresh CRT enumeration;
  * the exact centered-variance identity from those pair intersections;
  * F_2 against the independently enumerated pair tuples;
  * on the smallest case, F_3 against an independent triple-CRT enumeration.

Run:

    sage problems/3.2/research/scripts/q7697_pointwise_rainbow_verify.sage
"""

from sage.all import Integer, QQ, gcd, inverse_mod
from itertools import combinations
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
EXPERIMENT = HERE / "q7697_pointwise_rainbow_correlation.sage"
REFERENCE = HERE / "q7697_pointwise_rainbow_reference.json"

# Sage load() evaluates in this namespace.  Suppress the CLI entry point while
# importing the exact experiment functions.
_Q7697_LIBRARY_LOAD = True
load(str(EXPERIMENT))
del _Q7697_LIBRARY_LOAD


def crt2(a, m, b, n):
    """Canonical CRT representative in [0,mn), implemented independently."""
    a = Integer(a)
    b = Integer(b)
    m = Integer(m)
    n = Integer(n)
    if gcd(m, n) != 1:
        raise ValueError("CRT moduli must be coprime")
    step = ((b - a) * inverse_mod(m, n)) % n
    return (a + m*step) % (m*n)


def crt3(a, m, b, n, c, ell):
    first = crt2(a, m, b, n)
    return crt2(first, Integer(m)*Integer(n), c, ell)


def J_by_crt(state, p, q):
    M = state["M"]
    total = Integer(0)
    for r in state["zeros"][p]:
        for s in state["zeros"][q]:
            if crt2(r, p, s, q) < M:
                total += 1
    return total


def triple_count_by_crt(state):
    M = state["M"]
    total = Integer(0)
    primes = state["primes"]
    for p, q, ell in combinations(primes, 3):
        for r in state["zeros"][p]:
            for s in state["zeros"][q]:
                for t in state["zeros"][ell]:
                    if crt3(r, p, s, q, t, ell) < M:
                        total += 1
    return total


def verify_case(P, L, verify_triples):
    state = build_state(P, L)
    P = state["P"]
    M = state["M"]
    primes = state["primes"]
    A = state["A"]
    K = state["K"]

    # The early-representative interpretation requires the interval to be
    # shorter than every pair modulus in the block.
    if primes and not all(M < p*q for p, q in combinations(primes, 2)):
        raise AssertionError("test case left the early-CRT regime")

    pair_total = Integer(0)
    pair_delta_sum = QQ(0)
    for p, q in combinations(primes, 2):
        J_crt = J_by_crt(state, p, q)
        J_positions = Integer(len(set(state["positions"][p]).intersection(
            state["positions"][q]
        )))
        if J_crt != J_positions:
            raise AssertionError(
                "pair CRT mismatch at P=%s,L=%s,p=%s,q=%s" % (P, L, p, q)
            )
        pair_total += J_crt
        pair_delta_sum += QQ(J_crt) - QQ(A[p]*A[q], M)

    F2 = Integer(sum(k*(k-1) for k in K))
    if F2 != 2*pair_total:
        raise AssertionError("F2 does not equal twice the unordered pair count")

    mu = QQ(sum(K), M)
    variance_direct = sum(((QQ(k)-mu)**2 for k in K), QQ(0))
    diagonal = sum((QQ(A[p]) - QQ(A[p]*A[p], M) for p in primes), QQ(0))
    variance_crt = diagonal + 2*pair_delta_sum
    if variance_direct != variance_crt:
        raise AssertionError("variance identity failed")

    F3 = Integer(sum(k*(k-1)*(k-2) for k in K))
    triple_count = None
    if verify_triples:
        triple_count = triple_count_by_crt(state)
        if F3 != 6*triple_count:
            raise AssertionError("F3 does not equal six times unordered triples")

    return {
        "P": int(P),
        "L": int(L),
        "M": int(M),
        "pair_count": int(pair_total),
        "F2_ordered": int(F2),
        "F3_ordered": int(F3),
        "triple_count": None if triple_count is None else int(triple_count),
        "variance_num": str(variance_direct.numerator()),
        "variance_den": str(variance_direct.denominator()),
    }


def main_verify():
    reference = json.loads(REFERENCE.read_text(encoding="utf-8"))
    if reference.get("schema") != 1:
        raise AssertionError("unexpected reference schema")

    for p_text, expected in reference["known_zero_sets"].items():
        p = Integer(p_text)
        actual = list(apery_zero_set(p))
        if actual != expected:
            raise AssertionError(
                "known zero set mismatch at p=%s: %s != %s" % (p, actual, expected)
            )

    results = []
    for case in reference["cases"]:
        results.append(verify_case(
            case["P"], case["L"], case["verify_triples_by_crt"]
        ))

    print(json.dumps({"status": "PASS", "cases": results}, indent=2,
                     sort_keys=True))


if __name__ == "__main__":
    main_verify()
