#!/usr/bin/env python3
"""Audit the exact multiplicity ledger for Cartier state-prime descent.

For a first-cell block ``r=R,...,R+H-1``, write ``M=a*p+s``.  A prime
with

    p > R+H-1,  2*s < R,  p | b_a*b_s

divides every shell in the block.  This script checks the two exact
bookkeeping statements used in Section 63:

* all lower-residue-channel primes with fixed ``s`` divide
  ``gcd(M-s,b_s)``;
* all outer-digit-channel primes with fixed ``a`` divide ``rad(b_a)``;
* the quotient-one direct prefix and quotient-two diagonal families
  really occur among the state edges.

The ray and whole-shell congruences themselves are audited separately by
``q32_cartier_zero_segment_audit.py``.
"""

from collections import defaultdict
from functools import lru_cache
from math import comb, gcd, isqrt, prod


@lru_cache(maxsize=None)
def apery(n):
    return sum((comb(n, k) * comb(n + k, k)) ** 2 for k in range(n + 1))


@lru_cache(maxsize=None)
def apery_digit_mod(n, prime):
    assert 0 <= n < prime
    return sum(
        (comb(n, k) * comb(n + k, k)) ** 2
        for k in range(n + 1)
    ) % prime


def apery_mod(n, prime):
    """Apéry--Lucas evaluation from the base-prime digits."""

    out = 1
    while n:
        n, digit = divmod(n, prime)
        out = out * apery_digit_mod(digit, prime) % prime
    return out


def primes_up_to(limit):
    if limit < 2:
        return []
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[:2] = b"\x00\x00"
    for prime in range(2, isqrt(limit) + 1):
        if sieve[prime]:
            sieve[prime * prime : limit + 1 : prime] = b"\x00" * (
                (limit - prime * prime) // prime + 1
            )
    return [prime for prime in range(2, limit + 1) if sieve[prime]]


def zero_digit(n, prime):
    """Return one base-prime digit ``u`` for which ``prime | b_u``."""

    while n:
        n, digit = divmod(n, prime)
        if apery_digit_mod(digit, prime) == 0:
            return digit
    raise AssertionError((n, prime))


def state_edges(moment, start, width, primes):
    edges = []
    for prime in primes:
        if prime > moment:
            break
        quotient, residue = divmod(moment, prime)
        if not (
            prime > start + width - 1
            and 2 * residue < start
        ):
            continue
        zero_outer = apery_mod(quotient, prime) == 0
        zero_residue = apery_digit_mod(residue, prime) == 0
        if not (zero_outer or zero_residue):
            continue

        # Assign an overlap to the residue channel once, so every prime is
        # charged exactly once in the product ledger.
        if zero_residue:
            channel = "residue"
            child = residue
        else:
            channel = "outer"
            child = zero_digit(quotient, prime)
        assert apery_digit_mod(child, prime) == 0
        edges.append(
            (prime, quotient, residue, channel, child, prime + child)
        )
    return edges


def audit_multiplicity_ledger(limit=420):
    primes = primes_up_to(limit)
    edge_checks = 0
    residue_groups = 0
    outer_groups = 0
    direct_edges = 0
    diagonal_edges = 0

    for moment in range(40, limit + 1):
        width = max(1, round(moment ** (1 / 3)))
        start = moment // 2 - 2 * width
        if start <= 0 or start + width - 1 >= moment / 2:
            continue
        edges = state_edges(moment, start, width, primes)
        edge_set = {(edge[0], edge[1], edge[2]) for edge in edges}
        edge_checks += len(edges)

        by_residue = defaultdict(list)
        by_outer = defaultdict(list)
        for prime, quotient, residue, channel, child, child_row in edges:
            assert child_row == prime + child
            if channel == "residue":
                by_residue[residue].append(prime)
                if quotient >= 2:
                    assert 8 * child_row < 5 * moment + 8
            else:
                by_outer[quotient].append(prime)

        for residue, group in by_residue.items():
            radical = prod(group)
            assert (moment - residue) % radical == 0
            assert apery(residue) % radical == 0
            assert gcd(moment - residue, apery(residue)) % radical == 0
            residue_groups += 1

        for quotient, group in by_outer.items():
            radical = prod(group)
            assert apery(quotient) % radical == 0
            outer_groups += 1

        # Exact quotient-one fixed-point family.
        maximum_residue = (start - 1) // 2
        for residue in range(maximum_residue + 1):
            prime = moment - residue
            if prime not in primes or apery_digit_mod(residue, prime):
                continue
            assert (prime, 1, residue) in edge_set
            direct_edges += 1

        # Exact quotient-two diagonal family at the near-boundary scale.
        for prime in primes:
            if not (moment / 2 - width < prime <= moment / 2):
                continue
            residue = moment - 2 * prime
            if (
                apery_mod(2, prime)
                * apery_digit_mod(residue, prime)
                % prime
            ):
                continue
            assert (prime, 2, residue) in edge_set
            diagonal_edges += 1

    print("STATE_EDGE_CHECKS", edge_checks)
    print("RESIDUE_GROUP_CHECKS", residue_groups)
    print("OUTER_GROUP_CHECKS", outer_groups)
    print("DIRECT_FIXED_POINT_EDGES", direct_edges)
    print("QUOTIENT_TWO_DIAGONAL_EDGES", diagonal_edges)


def audit_displayed_examples():
    primes = primes_up_to(160)
    for moment in range(146, 151):
        quotient, residue = divmod(moment, 73)
        assert quotient == 2
        assert apery_mod(quotient, 73) == 0
        assert 2 * residue + 1 == 2 * (moment - 146) + 1

    quotient, residue = divmod(126, 61)
    assert (quotient, residue) == (2, 4)
    assert apery_digit_mod(residue, 61) == 0

    # Both examples occur in the state-edge scan for a suitable block.
    assert any(edge[0] == 73 for edge in state_edges(146, 1, 38, primes))
    assert any(edge[0] == 61 for edge in state_edges(126, 9, 38, primes))


def main():
    audit_displayed_examples()
    audit_multiplicity_ledger()
    print("PASS: Cartier state-prime descent ledger")


if __name__ == "__main__":
    main()
