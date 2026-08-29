#!/usr/bin/env python3
"""Finite sanity checks for the corrected Q32 RR68 average repair.

This script checks only:

* exact rational exponent arithmetic;
* the central-layer identity Z^2 = Z(Z-1) + Z;
* the literal RR68av high-corner competitors, target, and residual;
* the elementary congruence consequences of
      p*k == -rho (mod b), gcd(k,b)=gcd(rho,b)=1;
* bounded numbers of residue-class representatives in the physical p- and
  k-intervals; and
* the fact that zero-label weights are counted once in sum_k z_k.

It is not an asymptotic proof.  The asymptotic proof is the interval-counting
and moment argument in research/proofs/Q32_RR68_AVERAGE_REPAIR.md.
"""

from __future__ import annotations

from fractions import Fraction
from math import gcd


def ceil_div(a: int, b: int) -> int:
    if a < 0 or b <= 0:
        raise ValueError("ceil_div expects a >= 0 and b > 0")
    return (a + b - 1) // b


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


def check_exponent_ledger() -> None:
    Q = Fraction

    # Exact-k packet and fixed-theta W(theta).
    assert Q(-3) + Q(1) + Q(8, 5) == Q(-2, 5)

    # max W times the positive L1 ledger.
    assert Q(-2, 5) + Q(4, 5) == Q(2, 5)

    # Mean host term: outer N^-1, r^-1, L1 N^(4/5), and the central
    # factorial-plus-diagonal square moment N^(6/5).
    assert Q(-1) + Q(-1) + Q(4, 5) + Q(6, 5) == Q(0)

    # Weighted Cauchy on the same central layer.
    assert (Q(3, 5) + Q(6, 5)) / 2 == Q(9, 10)

    # Centered host term: outer N^-1, sqrt(r), weighted Cauchy, and
    # sqrt(S2)=N^(1/5)*sqrt(M/B).
    assert Q(-1) + Q(1, 2) + Q(9, 10) + Q(1, 5) == Q(3, 5)

    # AF1 plus the actual pointwise 2/3 zero bound.
    assert Q(2, 3) + Q(11, 5) == Q(43, 15)

    # Literal RR68av opposite high-corner ledger.
    crossing = (Q(13, 5) - Q(12, 5)) / (Q(1) + Q(1, 2))
    assert crossing == Q(2, 15)

    h_minus = Q(13, 5) - crossing
    h_plus = Q(12, 5) + crossing / 2
    assert h_minus == h_plus == Q(37, 15)

    # The maximum of a decreasing and an increasing affine function is
    # minimized at their crossing.  The finite grid is only a sanity check.
    for numerator in range(-120, 121):
        ell = Q(numerator, 60)
        competitor = max(Q(13, 5) - ell, Q(12, 5) + ell / 2)
        assert competitor >= Q(37, 15)

    # Correct fully-unconditional target and exact residual.
    for sigma in (Q(0), Q(1, 100), Q(1, 80), Q(1, 40)):
        target = Q(11, 5) - 2 * sigma
        residual = Q(37, 15) - target
        assert residual == Q(4, 15) + 2 * sigma
        assert residual > 0

    assert Q(4, 15) + 2 * Q(1, 40) == Q(19, 60)

    # The repair itself is sigma-free; sigma enters only the downstream target.
    sigma = Q(1, 40)

    # The proved cubic bound is much weaker than either staged sufficient
    # premise at the ridge endpoint.
    assert Q(43, 15) > Q(29, 15) - 6 * sigma
    assert Q(43, 15) > Q(2) - 6 * sigma


def check_central_factorial_plus_diagonal() -> None:
    """Check the exact local decomposition of the weighted square moment."""

    weights = [Fraction(1, 2), Fraction(2, 3), Fraction(5, 7), Fraction(3, 4)]
    zero_counts = [0, 1, 2, 5]

    # Exhaust every deletion mask to model the exact central host layer.
    for mask in range(1 << len(weights)):
        square = Fraction(0)
        factorial = Fraction(0)
        diagonal = Fraction(0)

        for i, (weight, z) in enumerate(zip(weights, zero_counts)):
            if not (mask & (1 << i)):
                continue
            square += weight * z * z
            factorial += weight * z * (z - 1)
            diagonal += weight * z

        assert square == factorial + diagonal

    # Exponent-level consequence: N^(6/5) + N^(3/5) is N^(6/5+o(1)).
    assert max(Fraction(6, 5), Fraction(3, 5)) == Fraction(6, 5)


def check_fixed_k_prime_class() -> None:
    """For fixed k, all matching p lie in one class modulo b."""

    for b in range(2, 31):
        for rho in range(1, b):
            if gcd(rho, b) != 1:
                continue
            for k in range(-3 * b, 3 * b + 1):
                if gcd(k, b) != 1:
                    continue

                residue = (-rho * pow(k % b, -1, b)) % b
                for n in range(5, 61):
                    lo, hi = n, 2 * n
                    matches = [
                        p for p in range(lo, hi + 1)
                        if (p * k + rho) % b == 0
                    ]
                    assert all(p % b == residue for p in matches)
                    assert all(gcd(p, b) == 1 for p in matches)

                    # A prime subset cannot be larger than the integer class.
                    prime_matches = [p for p in matches if is_prime(p)]
                    assert len(prime_matches) <= len(matches)

                    interval_size = hi - lo + 1
                    assert len(matches) <= ceil_div(interval_size, b)

                    # RR68c common-k consequence: every pair p,p' is in the
                    # same residue class modulo b.
                    for p in matches:
                        for p_prime in matches:
                            assert (p - p_prime) % b == 0


def check_fixed_prime_k_lifts() -> None:
    """For fixed p, the physical k interval contains O(1) lifts."""

    physical_constant = 3
    for b in range(2, 41):
        lo_k = -physical_constant * b
        hi_k = physical_constant * b
        k_interval_size = hi_k - lo_k + 1
        lift_bound = ceil_div(k_interval_size, b)

        for rho in range(1, b):
            if gcd(rho, b) != 1:
                continue
            for p in range(2, 81):
                solutions = [
                    k for k in range(lo_k, hi_k + 1)
                    if gcd(k, b) == 1 and (p * k + rho) % b == 0
                ]

                if gcd(p, b) != 1:
                    # The congruence and gcd(rho,b)=1 force this branch empty.
                    assert not solutions
                    continue

                residue = (-rho * pow(p % b, -1, b)) % b
                assert all(k % b == residue for k in solutions)
                assert len(solutions) <= lift_bound


def check_weighted_overlap_once() -> None:
    """Check sum_k z_k <= (max prime-to-k overlap) * sum_p Z(p)."""

    physical_constant = 2
    for n in range(7, 36):
        primes = [p for p in range(n, 2 * n + 1) if is_prime(p)]
        for b in range(2, min(16, n + 1)):
            physical_ks = range(-physical_constant * b,
                                physical_constant * b + 1)
            for rho in range(1, b):
                if gcd(rho, b) != 1:
                    continue

                packet: dict[int, list[int]] = {}
                for k in physical_ks:
                    if gcd(k, b) != 1:
                        continue
                    packet[k] = [
                        p for p in primes if (p * k + rho) % b == 0
                    ]

                # A deterministic positive stand-in for the already-contained
                # zero-label multiplicity Z(p).  Its arithmetic value is
                # irrelevant to this incidence identity.
                z_of_p = {p: 1 + (p % 5) for p in primes}

                sum_k_z_k = sum(
                    sum(z_of_p[p] for p in p_list)
                    for p_list in packet.values()
                )
                overlaps = {
                    p: sum(p in p_list for p_list in packet.values())
                    for p in primes
                }
                max_overlap = max(overlaps.values(), default=0)
                sum_p_z_p = sum(z_of_p.values())

                assert sum_k_z_k == sum(
                    z_of_p[p] * overlaps[p] for p in primes
                )
                assert sum_k_z_k <= max_overlap * sum_p_z_p

                # The physical interval has length O(b), so max_overlap is
                # bounded independently of n and b for fixed constant.
                expected_bound = ceil_div(2 * physical_constant * b + 1, b)
                assert max_overlap <= expected_bound


def main() -> None:
    check_exponent_ledger()
    check_central_factorial_plus_diagonal()
    check_fixed_k_prime_class()
    check_fixed_prime_k_lifts()
    check_weighted_overlap_once()

    print("PASS: RR68 core exponent ledger")
    print("PASS: central factorial-plus-diagonal square identity")
    print("PASS: literal RR68av high-corner, target, and residual")
    print("PASS: fixed-k source-prime residue class and shell bound")
    print("PASS: fixed-prime physical-k lift bound")
    print("PASS: zero labels counted once in the exact-k overlap")
    print("NOTE: finite computation is not the asymptotic proof")


if __name__ == "__main__":
    main()
