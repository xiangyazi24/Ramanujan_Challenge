#!/usr/bin/env python3
"""Exact audits for Cartier zero segments in the first-cell ray state.

Write ``M=a*p+s`` with ``0<=s<p``.  In the first cell and for ``r<p``,
the freshman's dream gives the exact digit factorization

    c_M((M-r)*kappa)
      = c_a(a*kappa) c_s((s-r)*kappa) (mod p).

For every nonzero first-cell ray ``kappa`` and every

    2*s < r < p,

the second factor is outside ``s*P``, so

    c_M((M-r)*kappa) = 0 (mod p).

Consequently, if ``p | b_a*b_s``, then the complete first-cell shell
``C_M(M-r)`` vanishes on that segment.  This strictly extends the
cyclic-rotation case ``s=0`` and explains the long state-gcd factors
73 near M=146 and 61 at M=126.
"""

from math import gcd

from q32_cartier_packet_audit import (
    apery,
    coefficient,
    polytope_points,
    primes_up_to,
    shell_batch,
)


RAYS = tuple(point for point in polytope_points(1) if point != (0, 0, 0))
assert len(RAYS) == 21


def audit_digit_factorization():
    checks = 0
    for prime in primes_up_to(23):
        if prime < 5:
            continue
        for quotient in range(1, 4):
            for residue in range(prime):
                moment = quotient * prime + residue
                upper = min(prime - 1, (moment - 1) // 2)
                for r in range(upper + 1):
                    node = moment - r
                    for ray in ((0, 0, 0),) + RAYS:
                        left = coefficient(
                            moment,
                            node * ray[0],
                            node * ray[1],
                            node * ray[2],
                            prime,
                        )
                        right = (
                            coefficient(
                                quotient,
                                quotient * ray[0],
                                quotient * ray[1],
                                quotient * ray[2],
                                prime,
                            )
                            * coefficient(
                                residue,
                                (residue - r) * ray[0],
                                (residue - r) * ray[1],
                                (residue - r) * ray[2],
                                prime,
                            )
                        ) % prime
                        assert left == right
                        checks += 1
    print("DIGIT_FACTORIZATION_CHECKS", checks)


def audit_ray_vanishing():
    checks = 0
    for prime in primes_up_to(23):
        if prime < 5:
            continue
        for quotient in range(1, 4):
            for residue in range(prime):
                moment = quotient * prime + residue
                upper = min(prime - 1, (moment - 1) // 2)
                for r in range(2 * residue + 1, upper + 1):
                    node = moment - r
                    for ray in RAYS:
                        assert (
                            coefficient(
                                moment,
                                node * ray[0],
                                node * ray[1],
                                node * ray[2],
                                prime,
                            )
                            == 0
                        )
                        checks += 1
    print("RAY_ZERO_CHECKS", checks)


def audit_shell_segments():
    checks = 0
    for prime in primes_up_to(43):
        if prime < 5:
            continue
        for quotient in range(1, 4):
            for residue in range(prime):
                if apery(quotient) * apery(residue) % prime:
                    continue
                moment = quotient * prime + residue
                upper = min(prime - 1, (moment - 1) // 2)
                rs = list(range(2 * residue + 1, upper + 1))
                if not rs:
                    continue
                shells = shell_batch(
                    moment,
                    (moment - r for r in rs),
                    prime,
                )
                for r in rs:
                    assert shells[moment - r] == 0
                    checks += 1
    print("SHELL_ZERO_CHECKS", checks)


def audit_top_target_state_exclusion():
    """A top-half target cannot also divide the augmented origin state.

    For ``M=p+s`` and ``p>5``, Apéry--Lucas gives
    ``b_M = b_1*b_s = 5*b_s (mod p)``.  If the candidate node ``p-1``
    is a target, then ``p | b_{s+1}``.  Two consecutive Apéry numbers
    below ``p`` cannot both vanish, so ``b_M`` is a unit modulo ``p``.
    """

    checks = 0
    hits = 0
    for prime in primes_up_to(199):
        if prime <= 5:
            continue
        values = [apery(index) % prime for index in range(prime + 1)]
        for residue in range(prime - 1):
            if values[residue + 1]:
                continue
            hits += 1
            assert values[residue]
            assert apery(prime + residue) % prime == 5 * values[residue] % prime
            checks += 1
    print("TOP_TARGET_STATE_EXCLUSIONS", checks, "TARGET_HITS", hits)


def state_window_gcd(moment, start, width):
    b = apery(moment)
    shells = shell_batch(
        moment,
        (moment - r for r in range(start, start + width)),
    )
    out = b
    for r in range(start, start + width):
        out = gcd(out, abs(shells[moment - r] - b))
    return out


def audit_observed_factors():
    # p=73, M=2p+s: the predicted start is 2s+1.
    for residue in range(5):
        moment = 2 * 73 + residue
        start = 2 * residue + 1
        gcd_value = state_window_gcd(moment, start, 38)
        assert gcd_value % 73 == 0
        if start > 0:
            shell_before = shell_batch(
                moment,
                [moment - (start - 1)],
                73,
            )
            # The theorem does not assert a converse, but these five
            # examples have the sharp endpoint seen in the state scan.
            assert shell_before[moment - (start - 1)] != 0
        print("STATE_FACTOR", moment, 73, "start", start, "gcd", gcd_value)

    # M=2*61+4 and 61 | b_4: the predicted start is 9.
    gcd_value = state_window_gcd(126, 9, 38)
    assert gcd_value % 61 == 0
    print("STATE_FACTOR", 126, 61, "start", 9, "gcd", gcd_value)


def main():
    audit_digit_factorization()
    audit_ray_vanishing()
    audit_shell_segments()
    audit_top_target_state_exclusion()
    audit_observed_factors()
    print("PASS: Cartier digit factorization and zero-segment theorem")


if __name__ == "__main__":
    main()
