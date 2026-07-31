#!/usr/bin/env python3
"""Independent audit of the corrected higher-Cartier shell formula.

For M = a*p+s and d = q*p-v, this checks

    C_M(d) == sum_{kappa in floor(a/q)P}
                  c_a(q*kappa)c_s(-v*kappa)       (mod p)

under the safety condition

    s + v*floor(M/d) < p.

It also checks the exact unrestricted Freshman's-dream convolution and
the first p >= 5 counterexample to the unsafe fixed-P formula.
"""

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


def scale(vector, scalar):
    return tuple(scalar * coordinate for coordinate in vector)


def subtract(left, right):
    return tuple(x - y for x, y in zip(left, right))


def in_polytope(vector, dilation):
    if dilation < 0:
        return False
    x, y, z = vector
    return (
        -dilation <= x <= dilation
        and -dilation <= y <= dilation
        and -dilation <= z <= dilation
        and x - y <= dilation
        and x - z <= dilation
    )


@lru_cache(maxsize=None)
def c(exponent, vector):
    if not in_polytope(vector, exponent):
        return 0
    return coefficient(exponent, *vector)


def corrected_packet(p, a, q, s, v):
    dilation = a // q
    return sum(
        c(a, scale(kappa, q)) * c(s, scale(kappa, -v))
        for kappa in polytope_points(dilation)
    ) % p


def full_freshman_convolution(p, a, q, s, v):
    """Exact mod-p shell after Freshman's dream, with every alias retained."""

    moment = a * p + s
    node = q * p - v
    shell_dilation = moment // node
    # beta = p*delta-v*kappa lies in sP, hence each delta coordinate
    # has absolute value at most floor((s+v*t)/p).
    radius = (s + v * shell_dilation) // p
    total = 0
    for kappa in polytope_points(shell_dilation):
        for delta in product(range(-radius, radius + 1), repeat=3):
            mu = subtract(scale(kappa, q), delta)
            beta = subtract(scale(delta, p), scale(kappa, v))
            if in_polytope(mu, a) and in_polytope(beta, s):
                total += c(a, mu) * c(s, beta)
    return total % p


def audit_support_domain():
    checks = 0
    for a in range(0, 13):
        for q in range(1, 13):
            m = a // q
            for kappa in polytope_points(12):
                assert in_polytope(scale(kappa, q), a) == in_polytope(
                    kappa, m
                )
                checks += 1
    return checks


def audit_general_formula():
    safe_checks = 0
    for p in primes_up_to(19):
        if p < 5:
            continue
        for a in range(1, 4):
            for s in range(p):
                moment = a * p + s
                nodes = []
                parameters = []
                for q in range(1, a + 2):
                    for v in range(0, min(4, q * p)):
                        node = q * p - v
                        if 1 <= node <= moment:
                            nodes.append(node)
                            parameters.append((q, v, node))
                values = shell_batch(moment, nodes, modulus=p)
                for q, v, node in parameters:
                    shell_dilation = moment // node
                    if s + v * shell_dilation < p:
                        # The safety inequality itself rules out every
                        # additional shell layer: t must equal floor(a/q).
                        assert shell_dilation == a // q
                        assert values[node] == corrected_packet(
                            p, a, q, s, v
                        )
                        safe_checks += 1

    # A few unrestricted checks retain every nonzero Cartier alias.  These
    # include both safe and unsafe cases and are compared to the original
    # shell, not to the corrected single-alias packet.
    unrestricted = (
        (5, 1, 1, 2, 2),
        (5, 2, 1, 4, 3),
        (7, 2, 2, 5, 2),
        (7, 3, 1, 6, 3),
        (11, 2, 1, 8, 2),
    )
    direct_checks = 0
    for p, a, q, s, v in unrestricted:
        moment = a * p + s
        node = q * p - v
        direct = shell_batch(moment, (node,), modulus=p)[node]
        assert direct == full_freshman_convolution(p, a, q, s, v)
        direct_checks += 1
    return direct_checks, safe_checks


def audit_boundary_failure():
    p, a, q, s = 5, 1, 1, 2
    moment = a * p + s
    values = shell_batch(moment, (3, 4, 5), modulus=p)
    assert (values[3], values[4], values[5]) == (3, 0, 0)
    fixed_p_pair = (
        corrected_packet(p, a, q, s, 2)
        + corrected_packet(p, a, q, s, 0)
    ) % p
    assert fixed_p_pair == 0
    assert (values[3] + values[5]) % p != fixed_p_pair


if __name__ == "__main__":
    domain_checks = audit_support_domain()
    direct_checks, safe_checks = audit_general_formula()
    audit_boundary_failure()
    print(
        "PASS:",
        domain_checks,
        "support-domain checks;",
        direct_checks,
        "direct shell/convolution checks;",
        safe_checks,
        "safe corrected-packet checks",
    )
