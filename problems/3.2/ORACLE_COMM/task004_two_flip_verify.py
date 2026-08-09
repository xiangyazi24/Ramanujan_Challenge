#!/usr/bin/env python3
"""Exact verification of the two-flip CRT reciprocity identity."""

from __future__ import annotations

import argparse
import cmath
import itertools
import math


def crt3(r: int, p: int, s: int, q: int, t: int, ell: int) -> int:
    modulus = p * q * ell
    value = 0
    for residue, prime in ((r, p), (s, q), (t, ell)):
        cofactor = modulus // prime
        value += residue * cofactor * pow(cofactor, -1, prime)
    return value % modulus


def two_flip_integer(r: int, p: int, s: int, q: int, t: int, ell: int) -> int:
    """Integer whose exponential divided by p q ell is the claimed RHS."""

    return (
        s
        + (r - s) * q * ell * pow(q * ell, -1, p)
        + (t - s) * p * q * pow(p * q, -1, ell)
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p", type=int, default=5)
    parser.add_argument("--q", type=int, default=7)
    parser.add_argument("--ell", type=int, default=11)
    parser.add_argument("--k-bound", type=int, default=25)
    args = parser.parse_args()

    p, q, ell = args.p, args.q, args.ell
    modulus = p * q * ell
    if any(math.gcd(a, b) != 1 for a, b in ((p, q), (p, ell), (q, ell))):
        raise ValueError("the three moduli must be pairwise coprime")

    exact_checks = 0
    max_phase_error = 0.0
    for r, s, t in itertools.product(range(p), range(q), range(ell)):
        m = crt3(r, p, s, q, t, ell)
        flipped = two_flip_integer(r, p, s, q, t, ell)
        if (flipped - m) % modulus:
            raise AssertionError((r, s, t, m, flipped, modulus))
        for k in range(-args.k_bound, args.k_bound + 1):
            # Exact equality in R/Z follows from the preceding congruence.
            if (k * (flipped - m)) % modulus:
                raise AssertionError((k, r, s, t, m, flipped))
            lhs = cmath.exp(2j * math.pi * k * m / modulus)
            rhs = (
                cmath.exp(
                    2j
                    * math.pi
                    * (k * (r - s) * pow(q * ell, -1, p) % p)
                    / p
                )
                * cmath.exp(
                    2j
                    * math.pi
                    * (k * (t - s) * pow(p * q, -1, ell) % ell)
                    / ell
                )
                * cmath.exp(2j * math.pi * k * s / modulus)
            )
            max_phase_error = max(max_phase_error, abs(lhs - rhs))
            exact_checks += 1

    print(f"moduli=({p},{q},{ell}), product={modulus}")
    print(f"exact checks={exact_checks}")
    print(f"max floating phase error={max_phase_error:.3e}")


if __name__ == "__main__":
    main()
