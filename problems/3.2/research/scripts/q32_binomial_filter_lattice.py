#!/usr/bin/env python3
"""Search the full allowed integer-valued Euler-filter affine lattice.

For q=3 coefficients C_d and sign eps, define

    D_k = sum_{d>=k} eps^d binom(d,k) C_d.

For m below the least candidate prime, every polynomial

    P(x)=1+sum_{1<=k<=m} a_k binom(x,k)

satisfies P(ap)=1 mod p at the four local degrees ap.  The attainable
certificate values are therefore exactly

    D_0 + gcd(D_1,...,D_m) Z.

This script computes the least nonzero centered residue for every allowed m.
It exhausts the whole affine lattice, not only a named orthogonal family.
"""

from __future__ import annotations

from math import comb, gcd, log

from q32_fixed_q_content import truncation_coefficients
from q32_strehl_gcd import franel_numbers, primes_up_to


INDICES = (40, 60, 80, 100, 120, 160, 200, 240, 300, 400, 600, 800)


def centered_nonzero(value: int, modulus: int) -> int:
    residue = value % modulus
    return min(residue, modulus - residue) if residue else modulus


def main() -> None:
    limit = max(INDICES)
    franel = franel_numbers(limit)
    primes = primes_up_to(limit)

    for n in INDICES:
        coefficients = truncation_coefficients(n, 3, franel)
        candidates = [
            prime for prime in primes if divmod(n, prime)[0] == 3
        ]
        maximum_degree = min(candidates) - 1
        outputs = []

        for sign in (-1, 1):
            base = sum(
                sign**degree * coefficients[degree]
                for degree in range(n + 1)
            )
            moment_gcd = 0
            best = (log(abs(base)) / n, 0, abs(base).bit_length())
            for order in range(1, maximum_degree + 1):
                moment = sum(
                    sign**degree
                    * comb(degree, order)
                    * coefficients[degree]
                    for degree in range(order, n + 1)
                )
                moment_gcd = gcd(moment_gcd, moment)
                least_value = centered_nonzero(base, abs(moment_gcd))
                record = (
                    log(least_value) / n if least_value > 1 else 0.0,
                    order,
                    least_value.bit_length(),
                )
                if record < best:
                    best = record
            outputs.append(best)

        print(
            f"n={n} mmax={maximum_degree} "
            f"minus={outputs[0][0]:.9f}@{outputs[0][1]} "
            f"plus={outputs[1][0]:.9f}@{outputs[1][1]}"
        )


if __name__ == "__main__":
    main()
