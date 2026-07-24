#!/usr/bin/env python3
"""Audit growing integer-polynomial Euler-derivative functionals.

For P in Z[x] with P(0)=1, the values

    P(theta) T(1),  P(theta) T(-1),  theta=c*d/dc,

remain selective on the q=3 candidate interval.  Positive moments are,
however, universally divisible by every candidate prime, so all attainable
values lie in one coarse residue class.  This script computes the exact
moment ideal and its least nonzero residue for representative n.
"""

from __future__ import annotations

from math import gcd, log, prod

from q32_fixed_q_content import truncation_coefficients
from q32_legendre_content import franel_numbers, primes_up_to


SAMPLES = (30, 65, 84, 100, 120, 140, 142, 150, 160)


def centered_nonzero_residue(value: int, modulus: int) -> int:
    residue = value % modulus
    result = min(residue, modulus - residue)
    return result or modulus


def moment_certificate(
    coefficients: list[int], sign: int
) -> tuple[int, int]:
    n = len(coefficients) - 1
    base = sum(
        sign**degree * coefficient
        for degree, coefficient in enumerate(coefficients)
    )
    powers = [1] * (n + 1)
    moment_ideal = 0
    for _ in range(1, n + 1):
        moment = 0
        for degree in range(1, n + 1):
            powers[degree] *= degree
            moment += (
                sign**degree * powers[degree] * coefficients[degree]
            )
        moment_ideal = gcd(moment_ideal, moment)
    return moment_ideal, centered_nonzero_residue(base, moment_ideal)


def main() -> None:
    limit = max(SAMPLES)
    franel = franel_numbers(limit)
    primes = primes_up_to(limit)
    for n in SAMPLES:
        coefficients = truncation_coefficients(n, 3, franel)
        candidate_primorial = prod(
            prime
            for prime in primes
            if 4 * prime > n and 3 * prime <= n
        )
        records = []
        for sign in (1, -1):
            ideal, residue = moment_certificate(coefficients, sign)
            assert ideal % candidate_primorial == 0
            records.append(
                (
                    sign,
                    log(residue) / n if residue > 1 else 0.0,
                    len(str(ideal // candidate_primorial)),
                )
            )
        print(
            f"n={n} candidate_rate="
            f"{log(candidate_primorial) / n if candidate_primorial > 1 else 0:.9f} "
            f"moment_records={records}"
        )


if __name__ == "__main__":
    main()
