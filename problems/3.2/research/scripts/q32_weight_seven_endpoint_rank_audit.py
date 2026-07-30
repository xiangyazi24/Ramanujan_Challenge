#!/usr/bin/env python3
"""Audit the next Apéry endpoint grade after the weight-five carrier.

Let

    Delta_p = A_(p-1)-1,
    h_p = (A_p-5+7 Delta_p)/p^5.

The proved mod-p^6 theorem says that the endpoint defects are

    A_(mp)-A_m       = C_m Delta_p + p^5 P_m h_p       (mod p^6),
    A_(mp-1)-A_(m-1) = F_m Delta_p + p^5 Q_m h_p       (mod p^6),

where P_m=-L_m/24 and Q_m=-M_m/24.  This script forms the next
normalized residues

    d_m = (...)/p^6 mod p,       f_m = (...)/p^6 mod p

using exact integer Apéry values and the full integer h_p.  For
p >= 11, p != 769, it checks the experimentally discovered rank-one
law

    d_m = R_m d_2,       f_m = S_m d_2                 (mod p),

where

    R_m = m^3[(60m^3-14m^2-51)A_m+(22m^2+3)A_(m-1)]
          / (288*769),

    S_m = m^3[-(22m^2+3)A_m+(60m^3+14m^2+51)A_(m-1)]
          / (288*769).

The fixed primes 7 and 769 are recorded separately.  At p=7 the
stable formula is false; at p=769 the normalization d_2 has a pole
and in fact d_2=0 while other coordinates need not vanish.
"""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from math import isqrt


def primes_at_most(limit: int) -> list[int]:
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[:2] = b"\x00\x00"
    for prime in range(2, isqrt(limit) + 1):
        if sieve[prime]:
            sieve[prime * prime : limit + 1 : prime] = b"\x00" * (
                (limit - prime * prime) // prime + 1
            )
    return [candidate for candidate in range(2, limit + 1) if sieve[candidate]]
from q32_weight_five_endpoint_sextic_target_audit import (
    apery_numbers,
    endpoint_coefficients,
    fraction_mod,
)


NORMALIZATION_PRIME = 769


def endpoint_rank_carriers(
    quotient: int, values: list[int]
) -> tuple[Fraction, Fraction]:
    m = quotient
    direct = Fraction(
        m**3
        * (
            (60 * m**3 - 14 * m**2 - 51) * values[m]
            + (22 * m**2 + 3) * values[m - 1]
        ),
        288 * NORMALIZATION_PRIME,
    )
    reflected = Fraction(
        m**3
        * (
            -(22 * m**2 + 3) * values[m]
            + (60 * m**3 + 14 * m**2 + 51) * values[m - 1]
        ),
        288 * NORMALIZATION_PRIME,
    )
    return direct, reflected


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prime-limit", type=int, default=1000)
    parser.add_argument("--quotient-limit", type=int, default=20)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    primes = [
        prime
        for prime in primes_at_most(args.prime_limit)
        if prime >= 7
    ]
    values = apery_numbers(max(primes) * args.quotient_limit)
    moments = {
        quotient: endpoint_coefficients(quotient)
        for quotient in range(1, args.quotient_limit + 1)
    }
    carriers = {
        quotient: endpoint_rank_carriers(quotient, values)
        for quotient in range(1, args.quotient_limit + 1)
    }

    assert carriers[1] == (Fraction(0), Fraction(0))
    assert carriers[2] == (
        Fraction(1),
        Fraction(-103, NORMALIZATION_PRIME),
    )

    checks: Counter[str] = Counter()
    exceptional_rows: dict[int, list[tuple[int, int]]] = {}

    for prime in primes:
        modulus5 = prime**5
        modulus6 = prime**6
        modulus7 = prime**7
        delta = values[prime - 1] - 1
        assert delta % prime**3 == 0

        direct_anchor_residual = values[prime] - 5 + 7 * delta
        assert direct_anchor_residual % modulus5 == 0
        h_scalar = direct_anchor_residual // modulus5

        direct_residues: list[int] = []
        reflected_residues: list[int] = []

        for quotient in range(1, args.quotient_limit + 1):
            direct_coefficient = (
                quotient**3
                * (
                    values[quotient - 1]
                    - 17 * values[quotient]
                )
                // 12
            )
            reflected_coefficient = (
                quotient**3
                * (
                    17 * values[quotient - 1]
                    - values[quotient]
                )
                // 12
            )
            direct_moment, reflected_moment = moments[quotient]
            direct_weight_five = -direct_moment / 24
            reflected_weight_five = -reflected_moment / 24

            direct_numerator = (
                values[quotient * prime]
                - values[quotient]
                - direct_coefficient * delta
                - modulus5
                * h_scalar
                * fraction_mod(direct_weight_five, modulus7)
            ) % modulus7
            reflected_numerator = (
                values[quotient * prime - 1]
                - values[quotient - 1]
                - reflected_coefficient * delta
                - modulus5
                * h_scalar
                * fraction_mod(reflected_weight_five, modulus7)
            ) % modulus7

            assert direct_numerator % modulus6 == 0
            assert reflected_numerator % modulus6 == 0
            direct_residues.append(
                direct_numerator // modulus6 % prime
            )
            reflected_residues.append(
                reflected_numerator // modulus6 % prime
            )
            checks["residual_integrality"] += 2

        if prime in (7, NORMALIZATION_PRIME):
            exceptional_rows[prime] = list(
                zip(direct_residues, reflected_residues)
            )
            continue

        base_scalar = direct_residues[1]
        for quotient in range(1, args.quotient_limit + 1):
            direct_carrier, reflected_carrier = carriers[quotient]
            assert direct_residues[quotient - 1] == (
                fraction_mod(direct_carrier, prime) * base_scalar
            ) % prime
            assert reflected_residues[quotient - 1] == (
                fraction_mod(reflected_carrier, prime) * base_scalar
            ) % prime
            checks["direct_rank_one"] += 1
            checks["reflected_rank_one"] += 1

    print(f"prime_limit={args.prime_limit}")
    print(f"quotient_limit={args.quotient_limit}")
    for key in sorted(checks):
        print(f"{key}_checks={checks[key]}")
    for prime in sorted(exceptional_rows):
        rows = exceptional_rows[prime]
        print(
            f"exceptional_p={prime}"
            f" d2={rows[1][0]} f2={rows[1][1]}"
            f" nonzero_coordinates="
            f"{sum(left != 0 or right != 0 for left, right in rows)}"
        )
    print("stable_prime_range=p>=11,p!=769")
    print("first_failure=None")
    print("failures=0")


if __name__ == "__main__":
    main()
