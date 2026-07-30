#!/usr/bin/env python3
"""Audit the conjectural eighth-precision endpoint rank-one law.

This is an exact modular audit, not a proof of the block-harmonic
identity.  It extends the proved precision-p^7 endpoint law by retaining
the full next digit of its scalar ``w_p`` and measuring the remaining
precision-p^8 defect.

For p >= 11, p not in {769, 22129}, define

    Delta = A_(p-1) - 1,
    H     = A_p - 5 + 7 Delta,

    w = (A_(2p) - 73 + 824 Delta - (752/5) H) / p^6,

    v = (A_(2p-1) - 5 - 8 Delta - (336/5) H
         + (103/769) p^6 w) / p^7.

The script checks, for every requested quotient m,

    A_(mp) - A_m
      = E_m Delta + P_m H + p^6 R_m w + p^7 C_m v  (mod p^8),

    A_(mp-1) - A_(m-1)
      = F_m Delta + Q_m H + p^6 S_m w + p^7 D_m v (mod p^8),

where E,F,P,Q,R,S are the proved lower-grade carriers and

    C_m = -m^3(P(m) A_m + Q(m) A_(m-1)) / N,
    D_m =  m^3(Q(m) A_m + P(-m) A_(m-1)) / N,

    P(m) = 3845 m^4 - 29268 m^3 + 36974 m^2 - 9112,
    Q(m) = 45371 m^4 - 58102 m^2 + 536,
    N    = 305911296 = 2^9 3^3 22129.

The normalization has C_1=D_1=C_2=0 and D_2=1.
"""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from math import isqrt

from q32_weight_five_endpoint_sextic_target_audit import (
    apery_numbers,
    endpoint_coefficients,
    fraction_mod,
)
from q32_weight_seven_endpoint_rank_audit import endpoint_rank_carriers


OLD_NORMALIZATION_PRIME = 769
NEW_NORMALIZATION_PRIME = 22129
NORMALIZATION = 305_911_296


def primes_at_most(limit: int) -> list[int]:
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[:2] = b"\x00\x00"
    for prime in range(2, isqrt(limit) + 1):
        if sieve[prime]:
            sieve[prime * prime : limit + 1 : prime] = b"\x00" * (
                (limit - prime * prime) // prime + 1
            )
    return [
        candidate
        for candidate in range(2, limit + 1)
        if sieve[candidate]
    ]


def weight_eight_carriers(
    quotient: int, values: list[int]
) -> tuple[Fraction, Fraction]:
    m = quotient
    direct_polynomial = (
        3845 * m**4
        - 29268 * m**3
        + 36974 * m**2
        - 9112
    )
    reflected_polynomial = (
        45371 * m**4 - 58102 * m**2 + 536
    )
    direct = Fraction(
        -m**3
        * (
            direct_polynomial * values[m]
            + reflected_polynomial * values[m - 1]
        ),
        NORMALIZATION,
    )
    reflected = Fraction(
        m**3
        * (
            reflected_polynomial * values[m]
            + (
                3845 * m**4
                + 29268 * m**3
                + 36974 * m**2
                - 9112
            )
            * values[m - 1]
        ),
        NORMALIZATION,
    )
    return direct, reflected


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prime-limit", type=int, default=500)
    parser.add_argument("--quotient-limit", type=int, default=12)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    primes = [
        prime
        for prime in primes_at_most(args.prime_limit)
        if prime >= 11
        and prime
        not in (OLD_NORMALIZATION_PRIME, NEW_NORMALIZATION_PRIME)
    ]
    values = apery_numbers(max(primes) * args.quotient_limit)
    moments = {
        quotient: endpoint_coefficients(quotient)
        for quotient in range(1, args.quotient_limit + 1)
    }
    old_carriers = {
        quotient: endpoint_rank_carriers(quotient, values)
        for quotient in range(1, args.quotient_limit + 1)
    }
    new_carriers = {
        quotient: weight_eight_carriers(quotient, values)
        for quotient in range(1, args.quotient_limit + 1)
    }

    assert new_carriers[1] == (Fraction(0), Fraction(0))
    assert new_carriers[2] == (Fraction(0), Fraction(1))

    checks: Counter[str] = Counter()
    zero_anchor_primes: list[int] = []

    for prime in primes:
        modulus6 = prime**6
        modulus7 = prime**7
        modulus8 = prime**8
        delta = values[prime - 1] - 1
        endpoint_h = values[prime] - 5 + 7 * delta
        assert delta % prime**3 == 0
        assert endpoint_h % prime**5 == 0

        endpoint_w = (
            Fraction(values[2 * prime] - 73 + 824 * delta)
            - Fraction(752, 5) * endpoint_h
        ) / modulus6
        endpoint_v = (
            Fraction(
                values[2 * prime - 1] - 5 - 8 * delta
            )
            - Fraction(336, 5) * endpoint_h
            + Fraction(103, OLD_NORMALIZATION_PRIME)
            * modulus6
            * endpoint_w
        ) / modulus7
        endpoint_v_residue = fraction_mod(endpoint_v, prime)
        if endpoint_v_residue == 0:
            zero_anchor_primes.append(prime)

        for quotient in range(1, args.quotient_limit + 1):
            direct_coefficient = Fraction(
                quotient**3
                * (
                    values[quotient - 1]
                    - 17 * values[quotient]
                ),
                12,
            )
            reflected_coefficient = Fraction(
                quotient**3
                * (
                    17 * values[quotient - 1]
                    - values[quotient]
                ),
                12,
            )
            direct_moment, reflected_moment = moments[quotient]
            direct_weight_five = -direct_moment / 24
            reflected_weight_five = -reflected_moment / 24
            direct_old, reflected_old = old_carriers[quotient]
            direct_new, reflected_new = new_carriers[quotient]

            direct_residual = (
                Fraction(values[quotient * prime] - values[quotient])
                - direct_coefficient * delta
                - direct_weight_five * endpoint_h
                - modulus6 * direct_old * endpoint_w
            )
            reflected_residual = (
                Fraction(
                    values[quotient * prime - 1]
                    - values[quotient - 1]
                )
                - reflected_coefficient * delta
                - reflected_weight_five * endpoint_h
                - modulus6 * reflected_old * endpoint_w
            )

            direct_modulus_residue = fraction_mod(
                direct_residual, modulus8
            )
            reflected_modulus_residue = fraction_mod(
                reflected_residual, modulus8
            )
            assert direct_modulus_residue % modulus7 == 0
            assert reflected_modulus_residue % modulus7 == 0
            checks["residual_integrality"] += 2

            direct_digit = (
                direct_modulus_residue // modulus7
            ) % prime
            reflected_digit = (
                reflected_modulus_residue // modulus7
            ) % prime
            assert direct_digit == (
                fraction_mod(direct_new, prime)
                * endpoint_v_residue
            ) % prime
            assert reflected_digit == (
                fraction_mod(reflected_new, prime)
                * endpoint_v_residue
            ) % prime
            checks["direct_rank_one"] += 1
            checks["reflected_rank_one"] += 1

    print(f"prime_limit={args.prime_limit}")
    print(f"quotient_limit={args.quotient_limit}")
    for name in sorted(checks):
        print(f"{name}_checks={checks[name]}")
    print(
        "normalization="
        f"{NORMALIZATION}=2^9*3^3*{NEW_NORMALIZATION_PRIME}"
    )
    print(
        "fixed_exceptions="
        f"{OLD_NORMALIZATION_PRIME},{NEW_NORMALIZATION_PRIME}"
    )
    print(
        "zero_anchor_primes="
        + (
            ",".join(map(str, zero_anchor_primes))
            if zero_anchor_primes
            else "none"
        )
    )
    print("status=exact_computational_audit_not_symbolic_proof")
    print("failures=0")


if __name__ == "__main__":
    main()
