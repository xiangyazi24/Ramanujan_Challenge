#!/usr/bin/env python3
"""Audit the direct and Strehl--Franel q=1 truncation certificates.

Write

    A_n = sum_k binom(n,k)^2 binom(n+k,k)^2
        = sum_k binom(n,k) binom(n+k,k) F_k,

where F_k=sum_i binom(k,i)^3 is the k-th Franel number.  If
J=floor((n-1)/3), both prefixes through J reduce to the folded A_j modulo
every top-half prime.  Their gcd therefore contains the complete q=1 bad
radical.  The Strehl prefix has exponential rate log(8), improving the
single-certificate rate log(16) of the direct prefix.
"""

from __future__ import annotations

from math import comb, factorial, gcd, log

from q32_newton import apery_numbers


LIMIT = 1200


def valuation(value: int, prime: int) -> int:
    exponent = 0
    while value % prime == 0:
        value //= prime
        exponent += 1
    return exponent


def primes_up_to(limit: int) -> list[int]:
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[:2] = b"\x00\x00"
    for prime in range(2, int(limit**0.5) + 1):
        if sieve[prime]:
            sieve[prime * prime : limit + 1 : prime] = b"\x00" * (
                (limit - prime * prime) // prime + 1
            )
    return [value for value, is_prime in enumerate(sieve) if is_prime]


def franel_numbers(limit: int) -> list[int]:
    return [
        sum(comb(index, part) ** 3 for part in range(index + 1))
        for index in range(limit + 1)
    ]


def certificates(n: int, franel: list[int]) -> tuple[int, int]:
    cutoff = (n - 1) // 3
    direct = 0
    strehl = 0
    for k in range(cutoff + 1):
        kernel = comb(n, k) * comb(n + k, k)
        direct += kernel * kernel
        strehl += kernel * franel[k]
    return direct, strehl


def main() -> None:
    apery = apery_numbers(LIMIT)
    franel = franel_numbers(LIMIT)
    primes = primes_up_to(LIMIT)
    records: list[tuple[int, int, float]] = []
    carrier_records: list[tuple[int, int, float]] = []

    # Independently verify the exact Strehl identity at modest indices.
    for n in range(81):
        total = sum(
            comb(n, k) * comb(n + k, k) * franel[k]
            for k in range(n + 1)
        )
        assert total == apery[n]

    for n in range(3, LIMIT + 1):
        direct, strehl = certificates(n, franel)
        common = gcd(direct, strehl)
        cutoff = (n - 1) // 3
        boundary_carrier = (
            comb(n, cutoff + 1) * comb(n + cutoff + 1, cutoff + 1)
        )
        scaled_difference = factorial(cutoff) ** 2 * (direct - strehl)
        assert scaled_difference % boundary_carrier == 0
        for prime in primes:
            if prime <= n // 2:
                continue
            if prime > n:
                break
            raw_index = n - prime
            folded_index = min(raw_index, prime - 1 - raw_index)
            assert direct % prime == apery[folded_index] % prime
            assert strehl % prime == apery[folded_index] % prime
            assert boundary_carrier % prime == 0
            if prime >= 7:
                expected_carrier_valuation = (
                    2 if 2 * n + 1 == 3 * prime else 1
                )
                assert (
                    valuation(boundary_carrier, prime)
                    == expected_carrier_valuation
                )
            assert (common % prime == 0) == (
                apery[folded_index] % prime == 0
            )
            assert (gcd(strehl, boundary_carrier) % prime == 0) == (
                apery[folded_index] % prime == 0
            )

        rate = log(common) / n if common > 1 else 0.0
        records.append((n, common, rate))
        carrier_common = gcd(strehl, boundary_carrier)
        carrier_rate = log(carrier_common) / n if carrier_common > 1 else 0.0
        carrier_records.append((n, carrier_common, carrier_rate))

    lower = 10
    while lower < LIMIT:
        upper = min(2 * lower, LIMIT)
        winner = max(
            (record for record in records if lower < record[0] <= upper),
            key=lambda record: record[2],
        )
        print(
            f"({lower},{upper}] max_log_gcd_over_n={winner[2]:.9f} "
            f"at_n={winner[0]} gcd={winner[1]}"
        )
        carrier_winner = max(
            (
                record
                for record in carrier_records
                if lower < record[0] <= upper
            ),
            key=lambda record: record[2],
        )
        print(
            f"({lower},{upper}] max_log_carrier_gcd_over_n="
            f"{carrier_winner[2]:.9f} at_n={carrier_winner[0]} "
            f"carrier_gcd={carrier_winner[1]}"
        )
        lower *= 2


if __name__ == "__main__":
    main()
