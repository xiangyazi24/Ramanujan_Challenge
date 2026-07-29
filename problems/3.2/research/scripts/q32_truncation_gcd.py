#!/usr/bin/env python3
"""Exact gcd audit for two truncations covering the full q=1 support.

For J=floor((n-1)/3), write the Apéry summand as

    T(n,k) = binom(n,k)^2 binom(n+k,k)^2.

For every q=1 candidate p, let j=min(n-p,2p-1-n).  Then the low prefix and
upper-half suffix satisfy L_n=A_j and H_n=4*A_j modulo p.  Hence every q=1
bad prime divides gcd(L_n,H_n).  The script verifies this and measures that
gcd.
"""

from __future__ import annotations

from math import comb, gcd, log, prod

from q32_newton import apery_numbers


LIMIT = 600


def primes_up_to(limit: int) -> list[int]:
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[:2] = b"\x00\x00"
    for prime in range(2, int(limit**0.5) + 1):
        if sieve[prime]:
            sieve[prime * prime : limit + 1 : prime] = b"\x00" * (
                (limit - prime * prime) // prime + 1
            )
    return [value for value, is_prime in enumerate(sieve) if is_prime]


def truncations(n: int) -> tuple[int, int]:
    cutoff = (n - 1) // 3
    low = sum(
        comb(n, k) ** 2 * comb(n + k, k) ** 2
        for k in range(cutoff + 1)
    )
    high = sum(
        comb(n, k) ** 2 * comb(n + k, k) ** 2
        for k in range((n + 1) // 2, n + 1)
    )
    return low, high


def main() -> None:
    apery = apery_numbers(LIMIT)
    primes = primes_up_to(LIMIT)
    records: list[tuple[int, int, float]] = []
    for n in range(3, LIMIT + 1):
        cutoff = (n - 1) // 3
        low, high = truncations(n)
        common = gcd(low, high)
        top_primes = [prime for prime in primes if n // 2 < prime <= n]
        top_primorial = prod(top_primes)
        assert (high - 4 * low) % (top_primorial * top_primorial) == 0

        for prime in top_primes:
            raw_index = n - prime
            folded_index = min(raw_index, prime - 1 - raw_index)
            assert 0 <= folded_index <= cutoff
            if apery[folded_index] % prime == 0:
                assert low % prime == 0
                assert high % prime == 0
                assert common % prime == 0

        rate = log(common) / n if common > 1 else 0.0
        records.append((n, common, rate))

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
        lower *= 2


if __name__ == "__main__":
    main()
