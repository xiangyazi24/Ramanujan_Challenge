#!/usr/bin/env python3
"""Fast audit of the full Newton--Padé determinantal divisor.

For M=floor(n/2), take the matrix with rows x=0,...,M-1,n and columns

    binom(x,k), A_x*binom(x,k),  0<=k<M.

Every non-boundary q=1 bad prime makes the last row equal an earlier row
modulo p, so it divides every maximal minor.  Unimodular Pascal elimination
shows that the gcd of all maximal minors is exactly

    gcd_{0<=j<M} L_j(n) * (A_n-A_j),

where L_j is the cardinal Lagrange coefficient for the nodes 0,...,M-1.
This formula avoids enumerating minors.
"""

from __future__ import annotations

from math import comb, gcd, log

from q32_newton import apery_numbers


LIMIT = 1200


def primes_up_to(limit: int) -> list[int]:
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[:2] = b"\x00\x00"
    for prime in range(2, int(limit**0.5) + 1):
        if sieve[prime]:
            sieve[prime * prime : limit + 1 : prime] = b"\x00" * (
                (limit - prime * prime) // prime + 1
            )
    return [value for value, is_prime in enumerate(sieve) if is_prime]


def determinantal_divisor(apery: list[int], n: int) -> int:
    half = n // 2
    if half == 0:
        return 1

    # L_0(n)=(-1)^(M-1) binom(n-1,M-1).
    cardinal = comb(n - 1, half - 1)
    if (half - 1) % 2:
        cardinal = -cardinal

    result = 0
    for index in range(half):
        difference = apery[n] - apery[index]
        if result:
            residue = (cardinal % result) * (difference % result) % result
            result = gcd(result, residue)
            if result == 1:
                break
        else:
            result = abs(cardinal * difference)

        if index + 1 < half:
            numerator = cardinal * (n - index) * (half - 1 - index)
            denominator = (index + 1) * (n - index - 1)
            assert numerator % denominator == 0
            cardinal = -(numerator // denominator)

    return result


def main() -> None:
    apery = apery_numbers(LIMIT)
    primes = primes_up_to(LIMIT)
    records = []
    for n in range(3, LIMIT + 1):
        common = determinantal_divisor(apery, n)
        if n >= 10:
            for prime in primes:
                if prime <= n // 2:
                    continue
                if prime > n:
                    break
                assert (common % prime == 0) == (apery[n] % prime == 0)
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
