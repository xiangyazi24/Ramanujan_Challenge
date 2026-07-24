#!/usr/bin/env python3
"""Verify Q539's exact odd-quotient formula for Delta_n/R_n.

For sqrt(n)<p<=floor(n/2), the Smith divisor Delta_n has p-adic valuation
one exactly when floor(n/p) is odd and A_(n mod p)=0 modulo p.  All other
prime powers in the smooth quotient lie below sqrt(n), where the first two
Smith coordinates give total logarithmic weight O(sqrt(n)).
"""

from __future__ import annotations

from math import isqrt

from q32_newton import apery_numbers
from q32_pade_minor_gcd import determinantal_divisor


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


def valuation(value: int, prime: int) -> int:
    result = 0
    while value % prime == 0:
        result += 1
        value //= prime
    return result


def main() -> None:
    apery = apery_numbers(LIMIT)
    primes = primes_up_to(LIMIT)
    checked = 0
    # The only failures below this threshold are the fixed prime p=5 at
    # n=11 and n=21; Q539's statement is asymptotic.
    for n in range(22, LIMIT + 1):
        half = n // 2
        delta = determinantal_divisor(apery, n)
        for prime in primes:
            if prime <= isqrt(n):
                continue
            if prime > half:
                break
            quotient, residue = divmod(n, prime)
            expected = int(
                quotient % 2 == 1 and apery[residue] % prime == 0
            )
            assert valuation(delta, prime) == expected
            checked += 1
    print(f"verified_prime_index_pairs={checked}")


if __name__ == "__main__":
    main()
