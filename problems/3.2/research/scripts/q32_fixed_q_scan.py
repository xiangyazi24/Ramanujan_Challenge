#!/usr/bin/env python3
"""Exact modular incidence scan for one fixed quotient slice.

For q=3 this records the primes p and residues r with n=q*p+r and
A_r=0 modulo p.  It is a reproducible sanity check for the first unresolved
odd slice in Q539, not a proof of a pointwise bound.
"""

from __future__ import annotations

from math import log


LIMIT = 20_000
QUOTIENT = 3


def primes_up_to(limit: int) -> tuple[bytearray, list[int]]:
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[:2] = b"\x00\x00"
    for prime in range(2, int(limit**0.5) + 1):
        if sieve[prime]:
            sieve[prime * prime : limit + 1 : prime] = b"\x00" * (
                (limit - prime * prime) // prime + 1
            )
    return sieve, [
        value for value, is_prime in enumerate(sieve) if is_prime
    ]


def apery_mod_prime(prime: int) -> list[int]:
    inverse = [0] * prime
    apery = [0] * prime
    inverse[1] = 1
    for k in range(2, prime):
        inverse[k] = prime - (prime // k) * inverse[prime % k] % prime
    apery[0] = 1
    apery[1] = 5 % prime
    for k in range(1, prime - 1):
        polynomial = (34 * k**3 + 51 * k**2 + 27 * k + 5) % prime
        apery[k + 1] = (
            (polynomial * apery[k] - k**3 * apery[k - 1])
            * pow(inverse[k + 1], 3, prime)
            % prime
        )
    return apery


def main() -> None:
    sieve, primes = primes_up_to(LIMIT)
    del sieve
    hits: list[list[tuple[int, int, int]]] = [
        [] for _ in range(LIMIT + 1)
    ]
    for prime in primes:
        if QUOTIENT * prime > LIMIT:
            break
        for residue, value in enumerate(apery_mod_prime(prime)):
            n = QUOTIENT * prime + residue
            if n <= LIMIT and value == 0:
                folded = min(residue, prime - 1 - residue)
                hits[n].append((prime, residue, folded))

    maximum = max(map(len, hits))
    print(f"maximum_hit_count={maximum}")
    for n, values in enumerate(hits):
        if len(values) == maximum and maximum:
            mass = sum(log(prime) for prime, _, _ in values)
            print(f"n={n} mass_over_n={mass/n:.12f} hits={values}")


if __name__ == "__main__":
    main()
