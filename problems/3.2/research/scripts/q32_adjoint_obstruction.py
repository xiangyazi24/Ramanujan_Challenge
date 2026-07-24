#!/usr/bin/env python3
"""Audit the lattice obstruction for interpolation residuals with g(n)=0.

For the Lagrange vector ell and arbitrary integer values z_j=g(j), the
constraint g(n)=0 is ell.z=0.  The values of

    -sum_j ell_j A_j z_j

on this saturated kernel form the ideal Theta_n Z, where Theta_n is the gcd
of all 2-by-2 minors

    ell_i * ell_j * (A_j-A_i).

Every non-boundary top-half prime divides every such minor, independently of
Apéry badness.  Thus Theta_n contains the full top-half primorial and no
nonzero g(n)=0 residual can have subexponential height.
"""

from __future__ import annotations

from math import comb, gcd, isqrt, log, prod

from q32_newton import apery_numbers
from q32_pade_minor_gcd import determinantal_divisor


LIMIT = 400


def primes_up_to(limit: int) -> list[int]:
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[:2] = b"\x00\x00"
    for prime in range(2, int(limit**0.5) + 1):
        if sieve[prime]:
            sieve[prime * prime : limit + 1 : prime] = b"\x00" * (
                (limit - prime * prime) // prime + 1
            )
    return [value for value, is_prime in enumerate(sieve) if is_prime]


def lagrange_vector(n: int) -> list[int]:
    half = n // 2
    return [
        (-1) ** (half - 1 - j)
        * comb(n, j)
        * comb(n - j - 1, half - 1 - j)
        for j in range(half)
    ]


def kernel_image_divisor(apery: list[int], n: int) -> int:
    ell = lagrange_vector(n)
    interpolated = sum(
        ell[j] * apery[j] for j in range(len(ell))
    )
    result = 0
    for j in range(len(ell)):
        result = gcd(
            result, ell[j] * (apery[j] - interpolated)
        )

    # At small indices, verify directly that this vertical-lattice formula
    # equals the gcd of every 2-by-2 minor.
    if n <= 40:
        pairwise = 0
        for i in range(len(ell)):
            for j in range(i + 1, len(ell)):
                pairwise = gcd(
                    pairwise,
                    ell[i] * ell[j] * (apery[j] - apery[i]),
                )
        assert result == pairwise
    return result


def main() -> None:
    apery = apery_numbers(LIMIT)
    primes = primes_up_to(LIMIT)
    for n in range(4, LIMIT + 1):
        half = n // 2
        ell = lagrange_vector(n)
        theta = kernel_image_divisor(apery, n)
        interpolated = sum(
            ell[j] * apery[j] for j in range(half)
        )
        delta = determinantal_divisor(apery, n)
        assert delta == gcd(theta, apery[n] - interpolated)
        universal_primes = [
            prime
            for prime in primes
            if half < prime <= n and not (n == 2 * prime - 1)
        ]
        universal_primorial = prod(universal_primes)
        assert theta % universal_primorial == 0
        for prime in universal_primes:
            value = theta
            valuation = 0
            while value % prime == 0:
                valuation += 1
                value //= prime
            assert valuation == 1

        # Q545's full large-prime classification: below the top interval,
        # only odd quotient Apéry-zero slices occur, apart from p=M.
        for prime in primes:
            # The classification uses A_0=1 != A_1=5 modulo p, so p=2,5
            # are fixed exceptions absorbed by the small-prime estimate.
            if prime <= max(5, isqrt(n)):
                continue
            if prime > n:
                break
            quotient, residue = divmod(n, prime)
            if quotient == 1:
                expected = n != 2 * prime - 1
            elif quotient % 2 == 1:
                expected = apery[residue] % prime == 0
            else:
                expected = prime == half
            value = theta
            valuation = 0
            while value % prime == 0:
                valuation += 1
                value //= prime
            assert valuation == int(expected)

    for n in (20, 40, 80, 160, 240, 320, 400):
        theta = kernel_image_divisor(apery, n)
        print(f"n={n} log_theta_over_n={log(theta)/n:.9f}")


if __name__ == "__main__":
    main()
