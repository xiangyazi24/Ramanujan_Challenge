#!/usr/bin/env python3
"""Verify the arbitrary-cutoff Lucas block identity.

Put

    n=q*p+r,  J=A*p+B,  0<=r,B<p.

For the Legendre--Euler transform T_(n,J), polynomial Lucas and Franel
Lucas give the exact coefficientwise congruence

    T_(q*p+r,A*p+B)(c)
      = A_r T_(q,A-1)(c^p)
        + K_(q,A)(c^p) g_A(c^p) T_(r,B)(c)          (mod p),

where

    K_(q,A)(c)=[y^A]Q_q(c+y).

The convention is T_(q,-1)=0, and a cutoff beyond the degree is saturated.
This identity exposes the final high-cutoff plateau in every quotient slice.
If r=p-1-j and p|A_j, its last-block onset is J=q*p+j.
"""

from __future__ import annotations

from math import comb

from q32_legendre_content import franel_numbers
from q32_newton import apery_numbers
from q32_q1_all_cutoff_profile import transform_coefficients_mod


PRIMES = (5, 7, 11, 13)
MAX_QUOTIENT = 3


def legendre_coefficient(n: int, index: int) -> int:
    return comb(n, index) * comb(n + index, index)


def polynomial_product(
    left: list[int], right: list[int], prime: int
) -> list[int]:
    result = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            result[i + j] = (result[i + j] + a * b) % prime
    return result


def kernel_polynomial(
    quotient: int, block: int, prime: int
) -> list[int]:
    """Coefficients of K_(quotient,block)(c)."""

    return [
        (
            legendre_coefficient(quotient, degree + block)
            * comb(degree + block, block)
        )
        % prime
        for degree in range(quotient - block + 1)
    ]


def shifted_franel_polynomial(
    block: int, prime: int, franel: list[int]
) -> list[int]:
    """Coefficients of g_block(c), in increasing degree order."""

    result = [0] * (block + 1)
    for index in range(block + 1):
        degree = block - index
        result[degree] = (
            comb(block, index) * (-1) ** degree * franel[index]
        ) % prime
    return result


def add_scaled_substitution(
    target: list[int],
    source: list[int],
    scale: int,
    dilation: int,
    prime: int,
) -> None:
    for degree, coefficient in enumerate(source):
        target[dilation * degree] = (
            target[dilation * degree] + scale * coefficient
        ) % prime


def predicted_block(
    quotient: int,
    residue: int,
    block: int,
    low_cutoff: int,
    prime: int,
    franel: list[int],
    apery: list[int],
) -> list[int]:
    n = quotient * prime + residue
    result = [0] * (n + 1)

    if block:
        completed_high_digit = transform_coefficients_mod(
            quotient, block - 1, prime, franel
        )
        add_scaled_substitution(
            result,
            completed_high_digit,
            apery[residue] % prime,
            prime,
            prime,
        )

    kernel = kernel_polynomial(quotient, block, prime)
    shifted = shifted_franel_polynomial(block, prime, franel)
    high_factor = polynomial_product(kernel, shifted, prime)
    low_factor = transform_coefficients_mod(
        residue, min(low_cutoff, residue), prime, franel
    )

    # Multiply high_factor(c^p) by low_factor(c) without materializing the
    # sparse substituted polynomial.
    product = [0] * (n + 1)
    for high_degree, high_coefficient in enumerate(high_factor):
        for low_degree, low_coefficient in enumerate(low_factor):
            product[prime * high_degree + low_degree] = (
                product[prime * high_degree + low_degree]
                + high_coefficient * low_coefficient
            ) % prime
    for degree, coefficient in enumerate(product):
        result[degree] = (result[degree] + coefficient) % prime
    return result


def main() -> None:
    limit = (MAX_QUOTIENT + 1) * max(PRIMES)
    franel = franel_numbers(limit)
    apery = apery_numbers(limit)
    identities = 0

    for prime in PRIMES:
        for quotient in range(MAX_QUOTIENT + 1):
            for residue in range(prime):
                n = quotient * prime + residue
                for cutoff in range(n + 1):
                    block, low_cutoff = divmod(cutoff, prime)
                    actual = transform_coefficients_mod(
                        n, cutoff, prime, franel
                    )
                    predicted = predicted_block(
                        quotient,
                        residue,
                        block,
                        low_cutoff,
                        prime,
                        franel,
                        apery,
                    )
                    assert actual == predicted, (
                        prime,
                        quotient,
                        residue,
                        cutoff,
                    )
                    identities += 1

    print(
        "arbitrary-cutoff Lucas block identity verified: "
        f"primes={PRIMES}, q<= {MAX_QUOTIENT}, identities={identities}"
    )


if __name__ == "__main__":
    main()
