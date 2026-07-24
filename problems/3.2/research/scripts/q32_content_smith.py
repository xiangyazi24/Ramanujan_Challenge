#!/usr/bin/env python3
"""Verify the integral Smith block behind fixed-q polynomial content.

For odd q, J=floor((n-q)/(2q+1)), and the coefficient matrix M taking
arbitrary input data x_0,...,x_J to the Legendre--Euler truncation, every
2-by-2 minor involving row zero is divisible by the squarefree product P of
the fixed-q candidate primes.  Unimodular row and column operations give

    M ~ diag(1, P*W).

For the Franel input this yields the exact content identity

    content(M F) = gcd(S, P*gcd_d(sum_i W[d,i] F_i)).

At every non-midpoint candidate p, the low triangular block of W is
invertible modulo p.  Hence a bad p occurs in the content to exponent one.
"""

from __future__ import annotations

from math import comb, gcd, prod

from q32_fixed_q_content import truncation_coefficients
from q32_newton import apery_numbers
from q32_strehl_gcd import franel_numbers, primes_up_to, valuation


LIMIT = 180


def kernel(n: int, index: int) -> int:
    return comb(n, index) * comb(n + index, index)


def coefficient_entry(n: int, cutoff: int, degree: int, index: int) -> int:
    if degree == 0:
        return kernel(n, index)
    if not (
        max(0, cutoff + 1 - degree)
        <= index
        <= min(cutoff, n - degree)
    ):
        return 0
    return (
        (-1) ** (cutoff - index)
        * kernel(n, index + degree)
        * comb(index + degree, index)
        * comb(degree - 1, cutoff - index)
    )


def candidate_primes(n: int, quotient: int, primes: list[int]) -> list[int]:
    return [
        prime
        for prime in primes
        if prime >= 7 and divmod(n, prime)[0] == quotient
    ]


def audit(quotient: int) -> None:
    assert quotient > 0 and quotient % 2 == 1
    apery = apery_numbers(LIMIT)
    franel = franel_numbers(LIMIT)
    primes = primes_up_to(LIMIT)

    for n in range(quotient, LIMIT + 1):
        cutoff = (n - quotient) // (2 * quotient + 1)
        if cutoff < 0:
            continue
        candidates = candidate_primes(n, quotient, primes)
        candidate_product = prod(candidates)
        matrix = [
            [
                coefficient_entry(n, cutoff, degree, index)
                for index in range(cutoff + 1)
            ]
            for degree in range(n + 1)
        ]
        assert matrix[0][0] == 1

        residual = [[0] * (cutoff + 1) for _ in range(n + 1)]
        for degree in range(1, n + 1):
            anchor = matrix[degree][0]
            for index in range(1, cutoff + 1):
                minor = (
                    matrix[degree][index]
                    - anchor * matrix[0][index]
                )
                assert minor % candidate_product == 0
                residual[degree][index] = minor // candidate_product

        coefficients = [
            sum(
                matrix[degree][index] * franel[index]
                for index in range(cutoff + 1)
            )
            for degree in range(n + 1)
        ]
        assert coefficients == truncation_coefficients(
            n, quotient, franel
        )
        content = 0
        for coefficient in coefficients:
            content = gcd(content, coefficient)
        strehl = coefficients[0]
        residual_gcd = 0
        for degree in range(1, n + 1):
            residual_value = sum(
                residual[degree][index] * franel[index]
                for index in range(1, cutoff + 1)
            )
            residual_gcd = gcd(residual_gcd, residual_value)
        assert content == gcd(
            strehl, candidate_product * residual_gcd
        )

        boundary = cutoff + 1
        carrier = kernel(n, boundary)
        for prime in candidates:
            _, raw_index = divmod(n, prime)
            folded = min(raw_index, prime - 1 - raw_index)
            is_midpoint = (
                2 * n + 1 == (2 * quotient + 1) * prime
            )
            expected_carrier_valuation = 2 if is_midpoint else 1
            assert (
                valuation(carrier, prime)
                == expected_carrier_valuation
            ), (
                quotient,
                n,
                prime,
                folded,
                valuation(carrier, prime),
                expected_carrier_valuation,
            )
            if not is_midpoint and cutoff:
                # The reverse-triangular diagonal entries are
                # +/- (carrier/P)*binom(boundary,index), all p-units.
                assert (carrier // candidate_product) % prime
                for index in range(1, cutoff + 1):
                    assert comb(boundary, index) % prime
            if (
                not is_midpoint
                and apery[folded] % prime == 0
            ):
                assert valuation(content, prime) == 1

    print(f"q={quotient}: Smith block verified through n={LIMIT}")


def main() -> None:
    audit(1)
    audit(3)


if __name__ == "__main__":
    main()
