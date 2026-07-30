#!/usr/bin/env python3
"""Exact audit of the mesoscopic Apéry-term radical reduction."""

from __future__ import annotations

import argparse
from fractions import Fraction
from math import comb, gcd, isqrt, lcm


def primes_upto(limit: int) -> list[int]:
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[:2] = b"\x00\x00"
    for prime in range(2, isqrt(limit) + 1):
        if sieve[prime]:
            sieve[prime * prime : limit + 1 : prime] = (
                b"\x00" * ((limit - prime * prime) // prime + 1)
            )
    return [prime for prime in range(2, limit + 1) if sieve[prime]]


def apery_values(limit: int) -> list[int]:
    values = [1, 5]
    for index in range(1, limit):
        polynomial = (
            34 * index**3 + 51 * index**2 + 27 * index + 5
        )
        numerator = (
            polynomial * values[index] - index**3 * values[index - 1]
        )
        denominator = (index + 1) ** 3
        assert numerator % denominator == 0
        values.append(numerator // denominator)
    return values[: limit + 1]


def finite_differences(values: list[int]) -> list[int]:
    row = values[:]
    coefficients = []
    while row:
        coefficients.append(row[0])
        row = [row[i + 1] - row[i] for i in range(len(row) - 1)]
    return coefficients


def evaluate_newton(coefficients: list[int], argument: int) -> int:
    return sum(
        coefficient * comb(argument - 1, degree)
        for degree, coefficient in enumerate(coefficients)
        if degree < argument
    )


def residue(value: Fraction, prime: int) -> int:
    assert value.denominator % prime
    return (
        value.numerator * pow(value.denominator, -1, prime)
    ) % prime


def four_scalar_rows(
    n: int,
    apery: list[int],
    denominator: int,
) -> tuple[int, int, int]:
    cutoff = (n - 1) // 2
    scaled_harmonic = [0] * (n + 1)
    running = 0
    for index in range(1, n + 1):
        scaled_harmonic[index] = running
        running += denominator // index

    x_zero = 0
    y_zero = 0
    g_kernel = [0] * (n + 1)
    for index in range(1, n + 1):
        x_zero += (
            (-1) ** (index + 1)
            * (denominator // index) ** 2
            * apery[n - index]
        )
        g_kernel[index] = (
            (-1) ** index
            * (denominator // index)
            * scaled_harmonic[index]
            * apery[n - index]
        )
        y_zero += g_kernel[index]

    inverse = [0] * (cutoff + 1)
    prescribed = []
    for index in range(1, cutoff + 1):
        inverse[index] = 1 - sum(
            inverse[shift] * apery[index - shift]
            for shift in range(1, index)
        )
        prescribed.append(
            (-1) ** (index + 1) * index * inverse[index]
        )
    newton = finite_differences(prescribed)
    y_star = sum(
        evaluate_newton(newton, index) * g_kernel[index]
        for index in range(1, n + 1)
    )
    return x_zero, y_zero, y_star


def outer_scalars(
    quotient: int,
    apery: list[int],
    harmonic: list[Fraction],
) -> tuple[Fraction, Fraction]:
    c_value = sum(
        Fraction(
            (-1) ** (index + 1) * apery[quotient - index],
            index**2,
        )
        for index in range(1, quotient + 1)
    )
    ell_value = sum(
        Fraction(
            (-1) ** index * apery[quotient - index],
            index,
        )
        * harmonic[index - 1]
        for index in range(1, quotient + 1)
    )
    return c_value, ell_value


def audit(maximum_n: int) -> None:
    primes = primes_upto(maximum_n)
    apery = apery_values(maximum_n)

    denominators = [1] * (maximum_n + 1)
    for n in range(1, maximum_n + 1):
        denominators[n] = lcm(denominators[n - 1], n)

    harmonic = [Fraction(0)]
    for n in range(1, maximum_n + 1):
        harmonic.append(harmonic[-1] + Fraction(1, n))
    outer = [
        None,
        *[
            outer_scalars(quotient, apery, harmonic)
            for quotient in range(1, maximum_n + 1)
        ],
    ]

    high_prime_checks = 0
    moving_checks = 0
    top_equivalences = 0
    apery_term_exceptions = 0
    four_scalar_exceptions = 0

    for n in range(10, maximum_n + 1):
        denominator = denominators[n]
        x_zero, y_zero, y_star = four_scalar_rows(
            n, apery, denominator
        )
        four_scalar = gcd(
            denominator**2,
            gcd(abs(x_zero), gcd(abs(y_zero), abs(y_star))),
        )

        for prime in primes:
            if not isqrt(n) < prime <= n:
                continue
            quotient, remainder = divmod(n, prime)
            assert quotient < prime
            moving = apery[remainder] % prime == 0
            in_apery_term = apery[n] % prime == 0
            in_four_scalar = four_scalar % prime == 0

            assert in_apery_term == (
                apery[quotient] * apery[remainder] % prime == 0
            )
            if moving:
                assert in_apery_term
                assert in_four_scalar
                moving_checks += 1
            if in_apery_term and not moving:
                assert apery[quotient] % prime == 0
                apery_term_exceptions += 1

            if prime > n // 2:
                assert in_apery_term == moving
                assert in_four_scalar == moving
                top_equivalences += 1
            else:
                c_value, ell_value = outer[quotient]
                unit = (denominator // prime) ** 2 % prime
                assert x_zero % prime == (
                    unit
                    * (apery[remainder] % prime)
                    * residue(c_value, prime)
                ) % prime
                assert y_zero % prime == (
                    unit
                    * (apery[remainder] % prime)
                    * residue(ell_value, prime)
                ) % prime
                if in_four_scalar and not moving:
                    assert residue(c_value, prime) == 0
                    assert residue(ell_value, prime) == 0
                    four_scalar_exceptions += 1

            high_prime_checks += 1

    print(
        {
            "maximum_n": maximum_n,
            "high_prime_checks": high_prime_checks,
            "moving_checks": moving_checks,
            "top_equivalences": top_equivalences,
            "apery_term_exceptions": apery_term_exceptions,
            "four_scalar_exceptions": four_scalar_exceptions,
            "failures": 0,
        }
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--maximum-n", type=int, default=180)
    audit(parser.parse_args().maximum_n)
