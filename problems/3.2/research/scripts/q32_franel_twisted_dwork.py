#!/usr/bin/env python3
"""Audit the Frobenius-twisted Dwork congruence for Franel periods.

Put

    P_j(x) = sum_m (-1)^m binom(j,m) F_m x^m
           = CT_u,v (1 - x Lambda(u,v))^j.

The coefficient-ring version of the Samol--van Straten congruence predicts

    P_(j+m p^r)(x) P_floor(j/p)(x^p)
      = P_j(x) P_(floor(j/p)+m p^(r-1))(x^p)  (mod p^r).

Taking an Euler logarithmic derivative and descending in r proves that
theta log P_j is congruence-preserving in the parameter j.  This script
checks both statements coefficientwise.
"""

from __future__ import annotations

from math import comb

from q32_strehl_gcd import franel_numbers


DEGREE_LIMIT = 90
PARAMETER_LIMIT = 90
PRIMES = (2, 3, 5, 7)


def period(parameter: int, franel: list[int]) -> list[int]:
    result = [0] * (DEGREE_LIMIT + 1)
    for degree in range(min(parameter, DEGREE_LIMIT) + 1):
        result[degree] = (
            (-1) ** degree
            * comb(parameter, degree)
            * franel[degree]
        )
    return result


def substitute_frobenius(polynomial: list[int], prime: int) -> list[int]:
    result = [0] * (DEGREE_LIMIT + 1)
    for degree in range(DEGREE_LIMIT // prime + 1):
        result[prime * degree] = polynomial[degree]
    return result


def multiply(left: list[int], right: list[int]) -> list[int]:
    result = [0] * (DEGREE_LIMIT + 1)
    for left_degree, left_coefficient in enumerate(left):
        if not left_coefficient:
            continue
        for right_degree in range(DEGREE_LIMIT - left_degree + 1):
            result[left_degree + right_degree] += (
                left_coefficient * right[right_degree]
            )
    return result


def logarithmic_derivative(polynomial: list[int]) -> list[int]:
    result = [0] * (DEGREE_LIMIT + 1)
    for degree in range(1, DEGREE_LIMIT + 1):
        result[degree] = (
            degree * polynomial[degree]
            - sum(
                polynomial[index] * result[degree - index]
                for index in range(1, degree + 1)
            )
        )
    return result


def main() -> None:
    franel = franel_numbers(DEGREE_LIMIT)
    periods = [
        period(parameter, franel)
        for parameter in range(PARAMETER_LIMIT + 1)
    ]
    logarithmic_derivatives = [
        logarithmic_derivative(polynomial) for polynomial in periods
    ]

    for prime in PRIMES:
        prime_power = prime
        exponent = 1
        while prime_power <= PARAMETER_LIMIT:
            for parameter in range(PARAMETER_LIMIT - prime_power + 1):
                lower = parameter // prime
                left = multiply(
                    periods[parameter + prime_power],
                    substitute_frobenius(periods[lower], prime),
                )
                right = multiply(
                    periods[parameter],
                    substitute_frobenius(
                        periods[lower + prime_power // prime], prime
                    ),
                )
                assert all(
                    (left[degree] - right[degree]) % prime_power == 0
                    for degree in range(DEGREE_LIMIT + 1)
                ), ("twisted Dwork", prime, exponent, parameter)

                assert all(
                    (
                        logarithmic_derivatives[parameter + prime_power][
                            degree
                        ]
                        - logarithmic_derivatives[parameter][degree]
                    )
                    % prime_power
                    == 0
                    for degree in range(DEGREE_LIMIT + 1)
                ), ("logarithmic derivative", prime, exponent, parameter)
            exponent += 1
            prime_power *= prime

    print(
        "verified twisted Dwork and parameter congruences through "
        f"degree {DEGREE_LIMIT}, parameter {PARAMETER_LIMIT}"
    )


if __name__ == "__main__":
    main()
