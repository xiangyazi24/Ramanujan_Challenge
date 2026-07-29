#!/usr/bin/env python3
"""Audit the congruence-preserving Franel logarithmic derivatives.

For

    Phi(alpha,x) = sum_m (-1)^m binom(alpha,m) F_m x^m

put H_n(alpha)=[x^n] theta log Phi(alpha,x).  The exact criterion for an
integer-valued polynomial H_n to preserve congruences on Z is

    lcm(1,...,k) | Delta^k H_n(0)       for every k.

That criterion also makes H_n' integer-valued, because
lcm(1,...,k) times the derivative of binom(alpha,k) is integer-valued.
This script checks the full triangular family of Newton coefficients.
"""

from __future__ import annotations

from math import comb, lcm

from q32_strehl_gcd import franel_numbers


LIMIT = 140


def logarithmic_derivative_coefficients(
    parameter: int, franel: list[int]
) -> list[int]:
    polynomial = [0] * (LIMIT + 1)
    for degree in range(parameter + 1):
        polynomial[degree] = (
            (-1) ** degree
            * comb(parameter, degree)
            * franel[degree]
        )
    quotient = [0] * (LIMIT + 1)
    for degree in range(1, LIMIT + 1):
        quotient[degree] = (
            degree * polynomial[degree]
            - sum(
                polynomial[index] * quotient[degree - index]
                for index in range(1, min(parameter, degree) + 1)
            )
        )
    return quotient


def main() -> None:
    franel = franel_numbers(LIMIT)
    values_by_parameter = [
        logarithmic_derivative_coefficients(parameter, franel)
        for parameter in range(LIMIT + 1)
    ]
    least_common_multiple = 1
    for newton_degree in range(LIMIT + 1):
        if newton_degree:
            least_common_multiple = lcm(
                least_common_multiple, newton_degree
            )
        for coefficient_degree in range(newton_degree, LIMIT + 1):
            differences = [
                values_by_parameter[parameter][coefficient_degree]
                for parameter in range(newton_degree + 1)
            ]
            while len(differences) > 1:
                differences = [
                    differences[index + 1] - differences[index]
                    for index in range(len(differences) - 1)
                ]
            assert differences[0] % least_common_multiple == 0, (
                coefficient_degree,
                newton_degree,
                differences[0],
                least_common_multiple,
            )
    print(
        "verified congruence-preserving Franel logarithmic derivatives "
        f"for 0 <= k <= n <= {LIMIT}"
    )


if __name__ == "__main__":
    main()
