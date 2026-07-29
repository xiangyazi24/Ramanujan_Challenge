#!/usr/bin/env python3
"""Exact quadratic-algebraic norm probe for the content polynomial.

If alpha satisfies alpha^2-s*alpha+t=0, evaluate T_n(alpha)=u+v*alpha
by Horner reduction.  Its rational norm is

    u^2+s*u*v+t*v^2.

The coefficient content Gamma_n squared divides this norm.  The scan tests
small irreducible monic quadratics as possible denominator-free substitutes
for the low raw height seen at rational centers in (-1,0).
"""

from __future__ import annotations

from math import isqrt, log

from q32_legendre_content import franel_numbers, truncation_coefficients


INDICES = (60, 120, 180)
SEARCH_RADIUS = 12


def is_square(value: int) -> bool:
    return value >= 0 and isqrt(value) ** 2 == value


def quadratic_norm(coefficients: list[int], trace: int, norm: int) -> int:
    scalar = 0
    alpha_coefficient = 0
    for coefficient in reversed(coefficients):
        scalar, alpha_coefficient = (
            coefficient - norm * alpha_coefficient,
            scalar + trace * alpha_coefficient,
        )
    return (
        scalar * scalar
        + trace * scalar * alpha_coefficient
        + norm * alpha_coefficient * alpha_coefficient
    )


def main() -> None:
    franel = franel_numbers(max(INDICES))
    quadratics = [
        (trace, norm)
        for trace in range(-SEARCH_RADIUS, SEARCH_RADIUS + 1)
        for norm in range(-SEARCH_RADIUS, SEARCH_RADIUS + 1)
        if not is_square(trace * trace - 4 * norm)
    ]
    for n in INDICES:
        coefficients = truncation_coefficients(n, franel)
        candidates = []
        for trace, norm in quadratics:
            value = quadratic_norm(coefficients, trace, norm)
            if value:
                candidates.append(
                    (log(abs(value)) / (2 * n), trace, norm)
                )
        candidates.sort()
        best = candidates[0]
        integer_rate = log(abs(coefficients[0])) / n
        print(
            f"n={n} best_quadratic_rate={best[0]:.9f} "
            f"polynomial=x^2-{best[1]}x+{best[2]} "
            f"integer_center_zero_rate={integer_rate:.9f}"
        )


if __name__ == "__main__":
    main()
