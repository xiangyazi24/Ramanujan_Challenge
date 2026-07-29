#!/usr/bin/env python3
"""Audit the linear-degree Krawtchouk filters proposed in ChatGPT Q575.

For the q=3 Legendre--Euler coefficient vector C_d and binary Krawtchouk
polynomials K_m(d;N), test both

    R^eps(n,N,m) = sum_d eps^d K_m(d;N) C_d.

When N is smaller than the least q=3 candidate prime, every candidate
satisfies p>N.  Hence
K_m(ap;N)=K_m(0;N)=binom(N,m) mod p, and

    R^- = -binom(N,m) A_j mod p,
    R^+ = 63 binom(N,m) A_j mod p.

The script verifies this selectivity and exhausts all 0<=m<=N<=J at a
representative set of n.  It is a diagnostic, not an asymptotic proof.
"""

from __future__ import annotations

from math import comb, log

from q32_fixed_q_content import truncation_coefficients
from q32_newton import apery_numbers
from q32_strehl_gcd import franel_numbers, primes_up_to


INDICES = (40, 60, 80, 100, 120, 160, 200, 240, 300, 400)


def filtered_value(
    coefficients: list[int], weights: list[int], sign: int
) -> int:
    return sum(
        sign**degree * weights[degree] * coefficients[degree]
        for degree in range(len(coefficients))
    )


def main() -> None:
    limit = max(INDICES)
    franel = franel_numbers(limit)
    apery = apery_numbers(limit)
    primes = primes_up_to(limit)

    for n in INDICES:
        coefficients = truncation_coefficients(n, 3, franel)
        cutoff = (n - 3) // 7
        candidates = [
            prime for prime in primes if divmod(n, prime)[0] == 3
        ]
        order_bound_limit = min(candidates) - 1
        baseline_minus = filtered_value(
            coefficients, [1] * (n + 1), -1
        )
        baseline_plus = filtered_value(
            coefficients, [1] * (n + 1), 1
        )
        best_minus = (
            log(abs(baseline_minus)) / n,
            0,
            0,
            baseline_minus.bit_length(),
        )
        best_plus = (
            log(abs(baseline_plus)) / n,
            0,
            0,
            baseline_plus.bit_length(),
        )

        for order_bound in range(1, order_bound_limit + 1):
            previous = [1] * (n + 1)
            current = [
                order_bound - 2 * degree for degree in range(n + 1)
            ]
            filters = [(0, previous), (1, current)]

            for order in range(1, order_bound):
                following: list[int] = []
                for degree in range(n + 1):
                    numerator = (
                        (order_bound - 2 * degree) * current[degree]
                        - (order_bound - order + 1) * previous[degree]
                    )
                    assert numerator % (order + 1) == 0
                    following.append(numerator // (order + 1))
                filters.append((order + 1, following))
                previous, current = current, following

            for order, weights in filters:
                minus_value = filtered_value(coefficients, weights, -1)
                plus_value = filtered_value(coefficients, weights, 1)
                assert weights[0] == comb(order_bound, order)
                for prime in candidates:
                    _, residue = divmod(n, prime)
                    folded = min(residue, prime - 1 - residue)
                    assert prime > order_bound
                    assert (
                        minus_value
                        + comb(order_bound, order) * apery[folded]
                    ) % prime == 0
                    assert (
                        plus_value
                        - 63 * comb(order_bound, order) * apery[folded]
                    ) % prime == 0
                if minus_value:
                    minus_record = (
                        log(abs(minus_value)) / n,
                        order_bound,
                        order,
                        minus_value.bit_length(),
                    )
                    if minus_record < best_minus:
                        best_minus = minus_record
                if plus_value:
                    plus_record = (
                        log(abs(plus_value)) / n,
                        order_bound,
                        order,
                        plus_value.bit_length(),
                    )
                    if plus_record < best_plus:
                        best_plus = plus_record

        print(
            f"n={n} J={cutoff} Nmax={order_bound_limit} "
            f"minus={best_minus[0]:.9f}@({best_minus[1]},"
            f"{best_minus[2]}) plus={best_plus[0]:.9f}@("
            f"{best_plus[1]},{best_plus[2]})"
        )


if __name__ == "__main__":
    main()
