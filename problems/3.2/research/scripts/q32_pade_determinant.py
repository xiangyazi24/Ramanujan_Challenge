#!/usr/bin/env python3
"""Exact Newton--Padé determinant audit for the q=1 Apéry support.

Rows are x=0,...,M-1,n.  Columns are binom(x,k), 0<=k<=a, followed by
A_x*binom(x,k), 0<=k<=b, where a+b=M-1.  A top-half bad prime makes the
last row congruent to the row x=n-p, except for the harmless odd boundary.

For each determinant we remove its entire M-smooth part.  This is more
favorable than removing only a known smooth factor, and directly tests
whether the remaining rough part could have subexponential height.
"""

from __future__ import annotations

from math import comb, gcd, log

from sympy import Matrix, primerange

from q32_newton import apery_numbers, evaluate_newton, forward_differences


SAMPLES = (20, 30, 40)


def determinant(apery: list[int], n: int, degree_a: int) -> int:
    half = n // 2
    degree_b = half - 1 - degree_a
    rows = []
    for value in [*range(half), n]:
        rows.append(
            [comb(value, k) for k in range(degree_a + 1)]
            + [
                apery[value] * comb(value, k)
                for k in range(degree_b + 1)
            ]
        )
    return abs(int(Matrix(rows).det(method="domain-ge")))


def rough_part(value: int, bound: int) -> int:
    result = value
    for prime in primerange(2, bound + 1):
        while result and result % prime == 0:
            result //= prime
    return result


def full_minor_gcd(apery: list[int], n: int) -> int:
    """Gcd of all maximal minors of the full binomial/Apéry matrix."""
    half = n // 2
    result = 0
    for degree in range(half):
        values = [apery[x] * comb(x, degree) for x in range(half)]
        coefficients = forward_differences(values)
        residual = (
            apery[n] * comb(n, degree)
            - evaluate_newton(coefficients, n)
        )
        result = gcd(result, abs(residual))
    return result


def main() -> None:
    apery = apery_numbers(max(SAMPLES))
    for n in SAMPLES:
        half = n // 2
        bad_primes = [
            prime
            for prime in primerange(half + 1, n + 1)
            if apery[n] % prime == 0 and n != 2 * prime - 1
        ]
        records = []
        family_gcd = 0
        endpoint = 0
        adjacent = 0
        for degree_a in range(half):
            degree_b = half - 1 - degree_a
            value = determinant(apery, n, degree_a)
            assert all(value % prime == 0 for prime in bad_primes)
            family_gcd = gcd(family_gcd, value)
            if degree_a == half - 1:
                endpoint = value
            elif degree_a == half - 2:
                adjacent = value
            rough = rough_part(value, half)
            records.append(
                (
                    log(rough) / n if rough > 1 else 0.0,
                    degree_a,
                    degree_b,
                    log(value) / n if value else float("-inf"),
                )
            )

        print(f"n={n}")
        two_gcd = gcd(endpoint, adjacent)
        assert family_gcd == two_gcd
        determinantal_divisor = full_minor_gcd(apery, n)
        assert family_gcd % determinantal_divisor == 0
        assert all(determinantal_divisor % prime == 0 for prime in bad_primes)
        print(
            f"  family_gcd={family_gcd} "
            f"log_family_gcd_over_n="
            f"{log(family_gcd) / n if family_gcd > 1 else 0.0:.9f}"
        )
        print(
            f"  full_minor_gcd={determinantal_divisor} "
            f"log_full_minor_gcd_over_n="
            f"{log(determinantal_divisor) / n if determinantal_divisor > 1 else 0.0:.9f}"
        )
        for rough_rate, degree_a, degree_b, full_rate in sorted(records)[:6]:
            print(
                f"  degrees=({degree_a},{degree_b}) "
                f"log_det_over_n={full_rate:.9f} "
                f"log_rough_over_n={rough_rate:.9f}"
            )


if __name__ == "__main__":
    main()
