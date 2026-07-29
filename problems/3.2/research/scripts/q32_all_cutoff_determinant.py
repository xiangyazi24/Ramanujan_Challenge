#!/usr/bin/env python3
"""Audit simultaneous-cutoff determinant amplification.

For every cutoff J, row subtraction changes the transform row T_J into

    h_J(c) = K_c(n,J) g_J(c).

A q=1 candidate prime with folded index j divides every coefficient of h_J
for J > j independently of whether the Apéry value A_j vanishes.  Hence an
all-cutoff determinant contains p^(H-j) universally; the bad condition can
add only one selective copy.  This script verifies that local obstruction
and measures the lowest/highest consecutive coefficient minors.
"""

from __future__ import annotations

from math import comb, log

from q32_strehl_gcd import franel_numbers, primes_up_to


LOWER_N = 10
UPPER_N = 40


def bareiss_determinant(matrix: list[list[int]]) -> int:
    work = [row[:] for row in matrix]
    size = len(work)
    if size == 0:
        return 1
    sign = 1
    previous = 1
    for pivot_index in range(size - 1):
        if work[pivot_index][pivot_index] == 0:
            replacement = next(
                (
                    row
                    for row in range(pivot_index + 1, size)
                    if work[row][pivot_index]
                ),
                None,
            )
            if replacement is None:
                return 0
            work[pivot_index], work[replacement] = (
                work[replacement],
                work[pivot_index],
            )
            sign = -sign
        pivot = work[pivot_index][pivot_index]
        for row in range(pivot_index + 1, size):
            for column in range(pivot_index + 1, size):
                numerator = (
                    work[row][column] * pivot
                    - work[row][pivot_index]
                    * work[pivot_index][column]
                )
                assert numerator % previous == 0
                work[row][column] = numerator // previous
        previous = pivot
    return sign * work[-1][-1]


def difference_rows(n: int, franel: list[int]) -> list[list[int]]:
    row_count = (n - 1) // 3 + 1
    rows: list[list[int]] = []
    for cutoff in range(row_count):
        kernel = [
            comb(n, degree)
            * comb(n + degree, degree)
            * comb(degree, cutoff)
            for degree in range(cutoff, n + 1)
        ]
        shifted_franel = [0] * (cutoff + 1)
        for index in range(cutoff + 1):
            shifted_franel[cutoff - index] = (
                (-1) ** (cutoff - index)
                * comb(cutoff, index)
                * franel[index]
            )
        coefficients = [0] * (n + 1)
        for left_degree, left in enumerate(kernel):
            for right_degree, right in enumerate(shifted_franel):
                coefficients[left_degree + right_degree] += left * right
        rows.append(coefficients)
    return rows


def main() -> None:
    franel = franel_numbers((UPPER_N - 1) // 3 + 1)
    primes = primes_up_to(UPPER_N)
    for n in range(LOWER_N, UPPER_N + 1):
        rows = difference_rows(n, franel)
        row_count = len(rows)

        for prime in primes:
            if not n / 2 < prime <= n:
                continue
            residue = n - prime
            folded = min(residue, prime - 1 - residue)
            assert folded < row_count
            for cutoff in range(folded + 1, row_count):
                assert all(
                    coefficient % prime == 0
                    for coefficient in rows[cutoff]
                ), (n, prime, folded, cutoff)

        low = abs(
            bareiss_determinant(
                [row[:row_count] for row in rows]
            )
        )
        high = abs(
            bareiss_determinant(
                [row[-row_count:] for row in rows]
            )
        )
        assert low and high
        print(
            f"n={n} rows={row_count} "
            f"low_rate={log(low)/(n*row_count):.9f} "
            f"high_rate={log(high)/(n*row_count):.9f}"
        )


if __name__ == "__main__":
    main()
