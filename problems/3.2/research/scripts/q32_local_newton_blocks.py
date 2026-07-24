#!/usr/bin/env python3
"""Audit local Newton divisors for close q=1 hits.

Fix a block j=u,...,u+H and interpolate A_(u+h) in the integer-valued
Newton basis.  Direct hits use the common extrapolation point x=n-u;
reflected hits use x=-(n+1+u).  If m is below every candidate prime in the
block, a hit prime divides every extrapolant of degree m>=H.  Their gcd is
therefore an exact local package for all hits in the block.

The construction is potentially useful for pair amplification, but bounding
the resulting gcd is a new local version of the original horizontal radical
problem.  This script verifies the divisibility and measures the package.
"""

from __future__ import annotations

from math import comb, gcd, log

from q32_newton import apery_numbers
from q32_strehl_gcd import primes_up_to


INDICES = (120, 240, 321, 400, 600)
BLOCK_LENGTHS = (5, 10, 20, 40)


def generalized_binomial(value: int, order: int) -> int:
    if value >= 0:
        return comb(value, order)
    return (-1) ** order * comb(-value + order - 1, order)


def newton_differences(values: list[int]) -> list[int]:
    differences = []
    row = values
    while row:
        differences.append(row[0])
        row = [
            row[index + 1] - row[index]
            for index in range(len(row) - 1)
        ]
    return differences


def block_divisor(
    apery: list[int],
    n: int,
    start: int,
    block_length: int,
    branch: str,
) -> tuple[int, int]:
    if branch == "direct":
        evaluation = n - start
        least_candidate = n - start - block_length
    elif branch == "reflected":
        evaluation = -(n + 1 + start)
        # Candidate primes increase with the folded index.  This lower bound
        # is deliberately conservative and handles the parity restriction.
        least_candidate = (n + 1 + start) // 2
    else:
        raise ValueError(branch)

    maximum_degree = min(
        len(apery) - 1 - start,
        least_candidate - 1,
    )
    assert block_length <= maximum_degree
    differences = newton_differences(
        apery[start : start + maximum_degree + 1]
    )

    extrapolant = sum(
        differences[order]
        * generalized_binomial(evaluation, order)
        for order in range(block_length + 1)
    )
    divisor = abs(extrapolant)
    for order in range(block_length + 1, maximum_degree + 1):
        extrapolant += (
            differences[order]
            * generalized_binomial(evaluation, order)
        )
        divisor = gcd(divisor, extrapolant)
        if divisor == 1:
            break
    return divisor, maximum_degree


def main() -> None:
    limit = max(INDICES)
    apery = apery_numbers(limit)
    prime_set = set(primes_up_to(limit))

    for n in INDICES:
        cutoff = (n - 1) // 3
        print(f"n={n}")
        for block_length in BLOCK_LENGTHS:
            if block_length > cutoff:
                continue
            for branch in ("direct", "reflected"):
                total_log = 0.0
                actual_log = 0.0
                maximum_rate = 0.0
                for start in range(
                    0, cutoff + 1, block_length + 1
                ):
                    local_length = min(
                        block_length, cutoff - start
                    )
                    divisor, _ = block_divisor(
                        apery,
                        n,
                        start,
                        local_length,
                        branch,
                    )
                    if divisor > 1:
                        local_log = log(divisor)
                        total_log += local_log
                        maximum_rate = max(
                            maximum_rate, local_log / n
                        )

                    for offset in range(local_length + 1):
                        index = start + offset
                        if branch == "direct":
                            prime = n - index
                            eligible = (
                                prime in prime_set
                                and prime > 2 * index
                            )
                        else:
                            doubled = n + 1 + index
                            prime = doubled // 2
                            eligible = (
                                doubled % 2 == 0
                                and prime in prime_set
                                and prime > 2 * index
                            )
                        if eligible and apery[index] % prime == 0:
                            assert divisor % prime == 0
                            actual_log += log(prime)

                print(
                    f"{branch} H={block_length} "
                    f"sum_log_divisor/n={total_log / n:.9f} "
                    f"actual_hit_log/n={actual_log / n:.9f} "
                    f"max_block_rate={maximum_rate:.9f}"
                )


if __name__ == "__main__":
    main()
