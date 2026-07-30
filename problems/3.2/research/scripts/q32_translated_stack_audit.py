#!/usr/bin/env python3
"""Exact translated-stack audits for the Section 56 hostile blocks.

For a core I=[D,D+N-1] and m>=2, put

    d0 = D-m+1,       L = N+m-2,
    K_m(I) = gcd(G_{d0+t,L}: 0<=t<m).

Every stencil contains I, so every target node in I survives.  Consecutive
subtraction gives the exact Euclidean normal form

    K_m = gcd(G_{d0,L},
              B_t Delta^(L+1) C_M(d0+t): 0<=t<m-1),
    B_t = binom(D+N+t,L).

The script verifies this identity, the finite-window Pascal normal form for
gcd(B_t), and the exact nuisance histories for the five hostile blocks.
Only Python's standard library is used.
"""

from argparse import ArgumentParser
from math import comb, gcd, isqrt, prod

from q32_cartier_packet_audit import shell_batch
from q32_newton_gcd_audit import (
    carrier_from_values,
    forward_difference,
)


ROWS = (
    (200, 128, 63, (139, 181), 47),
    (272, 180, 63, (191, 233), 1),
    (300, 180, 57, (191, 227), 1),
    (321, 168, 53, (179, 193, 211), 43),
    (755, 582, 161, (593, 733), 275),
)


def gcd_many(values):
    out = 0
    for value in values:
        out = gcd(out, abs(value))
    return out


def factor_small(number):
    """Exact trial-division factorization for the observed small nuisances."""

    number = abs(number)
    factors = {}
    divisor = 2
    while divisor * divisor <= number:
        while number % divisor == 0:
            factors[divisor] = factors.get(divisor, 0) + 1
            number //= divisor
        divisor = 3 if divisor == 2 else divisor + 2
    if number > 1:
        factors[number] = factors.get(number, 0) + 1
    return factors


def first_cell_cap(moment, D, N, requested):
    """Largest m whose entire translated union stays in M/2<d<=M."""

    left_cap = D - moment // 2
    right_cap = moment - D - N + 2
    return min(requested, left_cap, right_cap)


def prefactor_pascal_gcd(D, N, m):
    """Two equal presentations of gcd_t binom(D+N+t,N+m-2)."""

    length = N + m - 2
    translated = gcd_many(
        comb(D + N + shift, length)
        for shift in range(m - 1)
    )
    pascal = gcd_many(
        comb(D + N, length - shift)
        for shift in range(m - 1)
    )
    assert translated == pascal
    return translated


def analyze_row(n, D, N, targets, expected_m2, requested_m):
    moment = n - 1
    max_m = first_cell_cap(moment, D, N, requested_m)
    assert max_m >= 2

    lower = D - max_m + 1
    upper = D + N + max_m - 2
    values = shell_batch(moment, range(lower, upper + 1))
    assert all(moment // node == 1 for node in values)

    target_product = prod(targets)
    history = []
    first_clean = None
    for multiplicity in range(2, max_m + 1):
        length = N + multiplicity - 2
        d0 = D - multiplicity + 1
        carriers = [
            carrier_from_values(values, d0 + shift, length)
            for shift in range(multiplicity)
        ]
        stack_gcd = gcd_many(carriers)
        assert stack_gcd % target_product == 0
        nuisance = stack_gcd // target_product

        prefactors = [
            comb(D + N + shift, length)
            for shift in range(multiplicity - 1)
        ]
        high_differences = []
        for shift in range(multiplicity - 1):
            node = d0 + shift
            difference = forward_difference(
                [
                    values[index]
                    for index in range(node, node + length + 2)
                ],
                length + 1,
            )[0]
            high_differences.append(difference)
            assert carriers[shift] - carriers[shift + 1] == (
                (-1) ** (length + 1)
                * prefactors[shift]
                * difference
            )

        euclidean_gcd = gcd_many(
            [carriers[0]]
            + [
                prefactor * difference
                for prefactor, difference in zip(
                    prefactors, high_differences
                )
            ]
        )
        assert euclidean_gcd == stack_gcd

        prefactor_gcd = gcd_many(prefactors)
        assert prefactor_gcd == prefactor_pascal_gcd(
            D, N, multiplicity
        )
        nuisance_factors = factor_small(nuisance)
        classes = {
            prime: (
                "target-excess"
                if prime in targets
                else (
                    "common-prefactor"
                    if prefactor_gcd % prime == 0
                    else "high-difference"
                )
            )
            for prime in nuisance_factors
        }
        history.append(
            (
                multiplicity,
                nuisance,
                nuisance_factors,
                classes,
            )
        )
        if nuisance == 1 and first_clean is None:
            first_clean = multiplicity

    assert history[0][1] == expected_m2
    print(
        "STACK",
        n,
        "core",
        (D, N),
        "targets",
        targets,
        "max_m",
        max_m,
        "first_clean",
        first_clean,
    )
    for record in history:
        print(" ", record)
    return history


def audit_large_prefactor_classification():
    """Check the exact p>L common-multiple criterion on small rectangles."""

    for D in range(20, 45):
        for N in range(3, 9):
            for multiplicity in range(2, 7):
                length = N + multiplicity - 2
                prefactor_gcd = prefactor_pascal_gcd(
                    D, N, multiplicity
                )
                for prime in range(length + 1, D + N + 2):
                    if any(
                        prime % divisor == 0
                        for divisor in range(2, isqrt(prime) + 1)
                    ):
                        continue
                    common_multiple = any(
                        value % prime == 0
                        for value in range(D + 1, D + N + 1)
                    )
                    assert (
                        prefactor_gcd % prime == 0
                    ) == common_multiple


def main():
    parser = ArgumentParser()
    parser.add_argument("--max-m", type=int, default=12)
    args = parser.parse_args()

    audit_large_prefactor_classification()
    for row in ROWS:
        analyze_row(*row, args.max_m)
    print("PASS: translated-stack Euclidean/Smith nuisance audit")


if __name__ == "__main__":
    main()
