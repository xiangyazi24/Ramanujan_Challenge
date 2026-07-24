#!/usr/bin/env python3
"""Verify the local zero-absorber theorem for primitive Padé numerators.

For height H and numerator degree a, let (P_{H,a},Q_{H,a}) be the
primitive integral interpolation pair

    P_{H,a}(s) = A_s Q_{H,a}(s),  0 <= s <= H,
    deg P_{H,a} <= a,  deg Q_{H,a} <= H-a.

For a prime p>H, write Z={s<=H:p|A_s}.  Reduction over F_p gives

    P_{H,a} is the zero polynomial  iff  |Z|>a.

If |Z|<=a, its roots on the interpolation nodes consist of Z and at most
a-|Z| extra nodes; at a=|Z| they are exactly Z.  Consequently the common
candidate-prime support of degrees 0,...,A is the target support plus
only primes with |Z|>A.

The theorem is elementary.  This script checks primitive normalization
and the common-support corollary against the exact Padé implementation.
"""

from __future__ import annotations

import argparse

from q32_adjacent_pade_kappa import (
    apery_values,
    evaluate_newton,
    is_prime,
    newton_coefficients,
    primitive_pair,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("maximum_height", nargs="?", type=int, default=22)
    args = parser.parse_args()
    if args.maximum_height < 2:
        raise SystemExit("maximum_height must be at least 2")

    apery = apery_values(3 * args.maximum_height + 3)
    differences = newton_coefficients(apery)
    polynomial_checks = 0
    family_checks = 0

    for height in range(2, args.maximum_height + 1):
        n = 3 * height + 1
        numerators = [
            primitive_pair(
                height,
                height - numerator_degree,
                differences,
            )[0]
            for numerator_degree in range(height + 1)
        ]
        values_at_n = [
            evaluate_newton(numerator, n) for numerator in numerators
        ]

        for prime in range(height + 1, n + 1):
            if not is_prime(prime):
                continue
            zero_nodes = [
                node
                for node in range(height + 1)
                if apery[node] % prime == 0
            ]
            zero_count = len(zero_nodes)

            for numerator_degree, numerator in enumerate(numerators):
                zero_polynomial = all(
                    coefficient % prime == 0 for coefficient in numerator
                )
                assert zero_polynomial == (
                    zero_count > numerator_degree
                ), (
                    height,
                    prime,
                    numerator_degree,
                    zero_nodes,
                )
                if not zero_polynomial:
                    roots = [
                        node
                        for node in range(height + 1)
                        if evaluate_newton(numerator, node) % prime == 0
                    ]
                    assert set(zero_nodes) <= set(roots)
                    assert len(roots) <= numerator_degree
                    if numerator_degree == zero_count:
                        assert roots == zero_nodes
                polynomial_checks += 1

            candidate_node = n - prime
            if not 0 <= candidate_node <= height:
                continue
            common_so_far = True
            for cutoff, value in enumerate(values_at_n):
                common_so_far = common_so_far and value % prime == 0
                expected = (
                    candidate_node in zero_nodes or zero_count > cutoff
                )
                assert common_so_far == expected, (
                    height,
                    prime,
                    cutoff,
                    candidate_node,
                    zero_nodes,
                )
                family_checks += 1

    print(
        f"heights=2..{args.maximum_height} "
        f"polynomial_checks={polynomial_checks} "
        f"family_support_checks={family_checks}"
    )


if __name__ == "__main__":
    main()
