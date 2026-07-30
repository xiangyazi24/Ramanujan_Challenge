#!/usr/bin/env python3
"""Audit the precision-p^8 direct/reflected target elimination.

The proof uses only the exact shifted Apéry recurrence decompositions,
the already-proved precision-p^7 endpoint laws, and the mod-p
reflection of the distinguished Apéry solution.  A new endpoint digit
``v_p`` occurs, while the new target-dependent companion digit cancels
between the direct and reflected rows.
"""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from math import gcd

import sympy as sp

import q32_order_seven_target_audit as order_seven


OLD_ORDER = order_seven.ORDER
ORDER = 9
ENDPOINT_NORMALIZATION_PRIME = 769
TARGET_UNIT_PRIME = 18461


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prime-limit", type=int, default=1000)
    return parser.parse_args()


def endpoint_scalars(
    prime: int, values: list[int]
) -> tuple[int, int, Fraction, Fraction]:
    delta = values[prime - 1] - 1
    endpoint_h = values[prime] - 5 + 7 * delta
    endpoint_w = (
        Fraction(values[2 * prime] - 73 + 824 * delta)
        - Fraction(752, 5) * endpoint_h
    ) / prime**6
    endpoint_v = (
        Fraction(values[2 * prime - 1] - 5 - 8 * delta)
        - Fraction(336, 5) * endpoint_h
        + Fraction(103, ENDPOINT_NORMALIZATION_PRIME)
        * prime**6
        * endpoint_w
    ) / prime**7
    return delta, endpoint_h, endpoint_w, endpoint_v


def symbolic_checks() -> tuple[bool, bool]:
    x, h, w, v, j, p = sp.symbols("x h w v j p")
    direct = x * (1 - h / 5) - p**2 * h * j / 5
    reflected = (
        x
        * (
            1
            - sp.Rational(336, 25) * h
            + sp.Rational(103, 5 * ENDPOINT_NORMALIZATION_PRIME)
            * p**6
            * w
            - p**7 * v / 5
        )
        + sp.Rational(166144, 25) * p**2 * h * j
    )
    companion_elimination = sp.expand(
        166144 * direct
        + 5 * reflected
        - x
        * (
            166149
            - 33296 * h
            + sp.Rational(103, ENDPOINT_NORMALIZATION_PRIME)
            * p**6
            * w
            - p**7 * v
        )
    ) == 0

    correction = (
        33296 * h
        - sp.Rational(103, ENDPOINT_NORMALIZATION_PRIME)
        * p**6
        * w
        + p**7 * v
    )
    fixed_residual = sp.expand(
        (166144 + correction) * direct
        + 5 * reflected
        - 166149 * x
    )
    # The correction times (direct-x) has valuation at least ten:
    # this exact symbolic quotient is what is discarded modulo p^8.
    expected_residual = sp.expand(correction * (direct - x))
    fixed_elimination = fixed_residual == expected_residual
    return companion_elimination, fixed_elimination


def main() -> None:
    args = parse_args()
    primes = [
        prime
        for prime in order_seven.primes_at_most(args.prime_limit)
        if prime >= 11
        and prime != ENDPOINT_NORMALIZATION_PRIME
    ]
    maximum_prime = max(primes)
    values = order_seven.apery_numbers(2 * maximum_prime)

    # The shifted-solution routines read ORDER as a module global.
    order_seven.ORDER = ORDER
    try:
        direct_jets, companion_jets = (
            order_seven.shifted_fundamental_solutions(
                maximum_prime - 1
            )
        )
        j_jets = [
            order_seven.j_polynomial(companion_jets, index)
            for index in range(maximum_prime)
        ]
    finally:
        order_seven.ORDER = OLD_ORDER

    checks: Counter[str] = Counter()
    target_unit_exceptions: list[tuple[int, int]] = []

    for prime in primes:
        modulus8 = prime**8
        modulus9 = prime**9
        delta, endpoint_h, endpoint_w, endpoint_v = endpoint_scalars(
            prime, values
        )
        assert delta % prime**3 == 0
        assert endpoint_h % prime**5 == 0
        assert gcd(endpoint_w.denominator, prime) == 1
        assert gcd(endpoint_v.denominator, prime) == 1
        checks["endpoint_integrality"] += 2

        for index in range(prime):
            if values[index] % prime != 0:
                continue
            reflected_index = prime - 1 - index
            upper_index = prime + index
            assert 0 < index < prime - 1
            assert values[reflected_index] % prime == 0
            assert values[upper_index] % prime == 0
            x_value = values[upper_index] // prime % modulus8

            direct_u = order_seven.evaluate_mod(
                direct_jets[index], prime, modulus9
            )
            direct_j = order_seven.evaluate_mod(
                j_jets[index][:6], prime, modulus9
            )
            direct_numerator = (
                (5 - 7 * delta) * direct_u
                - prime**3 * (1 + delta) * direct_j
            ) % modulus9
            assert direct_numerator % prime == 0
            direct_residue = direct_numerator // prime % modulus8

            reflected_u = order_seven.evaluate_mod(
                direct_jets[reflected_index],
                -2 * prime,
                modulus9,
            )
            reflected_j = order_seven.evaluate_mod(
                j_jets[reflected_index][:6],
                -2 * prime,
                modulus9,
            )
            reflected_numerator = (
                (5 + 8 * delta) * reflected_u
                + 8
                * prime**3
                * (73 - 824 * delta)
                * reflected_j
            ) % modulus9
            assert reflected_numerator % prime == 0
            reflected_residue = (
                reflected_numerator // prime % modulus8
            )

            direct_companion = order_seven.fraction_mod(
                j_jets[index][0], prime
            )
            reflected_companion = order_seven.fraction_mod(
                j_jets[reflected_index][0], prime
            )
            assert direct_companion == reflected_companion
            checks["companion_reflection"] += 1

            direct_expected = (
                x_value
                * order_seven.fraction_mod(
                    Fraction(1) - Fraction(endpoint_h, 5),
                    modulus8,
                )
                - order_seven.fraction_mod(
                    Fraction(
                        prime**2
                        * endpoint_h
                        * direct_companion,
                        5,
                    ),
                    modulus8,
                )
            ) % modulus8
            assert direct_residue == direct_expected
            checks["direct_equation"] += 1

            reflected_factor = (
                Fraction(1)
                - Fraction(336, 25) * endpoint_h
                + Fraction(
                    103 * prime**6,
                    5 * ENDPOINT_NORMALIZATION_PRIME,
                )
                * endpoint_w
                - Fraction(prime**7, 5) * endpoint_v
            )
            reflected_expected = (
                x_value
                * order_seven.fraction_mod(
                    reflected_factor, modulus8
                )
                + order_seven.fraction_mod(
                    Fraction(
                        166144
                        * prime**2
                        * endpoint_h
                        * reflected_companion,
                        25,
                    ),
                    modulus8,
                )
            ) % modulus8
            assert reflected_residue == reflected_expected
            checks["reflected_equation"] += 1

            correction = (
                Fraction(33296 * endpoint_h)
                - Fraction(
                    103 * prime**6,
                    ENDPOINT_NORMALIZATION_PRIME,
                )
                * endpoint_w
                + Fraction(prime**7) * endpoint_v
            )
            fixed_left = (
                order_seven.fraction_mod(
                    Fraction(166144) + correction,
                    modulus8,
                )
                * direct_residue
                + 5 * reflected_residue
            ) % modulus8
            fixed_right = 166149 * x_value % modulus8
            assert fixed_left == fixed_right
            checks["fixed_target_law"] += 1

            small_w = order_seven.fraction_mod(
                endpoint_w, prime**2
            )
            small_v = order_seven.fraction_mod(endpoint_v, prime)
            small_correction = (
                Fraction(33296 * endpoint_h)
                - Fraction(
                    103 * prime**6 * small_w,
                    ENDPOINT_NORMALIZATION_PRIME,
                )
                + Fraction(prime**7 * small_v)
            )
            small_left = (
                order_seven.fraction_mod(
                    Fraction(166144) + small_correction,
                    modulus8,
                )
                * direct_residue
                + 5 * reflected_residue
            ) % modulus8
            assert small_left == fixed_left
            checks["small_endpoint_representatives"] += 1

            if prime == TARGET_UNIT_PRIME:
                target_unit_exceptions.append((prime, index))

    companion_symbolic, fixed_symbolic = symbolic_checks()
    assert companion_symbolic
    assert fixed_symbolic
    print(f"prime_limit={args.prime_limit}")
    print(
        "symbolic_companion_elimination="
        f"{int(companion_symbolic)}/1"
    )
    print(
        "symbolic_fixed_elimination="
        f"{int(fixed_symbolic)}/1"
    )
    for name in sorted(checks):
        print(f"{name}_checks={checks[name]}")
    print(
        "endpoint_normalization_exception="
        f"{ENDPOINT_NORMALIZATION_PRIME}"
    )
    print(
        "target_unit_exceptions="
        f"3,{ENDPOINT_NORMALIZATION_PRIME},"
        f"{TARGET_UNIT_PRIME}"
    )
    print(
        "observed_rows_at_target_unit_prime="
        f"{len(target_unit_exceptions)}"
    )
    print("failures=0")


if __name__ == "__main__":
    main()
