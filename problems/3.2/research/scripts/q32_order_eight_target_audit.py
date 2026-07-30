#!/usr/bin/env python3
"""Audit the denominator-free precision-p^8 target elimination.

The proof uses only the exact shifted Apéry recurrence decompositions,
the already-proved precision-p^7 endpoint laws, and the mod-p
reflection of the distinguished Apéry solution.  The target law uses
the single raw reflected endpoint residual

    F_p = A_(2p-1) - 5 - 8 Delta_p - (336/5) H_p.

This avoids the artificial p=769 exception introduced by splitting
F_p into the normalized coordinates w_p and v_p.
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
    direct_residual = (
        Fraction(values[2 * prime] - 73 + 824 * delta)
        - Fraction(752, 5) * endpoint_h
    )
    reflected_residual = (
        Fraction(values[2 * prime - 1] - 5 - 8 * delta)
        - Fraction(336, 5) * endpoint_h
    )
    return delta, endpoint_h, direct_residual, reflected_residual


def symbolic_checks() -> tuple[bool, bool, bool]:
    x, h, endpoint_f, j, p = sp.symbols("x h endpoint_f j p")
    direct = x * (1 - h / 5) - p**2 * h * j / 5
    reflected = (
        x * (1 - sp.Rational(336, 25) * h - endpoint_f / 5)
        + sp.Rational(166144, 25) * p**2 * h * j
    )
    companion_elimination = sp.expand(
        166144 * direct
        + 5 * reflected
        - x
        * (
            166149
            - 33296 * h
            - endpoint_f
        )
    ) == 0

    correction = 33296 * h + endpoint_f
    fixed_residual = sp.expand(
        (166144 + correction) * direct
        + 5 * reflected
        - 166149 * x
    )
    # The correction times (direct-x) has valuation at least ten:
    # this exact symbolic quotient is what is discarded modulo p^8.
    expected_residual = sp.expand(correction * (direct - x))
    fixed_elimination = fixed_residual == expected_residual

    endpoint_e, endpoint_v = sp.symbols("endpoint_e endpoint_v")
    split_f = (
        -sp.Rational(103, ENDPOINT_NORMALIZATION_PRIME) * endpoint_e
        + p**7 * endpoint_v
    )
    normalized_equivalence = sp.expand(
        correction.subs(endpoint_f, split_f)
        - (
            33296 * h
            - sp.Rational(103, ENDPOINT_NORMALIZATION_PRIME)
            * endpoint_e
            + p**7 * endpoint_v
        )
    ) == 0
    return (
        companion_elimination,
        fixed_elimination,
        normalized_equivalence,
    )


def main() -> None:
    args = parse_args()
    primes = [
        prime
        for prime in order_seven.primes_at_most(args.prime_limit)
        if prime >= 11
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
    normalization_prime_rows = 0

    for prime in primes:
        modulus8 = prime**8
        modulus9 = prime**9
        (
            delta,
            endpoint_h,
            direct_endpoint_residual,
            reflected_endpoint_residual,
        ) = endpoint_scalars(prime, values)
        assert delta % prime**3 == 0
        assert endpoint_h % prime**5 == 0
        endpoint_w = direct_endpoint_residual / prime**6
        reflected_digit = reflected_endpoint_residual / prime**6
        endpoint_bar_v = (
            ENDPOINT_NORMALIZATION_PRIME
            * reflected_endpoint_residual
            + 103 * direct_endpoint_residual
        ) / prime**7
        assert gcd(endpoint_w.denominator, prime) == 1
        assert gcd(reflected_digit.denominator, prime) == 1
        assert gcd(endpoint_bar_v.denominator, prime) == 1
        checks["endpoint_integrality"] += 3

        normalized_delta = delta // prime**3 % prime
        assert order_seven.fraction_mod(
            endpoint_w, prime
        ) == order_seven.fraction_mod(
            Fraction(24 * ENDPOINT_NORMALIZATION_PRIME, 5)
            * normalized_delta**2,
            prime,
        )
        assert order_seven.fraction_mod(
            reflected_digit, prime
        ) == order_seven.fraction_mod(
            -Fraction(24 * 103, 5) * normalized_delta**2,
            prime,
        )
        checks["endpoint_first_digit"] += 2

        endpoint_v: Fraction | None = None
        if prime != ENDPOINT_NORMALIZATION_PRIME:
            endpoint_v = (
                reflected_endpoint_residual
                + Fraction(103, ENDPOINT_NORMALIZATION_PRIME)
                * direct_endpoint_residual
            ) / prime**7
            assert gcd(endpoint_v.denominator, prime) == 1
            assert (
                reflected_endpoint_residual
                == -Fraction(103, ENDPOINT_NORMALIZATION_PRIME)
                * direct_endpoint_residual
                + prime**7 * endpoint_v
            )
            checks["normalized_v_integrality"] += 1

        companion_values = [
            j_jets[index][0] for index in range(prime)
        ]
        assert (
            values[prime - 2] * companion_values[prime - 1]
            - values[prime - 1] * companion_values[prime - 2]
            == Fraction(1, (prime - 1) ** 3)
        )
        checks["casoratian"] += 1
        for index in range(prime):
            universal_residual = (
                companion_values[prime - 1 - index]
                - companion_values[index]
                - companion_values[prime - 1] * values[index]
            )
            assert (
                order_seven.fraction_mod(universal_residual, prime)
                == 0
            )
            checks["universal_companion_reflection"] += 1

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

            forward_decomposition = (
                values[prime] * direct_u
                - prime**3 * values[prime - 1] * direct_j
            ) % modulus9
            reflected_decomposition = (
                values[2 * prime - 1] * reflected_u
                + 8
                * prime**3
                * values[2 * prime]
                * reflected_j
            ) % modulus9
            assert forward_decomposition == values[upper_index] % modulus9
            assert (
                reflected_decomposition
                == values[upper_index] % modulus9
            )
            checks["exact_shifted_decomposition"] += 2

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
                - reflected_endpoint_residual / 5
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
                + reflected_endpoint_residual
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

            integer_coefficient = (
                5 * values[2 * prime - 1]
                + 830695
                - 40 * delta
                + 166144 * endpoint_h
            )
            assert Fraction(integer_coefficient, 5) == (
                Fraction(166144) + correction
            )
            assert (
                integer_coefficient * direct_residue
                + 25 * reflected_residue
            ) % modulus8 == 830745 * x_value % modulus8
            checks["integer_fixed_target_law"] += 1

            small_reflected_digit = order_seven.fraction_mod(
                reflected_digit, prime**2
            )
            small_correction = (
                Fraction(33296 * endpoint_h)
                + prime**6 * small_reflected_digit
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

            if endpoint_v is not None:
                normalized_correction = (
                    Fraction(33296 * endpoint_h)
                    - Fraction(
                        103 * prime**6,
                        ENDPOINT_NORMALIZATION_PRIME,
                    )
                    * endpoint_w
                    + Fraction(prime**7) * endpoint_v
                )
                assert normalized_correction == correction
                checks["normalized_wv_equivalence"] += 1

            if prime == TARGET_UNIT_PRIME:
                target_unit_exceptions.append((prime, index))
            if prime == ENDPOINT_NORMALIZATION_PRIME:
                normalization_prime_rows += 1

    (
        companion_symbolic,
        fixed_symbolic,
        normalized_symbolic,
    ) = symbolic_checks()
    assert companion_symbolic
    assert fixed_symbolic
    assert normalized_symbolic
    print(f"prime_limit={args.prime_limit}")
    print(
        "symbolic_companion_elimination="
        f"{int(companion_symbolic)}/1"
    )
    print(
        "symbolic_fixed_elimination="
        f"{int(fixed_symbolic)}/1"
    )
    print(
        "symbolic_normalized_equivalence="
        f"{int(normalized_symbolic)}/1"
    )
    for name in sorted(checks):
        print(f"{name}_checks={checks[name]}")
    print(
        "wv_normalization_exception="
        f"{ENDPOINT_NORMALIZATION_PRIME}"
    )
    print(
        "fixed_inversion_exception="
        f"{TARGET_UNIT_PRIME}"
    )
    print(
        "observed_rows_at_target_unit_prime="
        f"{len(target_unit_exceptions)}"
    )
    print(
        "rows_recovered_at_wv_normalization_prime="
        f"{normalization_prime_rows}"
    )
    print("failures=0")


if __name__ == "__main__":
    main()
