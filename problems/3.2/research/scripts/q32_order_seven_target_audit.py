#!/usr/bin/env python3
"""Audit the order-seven target equations and the nonlinear W removal.

The original two-row linear audit treated W as independent.  H6
actually gives W=(24*769/5)(Delta/p^3)^2 mod p, so multiplying the
direct row by Delta^2 removes W and restores a fixed target law.
"""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from math import comb, gcd, isqrt

import sympy as sp

ORDER = 8


def primes_at_most(limit: int) -> list[int]:
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[:2] = b"\x00\x00"
    for prime in range(2, isqrt(limit) + 1):
        if sieve[prime]:
            sieve[prime * prime : limit + 1 : prime] = b"\x00" * (
                (limit - prime * prime) // prime + 1
            )
    return [
        candidate
        for candidate in range(2, limit + 1)
        if sieve[candidate]
    ]


def apery_numbers(limit: int) -> list[int]:
    values = [1, 5]
    for index in range(1, limit):
        numerator = (
            recurrence_polynomial(index) * values[index]
            - index**3 * values[index - 1]
        )
        denominator = (index + 1) ** 3
        assert numerator % denominator == 0
        values.append(numerator // denominator)
    return values


def recurrence_polynomial(index: int) -> int:
    return 34 * index**3 + 51 * index**2 + 27 * index + 5


def add(
    left: list[Fraction], right: list[Fraction]
) -> list[Fraction]:
    return [left[index] + right[index] for index in range(ORDER)]


def scale(
    polynomial: list[Fraction], scalar: int
) -> list[Fraction]:
    return [scalar * coefficient for coefficient in polynomial]


def multiply(
    left: list[Fraction], right: list[Fraction]
) -> list[Fraction]:
    result = [Fraction(0) for _ in range(ORDER)]
    for left_index, left_coefficient in enumerate(left):
        for right_index, right_coefficient in enumerate(right):
            degree = left_index + right_index
            if degree >= ORDER:
                break
            result[degree] += left_coefficient * right_coefficient
    return result


def inverse_shift_cube(shift: int) -> list[Fraction]:
    return [
        Fraction(
            (-1) ** degree * comb(degree + 2, 2),
            shift ** (degree + 3),
        )
        for degree in range(ORDER)
    ]


def shifted_fundamental_solutions(
    limit: int,
) -> tuple[list[list[Fraction]], list[list[Fraction]]]:
    zero = [Fraction(0) for _ in range(ORDER)]
    one = zero.copy()
    one[0] = Fraction(1)
    b_at_x = [
        Fraction(5),
        Fraction(27),
        Fraction(51),
        Fraction(34),
    ]
    direct = [one, multiply(b_at_x, inverse_shift_cube(1))]
    companion = [zero, one]

    for index in range(1, limit):
        shifted_b = [
            Fraction(recurrence_polynomial(index)),
            Fraction(102 * index**2 + 102 * index + 27),
            Fraction(102 * index + 51),
            Fraction(34),
        ]
        shifted_cube = [
            Fraction(index**3),
            Fraction(3 * index**2),
            Fraction(3 * index),
            Fraction(1),
        ]
        inverse_denominator = inverse_shift_cube(index + 1)
        for sequence in (direct, companion):
            numerator = add(
                multiply(shifted_b, sequence[index]),
                scale(
                    multiply(shifted_cube, sequence[index - 1]),
                    -1,
                ),
            )
            sequence.append(
                multiply(numerator, inverse_denominator)
            )
    return direct, companion


def fraction_mod(value: Fraction, modulus: int) -> int:
    assert gcd(value.denominator, modulus) == 1
    return (
        value.numerator
        * pow(value.denominator, -1, modulus)
        % modulus
    )


def evaluate_mod(
    polynomial: list[Fraction], argument: int, modulus: int
) -> int:
    return sum(
        fraction_mod(coefficient, modulus)
        * pow(argument, degree, modulus)
        for degree, coefficient in enumerate(polynomial)
    ) % modulus


def j_polynomial(
    companion: list[list[Fraction]], index: int
) -> list[Fraction]:
    result = [Fraction(0) for _ in range(ORDER)]
    for degree in range(ORDER):
        result[degree] = sum(
            companion[index][source_degree]
            * (-1) ** (degree - source_degree)
            * comb(degree - source_degree + 2, 2)
            for source_degree in range(degree + 1)
        )
    return result


def divide_residue(
    numerator: int, prime: int, modulus_before: int
) -> int:
    numerator %= modulus_before
    assert numerator % prime == 0
    return numerator // prime


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prime-limit", type=int, default=1000)
    return parser.parse_args()


def symbolic_correction_checks() -> tuple[bool, bool]:
    y, z = sp.symbols("y z")
    x_square = 5 * y + 2 * z
    endpoint_raw = 935 * x_square - 830 * y - 332 * z
    w_square = sp.expand(endpoint_raw - 769 * x_square) == 0

    x, h, delta_square = sp.symbols("x h delta_square")
    direct = x * (1 - h / 5)
    reflected = x * (
        1 - sp.Rational(336, 25) * h
        + sp.Rational(2472, 25) * delta_square
    )
    fixed = (
        (1680 + 2472 * delta_square) * direct
        - 25 * reflected
    )
    target_residual = sp.expand(fixed - 1655 * x)
    nonlinear_elimination = (
        target_residual
        == -sp.Rational(2472, 5) * delta_square * h * x
    )
    return w_square, nonlinear_elimination


def main() -> None:
    args = parse_args()
    primes = [
        prime
        for prime in primes_at_most(args.prime_limit)
        if prime >= 11 and prime != 769
    ]
    maximum_prime = max(primes)
    values = apery_numbers(2 * maximum_prime)
    direct_jets, companion_jets = shifted_fundamental_solutions(
        maximum_prime - 1
    )
    j_jets = [
        j_polynomial(companion_jets, index)
        for index in range(maximum_prime)
    ]
    checks: Counter[str] = Counter()

    for prime in primes:
        modulus6 = prime**6
        modulus7 = prime**7
        modulus8 = prime**8
        delta = values[prime - 1] - 1
        endpoint_h = values[prime] - 5 + 7 * delta
        assert delta % prime**3 == 0
        assert endpoint_h % prime**5 == 0

        endpoint_w = (
            Fraction(values[2 * prime] - 73 + 824 * delta)
            - Fraction(752, 5) * endpoint_h
        ) / modulus6
        assert gcd(endpoint_w.denominator, prime) == 1
        w_residue = fraction_mod(endpoint_w, prime)
        normalized_delta = delta // prime**3 % prime
        predicted_w = fraction_mod(
            Fraction(24 * 769, 5) * normalized_delta**2,
            prime,
        )
        assert w_residue == predicted_w
        checks["w_from_delta_square"] += 1

        reflected_w = (
            Fraction(
                values[2 * prime - 1] - 5 - 8 * delta
            )
            - Fraction(336, 5) * endpoint_h
        ) / modulus6
        assert fraction_mod(reflected_w, prime) == (
            -103 * pow(769, -1, prime) * w_residue
        ) % prime
        checks["endpoint_w"] += 1

        for index in range(prime):
            if values[index] % prime != 0:
                continue
            reflected_index = prime - 1 - index
            upper_index = prime + index
            assert values[reflected_index] % prime == 0
            assert values[upper_index] % prime == 0
            x_value = values[upper_index] // prime % modulus7

            direct_u = evaluate_mod(
                direct_jets[index], prime, modulus8
            )
            direct_j = evaluate_mod(
                j_jets[index][:5], prime, modulus8
            )
            direct_numerator = (
                (5 - 7 * delta) * direct_u
                - prime**3 * (1 + delta) * direct_j
            )
            direct_residue = divide_residue(
                direct_numerator, prime, modulus8
            )
            direct_factor = fraction_mod(
                Fraction(1) - Fraction(endpoint_h, 5),
                modulus7,
            )
            assert direct_residue == (
                x_value * direct_factor
            ) % modulus7
            checks["direct"] += 1

            reflected_u = evaluate_mod(
                direct_jets[reflected_index],
                -2 * prime,
                modulus8,
            )
            reflected_j = evaluate_mod(
                j_jets[reflected_index][:5],
                -2 * prime,
                modulus8,
            )
            reflected_numerator = (
                (5 + 8 * delta) * reflected_u
                + 8
                * prime**3
                * (73 - 824 * delta)
                * reflected_j
            )
            reflected_residue = divide_residue(
                reflected_numerator, prime, modulus8
            )
            reflected_factor = fraction_mod(
                Fraction(1)
                - Fraction(336, 25) * endpoint_h
                + Fraction(103, 5 * 769)
                * modulus6
                * w_residue,
                modulus7,
            )
            assert reflected_residue == (
                x_value * reflected_factor
            ) % modulus7
            checks["reflected"] += 1

            combined = (
                336 * direct_residue - 5 * reflected_residue
            ) % modulus7
            combined_factor = fraction_mod(
                Fraction(331)
                - Fraction(103, 769)
                * modulus6
                * w_residue,
                modulus7,
            )
            assert combined == (
                x_value * combined_factor
            ) % modulus7
            checks["combined"] += 1

            # H6 makes W a quadratic function of the already present
            # endpoint coordinate Delta.  Since Delta^2*H=0 mod p^7,
            # multiplying the direct row by Delta^2 cancels it:
            #
            # (1680+2472 Delta^2)D_7-25Z_7=1655x mod p^7.
            fixed_combined = (
                (1680 + 2472 * delta * delta) * direct_residue
                - 25 * reflected_residue
            ) % modulus7
            assert fixed_combined == 1655 * x_value % modulus7
            checks["fixed_combined"] += 1

            small_delta_square = (
                modulus6 * normalized_delta**2
            ) % modulus7
            assert small_delta_square == delta * delta % modulus7
            fixed_small = (
                (1680 + 2472 * small_delta_square)
                * direct_residue
                - 25 * reflected_residue
            ) % modulus7
            assert fixed_small == fixed_combined
            checks["small_coefficient"] += 1

    w_symbolic, elimination_symbolic = symbolic_correction_checks()
    assert w_symbolic
    assert elimination_symbolic
    print(f"prime_limit={args.prime_limit}")
    print(f"symbolic_w_square={int(w_symbolic)}/1")
    print(
        "symbolic_nonlinear_elimination="
        f"{int(elimination_symbolic)}/1"
    )
    for name in sorted(checks):
        print(f"{name}_checks={checks[name]}")
    print("normalization_exception=769")
    print("unit_exceptions=331,769")
    print("failures=0")


if __name__ == "__main__":
    main()
