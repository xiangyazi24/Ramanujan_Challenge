#!/usr/bin/env python3
"""Audit the weight-five endpoint law and sextic target elimination.

The script has two independent parts.

First it checks, for all requested endpoint quotients m,

    b_(mp)-b_m
      == E_m Delta_p + p^5 B_(p-5) L_m       (mod p^6),

    b_(mp-1)-b_(m-1)
      == F_m Delta_p + p^5 B_(p-5) M_m       (mod p^6),

where L_m and M_m are explicit finite binomial moments.

Second it constructs degree-six shifted recurrence jets.  At every
target p | b_r, n=p+r, s=p-1-r, it checks

    331 b_n/p == 336 D6_(p,r) - 5 Z6_(p,s)    (mod p^6)

and recovers the common cofactor modulo p^6 for p != 331.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from fractions import Fraction
from math import comb, gcd, isqrt


def primes_at_most(limit: int) -> list[int]:
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[:2] = b"\x00\x00"
    for prime in range(2, isqrt(limit) + 1):
        if sieve[prime]:
            sieve[prime * prime : limit + 1 : prime] = b"\x00" * (
                (limit - prime * prime) // prime + 1
            )
    return [candidate for candidate in range(2, limit + 1) if sieve[candidate]]


ORDER = 7


def recurrence_polynomial(index: int) -> int:
    return 34 * index**3 + 51 * index**2 + 27 * index + 5


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


def add(left: list[Fraction], right: list[Fraction]) -> list[Fraction]:
    return [left[index] + right[index] for index in range(ORDER)]


def scale(polynomial: list[Fraction], scalar: int) -> list[Fraction]:
    return [scalar * coefficient for coefficient in polynomial]


def multiply(
    left: list[Fraction], right: list[Fraction]
) -> list[Fraction]:
    result = [Fraction(0) for _ in range(ORDER)]
    for left_index, left_coefficient in enumerate(left):
        for right_index, right_coefficient in enumerate(right):
            if left_index + right_index >= ORDER:
                break
            result[left_index + right_index] += (
                left_coefficient * right_coefficient
            )
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
                scale(multiply(shifted_cube, sequence[index - 1]), -1),
            )
            sequence.append(multiply(numerator, inverse_denominator))

    return direct, companion


def fraction_mod(value: Fraction, modulus: int) -> int:
    assert gcd(value.denominator, modulus) == 1
    return value.numerator * pow(value.denominator, -1, modulus) % modulus


def evaluate_mod(
    polynomial: list[Fraction], argument: int, modulus: int
) -> int:
    return sum(
        fraction_mod(coefficient, modulus)
        * pow(argument, degree, modulus)
        for degree, coefficient in enumerate(polynomial)
    ) % modulus


def j_coefficient(
    companion: list[list[Fraction]], index: int, degree: int
) -> Fraction:
    """Coefficient of V_index(X)/(1+X)^3."""
    return sum(
        companion[index][source_degree]
        * (-1) ** (degree - source_degree)
        * comb(degree - source_degree + 2, 2)
        for source_degree in range(degree + 1)
    )


def apery_summand(upper: int, lower: int) -> int:
    return comb(upper, lower) ** 2 * comb(upper + lower, lower) ** 2


def endpoint_coefficients(
    quotient: int,
) -> tuple[Fraction, Fraction]:
    direct = Fraction(0)
    reflected = Fraction(0)
    m = quotient

    for j in range(m + 1):
        a = m - j
        anchor = Fraction(
            -4 * j * m**2 * (2 * j**2 + m**2 + 8),
            5,
        )
        direct_nonanchor = (
            Fraction(4)
            - Fraction(12 * a, 5)
            + Fraction(12 * a**2, 5)
            - Fraction(4 * m**2, 5)
            + 8 * m**3
        )
        direct += apery_summand(m, j) * (
            anchor + a**2 * direct_nonanchor
        )

        if j < m:
            c = m + j
            reflected_nonanchor = (
                Fraction(4)
                + Fraction(12 * c, 5)
                + Fraction(12 * c**2, 5)
                - Fraction(4 * m**2, 5)
                - 8 * m**3
            )
            reflected += apery_summand(m - 1, j) * (
                anchor + c**2 * reflected_nonanchor
            )

    return direct, reflected


def bernoulli_p_minus_5_mod(prime: int) -> int:
    """Compute B_(p-5) mod p from S_4 == (4/5)p B_(p-5)."""
    modulus = prime**2
    fourth_harmonic = sum(
        pow(index, -4, modulus) for index in range(1, prime)
    ) % modulus
    assert fourth_harmonic % prime == 0
    return (
        5
        * pow(4, -1, prime)
        * (fourth_harmonic // prime)
    ) % prime


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prime-limit", type=int, default=1000)
    parser.add_argument("--quotient-limit", type=int, default=12)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    primes = [
        prime
        for prime in primes_at_most(args.prime_limit)
        if prime >= 7
    ]
    maximum_prime = max(primes)
    maximum_index = max(
        args.quotient_limit * maximum_prime,
        2 * maximum_prime,
    )
    values = apery_numbers(maximum_index)
    direct_jets, companion_jets = shifted_fundamental_solutions(
        maximum_prime - 1
    )
    checks: Counter[str] = Counter()

    endpoint_moments = {
        quotient: endpoint_coefficients(quotient)
        for quotient in range(1, args.quotient_limit + 1)
    }
    assert endpoint_moments[1] == (Fraction(-24), Fraction(0))
    assert endpoint_moments[2] == (
        Fraction(-18048, 5),
        Fraction(-8064, 5),
    )

    target_rows: dict[int, list[int]] = defaultdict(list)
    local_combined: dict[tuple[int, int], int] = {}

    for prime in primes:
        modulus3 = prime**3
        modulus5 = prime**5
        modulus6 = prime**6
        modulus7 = prime**7
        delta = values[prime - 1] - 1
        assert delta % modulus3 == 0
        weight_five = bernoulli_p_minus_5_mod(prime)

        direct_new = values[prime] - 5 + 7 * delta
        reflected_new = values[2 * prime - 1] - 5 - 8 * delta
        assert direct_new % modulus5 == 0
        assert reflected_new % modulus5 == 0
        h_scalar = direct_new // modulus5 % prime
        reflected_scalar = reflected_new // modulus5 % prime
        assert h_scalar == -24 * weight_five % prime
        assert (
            5 * reflected_scalar - 336 * h_scalar
        ) % prime == 0
        checks["endpoint_scalar"] += 1

        for quotient in range(1, args.quotient_limit + 1):
            direct_coefficient = (
                quotient**3
                * (values[quotient - 1] - 17 * values[quotient])
                // 12
            )
            reflected_coefficient = (
                quotient**3
                * (17 * values[quotient - 1] - values[quotient])
                // 12
            )
            direct_moment, reflected_moment = endpoint_moments[quotient]

            direct_difference = (
                values[quotient * prime] - values[quotient]
            )
            reflected_difference = (
                values[quotient * prime - 1]
                - values[quotient - 1]
            )
            assert (
                direct_difference
                - direct_coefficient * delta
                - modulus5
                * weight_five
                * fraction_mod(direct_moment, modulus6)
            ) % modulus6 == 0
            assert (
                reflected_difference
                - reflected_coefficient * delta
                - modulus5
                * weight_five
                * fraction_mod(reflected_moment, modulus6)
            ) % modulus6 == 0
            checks["direct_endpoint"] += 1
            checks["reflected_endpoint"] += 1

        for index in range(prime):
            if values[index] % prime != 0:
                continue

            reflected_index = prime - 1 - index
            upper_index = prime + index
            assert values[reflected_index] % prime == 0
            assert values[upper_index] % prime == 0
            quotient = values[upper_index] // prime % modulus6
            target_rows[upper_index].append(prime)

            direct_j = sum(
                j_coefficient(companion_jets, index, degree)
                * prime**degree
                for degree in range(4)
            )
            direct_stripped = (
                (5 - 7 * delta)
                * evaluate_mod(
                    direct_jets[index], prime, modulus7
                )
                - prime**3 * (1 + delta)
                * fraction_mod(direct_j, modulus7)
            ) % modulus7
            assert direct_stripped % prime == 0
            normalized_direct = direct_stripped // prime % modulus6
            assert normalized_direct == (
                quotient
                * (
                    1
                    - prime**5
                    * h_scalar
                    * pow(5, -1, modulus6)
                )
            ) % modulus6
            checks["direct_target"] += 1

            reflected_j = sum(
                j_coefficient(
                    companion_jets, reflected_index, degree
                )
                * (-2 * prime) ** degree
                for degree in range(4)
            )
            reflected_stripped = (
                (5 + 8 * delta)
                * evaluate_mod(
                    direct_jets[reflected_index],
                    -2 * prime,
                    modulus7,
                )
                + 8 * prime**3 * (73 - 824 * delta)
                * fraction_mod(reflected_j, modulus7)
            ) % modulus7
            assert reflected_stripped % prime == 0
            normalized_reflected = (
                reflected_stripped // prime % modulus6
            )
            assert normalized_reflected == (
                quotient
                * (
                    1
                    - prime**5
                    * reflected_scalar
                    * pow(5, -1, modulus6)
                )
            ) % modulus6
            checks["reflected_target"] += 1

            combined = (
                336 * normalized_direct
                - 5 * normalized_reflected
            ) % modulus6
            assert combined == 331 * quotient % modulus6
            checks["sextic_elimination"] += 1
            if prime != 331:
                local_combined[(upper_index, prime)] = combined

    for upper_index, targets in target_rows.items():
        radical = 1
        for prime in targets:
            radical *= prime
        assert values[upper_index] % radical == 0
        common_quotient = values[upper_index] // radical

        for prime in targets:
            if prime == 331:
                continue
            modulus6 = prime**6
            combined = local_combined[(upper_index, prime)]
            radical_without_prime = radical // prime
            recovered = (
                combined
                * pow(
                    331 * radical_without_prime,
                    -1,
                    modulus6,
                )
            ) % modulus6
            assert recovered == common_quotient % modulus6
            checks["common_quotient_mod_p6"] += 1

    print(f"prime_limit={args.prime_limit}")
    print(f"quotient_limit={args.quotient_limit}")
    for key in sorted(checks):
        print(f"{key}_checks={checks[key]}")
    print("exceptional_determinant_prime=331")
    print("first_failure=None")
    print("failures=0")


if __name__ == "__main__":
    main()
