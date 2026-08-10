#!/usr/bin/env python3
"""Exact checks for gap_kernel_borel_gauges.tex."""

from fractions import Fraction
from math import gcd


def primes_through(limit: int):
    return [
        n
        for n in range(2, limit + 1)
        if all(n % d for d in range(2, int(n**0.5) + 1))
    ]


def p_value(index: int) -> int:
    return (2 * index + 1) * (17 * index * index + 17 * index + 5)


def apery_rows(length: int):
    apery = [Fraction(1), Fraction(5)]
    companion = [Fraction(0), Fraction(1)]
    for index in range(1, length - 1):
        for row in (apery, companion):
            row.append(
                Fraction(
                    p_value(index) * row[index]
                    - index**3 * row[index - 1],
                    (index + 1) ** 3,
                )
            )
    return apery, companion


def restart_values(start: int, length: int):
    previous = 0
    current = 1
    values = [current]
    for offset in range(length - 1):
        following = (
            p_value(start + offset) * current
            - (start + offset) ** 6 * previous
        )
        previous, current = current, following
        values.append(current)
    return values


def rising(start: int, length: int) -> int:
    result = 1
    for value in range(start + 1, start + length + 1):
        result *= value
    return result


def factorial(index: int) -> int:
    result = 1
    for value in range(2, index + 1):
        result *= value
    return result


def lcm_through(index: int) -> int:
    result = 1
    for value in range(2, index + 1):
        result = result * value // gcd(result, value)
    return result


def valuation_integer(value: int, prime: int) -> int:
    if value == 0:
        raise ValueError("valuation of zero")
    result = 0
    while value % prime == 0:
        value //= prime
        result += 1
    return result


def valuation_fraction(value: Fraction, prime: int) -> int:
    return (
        valuation_integer(abs(value.numerator), prime)
        - valuation_integer(value.denominator, prime)
    )


def reduce_fraction(value: Fraction, prime: int) -> int:
    assert value.denominator % prime
    return value.numerator * pow(value.denominator, prime - 2, prime) % prime


def main() -> None:
    apery, companion = apery_rows(90)
    recurrence_checks = 0
    gauge_checks = 0
    denominator_checks = 0
    integrality_checks = 0
    triangle_checks = 0

    for start in range(1, 11):
        values = restart_values(start, 12)
        for offset, current in enumerate(values):
            if offset:
                assert current == (
                    p_value(start + offset - 1) * values[offset - 1]
                    - (start + offset - 1) ** 6
                    * (0 if offset == 1 else values[offset - 2])
                )
            recurrence_checks += 1

            normalized = Fraction(current, rising(start, offset) ** 3)
            kernel = start**3 * (
                apery[start - 1] * companion[start + offset]
                - companion[start - 1] * apery[start + offset]
            )
            assert normalized == kernel
            gauges = []
            for allocation in range(4):
                denominator = (
                    factorial(offset) ** (3 - allocation)
                    * rising(start, offset) ** allocation
                )
                gauge = Fraction(current, denominator)
                gauges.append(gauge)
                binomial = Fraction(
                    rising(start, offset), factorial(offset)
                )
                assert gauge == binomial ** (3 - allocation) * normalized
                gauge_checks += 1
            binomial = Fraction(rising(start, offset), factorial(offset))
            for allocation in range(3):
                assert gauges[allocation + 1] * binomial == gauges[allocation]

            index = start + offset
            lcm_value = lcm_through(index)
            arithmetic_gauge = 6 * lcm_value**3 * normalized
            assert arithmetic_gauge.denominator == 1
            integrality_checks += 1

    for index in range(90):
        lcm_value = lcm_through(index)
        assert (6 * lcm_value**3 * companion[index]).denominator == 1
        assert apery[index].denominator == 1
        integrality_checks += 1

    for prime in primes_through(73):
        if prime < 7:
            continue
        for shift in range(min(9, prime)):
            scaled = prime**3 * companion[prime + shift]
            assert scaled.denominator % prime
            assert reduce_fraction(scaled, prime) == int(apery[shift]) % prime
            denominator_checks += 1
        assert valuation_fraction(companion[prime + 1], prime) == -3
        for allocation in range(4):
            gauge = (prime + 1) ** (3 - allocation) * companion[prime + 1]
            assert valuation_fraction(gauge, prime) == -3
            denominator_checks += 1

    for prime in primes_through(43):
        if prime < 7:
            continue
        for start in range(1, prime):
            values = restart_values(start, prime - start)
            for offset, current in enumerate(values):
                index = start + offset
                if index >= prime:
                    continue
                normalized = Fraction(current, rising(start, offset) ** 3)
                arithmetic_gauge = (
                    6 * lcm_through(index) ** 3 * normalized
                )
                assert arithmetic_gauge.denominator == 1
                assert (current % prime == 0) == (
                    arithmetic_gauge.numerator % prime == 0
                )
                triangle_checks += 1

    print("gap_kernel_borel_gauges_verify: PASS")
    print(f"recurrence_checks={recurrence_checks}")
    print(f"gauge_checks={gauge_checks}")
    print(f"denominator_checks={denominator_checks}")
    print(f"integrality_checks={integrality_checks}")
    print(f"triangle_checks={triangle_checks}")


if __name__ == "__main__":
    main()
