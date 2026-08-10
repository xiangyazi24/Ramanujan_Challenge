#!/usr/bin/env python3
"""Exact checks for gap_kernel_cartier_defect.tex."""

from fractions import Fraction


def primes_through(limit: int):
    return [
        value
        for value in range(2, limit + 1)
        if all(value % divisor for divisor in range(2, int(value**0.5) + 1))
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


def reduce_fraction(value: Fraction, prime: int) -> int:
    assert value.denominator % prime
    return value.numerator * pow(value.denominator, prime - 2, prime) % prime


def normalized_triangle(prime: int, start: int):
    maximum = prime - start - 1
    values = [1]
    previous = 0
    current = 1
    for offset in range(maximum):
        argument = start + offset
        following = (
            p_value(argument) * current - argument**3 * previous
        ) * pow((argument + 1) ** 3, prime - 2, prime) % prime
        values.append(following)
        previous, current = current, following
    return values


def operator_defect(prime: int, start: int, values):
    result = []
    for degree in range(len(values) + 2):
        current = values[degree] if degree < len(values) else 0
        previous = values[degree - 1] if 0 <= degree - 1 < len(values) else 0
        before_previous = (
            values[degree - 2] if 0 <= degree - 2 < len(values) else 0
        )
        result.append(
            (
                (start + degree) ** 3 * current
                - p_value(start + degree - 1) * previous
                + (start + degree - 1) ** 3 * before_previous
            )
            % prime
        )
    return result


def evaluate_polynomial(values, argument: int, prime: int) -> int:
    result = 0
    for coefficient in reversed(values):
        result = (result * argument + coefficient) % prime
    return result


def cartier(polynomial, prime: int):
    maximum = max(polynomial, default=-1)
    return {
        exponent // prime: coefficient % prime
        for exponent, coefficient in polynomial.items()
        if exponent % prime == 0 and coefficient % prime
    } if maximum >= 0 else {}


def main() -> None:
    apery, companion = apery_rows(80)
    defect_checks = 0
    boundary_checks = 0
    marked_fiber_checks = 0
    extraction_checks = 0
    degree_checks = 0
    bivariate_checks = 0

    for prime in primes_through(43):
        if prime < 5:
            continue
        apery_mod = [reduce_fraction(value, prime) for value in apery]
        companion_mod = [
            reduce_fraction(value, prime)
            for value in companion[:prime]
        ]
        assert apery_mod[prime - 1] == 1
        assert apery_mod[prime - 2] == 5 % prime
        assert companion_mod[prime - 1] == 0
        assert companion_mod[prime - 2] == 1

        aggregate_defect = {}
        expected_aggregate = {}
        for start in range(1, prime):
            values = normalized_triangle(prime, start)
            maximum = prime - start - 1
            defect = operator_defect(prime, start, values)
            expected = [0] * (maximum + 3)
            expected[0] = start**3 % prime
            expected[prime - start] = (
                -start**3 * apery_mod[start - 1]
            ) % prime
            assert defect == expected
            defect_checks += len(defect)

            next_value = start**3 * (
                apery[start - 1] * companion[prime]
                - companion[start - 1] * apery[prime]
            )
            scaled_next = prime**3 * next_value
            assert reduce_fraction(scaled_next, prime) == (
                start**3 * apery_mod[start - 1]
            ) % prime
            boundary_checks += 1

            if companion_mod[start - 1]:
                state = (
                    apery_mod[start - 1]
                    * pow(companion_mod[start - 1], prime - 2, prime)
                    % prime
                )
                marked = [
                    (apery_mod[index] - state * companion_mod[index]) % prime
                    for index in range(start - 1, prime)
                ]
                assert marked[0] == 0 and marked[1]
                inverse_next = pow(marked[1], prime - 2, prime)
                assert values == [
                    marked[offset + 1] * inverse_next % prime
                    for offset in range(len(values))
                ]
                marked_fiber_checks += len(values)

            for offset, coefficient in enumerate(values):
                mellin = sum(
                    evaluate_polynomial(values, point, prime)
                    * pow(point, prime - 1 - offset, prime)
                    for point in range(1, prime)
                ) % prime
                assert (-mellin) % prime == coefficient

                shifted = {
                    prime - offset + index: value
                    for index, value in enumerate(values)
                }
                expected_cartier = {1: coefficient} if coefficient else {}
                assert cartier(shifted, prime) == expected_cartier
                extraction_checks += 2

            if start <= prime - 2:
                actual_degree = max(
                    index for index, coefficient in enumerate(values)
                    if coefficient
                )
                assert actual_degree >= prime - start - 2
                assert values[-1] == -start**3 * companion_mod[start - 1] % prime
                assert values[-2] == start**3 * (
                    apery_mod[start - 1] - 5 * companion_mod[start - 1]
                ) % prime
                degree_checks += 3

            for degree, coefficient in enumerate(defect):
                if coefficient:
                    aggregate_defect[(start, degree)] = coefficient
            expected_aggregate[(start, 0)] = start**3 % prime
            boundary_coefficient = (
                -start**3 * apery_mod[start - 1]
            ) % prime
            if boundary_coefficient:
                expected_aggregate[(start, prime - start)] = boundary_coefficient

        assert aggregate_defect == expected_aggregate
        bivariate_checks += len(aggregate_defect)

    print("gap_kernel_cartier_defect_verify: PASS")
    print(f"defect_checks={defect_checks}")
    print(f"boundary_checks={boundary_checks}")
    print(f"marked_fiber_checks={marked_fiber_checks}")
    print(f"extraction_checks={extraction_checks}")
    print(f"degree_checks={degree_checks}")
    print(f"bivariate_checks={bivariate_checks}")


if __name__ == "__main__":
    main()
