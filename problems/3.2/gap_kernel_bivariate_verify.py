#!/usr/bin/env python3
"""Exact checks for gap_kernel_bivariate.tex."""

from fractions import Fraction
from math import prod


def P(n: int) -> int:
    return (2 * n + 1) * (17 * n * n + 17 * n + 5)


def companion_sequences(limit: int):
    b = [Fraction(1), Fraction(5)]
    c = [Fraction(0), Fraction(1)]
    for n in range(1, limit):
        denominator = (n + 1) ** 3
        b.append((P(n) * b[n] - n**3 * b[n - 1]) / denominator)
        c.append((P(n) * c[n] - n**3 * c[n - 1]) / denominator)
    return b, c


def normalized_restarts(a_max: int, m_max: int):
    values = {}
    for a in range(1, a_max + 1):
        previous, current = 0, 1
        for m in range(m_max + 1):
            denominator = prod((a + j) ** 3 for j in range(1, m + 1))
            values[a, m] = Fraction(current, denominator)
            previous, current = (
                current,
                P(a + m) * current - (a + m) ** 6 * previous,
            )
    return values


def valuation(number: int, prime: int) -> int:
    result = 0
    while number % prime == 0:
        number //= prime
        result += 1
    return result


def rational_valuation(value: Fraction, prime: int) -> int:
    return valuation(abs(value.numerator), prime) - valuation(value.denominator, prime)


def primes_through(limit: int):
    result = []
    for candidate in range(2, limit + 1):
        if all(candidate % q for q in range(2, int(candidate**0.5) + 1)):
            result.append(candidate)
    return result


def main() -> None:
    a_max = 10
    m_max = 10
    sequence_limit = 80
    b, c = companion_sequences(sequence_limit)
    values = normalized_restarts(a_max, m_max)

    casoratian_checks = 0
    pde_checks = 0

    # Build K(s/z,z) coefficientwise and then apply the nonnegative-z
    # projection.  In the CT formula the same retained term has u-exponent
    # -m and is paired with u^m from 1/(1-u).
    positive_part = {}
    constant_term = {}
    rejected_negative_terms = 0
    for r in range(a_max):
        for n in range(a_max + m_max + 1):
            a = r + 1
            m = n - r - 1
            coefficient = a**3 * (b[r] * c[n] - c[r] * b[n])
            if m < 0:
                rejected_negative_terms += 1
                continue
            positive_part[a, m] = positive_part.get((a, m), Fraction(0)) + coefficient
            constant_term[a, m] = constant_term.get((a, m), Fraction(0)) + coefficient

    positive_part_checks = 0

    for a in range(1, a_max + 1):
        for m in range(m_max + 1):
            expected = a**3 * (b[a - 1] * c[a + m] - c[a - 1] * b[a + m])
            assert values[a, m] == expected
            casoratian_checks += 1

            assert positive_part[a, m] == values[a, m]
            assert constant_term[a, m] == values[a, m]
            positive_part_checks += 1

            lhs = (a + m) ** 3 * values[a, m]
            if m >= 1:
                lhs -= P(a + m - 1) * values[a, m - 1]
            if m >= 2:
                lhs += (a + m - 1) ** 3 * values[a, m - 2]
            rhs = Fraction(a**3 if m == 0 else 0)
            assert lhs == rhs
            pde_checks += 1

    operator_checks = 0
    for n in range(sequence_limit):
        b_coeff = n**3 * b[n]
        c_coeff = n**3 * c[n]
        if n >= 1:
            b_coeff -= P(n - 1) * b[n - 1]
            c_coeff -= P(n - 1) * c[n - 1]
        if n >= 2:
            b_coeff += (n - 1) ** 3 * b[n - 2]
            c_coeff += (n - 1) ** 3 * c[n - 2]
        assert b_coeff == 0
        assert c_coeff == (1 if n == 1 else 0)
        operator_checks += 2

    denominator_checks = 0
    for prime in primes_through(73):
        if prime < 5:
            continue
        value = c[prime]
        assert rational_valuation(value, prime) == -3
        unit = value * prime**3
        assert unit.denominator % prime != 0
        residue = (unit.numerator * pow(unit.denominator, -1, prime)) % prime
        assert residue == 1
        denominator_checks += 1

    assert all(values[a, 0] == 1 for a in range(1, a_max + 1))
    assert all(values[1, m] == c[m + 1] for m in range(m_max + 1))

    print("gap_kernel_bivariate_verify: PASS")
    print(f"casoratian_checks={casoratian_checks}")
    print(f"positive_part_checks={positive_part_checks}")
    print(f"rejected_negative_terms={rejected_negative_terms}")
    print(f"pde_checks={pde_checks}")
    print(f"operator_checks={operator_checks}")
    print(f"prime_denominator_checks={denominator_checks}")


if __name__ == "__main__":
    main()
