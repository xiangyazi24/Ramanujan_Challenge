#!/usr/bin/env python3
"""Machine checks for the honest horizontal objects in the rank-two attack.

This script verifies six logically distinct facts.

1. Additive orthogonality gives the exact two-prime zero correlation.
2. Averaging that correlation over a full CRT period factors exactly.
3. A full *character-period* correlation of linear Mellin transforms has an
   exact gcd(p-1,q-1)-term formula (a genuine linearized power saving).
4. Interior Apéry coefficients themselves are Mellin coefficients, and the
   quadratic pullback gives the precise split-cover (quadratically twisted)
   projection.
5. The complex Euler/Kummer lift on that cover obeys the pointwise Weil bound,
   and its full two-prime character-period correlation collapses to gcd terms.
6. Exponentiating a Jacobi/character sum does not move the additive character
   inside the sum.  Thus the zero detector is not a Deligne sum merely because
   the coefficient being tested is one.

The last distinction is the obstruction that survives the explicit rank-two
formula.  All computations use only Python's standard library.
"""

from __future__ import annotations

import cmath
from math import gcd, lcm, pi, sqrt

from codex_hgk_coefficients import (
    apery_coefficients,
    branch_for_prime,
    branch_values_from_pullback,
    hypergeometric_hasse_coefficients,
    hypergeometric_parameters,
    legendre,
    polynomial_value,
)


def zero_set(prime: int) -> set[int]:
    return {
        index for index, value in enumerate(apery_coefficients(prime)) if value == 0
    }


def pair_count(prime: int, other: int, start: int, length: int) -> int:
    first_zeros = zero_set(prime)
    second_zeros = zero_set(other)
    return sum(
        index % prime in first_zeros and index % other in second_zeros
        for index in range(start, start + length)
    )


def additive_pair_count(prime: int, other: int, start: int, length: int) -> complex:
    first_values = apery_coefficients(prime)
    second_values = apery_coefficients(other)
    first_root = cmath.exp(2j * pi / prime)
    second_root = cmath.exp(2j * pi / other)
    result = 0j
    for first_mode in range(prime):
        for second_mode in range(other):
            result += sum(
                first_root ** (first_mode * first_values[index % prime])
                * second_root ** (second_mode * second_values[index % other])
                for index in range(start, start + length)
            )
    return result / (prime * other)


def verify_crt_averages() -> None:
    prime, other = 17, 19
    first_zeros = zero_set(prime)
    second_zeros = zero_set(other)
    assert first_zeros == {3, 13}
    assert second_zeros == {8, 10}

    period = prime * other
    complete = pair_count(prime, other, 0, period)
    assert complete == len(first_zeros) * len(second_zeros)

    length = 31
    sliding_total = sum(pair_count(prime, other, start, length) for start in range(period))
    assert sliding_total == length * len(first_zeros) * len(second_zeros)

    start = 23
    direct = pair_count(prime, other, start, length)
    orthogonal = additive_pair_count(prime, other, start, length)
    assert abs(orthogonal.imag) < 1e-8
    assert abs(orthogonal.real - direct) < 1e-8

    # Four-prime version: over the full CRT period every compatible zero
    # quadruple occurs exactly once.
    primes = (5, 11, 17, 19)
    four_period = 1
    for value in primes:
        four_period *= value
    expected = 1
    zero_sets = []
    for value in primes:
        zeros = zero_set(value)
        zero_sets.append(zeros)
        expected *= len(zeros)
    actual = sum(
        all(index % value in zeros for value, zeros in zip(primes, zero_sets))
        for index in range(four_period)
    )
    assert actual == expected

    print(
        "CRT averages: pair and four-prime complete periods factor exactly; "
        "sliding-interval mean and additive zero detector VERIFIED"
    )


def phi_mod(prime: int, value: int) -> int:
    """The rational pullback phi(x)=x(1-8x)/(1+x) over F_p."""

    assert (1 + value) % prime
    return value * (1 - 8 * value) * pow(1 + value, -1, prime) % prime


def cover_fibres(prime: int) -> dict[int, list[int]]:
    fibres = {value: [] for value in range(prime)}
    for x in range(prime):
        if (1 + x) % prime:
            fibres[phi_mod(prime, x)].append(x)
    return fibres


def verify_apery_mellin_and_split_cover() -> None:
    """Verify the shorter Mellin formula and its quadratic-cover projection.

    For 1 <= r <= p-2, ordinary and cyclic coefficient extraction agree, so

        b_r = -sum_{t != 0} A_p(t)t^{-r}.

    Summing after the degree-two pullback weights each t by
    1+eta(t^2-34t+1).  This is the exact split-cover projection that admits a
    direct complex Kummer lift; it is not silently identified with b_r.
    """

    for prime in (13, 29):
        apery = apery_coefficients(prime)
        hasse = hypergeometric_hasse_coefficients(prime)
        apery_values = [polynomial_value(apery, t, prime) for t in range(prime)]
        fibres = cover_fibres(prime)

        for t in range(prime):
            discriminant = (t * t - 34 * t + 1) % prime
            assert len(fibres[t]) == 1 + legendre(discriminant, prime)

        # The square-root branch disappears after squaring: on either branch
        # A_p(phi(x))=H_p(x)^2/(1+x)^(p-1)=H_p(x)^2 on F_p-points.
        for x in range(prime):
            if (1 + x) % prime:
                t = phi_mod(prime, x)
                assert apery_values[t] == polynomial_value(hasse, x, prime) ** 2 % prime

        split_values = [len(fibres[t]) * apery_values[t] % prime for t in range(prime)]
        for index in range(1, prime - 1):
            direct = -sum(
                apery_values[t] * pow(t, -index, prime)
                for t in range(1, prime)
            ) % prime
            assert direct == apery[index]

            split = -sum(
                split_values[t] * pow(t, -index, prime)
                for t in range(1, prime)
            ) % prime
            quadratic_twist = -sum(
                legendre(t * t - 34 * t + 1, prime)
                * apery_values[t]
                * pow(t, -index, prime)
                for t in range(1, prime)
            ) % prime
            assert split == (apery[index] + quadratic_twist) % prime

    print(
        "Apéry Mellin extraction and split quadratic-cover projection: "
        "p=13,29, every interior index VERIFIED"
    )


def primitive_root(prime: int) -> int:
    order = prime - 1
    factors = []
    remaining = order
    divisor = 2
    while divisor * divisor <= remaining:
        if remaining % divisor == 0:
            factors.append(divisor)
            while remaining % divisor == 0:
                remaining //= divisor
        divisor += 1
    if remaining > 1:
        factors.append(remaining)
    return next(
        candidate
        for candidate in range(2, prime)
        if all(pow(candidate, order // factor, prime) != 1 for factor in factors)
    )


def multiplicative_ordering(values: list[complex | int], prime: int) -> list[complex]:
    generator = primitive_root(prime)
    return [complex(values[pow(generator, exponent, prime)]) for exponent in range(prime - 1)]


def mellin_transform(values: list[complex], index: int) -> complex:
    order = len(values)
    root = cmath.exp(-2j * pi / order)
    return -sum(value * root ** (index * exponent) for exponent, value in enumerate(values))


def discrete_log_table(prime: int) -> list[int]:
    generator = primitive_root(prime)
    result = [-1] * prime
    value = 1
    for exponent in range(prime - 1):
        result[value] = exponent
        value = value * generator % prime
    return result


def complex_character(
    prime: int, exponent: int, value: int, logarithms: list[int]
) -> complex:
    value %= prime
    if value == 0:
        return 0j
    return cmath.exp(2j * pi * exponent * logarithms[value] / (prime - 1))


def euler_trace_on_cover(prime: int, x: int, logarithms: list[int]) -> complex | None:
    """Complex Kummer lift of H_p(x), away from the point z=infinity."""

    denominator = (1 - 2 * x) % prime
    if denominator == 0:
        return None
    z = 27 * x * x * pow(denominator, -3, prime) % prime
    first, second, _ = hypergeometric_parameters(prime)
    sign = -1 if (second + 1) & 1 else 1
    return sign * sum(
        complex_character(prime, first, y, logarithms)
        * complex_character(prime, second, 1 - y, logarithms)
        * complex_character(prime, first, 1 - z * y, logarithms)
        for y in range(prime)
    )


def euler_square_fibre_values(prime: int) -> tuple[list[complex], float]:
    """Push the squared Euler trace through the split part of the cover.

    The single x above z=infinity is punctured.  Its contribution is an
    explicit parity mode and is deliberately kept out of the uniform Weil
    estimate rather than hidden in it.
    """

    logarithms = discrete_log_table(prime)
    fibres = cover_fibres(prime)
    values = [0j] * prime
    maximum = 0.0
    for t, fibre in fibres.items():
        for x in fibre:
            trace = euler_trace_on_cover(prime, x, logarithms)
            if trace is None:
                continue
            maximum = max(maximum, abs(trace))
            values[t] += trace * trace
    return values, maximum


def verify_split_cover_trace_correlation() -> None:
    """Check the exact two-prime DFT collapse for the Kummer cover lift."""

    prime, other = 13, 29
    first_values, first_maximum = euler_square_fibre_values(prime)
    second_values, second_maximum = euler_square_fibre_values(other)

    # A rank-one Kummer sheaf on P^1 minus at most four points has H_c^1
    # dimension at most two here.  The exceptional z=0,1 fibres are smaller.
    assert first_maximum <= 2 * sqrt(prime) + 1e-8
    assert second_maximum <= 2 * sqrt(other) + 1e-8

    first = multiplicative_ordering(first_values, prime)
    second = multiplicative_ordering(second_values, other)
    first_order = prime - 1
    second_order = other - 1
    common = gcd(first_order, second_order)
    period = lcm(first_order, second_order)

    left = sum(
        mellin_transform(first, index)
        * mellin_transform(second, index).conjugate()
        for index in range(period)
    )
    right = period * sum(
        first[first_order * residue // common]
        * second[second_order * residue // common].conjugate()
        for residue in range(common)
    )
    assert abs(left - right) < 1e-7 * max(1.0, abs(right))

    # Each t-fibre has at most two points and each squared trace has size at
    # most 4p.  Hence the exact gcd-term expression is <=64 L*g*p*q.
    bound = 64 * period * common * prime * other
    assert abs(right) <= bound + 1e-7
    print(
        "split-cover Euler-square two-prime correlation: "
        f"Weil bound and {common}-term DFT collapse VERIFIED"
    )


def verify_linear_mellin_correlation() -> None:
    prime, other = 13, 29
    first_values = branch_values_from_pullback(prime, branch_for_prime(prime))
    second_values = branch_values_from_pullback(other, branch_for_prime(other))
    first = multiplicative_ordering(first_values, prime)
    second = multiplicative_ordering(second_values, other)

    first_order = prime - 1
    second_order = other - 1
    common = gcd(first_order, second_order)
    period = lcm(first_order, second_order)

    left = sum(
        mellin_transform(first, index)
        * mellin_transform(second, index).conjugate()
        for index in range(period)
    )
    right = period * sum(
        first[first_order * residue // common]
        * second[second_order * residue // common].conjugate()
        for residue in range(common)
    )
    assert abs(left - right) < 1e-7 * max(1.0, abs(right))

    print(
        "linear Mellin full-period identity: "
        f"p={prime}, q={other}, gcd(p-1,q-1)={common}, period={period}, VERIFIED"
    )


def verify_character_of_sum_obstruction() -> None:
    prime = 5
    root = cmath.exp(2j * pi / prime)
    values = (1, 2)
    character_after_sum = root ** (sum(values) % prime)
    sum_after_character = sum(root**value for value in values)
    assert abs(character_after_sum - sum_after_character) > 1e-3

    # This is the exact invalid interchange that would be needed to turn
    # e_p(a * (a Jacobi sum)) into a standard complete sum in its Jacobi
    # variables.
    print("nonlinear zero-detector interchange: explicit F_5 counterexample VERIFIED")


def verify_character_period_zero_correlation() -> None:
    prime, other = 17, 19
    first_period = prime - 1
    second_period = other - 1
    common = gcd(first_period, second_period)
    period = lcm(first_period, second_period)
    first_zeros = zero_set(prime) & set(range(first_period))
    second_zeros = zero_set(other) & set(range(second_period))

    direct = sum(
        index % first_period in first_zeros
        and index % second_period in second_zeros
        for index in range(period)
    )
    stratified = sum(
        sum(index % common == residue for index in first_zeros)
        * sum(index % common == residue for index in second_zeros)
        for residue in range(common)
    )
    assert direct == stratified
    print("zero events over lcm(p-1,q-1): gcd-stratified identity VERIFIED")


def main() -> None:
    verify_crt_averages()
    verify_apery_mellin_and_split_cover()
    verify_linear_mellin_correlation()
    verify_split_cover_trace_correlation()
    verify_character_of_sum_obstruction()
    verify_character_period_zero_correlation()


if __name__ == "__main__":
    main()
