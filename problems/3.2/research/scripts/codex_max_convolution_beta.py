#!/usr/bin/env python3
"""Verify the convolution shadow and the Bernoulli formula for beta_p."""

from fractions import Fraction as Q
from math import comb, isqrt


LIMIT = 500
ORDER = 60


def primes_below(n):
    sieve = bytearray(b"\x01") * n
    sieve[:2] = b"\x00\x00"
    for p in range(2, isqrt(n - 1) + 1):
        if sieve[p]:
            sieve[p * p : n : p] = b"\x00" * (((n - 1 - p * p) // p) + 1)
    return [p for p in range(2, n) if sieve[p]]


def apery(n):
    values = [1]
    for m in range(n - 1):
        previous = values[m - 1] if m else 0
        numerator = (
            (2 * m + 1) * (17 * m * m + 17 * m + 5) * values[m]
            - m**3 * previous
        )
        quotient, remainder = divmod(numerator, (m + 1) ** 3)
        assert remainder == 0
        values.append(quotient)
    return values


def square_root_series(coefficients):
    result = [Q(1)]
    for n in range(1, len(coefficients)):
        result.append(
            (Q(coefficients[n]) - sum(result[j] * result[n - j] for j in range(1, n))) / 2
        )
    return result


def divide_by_q(coefficients):
    result = []
    for n, value in enumerate(coefficients):
        result.append(Q(value) + (34 * result[n - 1] if n >= 1 else 0) - (result[n - 2] if n >= 2 else 0))
    return result


def convolution(values, n):
    return sum(values[j] * values[n - j] for j in range(n + 1))


def bernoulli_mod(n, p):
    values = [1]
    for m in range(1, n + 1):
        total = sum(comb(m + 1, k) * values[k] for k in range(m)) % p
        values.append(-total * pow(m + 1, -1, p) % p)
    return values[n]


def harmonic_mod(p, exponent, modulus):
    return sum(pow(k, -exponent, modulus) for k in range(1, p)) % modulus


def main():
    b = apery(max(LIMIT, ORDER))
    tau = square_root_series(b[:ORDER])
    sigma = square_root_series(divide_by_q(b[:ORDER]))
    sigma_square = [convolution(sigma, n) for n in range(ORDER)]
    for n in range(ORDER):
        assert convolution(tau, n) == b[n]
        corrected = sigma_square[n]
        if n >= 1:
            corrected -= 34 * sigma_square[n - 1]
        if n >= 2:
            corrected += sigma_square[n - 2]
        assert corrected == b[n]

    beta_checks = 0
    local_term_checks = 0
    harmonic_checks = 0
    for p in primes_below(LIMIT):
        if p < 5:
            continue
        bernoulli = bernoulli_mod(p - 3, p)
        beta = ((b[p] - 5) // (p**3)) % p
        assert (b[p] - 5) % (p**3) == 0
        assert beta == -14 * pow(3, -1, p) * bernoulli % p
        beta_checks += 1

        h1 = harmonic_mod(p, 1, p**3)
        h2 = harmonic_mod(p, 2, p**2)
        h3 = harmonic_mod(p, 3, p)
        assert h1 == -(p * p) * pow(3, -1, p**3) * bernoulli % (p**3)
        assert h2 == 2 * p * pow(3, -1, p**2) * bernoulli % (p**2)
        assert h3 == 0
        harmonic_checks += 3

        central = comb(2 * p - 1, p - 1) % (p**4)
        expected_central = (1 - 2 * p**3 * pow(3, -1, p**4) * bernoulli) % (p**4)
        assert central == expected_central

        if p < 100:
            interior = 0
            approximation = 0
            modulus = p**4
            for k in range(1, p):
                term = comb(p, k) ** 2 * comb(p + k, k) ** 2
                interior = (interior + term) % modulus
                invk = pow(k, -1, modulus)
                approximate_term = p * p * invk * invk * (1 + 2 * p * invk)
                assert term % modulus == approximate_term % modulus
                approximation = (approximation + approximate_term) % modulus
                local_term_checks += 1
            assert interior == approximation

    # Equal zero support does not determine convolution vanishing.
    modulus = 5
    left = [1, 1, 2]
    right = [1, 1, 1]
    assert all(value % modulus for value in left + right)
    assert sum(left[j] * left[2 - j] for j in range(3)) % modulus == 0
    assert sum(right[j] * right[2 - j] for j in range(3)) % modulus != 0

    print(f"tau/sigma convolution identities verified through n={ORDER - 1}")
    print(f"beta_p = -(14/3) B_(p-3) checked for {beta_checks} primes 5 <= p < {LIMIT}")
    print(f"harmonic/Glaisher ingredients checked: {harmonic_checks} congruences")
    print(f"interior summand expansion checked: {local_term_checks} terms")
    print("zero-support-only convolution law: DISPROVED by an explicit F_5 example")


if __name__ == "__main__":
    main()
