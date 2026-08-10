#!/usr/bin/env python3
"""Exact finite-field checks for toric_mellin_square.tex."""

from math import comb


def primes_through(limit: int):
    return [
        n
        for n in range(2, limit + 1)
        if all(n % d for d in range(2, int(n**0.5) + 1))
    ]


def chi(value: int, prime: int) -> int:
    value %= prime
    if value == 0:
        return 0
    return 1 if pow(value, (prime - 1) // 2, prime) == 1 else -1


def apery(index: int, prime: int) -> int:
    return sum(
        comb(index, k) ** 2 * comb(index + k, k) ** 2
        for k in range(index + 1)
    ) % prime


def legendre_coefficient(index: int, value: int, prime: int) -> int:
    return sum(
        comb(index, k) ** 2 * pow(value, k, prime)
        for k in range(index + 1)
    ) % prime


def inner_sum(index: int, value: int, prime: int) -> int:
    return sum(
        pow(r, index, prime)
        * chi((value - 1) ** 2 * r * r - 2 * (value + 1) * r + 1, prime)
        for r in range(1, prime)
    ) % prime


def lambda_value(x: int, y: int, z: int, prime: int) -> int:
    numerator = (
        (1 + x)
        * (1 + y)
        * (1 + z)
        * ((1 + y) * (1 + z) + x * y * z)
    ) % prime
    denominator = x * y * z % prime
    return numerator * pow(denominator, prime - 2, prime) % prime


def direct_fibers(prime: int):
    result = [0] * prime
    for x in range(1, prime):
        for y in range(1, prime):
            for z in range(1, prime):
                result[lambda_value(x, y, z, prime)] += 1
    return result


def quadratic_fiber(value: int, prime: int) -> int:
    result = (prime - 2) ** 2
    for t in range(1, prime):
        for w in range(1, prime):
            if w == 1 or t * w % prime == 1:
                continue
            correction = (
                t
                * (w - 1)
                * (1 - t * w)
                * pow(value * w % prime, prime - 2, prime)
            ) % prime
            trace = (t + 1 - correction) % prime
            result += chi(trace * trace - 4 * t, prime)
    return result


def polynomial_multiply(left, right, prime: int):
    result = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            result[i + j] = (result[i + j] + a * b) % prime
    return result


def main() -> None:
    fiber_checks = 0
    inner_checks = 0
    pairing_checks = 0
    duality_checks = 0
    tautology_checks = 0

    for prime in primes_through(19):
        if prime < 5:
            continue
        fibers = direct_fibers(prime)
        for value in range(1, prime):
            assert fibers[value] == quadratic_fiber(value, prime)
            fiber_checks += 1

    for prime in primes_through(43):
        if prime < 5:
            continue
        for index in range(1, prime - 1):
            expected = apery(index, prime)
            pairing = 0
            for value in range(1, prime):
                assert inner_sum(index, value, prime) == (
                    -legendre_coefficient(prime - 1 - index, value, prime)
                ) % prime
                inner_checks += 1
                pairing -= (
                    pow(value, index, prime)
                    * legendre_coefficient(index, value, prime)
                    * legendre_coefficient(prime - 1 - index, value, prime)
                )
            assert pairing % prime == expected
            pairing_checks += 1

            convolution = sum(
                comb(index, k) ** 2 * comb(prime - 1 - index, k) ** 2
                for k in range(index + 1)
            ) % prime
            assert convolution == expected
            tautology_checks += 1

        half = (prime - 1) // 2
        for index in range(1, half + 1):
            left = [comb(prime - 1 - index, k) ** 2 % prime
                    for k in range(prime - index)]
            right = [comb(index, k) ** 2 % prime for k in range(index + 1)]
            factor = [comb(prime - 1 - 2 * index, k) * (-1) ** k % prime
                      for k in range(prime - 2 * index)]
            product = polynomial_multiply(right, factor, prime)
            assert left == product
            duality_checks += 1

    print("toric_mellin_square_verify: PASS")
    print(f"fiber_checks={fiber_checks}")
    print(f"inner_checks={inner_checks}")
    print(f"pairing_checks={pairing_checks}")
    print(f"duality_checks={duality_checks}")
    print(f"tautology_checks={tautology_checks}")


if __name__ == "__main__":
    main()
