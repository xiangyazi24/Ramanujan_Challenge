#!/usr/bin/env python3
"""Exact checks for apery_mellin_carry_polygon.tex."""

from collections import Counter
from fractions import Fraction
from math import comb, isqrt


def primes_up_to(limit: int) -> list[int]:
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[:2] = b"\x00\x00"
    for divisor in range(2, isqrt(limit) + 1):
        if sieve[divisor]:
            start = divisor * divisor
            sieve[start : limit + 1 : divisor] = b"\x00" * (
                (limit - start) // divisor + 1
            )
    return [number for number, flag in enumerate(sieve) if flag]


def morita_gamma_integer(argument: int, prime: int) -> int:
    product = 1
    for factor in range(1, argument):
        if factor % prime:
            product *= factor
    return product if argument % 2 == 0 else -product


def unit_fraction(j: int, k: int, prime: int) -> Fraction:
    return Fraction(
        morita_gamma_integer(j + k + 1, prime),
        morita_gamma_integer(k + 1, prime) ** 2
        * morita_gamma_integer(j - k + 1, prime),
    )


def fraction_mod(value: Fraction, prime: int) -> int:
    numerator = value.numerator % prime
    denominator = value.denominator % prime
    assert denominator
    return numerator * pow(denominator, -1, prime) % prime


def valuation(value: int, prime: int) -> int:
    order = 0
    while value and value % prime == 0:
        value //= prime
        order += 1
    return order


def apery_mod(index: int, prime: int) -> int:
    return sum(
        (comb(index, k) * comb(index + k, k)) ** 2
        for k in range(index + 1)
    ) % prime


def recurrence_zero_set(prime: int) -> tuple[int, ...]:
    previous, current = 1 % prime, 5 % prime
    zeros: list[int] = []
    if previous == 0:
        zeros.append(0)
    if prime > 1 and current == 0:
        zeros.append(1)
    for n in range(1, prime - 1):
        polynomial = (34 * n**3 + 51 * n**2 + 27 * n + 5) % prime
        following = (
            polynomial * current - pow(n, 6, prime) * previous
        ) % prime
        previous, current = current, following
        if current == 0:
            zeros.append(n + 1)
    return tuple(zeros)


def real_sqrt_window_max(zeros: tuple[int, ...], prime: int) -> int:
    maximum = 0
    for left in range(len(zeros)):
        for right in range(left, len(zeros)):
            if (zeros[right] - zeros[left]) ** 2 <= prime:
                maximum = max(maximum, right - left + 1)
    return maximum


gamma_checks = 0
polygon_checks = 0
residue_checks = 0
for prime in primes_up_to(43):
    if prime == 2:
        continue
    for j in range(prime):
        valuations: list[int] = []
        theta = 0
        for k in range(j + 1):
            integer_factor = comb(j, k) * comb(j + k, k)
            carry = int(j + k >= prime)
            unit = unit_fraction(j, k, prime)
            assert Fraction(integer_factor) == prime**carry * unit
            assert Fraction(integer_factor**2) == prime ** (2 * carry) * unit**2
            assert valuation(integer_factor**2, prime) == 2 * carry
            gamma_checks += 3
            valuations.append(2 * carry)
            if not carry:
                theta += fraction_mod(unit**2, prime)
        unit_count = min(j, prime - 1 - j) + 1
        high_count = max(0, 2 * j + 1 - prime)
        assert sorted(valuations) == [0] * unit_count + [2] * high_count
        assert len(valuations) == unit_count + high_count
        if 1 <= j <= prime - 2:
            assert unit_count >= 2
        assert theta % prime == apery_mod(j, prime)
        polygon_checks += 3
        residue_checks += 1


patterns: dict[int, list[int]] = {}
for prime in (11, 13):
    patterns[prime] = [
        (comb(5, k) * comb(5 + k, k)) ** 2 % prime
        for k in range(6)
    ]
assert patterns[11] == [1, 9, 1, 1, 9, 1]
assert patterns[13] == [1, 3, 4, 1, 10, 12]
assert sum(patterns[11]) % 11 == 0
assert sum(patterns[13]) % 13 == 5
b5 = sum((comb(5, k) * comb(5 + k, k)) ** 2 for k in range(6))
assert b5 == 819005
assert valuation(b5, 11) == 1
assert b5 % 13 == 5
counterexample_checks = 7


rows: list[tuple[int, tuple[int, ...]]] = []
for prime in primes_up_to(200):
    recurrence_zeros = recurrence_zero_set(prime)
    direct_zeros = tuple(
        index for index in range(prime) if apery_mod(index, prime) == 0
    )
    assert recurrence_zeros == direct_zeros
    if prime > 2:
        assert all(prime - 1 - index in recurrence_zeros for index in recurrence_zeros)
        assert all(index + 1 not in recurrence_zeros for index in recurrence_zeros)
    rows.append((prime, recurrence_zeros))

histogram = Counter(len(zeros) for _, zeros in rows)
assert histogram == Counter({0: 27, 2: 16, 4: 2, 1: 1})
assert max(len(zeros) for _, zeros in rows) == 4
assert [prime for prime, zeros in rows if len(zeros) == 4] == [181, 193]
local_maximum = max(real_sqrt_window_max(zeros, prime) for prime, zeros in rows)
assert local_maximum == 2
assert [
    prime
    for prime, zeros in rows
    if real_sqrt_window_max(zeros, prime) == local_maximum
] == [5, 19, 37]
audit_checks = len(rows) * 3 + 5


print("apery_mellin_carry_polygon_verify: PASS")
print(f"gamma_checks={gamma_checks}")
print(f"polygon_checks={polygon_checks}")
print(f"residue_checks={residue_checks}")
print(f"counterexample_checks={counterexample_checks}")
print(f"audit_checks={audit_checks}")
