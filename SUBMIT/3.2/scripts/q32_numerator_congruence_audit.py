#!/usr/bin/env python3
"""p^3 a_n = (6/5) b_n (mod p) for 7 <= p, n/2 < p <= n.

This makes the denominator-defect law a theorem in the direction the reduction needs:
if p does not divide b_n then p^3 a_n is a p-unit, so v_p(D_n) = 3 and e_p(n) = 0.
Hence only primes dividing b_n can contribute to G_n = gcd(d_n a_n, d_n b_n).

The constant 6/5 = a_1/b_1 is the Apery-Lucas rank-one structure for the vector solution:
modulo p the pair (a_{p+r}, b_{p+r}) is proportional to (a_1, b_1) = (6,5).

Audit below: every pair (n,p) with 7 <= p, n/2 < p <= n, 20 <= n < 260.
"""

from fractions import Fraction
import sys


def sequences(limit):
    a = {0: Fraction(0), 1: Fraction(6)}
    b = {0: Fraction(1), 1: Fraction(5)}
    for n in range(1, limit):
        lead = Fraction(34 * n ** 3 + 51 * n ** 2 + 27 * n + 5, (n + 1) ** 3)
        tail = Fraction(n ** 3, (n + 1) ** 3)
        a[n + 1] = lead * a[n] - tail * a[n - 1]
        b[n + 1] = lead * b[n] - tail * b[n - 1]
    return a, b


def primes_up_to(limit):
    sieve = bytearray([1]) * (limit + 1)
    sieve[0:2] = b"\x00\x00"
    for i in range(2, int(limit ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i::i] = bytearray(len(sieve[i * i::i]))
    return [i for i in range(limit + 1) if sieve[i]]


def audit(limit=260):
    a, b = sequences(limit)
    checks = 0
    divisible = 0
    for n in range(20, limit):
        for prime in primes_up_to(n):
            if prime < 7 or prime * 2 <= n:
                continue
            scaled = a[n] * prime ** 3
            assert scaled.denominator % prime, (n, prime)
            left = scaled.numerator * pow(scaled.denominator, -1, prime) % prime
            right = 6 * pow(5, -1, prime) % prime * (int(b[n]) % prime) % prime
            assert left == right, (n, prime, left, right)
            checks += 1
            if int(b[n]) % prime == 0:
                divisible += 1
    return checks, divisible


if __name__ == "__main__":
    checks, divisible = audit()
    print("NUMERATOR_CONGRUENCE_CHECKS", checks)
    print("OF_WHICH_p_DIVIDES_b_n", divisible)
    print("Q32_NUMERATOR_CONGRUENCE_AUDIT=PASS")
