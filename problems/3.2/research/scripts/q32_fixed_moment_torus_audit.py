#!/usr/bin/env python3
"""Audit the fixed-moment finite-torus identity for Apéry numbers.

For n = a*p + r and M = n - a = a*(p-1) + r, the script checks

    b_r = -sum_{(x,y,z) in (F_p^*)^3} Lambda(x,y,z)^M  (mod p)

for every interior residue 1 <= r <= p-2.  It also checks the constant
term model for small exponents.
"""

from collections import Counter, defaultdict
from math import comb


def primes_up_to(limit):
    sieve = [True] * (limit + 1)
    sieve[:2] = [False, False]
    for p in range(2, int(limit**0.5) + 1):
        if sieve[p]:
            sieve[p * p : limit + 1 : p] = [False] * (
                (limit - p * p) // p + 1
            )
    return [p for p, is_prime in enumerate(sieve) if is_prime]


def apery(n):
    return sum(comb(n, k) ** 2 * comb(n + k, k) ** 2 for k in range(n + 1))


def multiply(left, right):
    out = defaultdict(int)
    for u, a in left.items():
        for v, b in right.items():
            out[tuple(u[i] + v[i] for i in range(3))] += a * b
    return dict(out)


def lambda_polynomial():
    one = {(0, 0, 0): 1}
    x = {(0, 0, 0): 1, (1, 0, 0): 1}
    y = {(0, 0, 0): 1, (0, 1, 0): 1}
    z = {(0, 0, 0): 1, (0, 0, 1): 1}
    yz = multiply(y, z)
    bracket = dict(yz)
    bracket[(1, 1, 1)] = bracket.get((1, 1, 1), 0) + 1
    numerator = one
    for factor in (x, y, z, bracket):
        numerator = multiply(numerator, factor)
    return {
        (u[0] - 1, u[1] - 1, u[2] - 1): coefficient
        for u, coefficient in numerator.items()
    }


LAMBDA = lambda_polynomial()


def constant_terms(limit):
    power = {(0, 0, 0): 1}
    values = []
    for n in range(limit + 1):
        values.append(power.get((0, 0, 0), 0))
        power = multiply(power, LAMBDA)
    return values


def lambda_value(x, y, z, p):
    numerator = (
        (1 + x)
        * (1 + y)
        * (1 + z)
        * ((1 + y) * (1 + z) + x * y * z)
    )
    denominator = x * y * z
    return numerator * pow(denominator, -1, p) % p


def value_distribution(p):
    counts = Counter()
    for x in range(1, p):
        for y in range(1, p):
            for z in range(1, p):
                counts[lambda_value(x, y, z, p)] += 1
    return counts


def torus_moment(counts, exponent, p):
    return sum(count * pow(value, exponent, p) for value, count in counts.items()) % p


def main():
    assert sum(LAMBDA.values()) == 40
    assert all(-1 <= coordinate <= 1 for u in LAMBDA for coordinate in u)

    ct_values = constant_terms(8)
    for n, value in enumerate(ct_values):
        assert value == apery(n), (n, value, apery(n))

    checks = 0
    for p in primes_up_to(43):
        counts = value_distribution(p)
        for a in range(1, 6):
            for r in range(1, p - 1):
                moment = torus_moment(counts, a * (p - 1) + r, p)
                assert (-moment - apery(r)) % p == 0, (p, a, r, moment)
                checks += 1

    print(
        "PASS:",
        f"CT Lambda^n=b_n for n<=8; fixed-moment identity in {checks} cases",
    )


if __name__ == "__main__":
    main()
