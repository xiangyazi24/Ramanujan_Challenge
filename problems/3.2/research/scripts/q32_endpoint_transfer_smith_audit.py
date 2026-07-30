#!/usr/bin/env python3
"""Exact local Smith audit for the endpoint transfer in Q5615.

For primes p < q < 2p, let T carry the canonical Apéry state from
2p-1 to 2q-1:

    y_(2q-1) = T y_(2p-1).

The entries are computed as exact fractions.  For a 2 by 2 rational
matrix, the two local Smith exponents are determined by the minimum
valuation of its entries and the valuation of its determinant.

The short-gap range q-p < p/2 is the one relevant to the proposed
pure-cross carrier.  There the expected exponents are

    at p: (-3, 3),     at q: (0, 0).

The full range has two endpoint effects:

* q-p = (p+1)/2 gives p-exponents (-3, 0), because 2q-1=3p;
* q-p = p-1 gives q-exponents (0, 3), because 2p-1=q.

The second case cannot occur for an actual top-half target pair with
q <= n < 2p and q > 5, since it forces q=n and b_q = 5 (mod q).
"""

from __future__ import annotations

from fractions import Fraction


LIMIT = 100


def primes_below(limit: int) -> list[int]:
    sieve = bytearray(b"\x01") * limit
    sieve[:2] = b"\x00\x00"
    for prime in range(2, int(limit**0.5) + 1):
        if sieve[prime]:
            sieve[prime * prime : limit : prime] = b"\x00" * (
                (limit - 1 - prime * prime) // prime + 1
            )
    return [value for value, flag in enumerate(sieve) if flag]


def apery_coefficient(index: int) -> int:
    return (
        34 * index**3
        + 51 * index**2
        + 27 * index
        + 5
    )


def multiply(
    left: list[list[Fraction]], right: list[list[Fraction]]
) -> list[list[Fraction]]:
    return [
        [
            sum(left[row][mid] * right[mid][column] for mid in range(2))
            for column in range(2)
        ]
        for row in range(2)
    ]


def step(index: int) -> list[list[Fraction]]:
    denominator = (index + 1) ** 3
    return [
        [
            Fraction(apery_coefficient(index), denominator),
            Fraction(-(index**3), denominator),
        ],
        [Fraction(1), Fraction(0)],
    ]


def endpoint_transfer(prime: int, upper_prime: int) -> list[list[Fraction]]:
    result = [
        [Fraction(1), Fraction(0)],
        [Fraction(0), Fraction(1)],
    ]
    for index in range(2 * prime - 1, 2 * upper_prime - 1):
        result = multiply(step(index), result)
    return result


def integer_valuation(value: int, prime: int) -> int:
    if value == 0:
        return 10**9
    result = 0
    while value % prime == 0:
        value //= prime
        result += 1
    return result


def valuation(value: Fraction, prime: int) -> int:
    return integer_valuation(value.numerator, prime) - integer_valuation(
        value.denominator, prime
    )


def smith_exponents(
    matrix: list[list[Fraction]], prime: int
) -> tuple[int, int]:
    first = min(
        valuation(entry, prime) for row in matrix for entry in row
    )
    determinant = (
        matrix[0][0] * matrix[1][1]
        - matrix[0][1] * matrix[1][0]
    )
    return first, valuation(determinant, prime) - first


def expected_exponents(
    prime: int, upper_prime: int
) -> tuple[tuple[int, int], tuple[int, int]]:
    gap = upper_prime - prime
    at_prime = (-3, 0) if 2 * gap == prime + 1 else (-3, 3)
    at_upper = (0, 3) if gap == prime - 1 else (0, 0)
    return at_prime, at_upper


def gap_continuant(start: int, gap: int) -> int:
    """Return the normalized Apéry Casoratian N_gap(start)."""

    if gap == 0:
        return 0
    previous, current = 0, 1
    for offset in range(1, gap):
        index = start + offset
        previous, current = (
            current,
            apery_coefficient(index) * current - index**6 * previous,
        )
    return current


def cube_product(lower: int, upper: int) -> int:
    result = 1
    for value in range(lower, upper + 1):
        result *= value**3
    return result


def cross_casoratian(first: int, second: int) -> Fraction:
    """Return a_i b_j-b_i a_j in the normalization used in the note."""

    if first == second:
        return Fraction(0)
    if first > second:
        return -cross_casoratian(second, first)
    return Fraction(
        6 * gap_continuant(first, second - first),
        cube_product(first + 1, second),
    )


def check_nested_example() -> None:
    """Audit the correctly oriented n=321 nested reflection example."""

    lower, inner_left, inner_right, upper = 36, 64, 128, 142
    prime, lower_prime = 193, 179

    target_inner = gap_continuant(inner_left, inner_right - inner_left)
    target_outer = gap_continuant(lower, upper - lower)
    assert target_inner % prime == 0
    assert target_inner % lower_prime != 0
    assert target_outer % lower_prime == 0
    assert target_outer % prime != 0

    for first, second in (
        (lower, inner_left),
        (inner_right, upper),
        (lower, inner_right),
        (inner_left, upper),
    ):
        value = gap_continuant(first, second - first)
        assert value % prime != 0
        assert value % lower_prime != 0

    pluecker = (
        cross_casoratian(lower, inner_left)
        * cross_casoratian(inner_right, upper)
        - cross_casoratian(lower, inner_right)
        * cross_casoratian(inner_left, upper)
        + cross_casoratian(lower, upper)
        * cross_casoratian(inner_left, inner_right)
    )
    assert pluecker == 0


def main() -> None:
    primes = primes_below(2 * LIMIT)
    checks = 0
    short_checks = 0
    histogram: dict[
        tuple[tuple[int, int], tuple[int, int]], int
    ] = {}

    for prime in (value for value in primes if 5 <= value < LIMIT):
        for upper_prime in (
            value for value in primes if prime < value < 2 * prime
        ):
            transfer = endpoint_transfer(prime, upper_prime)
            observed = (
                smith_exponents(transfer, prime),
                smith_exponents(transfer, upper_prime),
            )
            expected = expected_exponents(prime, upper_prime)
            assert observed == expected, (
                prime,
                upper_prime,
                observed,
                expected,
            )
            if 2 * (upper_prime - prime) < prime:
                assert observed == ((-3, 3), (0, 0))
                short_checks += 1
            histogram[observed] = histogram.get(observed, 0) + 1
            checks += 1

    print(f"PASS: {checks} exact prime-pair transfers")
    print(f"PASS: {short_checks} short-gap transfers have p=(-3,3), q=(0,0)")
    for exponents, count in sorted(histogram.items()):
        print(f"{exponents}: {count}")
    check_nested_example()
    print("PASS: corrected n=321 nested Pluecker and local-unit audit")


if __name__ == "__main__":
    main()
