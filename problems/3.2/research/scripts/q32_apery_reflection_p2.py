#!/usr/bin/env python3
"""Verify the exact Apéry reflection congruence modulo p^2.

For

    A_j = sum_{k=0}^j binom(j,k)^2 binom(j+k,k)^2

and an odd prime p with 0 <= j <= (p-1)/2, expansion of the two
binomial factors gives

    A_(p-1-j) = A_j - 2*p*W_j  (mod p^2),

where

    W_j = sum_{k=0}^j
            binom(j,k)^2 binom(j+k,k)^2
            (H_(j+k)-H_(j-k))                 (mod p).

All harmonic denominators are nonzero modulo p.  Terms k>j in
A_(p-1-j) contain a squared factor p and vanish modulo p^2.

Thus p|A_j does not force the normalized reflection difference to vanish:

    (A_(p-1-j)-A_j)/p = -2 W_j                (mod p).

The script verifies the formula termwise for every prime p<=500 and records
explicit bad-zero counterexamples to the hoped-for extra p^2 divisibility.
It then uses the differentiated Apéry recurrence to scan through p<=5000.
Every bad lower-half zero in that range is simple as a zero of the truncated
index polynomial; this is computational evidence, not an all-prime theorem.
"""

from __future__ import annotations

from math import comb

from q32_strehl_gcd import primes_up_to


TERM_LIMIT = 500
RECURRENCE_LIMIT = 5_000


def apery_mod_prime_square(prime: int) -> list[int]:
    modulus = prime * prime
    values = [0] * prime
    values[0] = 1
    if prime > 1:
        values[1] = 5
    for index in range(1, prime - 1):
        polynomial = (
            34 * index**3
            + 51 * index**2
            + 27 * index
            + 5
        )
        values[index + 1] = (
            (
                polynomial * values[index]
                - index**3 * values[index - 1]
            )
            * pow(index + 1, -3, modulus)
        ) % modulus
    return values


def harmonic_correction(
    index: int, prime: int, harmonic: list[int]
) -> int:
    result = 0
    for k in range(1, index + 1):
        harmonic_difference = (
            harmonic[index + k] - harmonic[index - k]
        ) % prime
        term = (
            comb(index, k) ** 2
            * comb(index + k, k) ** 2
        ) % prime
        result = (result + term * harmonic_difference) % prime
    return result


def recurrence_polynomial(index: int) -> int:
    return (
        34 * index**3
        + 51 * index**2
        + 27 * index
        + 5
    )


def recurrence_polynomial_derivative(index: int) -> int:
    return 102 * index**2 + 102 * index + 27


def index_derivatives_mod_prime(
    prime: int, values: list[int]
) -> list[int]:
    """Differentiate the Apéry recurrence in its index parameter."""

    derivatives = [0] * prime
    derivatives[0] = 0
    derivatives[1] = 12 % prime
    for index in range(1, prime - 1):
        rhs = (
            recurrence_polynomial_derivative(index)
            * (values[index] % prime)
            + recurrence_polynomial(index) * derivatives[index]
            - 3 * index**2 * (values[index - 1] % prime)
            - index**3 * derivatives[index - 1]
            - 3 * (index + 1) ** 2 * (values[index + 1] % prime)
        ) % prime
        derivatives[index + 1] = (
            rhs * pow((index + 1) ** 3, -1, prime)
        ) % prime
    return derivatives


def main() -> None:
    counterexamples: list[tuple[int, int, int]] = []
    for prime in primes_up_to(TERM_LIMIT):
        if prime == 2:
            continue
        modulus = prime * prime
        values = apery_mod_prime_square(prime)
        derivatives = index_derivatives_mod_prime(prime, values)
        inverse = [0] * prime
        inverse[1] = 1
        for value in range(2, prime):
            inverse[value] = (
                prime
                - (prime // value) * inverse[prime % value] % prime
            )
        harmonic = [0] * prime
        for value in range(1, prime):
            harmonic[value] = (
                harmonic[value - 1] + inverse[value]
            ) % prime
        # Exclude the central fixed point.  Its index derivative can vanish
        # for the formal reflection-symmetry reason and is not part of the
        # two-branch q=1 range audited here.
        for index in range((prime - 1) // 2):
            value = values[index]
            reflected = values[prime - 1 - index]
            correction = harmonic_correction(index, prime, harmonic)
            assert derivatives[index] == 2 * correction % prime
            assert (
                reflected - value + 2 * prime * correction
            ) % modulus == 0, (prime, index)

            if value % prime == 0:
                normalized_difference = (
                    (reflected - value) // prime
                ) % prime
                assert normalized_difference == -2 * correction % prime
                if normalized_difference:
                    counterexamples.append(
                        (prime, index, normalized_difference)
                    )

    expected = {
        (5, 1, 3),
        (17, 3, 16),
        (19, 8, 17),
        (73, 2, 9),
    }
    assert expected.issubset(set(counterexamples))

    zero_count = 0
    zero_primes: set[int] = set()
    multiple_roots: list[tuple[int, int]] = []
    reflection_coincidences: list[tuple[int, int]] = []
    lower_prime_squares: list[tuple[int, int]] = []
    reflected_prime_squares: list[tuple[int, int]] = []
    for prime in primes_up_to(RECURRENCE_LIMIT):
        if prime == 2:
            continue
        modulus = prime * prime
        values = apery_mod_prime_square(prime)
        derivatives = index_derivatives_mod_prime(prime, values)
        for index in range((prime - 1) // 2):
            reflected_index = prime - 1 - index
            assert values[reflected_index] == (
                values[index] - prime * derivatives[index]
            ) % modulus
            if values[index] % prime:
                continue
            zero_count += 1
            zero_primes.add(prime)
            if derivatives[index] == 0:
                multiple_roots.append((prime, index))
            if (
                values[reflected_index] - values[index]
            ) % modulus == 0:
                reflection_coincidences.append((prime, index))
            if values[index] % modulus == 0:
                lower_prime_squares.append((prime, index))
            if values[reflected_index] % modulus == 0:
                reflected_prime_squares.append((prime, index))

    assert not multiple_roots
    assert not reflection_coincidences
    assert lower_prime_squares == [(17, 3)]
    assert not reflected_prime_squares
    print(
        "Apéry reflection mod p^2 verified through "
        f"p<={TERM_LIMIT}; selective counterexamples={len(counterexamples)}"
    )
    print("first counterexamples:", counterexamples[:12])
    print(
        f"recurrence scan through p<={RECURRENCE_LIMIT}: "
        f"lower_half_zeros={zero_count}, zero_primes={len(zero_primes)}, "
        f"multiple_roots={multiple_roots}, "
        f"reflection_p2_coincidences={reflection_coincidences}, "
        f"lower_p2={lower_prime_squares}, "
        f"reflected_p2={reflected_prime_squares}"
    )


if __name__ == "__main__":
    main()
