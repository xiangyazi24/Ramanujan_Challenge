#!/usr/bin/env python3
"""Audit the exact all-cutoff profile in the q=1 prime slice.

For a cutoff J, let T_(n,J)(c) be the Legendre--Euler transform.  Write

    n = p+r,  0 <= r < p,  j=min(r,p-1-r),

where p is an odd prime.  Put x=c^p.  Lucas factorization of the Legendre
kernel and of the shifted Franel moments gives

    T_(n,J) = (1+2x) T_(r,J)                         (J<p),

and, for J=p+B,

    T_(n,p+B)
      = A_r(1+2x)+2(2-x)T_(r,B)                    (mod p).

Modulo p the effective degree of the Legendre kernel Q_r is j.  Therefore
T_(r,J)=A_r once J>=j.  If J<j, the degree-j coefficient of T_(r,J) is

    (-1)^J [t^j]Q_r(t) binom(j-1,J),

which is nonzero modulo p.  Since A_r=A_j modulo p, a bad q=1 prime
p|A_j has the exact cutoff-content profile

    p | content(T_(n,J))
      iff J belongs to [j,p-1] union [p+j,n].

Thus each bad prime has two divisibility plateaux separated by the exact
gap [p,p+j-1].  For a direct prime r=j the high plateau consists only of
J=n; for a reflected prime r=p-1-j it has length p-2j.
"""

from __future__ import annotations

from math import comb

from q32_legendre_content import franel_numbers, primes_up_to
from q32_newton import apery_numbers


PRIME_LIMIT = 31


def transform_coefficients_mod(
    n: int, cutoff: int, prime: int, franel: list[int]
) -> list[int]:
    """Return the monomial coefficients of T_(n,cutoff) modulo prime."""

    coefficients = [0] * (n + 1)
    coefficients[0] = sum(
        comb(n, index)
        * comb(n + index, index)
        * franel[index]
        for index in range(cutoff + 1)
    ) % prime

    for degree in range(1, n + 1):
        lower = max(0, cutoff + 1 - degree)
        upper = min(cutoff, n - degree)
        coefficients[degree] = sum(
            (-1) ** (cutoff - index)
            * comb(n, index + degree)
            * comb(n + index + degree, index + degree)
            * comb(index + degree, index)
            * comb(degree - 1, cutoff - index)
            * franel[index]
            for index in range(lower, upper + 1)
        ) % prime
    return coefficients


def add_scaled_shift(
    target: list[int],
    source: list[int],
    scale: int,
    shift: int,
    prime: int,
) -> None:
    for degree, coefficient in enumerate(source):
        target[degree + shift] = (
            target[degree + shift] + scale * coefficient
        ) % prime


def predicted_coefficients(
    n: int,
    cutoff: int,
    prime: int,
    franel: list[int],
    apery: list[int],
) -> list[int]:
    residue = n - prime
    result = [0] * (n + 1)

    if cutoff < prime:
        lower = transform_coefficients_mod(
            residue, min(cutoff, residue), prime, franel
        )
        add_scaled_shift(result, lower, 1, 0, prime)
        add_scaled_shift(result, lower, 2, prime, prime)
        return result

    lower_cutoff = cutoff - prime
    lower = transform_coefficients_mod(
        residue, min(lower_cutoff, residue), prime, franel
    )
    folded_value = apery[residue] % prime
    result[0] = folded_value
    result[prime] = 2 * folded_value % prime
    add_scaled_shift(result, lower, 4, 0, prime)
    add_scaled_shift(result, lower, -2, prime, prime)
    return result


def content_is_zero(coefficients: list[int], prime: int) -> bool:
    return all(coefficient % prime == 0 for coefficient in coefficients)


def main() -> None:
    primes = [prime for prime in primes_up_to(PRIME_LIMIT) if prime > 2]
    limit = 2 * PRIME_LIMIT
    franel = franel_numbers(limit)
    apery = apery_numbers(limit)
    identity_count = 0
    selective_count = 0

    for prime in primes:
        for folded in range((prime - 1) // 2 + 1):
            for residue in {folded, prime - 1 - folded}:
                n = prime + residue
                for cutoff in range(n + 1):
                    actual = transform_coefficients_mod(
                        n, cutoff, prime, franel
                    )
                    predicted = predicted_coefficients(
                        n, cutoff, prime, franel, apery
                    )
                    assert actual == predicted, (
                        prime,
                        residue,
                        cutoff,
                    )
                    identity_count += 1

                    if apery[folded] % prime:
                        continue
                    expected = (
                        folded <= cutoff < prime
                        or prime + folded <= cutoff <= n
                    )
                    assert content_is_zero(actual, prime) == expected, (
                        prime,
                        residue,
                        folded,
                        cutoff,
                    )
                    selective_count += 1

                    if cutoff < folded:
                        effective = transform_coefficients_mod(
                            residue, cutoff, prime, franel
                        )
                        assert effective[folded] % prime
                    elif prime <= cutoff < prime + folded:
                        effective = transform_coefficients_mod(
                            residue, cutoff - prime, prime, franel
                        )
                        assert effective[folded] % prime

    print(
        "q=1 all-cutoff identities verified for odd primes through "
        f"{PRIME_LIMIT}: identities={identity_count}, "
        f"selective_bad_cases={selective_count}"
    )


if __name__ == "__main__":
    main()
