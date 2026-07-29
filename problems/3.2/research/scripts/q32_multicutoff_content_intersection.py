#!/usr/bin/env python3
"""Audit the intersection of all safe q=1 cutoff-content ideals.

Let Gamma_(n,J) be the content of the Legendre--Euler transform T_(n,J).
The Franel tail theorem proves Gamma_(n,J) | A_n for every cutoff J.
Put

    H = floor((n-1)/3),  M = floor(n/2),
    I_n = gcd_{H <= J <= M} Gamma_(n,J).

Every q=1 bad prime divides I_n: its folded index is at most H and every
cutoff in this interval is still below the prime.

There is an exact fast reduction.  Consecutive cutoffs satisfy

    T_(n,J) - T_(n,J-1) = K_(n,J) g_J,

where g_J is primitive because its top coefficient is +/-F_0=+/-1.
Gauss's lemma therefore gives

    content(T_(n,J)-T_(n,J-1)) = content(K_(n,J)) = kappa_(n,J),

and elementary ideal arithmetic gives

    I_n = gcd(Gamma_(n,H), kappa_(n,H+1), ..., kappa_(n,M)).

Here

    kappa_(n,J)
      = gcd_{J <= k <= n}
          binom(n,k) binom(n+k,k) binom(k,J).

This avoids recomputing every transformed polynomial.  The construction is
a clean target package, not yet a height proof: bounding log(I_n)=o(n)
would already prove the q=1 slice.
"""

from __future__ import annotations

from math import comb, gcd

from q32_fixed_q_content import truncation_content
from q32_newton import apery_numbers
from q32_strehl_gcd import franel_numbers, primes_up_to


SAMPLES = (20, 27, 30, 39, 60, 90, 120, 180, 240, 321)
DIRECT_SAMPLES = (30, 56, 142, 180, 200, 321, 394, 400)


def kappa_against(n: int, cutoff: int, modulus: int) -> int:
    """Return gcd(modulus,kappa_(n,cutoff)) without forming kappa."""

    result = modulus
    for index in range(cutoff, n + 1):
        coefficient = (
            comb(n, index)
            * comb(n + index, index)
            * comb(index, cutoff)
        )
        result = gcd(result, coefficient % result)
        if result == 1:
            break
    return result


def multicutoff_intersection(n: int, franel: list[int]) -> int:
    lower = (n - 1) // 3
    upper = n // 2
    result = truncation_content(n, 1, franel)
    for cutoff in range(lower + 1, upper + 1):
        result = kappa_against(n, cutoff, result)
        if result == 1:
            break
    return result


def direct_multicutoff_intersection(n: int, franel: list[int]) -> int:
    """Intersect every cutoff in the full common direct-safe interval.

    The interval is

        H <= J <= n-H-1 = floor(2*n/3).

    A direct q=1 candidate has p=n-j with j<=H, hence every displayed
    cutoff satisfies j<=J<p.
    """

    lower = (n - 1) // 3
    upper = n - lower - 1
    result = truncation_content(n, 1, franel)
    for cutoff in range(lower + 1, upper + 1):
        result = kappa_against(n, cutoff, result)
        if result == 1:
            break
    return result


def q1_target_radical(
    n: int, apery: list[int], primes: list[int]
) -> int:
    result = 1
    for prime in primes:
        if not n / 2 < prime <= n:
            continue
        folded = min(n - prime, 2 * prime - 1 - n)
        if apery[folded] % prime == 0:
            result *= prime
    return result


def direct_q1_target_radical(
    n: int, apery: list[int], primes: list[int]
) -> int:
    result = 1
    lower = (n - 1) // 3
    for prime in primes:
        if not n / 2 < prime <= n:
            continue
        index = n - prime
        if index <= lower and apery[index] % prime == 0:
            result *= prime
    return result


def large_prime_survives_direct_kappas(n: int, prime: int) -> bool:
    """Test the Lucas support condition when prime**2 > 2*n.

    In this range, if n=q*p+r then q<=(p-1)/2.  The nonzero exponent
    support of Q_n modulo p is therefore the union

        [b*p, b*p+j],  0<=b<=q,

    where j=min(r,p-1-r).  The prime divides every kappa in the central
    third exactly when that support misses the central third.
    """

    assert prime * prime > 2 * n
    quotient, residue = divmod(n, prime)
    folded = min(residue, prime - 1 - residue)
    lower = (n + 2) // 3
    upper = (2 * n) // 3
    return not any(
        multiple * prime <= upper
        and multiple * prime + folded >= lower
        for multiple in range(quotient + 1)
    )


def main() -> None:
    limit = max(max(SAMPLES), max(DIRECT_SAMPLES))
    franel = franel_numbers(limit)
    apery = apery_numbers(limit)
    primes = primes_up_to(limit)

    expected = {
        20: 17,
        27: 19,
        30: 85,
        39: 31,
        60: 1,
        90: 1,
        120: 1,
        180: 55,
        240: 5,
        321: 7_289_417,
    }
    for n in SAMPLES:
        intersection = multicutoff_intersection(n, franel)
        target = q1_target_radical(n, apery, primes)
        assert intersection == expected[n], (n, intersection)
        assert intersection % target == 0, (n, intersection, target)
        print(
            f"n={n} intersection={intersection} "
            f"target={target} nuisance={intersection // target}"
        )

    # The three-hit example is purified exactly.
    assert expected[321] == 179 * 193 * 211

    direct_expected = {
        30: 5,
        56: 61,
        142: 145_885,
        180: 55,
        200: 25_159,
        321: 1,
        394: 85,
        400: 1,
    }
    for n in DIRECT_SAMPLES:
        intersection = direct_multicutoff_intersection(n, franel)
        target = direct_q1_target_radical(n, apery, primes)
        assert intersection == direct_expected[n], (n, intersection)
        assert intersection % target == 0, (n, intersection, target)
        print(
            f"direct n={n} intersection={intersection} "
            f"target={target} nuisance={intersection // target}"
        )

    # Exact large-prime support classification for the difference ideals:
    # among primes p<=n with p^2>2*n, the only survivors are direct q=1
    # candidates.  This verifies the floor-sensitive theorem on a broad
    # finite range; the proof is the central-third integer-multiple lemma.
    check_limit = 2_000
    check_primes = primes_up_to(check_limit)
    for n in range(5, check_limit + 1):
        lower = (n - 1) // 3
        for prime in check_primes:
            if prime * prime <= 2 * n:
                continue
            if prime > n:
                break
            survives = large_prime_survives_direct_kappas(n, prime)
            quotient, residue = divmod(n, prime)
            expected_survival = quotient == 1 and residue <= lower
            assert survives == expected_survival, (
                n,
                prime,
                survives,
                expected_survival,
            )


if __name__ == "__main__":
    main()
