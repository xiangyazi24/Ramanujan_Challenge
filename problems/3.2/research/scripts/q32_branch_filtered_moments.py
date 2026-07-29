#!/usr/bin/env python3
"""Audit the q=1 branch-filtered safe Taylor-moment certificates.

For the Legendre--Euler coefficient vector C_d, put

    D_k(eps) = sum_d eps^d binom(d,k) C_d,
    g_m(eps) = gcd(D_1(eps), ..., D_m(eps)).

The affine values obtained from integer-valued filters of degree at most m
are D_0(eps) + g_m(eps) Z.  Their least nonzero absolute value is the
centered residue mu_m(eps).  We intersect this value with the two exact
q=1 factorial carriers

    Bminus = binom(n,K),       Bplus = binom(n+K,K),
    K = floor((n-1)/3) + 1.

The safe cutoff for a branch is one below the least top-half prime in that
carrier.  This script verifies that the top-half support of
gcd(Bbranch, mu_m(eps)) is exactly the corresponding bad-prime support and
prints the remaining (usually small) nuisance factor.
"""

from __future__ import annotations

from math import comb, gcd, log

from q32_fixed_q_content import truncation_coefficients
from q32_newton import apery_numbers
from q32_strehl_gcd import franel_numbers, primes_up_to


INDICES = (400, 600, 800)


def centered_nonzero(value: int, modulus: int) -> int:
    residue = value % modulus
    return min(residue, modulus - residue) if residue else modulus


def safe_residue(coefficients: list[int], sign: int, degree: int) -> int:
    base = sum(
        sign**index * coefficient
        for index, coefficient in enumerate(coefficients)
    )
    moment_gcd = 0
    for order in range(1, degree + 1):
        moment = sum(
            sign**index * comb(index, order) * coefficients[index]
            for index in range(order, len(coefficients))
        )
        moment_gcd = gcd(moment_gcd, moment)
    return centered_nonzero(base, abs(moment_gcd))


def factor_over(value: int, primes: list[int]) -> list[tuple[int, int]]:
    factors = []
    remaining = value
    for prime in primes:
        if prime * prime > remaining:
            break
        exponent = 0
        while remaining % prime == 0:
            remaining //= prime
            exponent += 1
        if exponent:
            factors.append((prime, exponent))
    if remaining > 1:
        factors.append((remaining, 1))
    return factors


def main() -> None:
    limit = max(INDICES)
    primes = primes_up_to(limit)
    franel = franel_numbers(limit)
    apery = apery_numbers(limit // 3 + 2)

    for n in INDICES:
        cutoff = (n - 1) // 3
        boundary = cutoff + 1
        coefficients = truncation_coefficients(n, 1, franel)
        carriers = (
            ("direct", comb(n, boundary)),
            ("reflected", comb(n + boundary, boundary)),
        )
        branch_values = []

        for branch, carrier in carriers:
            support = [
                prime
                for prime in primes
                if n / 2 < prime <= n and carrier % prime == 0
            ]
            degree = min(support) - 1 if support else 0
            sign_values = []
            for sign in (-1, 1):
                residue = (
                    safe_residue(coefficients, sign, degree)
                    if support
                    else 1
                )
                branch_gcd = gcd(carrier, residue)
                bad_support = []
                for prime in support:
                    raw = n - prime
                    folded = min(raw, prime - 1 - raw)
                    if apery[folded] % prime == 0:
                        bad_support.append(prime)
                    assert (branch_gcd % prime == 0) == (
                        apery[folded] % prime == 0
                    )
                bad_radical = 1
                for prime in bad_support:
                    bad_radical *= prime
                nuisance = branch_gcd // gcd(branch_gcd, bad_radical)
                sign_values.append(
                    (
                        log(branch_gcd) / n if branch_gcd > 1 else 0.0,
                        sign,
                        branch_gcd,
                        bad_radical,
                        nuisance,
                    )
                )

            best = min(sign_values)
            branch_values.append(best[2])
            print(
                f"n={n} branch={branch} m={degree} sign={best[1]:+d} "
                f"rate={best[0]:.12f} H={best[2]} "
                f"bad={best[3]} nuisance={best[4]} "
                f"factors={factor_over(best[2], primes)}"
            )

        combined = branch_values[0] * branch_values[1]
        print(
            f"n={n} product_rate="
            f"{(log(combined) / n if combined > 1 else 0.0):.12f}"
        )


if __name__ == "__main__":
    main()
