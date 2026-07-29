#!/usr/bin/env python3
"""Round 4, Task B: eta-product coefficients and the Apéry center.

We verify, for every odd prime p <= LIMIT, the congruence

    b_{(p-1)/2} == a(p) (mod p),

where

    eta(2z)^4 eta(4z)^4 = sum_{n>=0} a(n) q^n

and

    b_n = sum_{k=0}^n binom(n,k)^2 binom(n+k,k)^2.

The eta coefficients are computed independently in two exact ways:

1. a logarithmic-derivative recurrence using divisor sums;
2. Euler's pentagonal theorem followed by sparse convolutions.

The Apéry center is also computed independently from its three-term
recurrence and from its defining binomial sum modulo p.  The default report
is written to /tmp/round4_task_b.txt.

The theorem is an odd-prime statement.  At p=2, a(2)=0, but (p-1)/2 is not
an integer, so the center and the Beukers congruence are not defined.  The
report records this boundary case explicitly and excludes p=2 from the list
of non-ordinary primes.
"""

from __future__ import annotations

import argparse
import hashlib
from math import isqrt
from pathlib import Path
from typing import List, Sequence, Tuple


DEFAULT_LIMIT = 10_000
DEFAULT_OUTPUT = Path("/tmp/round4_task_b.txt")


def primes_up_to(limit: int) -> List[int]:
    """Return all primes <= limit by an exact Eratosthenes sieve."""
    if limit < 2:
        return []
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[0:2] = b"\x00\x00"
    for p in range(2, isqrt(limit) + 1):
        if sieve[p]:
            start = p * p
            sieve[start : limit + 1 : p] = b"\x00" * (
                (limit - start) // p + 1
            )
    return [n for n in range(2, limit + 1) if sieve[n]]


def sigma_one_up_to(limit: int) -> List[int]:
    """Return sigma_1(n) for 0 <= n <= limit."""
    sigma = [0] * (limit + 1)
    for divisor in range(1, limit + 1):
        for multiple in range(divisor, limit + 1, divisor):
            sigma[multiple] += divisor
    return sigma


def eta_reduced_coefficients_log_derivative(q_limit: int) -> List[int]:
    """Compute c_m in eta(2z)^4 eta(4z)^4 = q sum c_m q^(2m).

    Put x=q^2 and

        F(x) = product_{n>=1} (1-x^n)^4 (1-x^(2n))^4
             = sum_{m>=0} c_m x^m.

    The logarithmic derivative is

        x F'(x)/F(x)
          = sum_{k>=1} L_k x^k,
        L_k = -4 sigma_1(k) - 8 [2|k] sigma_1(k/2).

    Hence m*c_m = sum_{k=1}^m L_k*c_{m-k}.  Every division below is
    checked to be exact.
    """
    reduced_limit = (q_limit - 1) // 2
    sigma = sigma_one_up_to(reduced_limit)
    logarithmic_derivative = [0] * (reduced_limit + 1)
    for k in range(1, reduced_limit + 1):
        value = -4 * sigma[k]
        if k % 2 == 0:
            value -= 8 * sigma[k // 2]
        logarithmic_derivative[k] = value

    coefficients = [0] * (reduced_limit + 1)
    coefficients[0] = 1
    for m in range(1, reduced_limit + 1):
        numerator = sum(
            logarithmic_derivative[k] * coefficients[m - k]
            for k in range(1, m + 1)
        )
        quotient, remainder = divmod(numerator, m)
        if remainder != 0:
            raise AssertionError(
                f"non-integral logarithmic-derivative coefficient at m={m}"
            )
        coefficients[m] = quotient
    return coefficients


def euler_pentagonal_terms(limit: int) -> List[Tuple[int, int]]:
    """Sparse terms (exponent, coefficient) of prod(1-x^n)."""
    terms: List[Tuple[int, int]] = [(0, 1)]
    k = 1
    while True:
        lower = k * (3 * k - 1) // 2
        if lower > limit:
            break
        sign = -1 if k % 2 else 1
        terms.append((lower, sign))
        upper = k * (3 * k + 1) // 2
        if upper <= limit:
            terms.append((upper, sign))
        k += 1
    terms.sort()
    return terms


def multiply_by_sparse_euler(
    polynomial: Sequence[int],
    pentagonal_terms: Sequence[Tuple[int, int]],
    scale: int,
    limit: int,
) -> List[int]:
    """Multiply by prod(1-x^(scale*n)), truncated through x^limit."""
    result = [0] * (limit + 1)
    for exponent, sign in pentagonal_terms:
        shift = scale * exponent
        if shift > limit:
            break
        stop = limit - shift
        if sign == 1:
            for index in range(stop + 1):
                result[index + shift] += polynomial[index]
        else:
            for index in range(stop + 1):
                result[index + shift] -= polynomial[index]
    return result


def eta_reduced_coefficients_pentagonal(q_limit: int) -> List[int]:
    """Independently compute the reduced eta coefficients by Euler products."""
    reduced_limit = (q_limit - 1) // 2
    pentagonal_terms = euler_pentagonal_terms(reduced_limit)
    coefficients = [0] * (reduced_limit + 1)
    coefficients[0] = 1

    # F(x) = E(x)^4 E(x^2)^4, with E(x)=prod_{n>=1}(1-x^n).
    for _ in range(4):
        coefficients = multiply_by_sparse_euler(
            coefficients, pentagonal_terms, 1, reduced_limit
        )
    for _ in range(4):
        coefficients = multiply_by_sparse_euler(
            coefficients, pentagonal_terms, 2, reduced_limit
        )
    return coefficients


def eta_q_coefficients(reduced: Sequence[int], q_limit: int) -> List[int]:
    """Expand q*F(q^2) as the full list a(0),...,a(q_limit)."""
    coefficients = [0] * (q_limit + 1)
    for m, coefficient in enumerate(reduced):
        exponent = 2 * m + 1
        if exponent <= q_limit:
            coefficients[exponent] = coefficient
    return coefficients


def inverses_up_to(limit: int, p: int) -> List[int]:
    """Return k^(-1) mod p for 1 <= k <= limit < p."""
    inverses = [0] * (limit + 1)
    if limit >= 1:
        inverses[1] = 1
    for k in range(2, limit + 1):
        inverses[k] = (-((p // k) * inverses[p % k])) % p
    return inverses


def apery_center_recurrence(p: int, inverses: Sequence[int]) -> int:
    """Compute b_{(p-1)/2} mod p from the Apéry recurrence."""
    center = (p - 1) // 2
    if center == 0:
        return 1 % p
    previous = 1 % p
    current = 5 % p
    for n in range(1, center):
        n2 = n * n
        n3 = n2 * n
        middle = 34 * n3 + 51 * n2 + 27 * n + 5
        numerator = (middle * current - n3 * previous) % p
        inverse_cube = pow(inverses[n + 1], 3, p)
        following = numerator * inverse_cube % p
        previous, current = current, following
    return current


def apery_center_binomial(p: int, inverses: Sequence[int]) -> int:
    """Compute the same center directly from the defining binomial sum."""
    n = (p - 1) // 2
    binomial_n_k = 1
    binomial_n_plus_k_k = 1
    total = 1
    for k in range(1, n + 1):
        inverse_k = inverses[k]
        binomial_n_k = binomial_n_k * (n - k + 1) * inverse_k % p
        binomial_n_plus_k_k = (
            binomial_n_plus_k_k * (n + k) * inverse_k % p
        )
        term = (
            binomial_n_k
            * binomial_n_k
            * binomial_n_plus_k_k
            * binomial_n_plus_k_k
        ) % p
        total = (total + term) % p
    return total


def coefficient_digest(coefficients: Sequence[int]) -> str:
    """Stable digest of the exact decimal coefficient vector."""
    payload = ",".join(str(value) for value in coefficients).encode("ascii")
    return hashlib.sha256(payload).hexdigest()


def format_q_sample(coefficients: Sequence[int], last_exponent: int) -> str:
    """Format the nonzero terms through q^last_exponent."""
    terms = []
    for exponent in range(last_exponent + 1):
        coefficient = coefficients[exponent]
        if coefficient:
            terms.append(f"a({exponent})={coefficient}")
    return ", ".join(terms)


def build_report(limit: int) -> str:
    if limit < 3:
        raise ValueError("LIMIT must be at least 3")

    reduced_log = eta_reduced_coefficients_log_derivative(limit)
    reduced_pentagonal = eta_reduced_coefficients_pentagonal(limit)
    if reduced_log != reduced_pentagonal:
        mismatch = next(
            index
            for index, pair in enumerate(zip(reduced_log, reduced_pentagonal))
            if pair[0] != pair[1]
        )
        raise AssertionError(
            "eta algorithms disagree at "
            f"q^{2 * mismatch + 1}: {reduced_log[mismatch]} != "
            f"{reduced_pentagonal[mismatch]}"
        )

    eta_coefficients = eta_q_coefficients(reduced_log, limit)
    primes = primes_up_to(limit)
    odd_primes = [p for p in primes if p != 2]

    rows = []
    apery_mismatches = []
    beukers_mismatches = []
    nonordinary = []
    for p in odd_primes:
        center = (p - 1) // 2
        inverses = inverses_up_to(center, p)
        b_recurrence = apery_center_recurrence(p, inverses)
        b_binomial = apery_center_binomial(p, inverses)
        a_p = eta_coefficients[p]
        a_mod_p = a_p % p
        if b_recurrence != b_binomial:
            apery_mismatches.append((p, b_recurrence, b_binomial))
        if b_recurrence != a_mod_p:
            beukers_mismatches.append((p, a_mod_p, b_recurrence))
        if a_mod_p == 0:
            nonordinary.append((p, a_p, center, b_recurrence))
        rows.append((p, a_p, a_mod_p, center, b_recurrence, b_binomial))

    if apery_mismatches:
        raise AssertionError(f"Apéry cross-check failures: {apery_mismatches[:5]}")
    if beukers_mismatches:
        raise AssertionError(f"Beukers congruence failures: {beukers_mismatches[:5]}")

    lines = [
        "ROUND 4 -- TASK B: CENTER-MODULAR VERIFICATION",
        "=" * 72,
        f"limit: {limit}",
        "",
        "Definitions checked against problems/3.2/proof.tex:",
        "  b_0=1, b_1=5, and",
        "  (n+1)^3 b_(n+1) = (34n^3+51n^2+27n+5)b_n - n^3 b_(n-1).",
        "  Equivalently, b_n=sum_k binom(n,k)^2 binom(n+k,k)^2.",
        "",
        "Eta normalization:",
        "  eta(2z)^4 eta(4z)^4",
        "    = q * product_(n>=1) (1-q^(2n))^4 (1-q^(4n))^4",
        "    = q * F(q^2),",
        "  so a(2m)=0 and a(2m+1)=[x^m]F(x).",
        "",
        "Independent exact cross-checks:",
        "  [PASS] logarithmic-derivative/divisor-sum eta expansion equals",
        "         the Euler-pentagonal sparse-product expansion at every",
        f"         coefficient a(n), 0 <= n <= {limit}.",
        f"  [PASS] Apéry recurrence equals the defining binomial sum at the",
        f"         center for all {len(odd_primes)} odd primes <= {limit}.",
        f"  [PASS] Beukers congruence holds for all {len(odd_primes)} odd",
        f"         primes <= {limit}; number of failures = 0.",
        "",
        "Coefficient audit:",
        f"  number of a(n) values (including a(0)): {len(eta_coefficients)}",
        f"  SHA256 of comma-separated [a(0),...,a({limit})]:",
        f"    {coefficient_digest(eta_coefficients)}",
        f"  initial q-expansion data: {format_q_sample(eta_coefficients, 31)}",
        "",
        "Boundary p=2:",
        "  a(2)=0, hence the literal divisibility 2|a(2) is true, but",
        "  (p-1)/2=1/2 is not an integer.  Thus the center and Beukers",
        "  congruence are undefined at p=2.  As in the stated expected",
        "  answer, 'p|a(p)' below means odd (good) primes.",
        "",
        "All odd primes p <= limit with p | a(p):",
    ]

    for p, a_p, center, b_center in nonordinary:
        lines.append(
            f"  p={p}: a(p)={a_p}, a(p)/p={a_p // p}, "
            f"center={(p - 1) // 2}, b_center mod p={b_center}"
        )
    lines.extend(
        [
            f"  count = {len(nonordinary)}",
            f"  primes = {[entry[0] for entry in nonordinary]}",
            "",
            "Complete odd-prime congruence table:",
            "  columns: p, a(p), a(p) mod p, center, "
            "b_center(recurrence), b_center(binomial), status, p|a(p)",
        ]
    )
    for p, a_p, a_mod_p, center, b_recurrence, b_binomial in rows:
        lines.append(
            f"  {p:5d} {a_p:10d} {a_mod_p:5d} {center:5d} "
            f"{b_recurrence:5d} {b_binomial:5d} PASS "
            f"{'YES' if a_mod_p == 0 else 'no'}"
        )

    lines.extend(
        [
            "",
            "Conclusion:",
            f"  Among odd primes p <= {limit}, exactly "
            f"{[entry[0] for entry in nonordinary]} divide a(p).",
            "  Both have b_((p-1)/2)=a(p)=0 (mod p).",
        ]
    )
    return "\n".join(lines) + "\n"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--limit", type=int, default=DEFAULT_LIMIT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    report = build_report(args.limit)
    args.output.write_text(report, encoding="utf-8")
    print(report, end="")
    print(f"Report written to {args.output}")


if __name__ == "__main__":
    main()
