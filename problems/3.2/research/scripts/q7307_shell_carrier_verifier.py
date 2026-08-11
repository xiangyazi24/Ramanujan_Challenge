#!/usr/bin/env python3
"""Standalone exact verifier for the shell carrier proposed in Q7299.

This is a finite mechanical audit.  It does not infer asymptotics and it
does not promote any tested identity to a theorem.

For Apéry numbers b_n, Franel numbers F_a, and

    L(n,k) = C(n,k) C(n+k,k),
    Q_n(t) = sum_k L(n,k) t^k,
    K_{n,a}(c) = [y^a] Q_n(c+y),
    g_a(c) = sum_u C(a,u) (-c)^(a-u) F_u,
    T_{n,J}(c) = sum_{a=0}^J K_{n,a}(c) g_a(c),
    Gamma_{n,J} = content_c T_{n,J}(c),

Q7299 defines, for X <= m < X^2,

    E_X = C(m,X) C(m+X,X),
    I_X = gcd(Gamma_{m,X-1}, Gamma_{m,X}),
    C_X = rad(gcd(I_X,E_X)).

The default run performs four independent finite sweeps:

1. Direct formula / recurrence / Strehl cross-checks for b_n and F_n.
2. Explicit polynomial construction versus the closed coefficient formula.
3. Every cutoff Gamma_{n,J} through n=96, including Gamma|b_n and the
   complete local cutoff support test for odd primes p<=197 with n<p^2.
4. Every shell pair 3<=X<=24, X<=m<X^2, plus wider endpoint and
   valuation stress sweeps through X=128 (K=X) and X=40 (every K<p).

Only the Python standard library is used.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import struct
import sys
import time
import traceback
from collections import Counter
from math import comb, gcd
from pathlib import Path
from typing import Any, Iterable


AUDITED_MAIN_COMMIT = "8402484f3765c0a2a1e979bb724478096598c35e"
Q7299_COMMIT = "f0b1b9cafba97cc7cf64c47b4e975e98187c7d75"
Q7299_BLOB = "18e418c5a3bcfcd24cd6a0eef1c4b8dbb8b1c531"
FRANEL_TAIL_PROOF_BLOB = "381f76cb9619f4eea5f37a65b238c55e11069e3b"
ALL_CUTOFF_PROOF_BLOB = "10e1b75e3e0edde437acc89c6ec08317179076ca"
FIXED_Q_SCRIPT_BLOB = "2d537ff8b80fa119f16c2105228795b40966b798"


def canonical_json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


class Audit:
    def __init__(self) -> None:
        self.failures: list[dict[str, Any]] = []
        self.counters: Counter[str] = Counter()

    def count(self, name: str, amount: int = 1) -> None:
        self.counters[name] += amount

    def fail(self, kind: str, **data: Any) -> None:
        record = {"kind": kind, **data}
        self.failures.append(record)
        print("FAILURE " + canonical_json(record), flush=True)

    def check(self, condition: bool, kind: str, **data: Any) -> bool:
        self.count("assertions")
        if not condition:
            self.fail(kind, **data)
            return False
        return True


class JsonDigest:
    def __init__(self) -> None:
        self._hash = hashlib.sha256()
        self.count = 0

    def add(self, record: Any) -> None:
        self._hash.update(canonical_json(record).encode("utf-8"))
        self._hash.update(b"\n")
        self.count += 1

    def hexdigest(self) -> str:
        return self._hash.hexdigest()


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def primes_up_to(limit: int) -> list[int]:
    if limit < 2:
        return []
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[0:2] = b"\x00\x00"
    stop = int(limit**0.5)
    for prime in range(2, stop + 1):
        if sieve[prime]:
            start = prime * prime
            count = (limit - start) // prime + 1
            sieve[start : limit + 1 : prime] = b"\x00" * count
    return [n for n in range(2, limit + 1) if sieve[n]]


def vp_factorial(n: int, prime: int) -> int:
    total = 0
    while n:
        n //= prime
        total += n
    return total


def vp_binomial(n: int, k: int, prime: int) -> int:
    if not 0 <= k <= n:
        raise ValueError((n, k, prime))
    return (
        vp_factorial(n, prime)
        - vp_factorial(k, prime)
        - vp_factorial(n - k, prime)
    )


def franel_numbers(limit: int) -> list[int]:
    return [sum(comb(n, k) ** 3 for k in range(n + 1)) for n in range(limit + 1)]


def apery_direct(n: int) -> int:
    return sum(comb(n, k) ** 2 * comb(n + k, k) ** 2 for k in range(n + 1))


def apery_recurrence(limit: int, audit: Audit) -> list[int]:
    if limit == 0:
        return [1]
    values = [1, 5]
    for n in range(1, limit):
        numerator = (
            (34 * n**3 + 51 * n**2 + 27 * n + 5) * values[n]
            - n**3 * values[n - 1]
        )
        denominator = (n + 1) ** 3
        audit.check(
            numerator % denominator == 0,
            "apery_recurrence_nonintegral",
            n=n,
            numerator=numerator,
            denominator=denominator,
        )
        values.append(numerator // denominator)
        audit.count("apery_recurrence_steps")
    return values


def legendre_row(n: int) -> list[int]:
    return [comb(n, k) * comb(n + k, k) for k in range(n + 1)]


def strehl_value(row: list[int], franel: list[int], cutoff: int) -> int:
    return sum(row[k] * franel[k] for k in range(cutoff + 1))


def fast_cutoff_coefficients(
    n: int, cutoff: int, row: list[int], franel: list[int]
) -> list[int]:
    """Closed integer coefficient formula used independently in the audit."""
    coefficients = [strehl_value(row, franel, cutoff)] + [0] * n
    for degree in range(1, n + 1):
        lower = max(0, cutoff + 1 - degree)
        upper = min(cutoff, n - degree)
        total = 0
        for index in range(lower, upper + 1):
            term = (
                row[index + degree]
                * comb(index + degree, index)
                * comb(degree - 1, cutoff - index)
                * franel[index]
            )
            total += -term if (cutoff - index) & 1 else term
        coefficients[degree] = total
    return coefficients


def gamma_content(
    n: int, cutoff: int, row: list[int], franel: list[int]
) -> tuple[int, int]:
    """Return (Gamma_{n,J}, constant coefficient) by exact incremental gcd."""
    constant = strehl_value(row, franel, cutoff)
    result = abs(constant)
    for degree in range(1, n + 1):
        if result == 1:
            break
        lower = max(0, cutoff + 1 - degree)
        upper = min(cutoff, n - degree)
        residue = 0
        for index in range(lower, upper + 1):
            term = (
                row[index + degree]
                * comb(index + degree, index)
                * comb(degree - 1, cutoff - index)
                * franel[index]
            )
            if (cutoff - index) & 1:
                residue -= term % result
            else:
                residue += term % result
            residue %= result
        result = gcd(result, residue)
    return result, constant


def explicit_cutoff_coefficients(
    n: int, cutoff: int, row: list[int], franel: list[int]
) -> list[int]:
    """Construct K_{n,a}, g_a and T_{n,J} literally by polynomial products."""
    result = [0] * (n + 1)
    for a in range(cutoff + 1):
        # K_{n,a}(c) = sum_{d=0}^{n-a} L(n,a+d) C(a+d,a) c^d.
        k_coefficients = [
            row[a + degree] * comb(a + degree, a)
            for degree in range(n - a + 1)
        ]
        # [c^e] g_a(c) = (-1)^e C(a,e) F_{a-e}.
        g_coefficients = [
            (-1 if exponent & 1 else 1)
            * comb(a, exponent)
            * franel[a - exponent]
            for exponent in range(a + 1)
        ]
        for left_degree, left in enumerate(k_coefficients):
            for right_degree, right in enumerate(g_coefficients):
                result[left_degree + right_degree] += left * right
    return result


def content(coefficients: Iterable[int]) -> int:
    result = 0
    for coefficient in coefficients:
        result = gcd(result, coefficient)
    return abs(result)


def radical_factorization(
    value: int, primes: list[int], factor_bound: int, audit: Audit, **context: Any
) -> tuple[int, list[int]]:
    """Factor a positive divisor of E_X, whose prime factors are <= factor_bound."""
    if value < 1:
        audit.fail("nonpositive_radical_input", value=value, **context)
        return 0, []
    remaining = value
    radical = 1
    factors: list[int] = []
    for prime in primes:
        if prime > factor_bound or prime * prime > remaining:
            break
        if remaining % prime == 0:
            factors.append(prime)
            radical *= prime
            while remaining % prime == 0:
                remaining //= prime
    if remaining > 1:
        audit.check(
            remaining <= factor_bound,
            "radical_residual_exceeds_endpoint_bound",
            residual=remaining,
            factor_bound=factor_bound,
            **context,
        )
        factors.append(remaining)
        radical *= remaining
    return radical, factors


def sequence_hash(values: list[int]) -> str:
    digest = hashlib.sha256()
    for index, value in enumerate(values):
        digest.update(f"{index}:{value}\n".encode("ascii"))
    return digest.hexdigest()


def shell_primes(primes: list[int], x_value: int) -> list[int]:
    return [prime for prime in primes if x_value < prime <= 2 * x_value]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--explicit-n-max", type=int, default=16)
    parser.add_argument("--direct-apery-n-max", type=int, default=96)
    parser.add_argument("--all-cutoff-n-max", type=int, default=96)
    parser.add_argument("--all-cutoff-prime-max", type=int, default=197)
    parser.add_argument("--shell-x-max", type=int, default=24)
    parser.add_argument("--endpoint-all-k-x-max", type=int, default=40)
    parser.add_argument("--endpoint-x-max", type=int, default=128)
    return parser.parse_args()


def run(args: argparse.Namespace) -> int:
    started = time.perf_counter()
    audit = Audit()

    audit.check(args.explicit_n_max >= 0, "bad_configuration", field="explicit_n_max")
    audit.check(args.direct_apery_n_max >= 0, "bad_configuration", field="direct_apery_n_max")
    audit.check(args.all_cutoff_n_max >= 0, "bad_configuration", field="all_cutoff_n_max")
    audit.check(args.all_cutoff_prime_max >= 3, "bad_configuration", field="all_cutoff_prime_max")
    audit.check(args.shell_x_max >= 3, "bad_configuration", field="shell_x_max")
    audit.check(
        args.endpoint_all_k_x_max >= 3,
        "bad_configuration",
        field="endpoint_all_k_x_max",
    )
    audit.check(args.endpoint_x_max >= 3, "bad_configuration", field="endpoint_x_max")

    config = {
        "explicit_n": [0, args.explicit_n_max],
        "direct_apery_n": [0, args.direct_apery_n_max],
        "all_cutoff_n": [0, args.all_cutoff_n_max],
        "all_cutoff_odd_primes": [3, args.all_cutoff_prime_max],
        "exact_shell_X": [3, args.shell_x_max],
        "exact_shell_m_rule": "X <= m < X^2",
        "all_K_endpoint_X": [3, args.endpoint_all_k_x_max],
        "all_K_rule": "0 <= K < p for every X<p<=2X and X<=m<X^2",
        "wide_endpoint_X": [3, args.endpoint_x_max],
        "wide_endpoint_K": "K=X",
    }

    script_path = Path(__file__)
    script_sha256 = sha256_bytes(script_path.read_bytes())
    references = {
        "audited_main_commit": AUDITED_MAIN_COMMIT,
        "q7299_commit": Q7299_COMMIT,
        "q7299_blob": Q7299_BLOB,
        "franel_tail_proof_blob": FRANEL_TAIL_PROOF_BLOB,
        "all_cutoff_proof_blob": ALL_CUTOFF_PROOF_BLOB,
        "fixed_q_script_blob": FIXED_Q_SCRIPT_BLOB,
        "audit_input_commit": os.environ.get("GITHUB_SHA", "not-running-in-github-actions"),
        "github_run_id": os.environ.get("GITHUB_RUN_ID", "not-running-in-github-actions"),
        "verifier_sha256": script_sha256,
    }

    print("Q7307_EXACT_SHELL_CARRIER_AUDIT")
    print("REFERENCES " + canonical_json(references))
    print("CONFIG " + canonical_json(config))

    sequence_limit = max(
        args.explicit_n_max,
        args.direct_apery_n_max,
        args.all_cutoff_n_max,
        args.shell_x_max**2 - 1,
    )
    factor_bound = args.shell_x_max**2 + args.shell_x_max
    prime_limit = max(
        factor_bound,
        args.all_cutoff_prime_max,
        2 * args.endpoint_x_max,
        2 * args.endpoint_all_k_x_max,
    )
    primes = primes_up_to(prime_limit)
    audit.count("sieved_primes", len(primes))

    franel = franel_numbers(sequence_limit)
    apery_rec = apery_recurrence(sequence_limit, audit)
    apery_strehl: list[int] = []
    for n in range(sequence_limit + 1):
        value = strehl_value(legendre_row(n), franel, n)
        apery_strehl.append(value)
        audit.check(
            value == apery_rec[n],
            "apery_strehl_recurrence_mismatch",
            n=n,
            strehl=value,
            recurrence=apery_rec[n],
        )
        audit.count("apery_strehl_full_checks")

    known_franel = [1, 2, 10, 56, 346, 2252]
    known_apery = [1, 5, 73, 1445, 33001, 819005]
    audit.check(
        franel[: len(known_franel)] == known_franel,
        "known_franel_prefix_mismatch",
        actual=franel[: len(known_franel)],
        expected=known_franel,
    )
    audit.check(
        apery_rec[: len(known_apery)] == known_apery,
        "known_apery_prefix_mismatch",
        actual=apery_rec[: len(known_apery)],
        expected=known_apery,
    )

    for n in range(args.direct_apery_n_max + 1):
        direct = apery_direct(n)
        audit.check(
            direct == apery_rec[n],
            "apery_direct_recurrence_mismatch",
            n=n,
            direct=direct,
            recurrence=apery_rec[n],
        )
        audit.count("apery_direct_checks")

    sequence_digests = {
        "franel_0_to_limit": sequence_hash(franel),
        "apery_recurrence_0_to_limit": sequence_hash(apery_rec),
        "apery_strehl_0_to_limit": sequence_hash(apery_strehl),
    }
    print(
        "SEQUENCES "
        + canonical_json(
            {
                "limit": sequence_limit,
                "digests": sequence_digests,
                "direct_apery_limit": args.direct_apery_n_max,
            }
        )
    )

    explicit_digest = JsonDigest()
    for n in range(args.explicit_n_max + 1):
        row = legendre_row(n)
        for cutoff in range(n + 1):
            explicit = explicit_cutoff_coefficients(n, cutoff, row, franel)
            fast = fast_cutoff_coefficients(n, cutoff, row, franel)
            audit.check(
                explicit == fast,
                "explicit_fast_coefficient_mismatch",
                n=n,
                cutoff=cutoff,
                explicit=explicit,
                fast=fast,
            )
            gamma, constant = gamma_content(n, cutoff, row, franel)
            audit.check(
                gamma == content(explicit),
                "explicit_fast_content_mismatch",
                n=n,
                cutoff=cutoff,
                explicit_content=content(explicit),
                fast_content=gamma,
            )
            audit.check(
                constant == explicit[0],
                "strehl_prefix_constant_mismatch",
                n=n,
                cutoff=cutoff,
                constant=constant,
                explicit_constant=explicit[0],
            )
            explicit_digest.add(
                {"n": n, "J": cutoff, "coefficients": explicit, "Gamma": gamma}
            )
            audit.count("explicit_cutoff_polynomials")
            audit.count("explicit_coefficients_compared", n + 1)
        full = explicit_cutoff_coefficients(n, n, row, franel)
        audit.check(
            full == [apery_rec[n]] + [0] * n,
            "explicit_full_strehl_transform_mismatch",
            n=n,
            coefficients=full,
            expected_constant=apery_rec[n],
        )

    print(
        "EXPLICIT_POLYNOMIALS "
        + canonical_json(
            {
                "n_range": [0, args.explicit_n_max],
                "records": explicit_digest.count,
                "sha256": explicit_digest.hexdigest(),
            }
        )
    )

    all_cutoff_digest = JsonDigest()
    cutoff_primes = [p for p in primes if 3 <= p <= args.all_cutoff_prime_max]
    for n in range(args.all_cutoff_n_max + 1):
        row = legendre_row(n)
        running_prefix = 0
        for cutoff in range(n + 1):
            running_prefix += row[cutoff] * franel[cutoff]
            gamma, constant = gamma_content(n, cutoff, row, franel)
            audit.check(
                gamma > 0,
                "nonpositive_gamma",
                n=n,
                cutoff=cutoff,
                gamma=gamma,
            )
            audit.check(
                constant == running_prefix,
                "running_strehl_prefix_mismatch",
                n=n,
                cutoff=cutoff,
                constant=constant,
                running_prefix=running_prefix,
            )
            audit.check(
                apery_rec[n] % gamma == 0,
                "gamma_does_not_divide_apery",
                n=n,
                cutoff=cutoff,
                gamma=gamma,
                apery=apery_rec[n],
            )
            audit.count("gamma_divisibility_checks")
            local_checks = 0
            for prime in cutoff_primes:
                if n >= prime * prime:
                    continue
                residue = n % prime
                folded = min(residue, prime - 1 - residue)
                expected = cutoff >= folded and apery_rec[folded] % prime == 0
                actual = gamma % prime == 0
                audit.check(
                    actual == expected,
                    "all_cutoff_local_support_mismatch",
                    n=n,
                    cutoff=cutoff,
                    prime=prime,
                    quotient=n // prime,
                    residue=residue,
                    folded=folded,
                    gamma=gamma,
                    actual=actual,
                    expected=expected,
                    folded_apery_mod_p=apery_rec[folded] % prime,
                )
                local_checks += 1
            audit.count("all_cutoff_local_prime_checks", local_checks)
            if cutoff == n:
                audit.check(
                    gamma == apery_rec[n],
                    "full_cutoff_content_not_apery",
                    n=n,
                    gamma=gamma,
                    apery=apery_rec[n],
                )
            all_cutoff_digest.add(
                {
                    "n": n,
                    "J": cutoff,
                    "Gamma": gamma,
                    "constant": constant,
                    "local_checks": local_checks,
                }
            )

    print(
        "ALL_CUTOFFS "
        + canonical_json(
            {
                "n_range": [0, args.all_cutoff_n_max],
                "odd_prime_range": [3, args.all_cutoff_prime_max],
                "condition": "n < p^2",
                "records": all_cutoff_digest.count,
                "sha256": all_cutoff_digest.hexdigest(),
            }
        )
    )

    shell_digest = JsonDigest()
    gamma_digest = JsonDigest()
    shell_valuation_histogram: Counter[int] = Counter()
    shell_factor_histogram: Counter[int] = Counter()
    shell_nontrivial_samples: list[dict[str, Any]] = []
    shell_hit_count = 0
    shell_target_rows = 0
    max_gamma_bits = 0
    max_c_bits = 0

    for x_value in range(3, args.shell_x_max + 1):
        x_primes = shell_primes(primes, x_value)
        audit.count("exact_shell_X_values")
        for m in range(x_value, x_value * x_value):
            row = legendre_row(m)
            gamma_left, constant_left = gamma_content(
                m, x_value - 1, row, franel
            )
            gamma_right, constant_right = gamma_content(m, x_value, row, franel)
            prefix_left = strehl_value(row, franel, x_value - 1)
            prefix_right = prefix_left + row[x_value] * franel[x_value]
            audit.check(
                constant_left == prefix_left,
                "shell_left_prefix_mismatch",
                X=x_value,
                m=m,
                constant=constant_left,
                prefix=prefix_left,
            )
            audit.check(
                constant_right == prefix_right,
                "shell_right_prefix_mismatch",
                X=x_value,
                m=m,
                constant=constant_right,
                prefix=prefix_right,
            )
            audit.check(
                apery_rec[m] % gamma_left == 0,
                "shell_left_gamma_does_not_divide_apery",
                X=x_value,
                m=m,
                gamma=gamma_left,
                apery=apery_rec[m],
            )
            audit.check(
                apery_rec[m] % gamma_right == 0,
                "shell_right_gamma_does_not_divide_apery",
                X=x_value,
                m=m,
                gamma=gamma_right,
                apery=apery_rec[m],
            )
            audit.count("exact_shell_gamma_divisibility_checks", 2)

            intersection = gcd(gamma_left, gamma_right)
            endpoint = comb(m, x_value) * comb(m + x_value, x_value)
            gcd_ie = gcd(intersection, endpoint)
            carrier, carrier_factors = radical_factorization(
                gcd_ie,
                primes,
                m + x_value,
                audit,
                X=x_value,
                m=m,
            )
            audit.check(
                apery_rec[m] % intersection == 0,
                "shell_intersection_does_not_divide_apery",
                X=x_value,
                m=m,
                intersection=intersection,
                apery=apery_rec[m],
            )
            audit.check(
                apery_rec[m] % carrier == 0,
                "shell_carrier_does_not_divide_apery",
                X=x_value,
                m=m,
                carrier=carrier,
                apery=apery_rec[m],
            )
            audit.check(
                gcd_ie % carrier == 0,
                "shell_radical_not_divisor_of_gcd",
                X=x_value,
                m=m,
                gcd_ie=gcd_ie,
                carrier=carrier,
            )
            audit.check(
                carrier == (1 if not carrier_factors else __import__("math").prod(carrier_factors)),
                "shell_radical_factor_product_mismatch",
                X=x_value,
                m=m,
                carrier=carrier,
                factors=carrier_factors,
            )
            for factor in carrier_factors:
                audit.check(
                    carrier % (factor * factor) != 0,
                    "shell_carrier_not_squarefree",
                    X=x_value,
                    m=m,
                    carrier=carrier,
                    factor=factor,
                )
                shell_factor_histogram[factor] += 1

            target_primes: list[int] = []
            prime_records: list[dict[str, Any]] = []
            for prime in x_primes:
                quotient, residue = divmod(m, prime)
                folded = min(residue, prime - 1 - residue)
                direct_bad = apery_rec[residue] % prime == 0
                folded_bad = apery_rec[folded] % prime == 0
                if direct_bad:
                    target_primes.append(prime)
                    shell_hit_count += 1
                audit.check(
                    m < prime * prime,
                    "shell_not_two_digit",
                    X=x_value,
                    m=m,
                    prime=prime,
                )
                audit.check(
                    quotient < prime - 1,
                    "shell_endpoint_quotient_boundary_failure",
                    X=x_value,
                    m=m,
                    prime=prime,
                    quotient=quotient,
                )
                audit.check(
                    folded <= x_value - 1,
                    "shell_folded_index_exceeds_left_cutoff",
                    X=x_value,
                    m=m,
                    prime=prime,
                    residue=residue,
                    folded=folded,
                )
                audit.check(
                    direct_bad == folded_bad,
                    "shell_reflection_zero_mismatch",
                    X=x_value,
                    m=m,
                    prime=prime,
                    residue=residue,
                    folded=folded,
                    direct_mod=apery_rec[residue] % prime,
                    folded_mod=apery_rec[folded] % prime,
                )
                left_support = gamma_left % prime == 0
                right_support = gamma_right % prime == 0
                audit.check(
                    left_support == direct_bad,
                    "shell_left_gamma_support_mismatch",
                    X=x_value,
                    m=m,
                    prime=prime,
                    gamma=gamma_left,
                    residue=residue,
                    expected=direct_bad,
                    actual=left_support,
                )
                audit.check(
                    right_support == direct_bad,
                    "shell_right_gamma_support_mismatch",
                    X=x_value,
                    m=m,
                    prime=prime,
                    gamma=gamma_right,
                    residue=residue,
                    expected=direct_bad,
                    actual=right_support,
                )

                endpoint_valuation = vp_binomial(m, x_value, prime) + vp_binomial(
                    m + x_value, x_value, prime
                )
                predicted_valuation = int(residue < x_value) + int(
                    prime - residue <= x_value
                )
                shell_valuation_histogram[endpoint_valuation] += 1
                audit.check(
                    endpoint_valuation == predicted_valuation,
                    "shell_endpoint_valuation_mismatch",
                    X=x_value,
                    m=m,
                    prime=prime,
                    quotient=quotient,
                    residue=residue,
                    exact=endpoint_valuation,
                    predicted=predicted_valuation,
                )
                audit.check(
                    endpoint_valuation in (1, 2),
                    "shell_endpoint_does_not_capture_prime",
                    X=x_value,
                    m=m,
                    prime=prime,
                    valuation=endpoint_valuation,
                )
                audit.check(
                    (endpoint % prime == 0) == (endpoint_valuation > 0),
                    "shell_endpoint_integer_valuation_disagreement",
                    X=x_value,
                    m=m,
                    prime=prime,
                    endpoint_mod=endpoint % prime,
                    valuation=endpoint_valuation,
                )
                carrier_support = carrier % prime == 0
                audit.check(
                    carrier_support == direct_bad,
                    "shell_carrier_support_mismatch",
                    X=x_value,
                    m=m,
                    prime=prime,
                    residue=residue,
                    carrier=carrier,
                    carrier_support=carrier_support,
                    apery_mod=apery_rec[residue] % prime,
                    expected=direct_bad,
                )
                audit.count("exact_shell_prime_equivalences")
                prime_records.append(
                    {
                        "p": prime,
                        "q": quotient,
                        "r": residue,
                        "j": folded,
                        "bad": direct_bad,
                        "endpoint_v": endpoint_valuation,
                    }
                )

            if target_primes:
                shell_target_rows += 1
            nuisance = [
                factor
                for factor in carrier_factors
                if not (x_value < factor <= 2 * x_value)
            ]
            if carrier > 1 and len(shell_nontrivial_samples) < 25:
                shell_nontrivial_samples.append(
                    {
                        "X": x_value,
                        "m": m,
                        "C_X": carrier,
                        "factors": carrier_factors,
                        "target_primes": target_primes,
                        "nuisance_factors": nuisance,
                    }
                )
            max_gamma_bits = max(
                max_gamma_bits, gamma_left.bit_length(), gamma_right.bit_length()
            )
            max_c_bits = max(max_c_bits, carrier.bit_length())
            gamma_digest.add(
                {
                    "X": x_value,
                    "m": m,
                    "Gamma_X_minus_1": gamma_left,
                    "Gamma_X": gamma_right,
                    "I_X": intersection,
                }
            )
            shell_digest.add(
                {
                    "X": x_value,
                    "m": m,
                    "E_X": endpoint,
                    "I_X": intersection,
                    "C_X": carrier,
                    "carrier_factors": carrier_factors,
                    "target_primes": target_primes,
                    "prime_records": prime_records,
                }
            )
            audit.count("exact_shell_pairs")

    print(
        "EXACT_SHELL "
        + canonical_json(
            {
                "X_range": [3, args.shell_x_max],
                "m_rule": "X <= m < X^2",
                "pairs": shell_digest.count,
                "prime_equivalences": audit.counters[
                    "exact_shell_prime_equivalences"
                ],
                "target_hits": shell_hit_count,
                "rows_with_target_hits": shell_target_rows,
                "endpoint_valuation_histogram": dict(
                    sorted(shell_valuation_histogram.items())
                ),
                "gamma_records_sha256": gamma_digest.hexdigest(),
                "shell_records_sha256": shell_digest.hexdigest(),
                "max_gamma_bits": max_gamma_bits,
                "max_C_X_bits": max_c_bits,
                "nontrivial_C_X_samples": shell_nontrivial_samples,
            }
        )
    )

    wide_endpoint_hash = hashlib.sha256()
    wide_endpoint_histogram: Counter[int] = Counter()
    wide_endpoint_records = 0
    for x_value in range(3, args.endpoint_x_max + 1):
        for prime in shell_primes(primes, x_value):
            for m in range(x_value, x_value * x_value):
                quotient, residue = divmod(m, prime)
                folded = min(residue, prime - 1 - residue)
                exact = vp_binomial(m, x_value, prime) + vp_binomial(
                    m + x_value, x_value, prime
                )
                predicted = int(residue < x_value) + int(
                    prime - residue <= x_value
                )
                audit.check(
                    quotient < prime - 1,
                    "wide_endpoint_quotient_boundary_failure",
                    X=x_value,
                    m=m,
                    prime=prime,
                    quotient=quotient,
                )
                audit.check(
                    m + x_value < prime * prime,
                    "wide_endpoint_sum_reaches_p_squared",
                    X=x_value,
                    m=m,
                    prime=prime,
                )
                audit.check(
                    folded <= x_value - 1,
                    "wide_endpoint_folded_cutoff_failure",
                    X=x_value,
                    m=m,
                    prime=prime,
                    folded=folded,
                )
                audit.check(
                    exact == predicted,
                    "wide_endpoint_valuation_mismatch",
                    X=x_value,
                    m=m,
                    prime=prime,
                    residue=residue,
                    exact=exact,
                    predicted=predicted,
                )
                audit.check(
                    exact in (1, 2),
                    "wide_endpoint_support_failure",
                    X=x_value,
                    m=m,
                    prime=prime,
                    exact=exact,
                )
                wide_endpoint_histogram[exact] += 1
                wide_endpoint_hash.update(
                    struct.pack(">IIIII", x_value, m, prime, residue, exact)
                )
                wide_endpoint_records += 1
    audit.count("wide_endpoint_records", wide_endpoint_records)
    print(
        "WIDE_ENDPOINT "
        + canonical_json(
            {
                "X_range": [3, args.endpoint_x_max],
                "K": "X",
                "m_rule": "X <= m < X^2",
                "records": wide_endpoint_records,
                "valuation_histogram": dict(sorted(wide_endpoint_histogram.items())),
                "sha256": wide_endpoint_hash.hexdigest(),
            }
        )
    )

    all_k_hash = hashlib.sha256()
    all_k_histogram: Counter[int] = Counter()
    all_k_records = 0
    for x_value in range(3, args.endpoint_all_k_x_max + 1):
        for prime in shell_primes(primes, x_value):
            for m in range(x_value, x_value * x_value):
                quotient, residue = divmod(m, prime)
                folded = min(residue, prime - 1 - residue)
                audit.check(
                    quotient < prime - 1,
                    "all_K_quotient_boundary_failure",
                    X=x_value,
                    m=m,
                    prime=prime,
                    quotient=quotient,
                )
                for cutoff in range(prime):
                    exact = vp_binomial(m, cutoff, prime) + vp_binomial(
                        m + cutoff, cutoff, prime
                    )
                    predicted = int(residue < cutoff) + int(
                        prime - residue <= cutoff
                    )
                    audit.check(
                        exact == predicted,
                        "all_K_endpoint_valuation_mismatch",
                        X=x_value,
                        m=m,
                        prime=prime,
                        K=cutoff,
                        quotient=quotient,
                        residue=residue,
                        exact=exact,
                        predicted=predicted,
                    )
                    audit.check(
                        (exact > 0) == (cutoff > folded),
                        "all_K_endpoint_cutoff_boundary_mismatch",
                        X=x_value,
                        m=m,
                        prime=prime,
                        K=cutoff,
                        residue=residue,
                        folded=folded,
                        exact=exact,
                        expected_support=cutoff > folded,
                    )
                    audit.check(
                        exact in (0, 1, 2),
                        "all_K_endpoint_valuation_out_of_range",
                        X=x_value,
                        m=m,
                        prime=prime,
                        K=cutoff,
                        exact=exact,
                    )
                    all_k_histogram[exact] += 1
                    all_k_hash.update(
                        struct.pack(">IIIIII", x_value, m, prime, cutoff, residue, exact)
                    )
                    all_k_records += 1
    audit.count("all_K_endpoint_records", all_k_records)
    print(
        "ALL_K_ENDPOINT "
        + canonical_json(
            {
                "X_range": [3, args.endpoint_all_k_x_max],
                "m_rule": "X <= m < X^2",
                "K_rule": "0 <= K < p",
                "records": all_k_records,
                "valuation_histogram": dict(sorted(all_k_histogram.items())),
                "sha256": all_k_hash.hexdigest(),
            }
        )
    )

    deterministic_summary = {
        "verdict": "PASS" if not audit.failures else "FAIL",
        "references": references,
        "config": config,
        "sequence_limit": sequence_limit,
        "sequence_digests": sequence_digests,
        "explicit_polynomials_sha256": explicit_digest.hexdigest(),
        "all_cutoff_records_sha256": all_cutoff_digest.hexdigest(),
        "gamma_records_sha256": gamma_digest.hexdigest(),
        "shell_records_sha256": shell_digest.hexdigest(),
        "wide_endpoint_sha256": wide_endpoint_hash.hexdigest(),
        "all_K_endpoint_sha256": all_k_hash.hexdigest(),
        "shell_endpoint_valuation_histogram": dict(
            sorted(shell_valuation_histogram.items())
        ),
        "wide_endpoint_valuation_histogram": dict(
            sorted(wide_endpoint_histogram.items())
        ),
        "all_K_endpoint_valuation_histogram": dict(
            sorted(all_k_histogram.items())
        ),
        "counters": dict(sorted(audit.counters.items())),
        "failure_count": len(audit.failures),
        "failures": audit.failures,
    }
    summary_text = canonical_json(deterministic_summary)
    summary_sha256 = sha256_bytes(summary_text.encode("utf-8"))
    elapsed = time.perf_counter() - started

    print("FAILURE_COUNT " + str(len(audit.failures)))
    print("SUMMARY_JSON " + summary_text)
    print("SUMMARY_SHA256 " + summary_sha256)
    print(f"ELAPSED_SECONDS {elapsed:.6f}")
    print("MECHANICAL_SCOPE_ONLY no_asymptotics_no_proof_inference")
    return 0 if not audit.failures else 1


def main() -> int:
    try:
        return run(parse_args())
    except Exception as exc:  # A fatal implementation/runtime failure is itself reported.
        print(
            "FATAL_EXCEPTION "
            + canonical_json(
                {
                    "type": type(exc).__name__,
                    "message": str(exc),
                    "traceback": traceback.format_exc(),
                }
            ),
            flush=True,
        )
        return 2


if __name__ == "__main__":
    sys.exit(main())
