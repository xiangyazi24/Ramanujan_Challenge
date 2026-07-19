#!/usr/bin/env python3
"""Exact exploration of the separated Apéry gap resultants.

The script performs the finite computation requested in
``CODEX_SPEC_meso.md``.  It uses exact integer polynomials throughout.
Numerical output is evidence and reconnaissance, not a proof of any
uniform estimate.

Run from the repository root with

    python3 problems/3.2/meso_explore.py

SageMath 10.0 or later is sufficient.  If the file is launched with ordinary
Python it re-executes itself under ``sage -python``.  The default run computes every
``R_{d,r}`` with ``2 <= d,r`` and ``d+r <= 24``, removes all prime factors
up to 5000, and computes the exact affine-root statistics for every prime
up to 5000.  Full factorisation of the largest (2246-digit) resultants is
deliberately not attempted; the output gives a length-and-hash fingerprint
of the remaining cofactor instead of silently treating it as prime.
"""

from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import dataclass
import hashlib
from math import sqrt
import os
from pathlib import Path
import shutil
import sys
import time
from typing import Iterable

SAGE_CACHE = "/tmp/ramanujan-meso-sage"
os.environ.setdefault("DOT_SAGE", SAGE_CACHE)
Path(os.environ["DOT_SAGE"]).mkdir(parents=True, exist_ok=True)

try:
    from sage.all import GF, PolynomialRing, ZZ  # type: ignore[import-not-found]
except ModuleNotFoundError:
    sage = shutil.which("sage")
    if sage is None:
        raise SystemExit("SageMath is required (the `sage` executable was not found)")
    environment = os.environ.copy()
    environment["DOT_SAGE"] = SAGE_CACHE
    os.execvpe(sage, [sage, "-python", __file__, *sys.argv[1:]], environment)


POLYNOMIAL_RING = PolynomialRing(ZZ, "x")
X = POLYNOMIAL_RING.gen()


def apery_coefficient(argument):
    """The cubic P(argument) in the cleared Apéry recurrence."""

    return (2 * argument + 1) * (
        17 * argument**2 + 17 * argument + 5
    )


def build_gap_polynomials(max_h: int) -> list:
    """Return ``N_0,...,N_max_h`` in ``ZZ[x]``."""

    if max_h < 1:
        raise ValueError("max_h must be positive")
    polynomials = [POLYNOMIAL_RING.zero(), POLYNOMIAL_RING.one()]
    for h in range(1, max_h):
        next_polynomial = POLYNOMIAL_RING(
            apery_coefficient(X + h) * polynomials[h]
            - (X + h) ** 6 * polynomials[h - 1]
        )
        polynomials.append(next_polynomial)
    return polynomials


def ordered_pairs(height: int) -> list[tuple[int, int]]:
    return [
        (d, r)
        for d in range(2, height - 1)
        for r in range(2, height - d + 1)
    ]


def canonical_pairs(height: int) -> list[tuple[int, int]]:
    return [(d, r) for d, r in ordered_pairs(height) if d <= r]


def shifted(poly, amount: int):
    return POLYNOMIAL_RING(poly(X + amount))


def separated_resultant(
    n_polynomials: list, d: int, r: int
) -> int:
    """Return ``abs(Res_x(N_d(x),N_r(x+d)))`` exactly."""

    value = n_polynomials[d].resultant(shifted(n_polynomials[r], d))
    return abs(int(value))


def exact_resultants(
    n_polynomials: list, height: int
) -> dict[tuple[int, int], int]:
    return {
        pair: separated_resultant(n_polynomials, *pair)
        for pair in canonical_pairs(height)
    }


def lookup_resultant(
    resultants: dict[tuple[int, int], int], d: int, r: int
) -> int:
    return resultants[(min(d, r), max(d, r))]


def prime_sieve(limit: int) -> list[int]:
    if limit < 2:
        return []
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[0:2] = b"\x00\x00"
    for prime in range(2, int(limit**0.5) + 1):
        if sieve[prime]:
            start = prime * prime
            sieve[start : limit + 1 : prime] = b"\x00" * (
                (limit - start) // prime + 1
            )
    return [number for number, flag in enumerate(sieve) if flag]


def trial_factor(
    value: int, primes: Iterable[int]
) -> tuple[tuple[tuple[int, int], ...], int]:
    """Remove the supplied primes and return factors plus exact cofactor."""

    cofactor = abs(value)
    factors: list[tuple[int, int]] = []
    for prime in primes:
        if cofactor % prime:
            continue
        exponent = 0
        while cofactor % prime == 0:
            cofactor //= prime
            exponent += 1
        factors.append((prime, exponent))
    return tuple(factors), cofactor


def apery_values(max_index: int) -> list[int]:
    values = [1, 5]
    for n in range(1, max_index):
        numerator = int(apery_coefficient(n)) * values[n] - n**3 * values[n - 1]
        denominator = (n + 1) ** 3
        assert numerator % denominator == 0
        values.append(numerator // denominator)
    return values[: max_index + 1]


def pell_values(max_index: int) -> list[int]:
    values = [0, 1]
    for _ in range(1, max_index):
        values.append(34 * values[-1] - values[-2])
    return values[: max_index + 1]


def exact_digest(resultants: dict[tuple[int, int], int]) -> str:
    payload = "\n".join(
        f"{d},{r},{resultants[d, r]}" for d, r in sorted(resultants)
    )
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def cofactor_tag(cofactor: int) -> str:
    if cofactor == 1:
        return "1"
    decimal = str(cofactor)
    digest = hashlib.sha256(decimal.encode("ascii")).hexdigest()[:16]
    return f"C{len(decimal)}[{digest}]"


def format_factorisation(
    factors: tuple[tuple[int, int], ...], cofactor: int
) -> str:
    pieces = [str(p) if exponent == 1 else f"{p}^{exponent}" for p, exponent in factors]
    if cofactor != 1:
        pieces.append(cofactor_tag(cofactor))
    return " * ".join(pieces) if pieces else "1"


def polynomial_data_mod_prime(
    n_polynomials: list, max_h: int, prime: int
) -> tuple[list, list[set[int]]]:
    """Reduce the ``N_h`` and compute their F_p-roots exactly."""

    finite_field = GF(prime)
    finite_ring = PolynomialRing(finite_field, "x")
    reductions = [finite_ring.zero() for _ in range(max_h + 1)]
    roots = [set() for _ in range(max_h + 1)]
    for h in range(2, max_h + 1):
        reduced = finite_ring(n_polynomials[h])
        reductions[h] = reduced
        if reduced.is_zero():
            roots[h] = set(range(prime))
        else:
            roots[h] = {int(root) for root in reduced.roots(multiplicities=False)}
    return reductions, roots


@dataclass(frozen=True)
class PrimeStatistics:
    prime: int
    resultant_support: int
    algebraic_support: int
    fp_support: int
    energy: int
    max_multiplicity: int

    @property
    def infinity_support(self) -> int:
        return self.resultant_support - self.algebraic_support

    @property
    def extension_only_support(self) -> int:
        return self.algebraic_support - self.fp_support


def prime_statistics(
    height: int,
    prime_limit: int,
    n_polynomials: list,
    resultants: dict[tuple[int, int], int],
) -> list[PrimeStatistics]:
    pairs = ordered_pairs(height)
    max_h = height - 2
    output: list[PrimeStatistics] = []
    for prime in prime_sieve(prime_limit):
        reductions, roots = polynomial_data_mod_prime(
            n_polynomials, max_h, prime
        )
        resultant_support = 0
        algebraic_support = 0
        fp_support = 0
        energy = 0
        max_multiplicity = 0
        for d, r in pairs:
            if lookup_resultant(resultants, d, r) % prime == 0:
                resultant_support += 1
            first = reductions[d]
            second = reductions[r](reductions[r].parent().gen() + d)
            if first.is_zero() and second.is_zero():
                has_algebraic_root = True
            elif first.is_zero():
                has_algebraic_root = second.degree() > 0
            elif second.is_zero():
                has_algebraic_root = first.degree() > 0
            else:
                has_algebraic_root = first.gcd(second).degree() > 0
            if has_algebraic_root:
                algebraic_support += 1
            multiplicity = sum(
                1 for x_value in roots[d] if (x_value + d) % prime in roots[r]
            )
            if multiplicity:
                fp_support += 1
                energy += multiplicity
                max_multiplicity = max(max_multiplicity, multiplicity)
        output.append(
            PrimeStatistics(
                prime,
                resultant_support,
                algebraic_support,
                fp_support,
                energy,
                max_multiplicity,
            )
        )
    return output


def distribution(records: Iterable[PrimeStatistics], field: str) -> str:
    counts = Counter(getattr(record, field) for record in records)
    return "{" + ", ".join(f"{key}:{counts[key]}" for key in sorted(counts)) + "}"


def print_resultant_summary(
    height: int,
    factor_bound: int,
    resultants: dict[tuple[int, int], int],
) -> None:
    primes = prime_sieve(factor_bound)
    print("EXACT RESULTANTS")
    print(f"H={height}; symmetric representatives={len(resultants)}")
    print(f"ordered pairs={len(ordered_pairs(height))}")
    print(f"sha256={exact_digest(resultants)}")
    largest_pair = max(resultants, key=lambda pair: len(str(resultants[pair])))
    print(
        f"largest=R_{largest_pair[0]},{largest_pair[1]}; "
        f"digits={len(str(resultants[largest_pair]))}"
    )
    print(f"certified trial-factor bound={factor_bound}")
    print("d r digits factorisation-below-bound-and-cofactor-fingerprint")
    for d, r in sorted(resultants, key=lambda pair: (sum(pair), pair)):
        factors, cofactor = trial_factor(resultants[d, r], primes)
        print(
            f"{d:2d} {r:2d} {len(str(resultants[d, r])):4d} "
            f"{format_factorisation(factors, cofactor)}"
        )


def print_sequence_overlap(
    height: int,
    factor_bound: int,
    resultants: dict[tuple[int, int], int],
) -> None:
    primes = prime_sieve(factor_bound)
    apery = apery_values(height - 2)
    pell = pell_values(height - 2)
    apery_support = {p for p in primes if any(value % p == 0 for value in apery[1:])}
    pell_support = {p for p in primes if any(value % p == 0 for value in pell[1:])}
    resultant_support = {
        p for p in primes if any(value % p == 0 for value in resultants.values())
    }
    print("\nSEQUENCE-FACTOR OVERLAP")
    print(f"resultant primes <= {factor_bound}: {len(resultant_support)}")
    print(f"also divide an Apéry b_j, j<={height-2}: {len(resultant_support & apery_support)}")
    print(f"also divide a Pell ell_j, j<={height-2}: {len(resultant_support & pell_support)}")
    unexplained = sorted(resultant_support - apery_support - pell_support)
    print(f"in neither known sequence: {len(unexplained)}")
    print("first unexplained primes:", " ".join(map(str, unexplained[:40])))


def print_prime_summary(
    height: int, prime_limit: int, records: list[PrimeStatistics]
) -> None:
    mesoscopic = [record for record in records if record.prime >= height**2]
    print("\nPRIME STATISTICS")
    print(f"primes <= {prime_limit}: {len(records)}")
    print(
        f"mesoscopic primes {height**2} <= p <= {prime_limit}: "
        f"{len(mesoscopic)}"
    )
    for name in (
        "resultant_support",
        "algebraic_support",
        "fp_support",
        "energy",
    ):
        print(f"{name}: {distribution(mesoscopic, name)}")
    if mesoscopic:
        max_v = max(mesoscopic, key=lambda record: record.resultant_support)
        max_e = max(mesoscopic, key=lambda record: record.energy)
        print(
            f"max V={max_v.resultant_support} at p={max_v.prime}; "
            f"V/H={max_v.resultant_support/height:.6f}"
        )
        print(
            f"max E={max_e.energy} at p={max_e.prime}; "
            f"E/H^(3/2)={max_e.energy/(height*sqrt(height)):.6f}"
        )
        false_only = [
            record.prime
            for record in mesoscopic
            if record.resultant_support and not record.fp_support
        ]
        print("V>0 but F_p support=0:", " ".join(map(str, false_only)))
    all_max_v = max(records, key=lambda record: record.resultant_support)
    print(
        f"all-prime maximum V={all_max_v.resultant_support} "
        f"at p={all_max_v.prime}"
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--height", type=int, default=24)
    parser.add_argument("--prime-limit", type=int, default=5000)
    parser.add_argument("--factor-bound", type=int, default=5000)
    parser.add_argument(
        "--skip-factor-table",
        action="store_true",
        help="omit the 121-row certified partial-factor table",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.height < 4:
        raise SystemExit("height must be at least 4")
    started = time.monotonic()
    n_polynomials = build_gap_polynomials(args.height - 2)
    resultants = exact_resultants(n_polynomials, args.height)
    print(f"exact-resultant seconds={time.monotonic()-started:.3f}")
    if not args.skip_factor_table:
        print_resultant_summary(args.height, args.factor_bound, resultants)
    print_sequence_overlap(args.height, args.factor_bound, resultants)
    records = prime_statistics(
        args.height, args.prime_limit, n_polynomials, resultants
    )
    print_prime_summary(args.height, args.prime_limit, records)
    print(f"total seconds={time.monotonic()-started:.3f}")


if __name__ == "__main__":
    main()
