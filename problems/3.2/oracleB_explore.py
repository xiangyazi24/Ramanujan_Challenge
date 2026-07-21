#!/usr/bin/env python3
"""Finite exploration for the arithmetic-oracle audit.

The script deliberately keeps two different notions of zero separate:

* ``Z_p = {j : b_j == 0 (mod p)}`` is the coefficient/Mellin zero-set;
* roots of ``H_p(t) = sum_j b_j t^j`` are evaluation (fiber) zeros.

The ordinary Python phase recomputes every coefficient zero-set for
``p <= 20000`` from the division-free recurrence and performs the requested
distribution and fixed-anchor tests.  The optional Sage worker verifies the
polynomial identity

    H_p(t) = (t^2 - 34t + 1)^epsilon_p B_p(t)^2

exactly over ``F_p``.  In ``--hasse-mode all`` it extracts every linear/
quadratic factor-degree type using Frobenius gcds, but only asks Sage for the
explicit factor polynomials at a sparse sample: complete factorization at
every prime would answer no additional logical question and is much costlier.

All statistics are diagnostics, never inputs to the propositions in
``oracleB_result.tex``.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from dataclasses import dataclass
from hashlib import sha256
import json
from math import comb, exp, floor, isqrt, log, nextafter, sqrt
from pathlib import Path
import random
import shutil
import struct
import subprocess
import sys
import tempfile
from typing import Callable, Iterable, Sequence


DEFAULT_LIMIT = 20_000
DEFAULT_SEED = 32_032
SMALL_MODULI = (3, 5, 7, 11, 13)
FIXED_ANCHORS = (0, 1, 2, 3, 5, 10, 20, 50, 100, 1000)
CHECKPOINTS = (100, 500, 1000, 5000, 10_000, 20_000)
FACTOR_SAMPLES = (
    5,
    7,
    11,
    13,
    17,
    19,
    23,
    29,
    47,
    73,
    101,
    251,
    503,
    1009,
    2003,
    5003,
    10_007,
    19_997,
)
EXACT_SQUARE_CLASSES = {1, 5, 7, 11}
CORRECTED_SQUARE_CLASSES = {13, 17, 19, 23}


@dataclass(frozen=True)
class ZeroRecord:
    prime: int
    zeros: tuple[int, ...]

    @property
    def first_zero(self) -> int | None:
        return self.zeros[0] if self.zeros else None

    @property
    def center_hit(self) -> bool:
        return (self.prime - 1) // 2 in self.zeros

    @property
    def pair_count(self) -> int:
        return (len(self.zeros) - int(self.center_hit)) // 2

    @property
    def gap(self) -> int | None:
        if not self.zeros:
            return None
        return self.prime - 1 - 2 * self.zeros[0]


def primes_up_to(limit: int) -> list[int]:
    """Return all primes at most ``limit``."""

    if limit < 2:
        return []
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[:2] = b"\x00\x00"
    for divisor in range(2, isqrt(limit) + 1):
        if not sieve[divisor]:
            continue
        start = divisor * divisor
        sieve[start : limit + 1 : divisor] = b"\x00" * (
            (limit - start) // divisor + 1
        )
    return [number for number, flag in enumerate(sieve) if flag]


def apery_zero_set(prime: int) -> tuple[int, ...]:
    """Compute ``Z_p`` with ``A_n=(n!)^3 b_n`` and no modular division.

    For ``n < p``, multiplication by ``(n!)^3`` does not alter vanishing.
    The recurrence is

      A_{n+1}=P(n)A_n-n^6A_{n-1},
      P(n)=34n^3+51n^2+27n+5.
    """

    if prime < 3:
        raise ValueError("the recurrence audit is used only for odd primes")
    previous = 1 % prime
    current = 5 % prime
    zeros: list[int] = []
    if previous == 0:
        zeros.append(0)
    if current == 0:
        zeros.append(1)
    for index in range(1, prime - 1):
        coefficient = (
            34 * index**3 + 51 * index**2 + 27 * index + 5
        ) % prime
        following = (
            coefficient * current
            - pow(index, 6, prime) * previous
        ) % prime
        previous, current = current, following
        if current == 0:
            zeros.append(index + 1)
    return tuple(zeros)


def apery_integers(max_index: int) -> list[int]:
    """Return the exact positive Apéry integers through ``b_max_index``."""

    values = [1]
    if max_index == 0:
        return values
    values.append(5)
    for index in range(1, max_index):
        coefficient = 34 * index**3 + 51 * index**2 + 27 * index + 5
        numerator = coefficient * values[index] - index**3 * values[index - 1]
        denominator = (index + 1) ** 3
        if numerator % denominator:
            raise ArithmeticError("the integral Apéry recurrence did not divide")
        values.append(numerator // denominator)
    return values


def load_binary_zero_sets(path: Path, limit: int) -> tuple[dict[int, tuple[int, ...]], str]:
    """Read the repository's headerless little-endian ``(p,r)`` data."""

    raw = path.read_bytes()
    if len(raw) % 8:
        raise ArithmeticError("data_zp_pairs.bin has a partial record")
    grouped: dict[int, list[int]] = defaultdict(list)
    previous = (-1, -1)
    for prime, residue in struct.iter_unpack("<II", raw):
        if (prime, residue) <= previous:
            raise ArithmeticError("binary zero records are not strictly sorted")
        previous = (prime, residue)
        if prime <= limit:
            grouped[prime].append(residue)
    return (
        {prime: tuple(residues) for prime, residues in grouped.items()},
        sha256(raw).hexdigest(),
    )


def validate_zero_records(records: Sequence[ZeroRecord]) -> None:
    """Check the exact structural identities used by the exploration."""

    for record in records:
        prime = record.prime
        zeros = record.zeros
        if any(not 0 <= residue < prime for residue in zeros):
            raise ArithmeticError(f"out-of-range zero at p={prime}")
        if tuple(sorted(set(zeros))) != zeros:
            raise ArithmeticError(f"duplicate/unsorted zero-set at p={prime}")
        zero_set = set(zeros)
        if any(prime - 1 - residue not in zero_set for residue in zeros):
            raise ArithmeticError(f"reflection failed at p={prime}")
        if any((residue + 1) % prime in zero_set for residue in zeros):
            raise ArithmeticError(f"consecutive zeros at p={prime}")


def quantiles(values: Sequence[float], probabilities: Iterable[float]) -> dict[str, float]:
    """Linear-interpolation empirical quantiles."""

    ordered = sorted(values)
    if not ordered:
        return {}
    answer: dict[str, float] = {}
    for probability in probabilities:
        position = probability * (len(ordered) - 1)
        lower = floor(position)
        upper = min(lower + 1, len(ordered) - 1)
        weight = position - lower
        answer[f"{probability:.2f}"] = (
            (1.0 - weight) * ordered[lower] + weight * ordered[upper]
        )
    return answer


def uniform_ks(values: Sequence[float]) -> tuple[float, float]:
    """One-sample KS statistic and its standard continuous-null p-value.

    The p-value is only reported for the doublet diagnostic.  The exact null
    is a fine discrete grid depending on p, so the asymptotic p-value should
    not be read as an arithmetic theorem.
    """

    ordered = sorted(values)
    count = len(ordered)
    if count == 0:
        return 0.0, 1.0
    d_plus = max((index + 1) / count - value for index, value in enumerate(ordered))
    d_minus = max(value - index / count for index, value in enumerate(ordered))
    statistic = max(d_plus, d_minus)
    scaled = (sqrt(count) + 0.12 + 0.11 / sqrt(count)) * statistic
    tail = 0.0
    for order in range(1, 200):
        term = 2.0 * (-1.0 if order % 2 == 0 else 1.0) * exp(
            -2.0 * order * order * scaled * scaled
        )
        tail += term
        if abs(term) < 1e-14:
            break
    return statistic, min(1.0, max(0.0, tail))


def reflected_min_cdf(records: Sequence[ZeroRecord], argument: float) -> float:
    """CDF of ``2 min(Z_p)/(p-1)`` under the fixed-cardinality null.

    For each p, its actual number of noncentral reflection pairs is retained,
    but those pairs are sampled uniformly without replacement.  A center hit
    is also retained.  The returned CDF is the average of these non-identical
    per-prime null CDFs.
    """

    if not records:
        return 0.0
    total = 0.0
    for record in records:
        prime = record.prime
        pair_count = record.pair_count
        representatives = (prime - 3) // 2
        if pair_count == 0:
            # The only possible active record here is a center singleton.
            total += float(argument >= 1.0)
            continue
        threshold = floor(argument * (prime - 1) / 2 + 1e-14)
        threshold = max(0, min(representatives, threshold))
        if threshold == 0:
            continue
        if representatives - threshold < pair_count:
            total += 1.0
            continue
        survival = 1.0
        for offset in range(pair_count):
            survival *= (
                (representatives - threshold - offset)
                / (representatives - offset)
            )
        total += 1.0 - survival
    return total / len(records)


def ks_against_cdf(values: Sequence[float], cdf: Callable[[float], float]) -> float:
    """KS distance from a possibly discrete right-continuous CDF."""

    ordered = sorted(values)
    count = len(ordered)
    if count == 0:
        return 0.0
    distance = 0.0
    for index, value in enumerate(ordered):
        distance = max(
            distance,
            abs((index + 1) / count - cdf(value)),
            abs(index / count - cdf(nextafter(value, float("-inf")))),
        )
    return distance


def chi_square_statistic(row_labels: Sequence[int], column_labels: Sequence[int]) -> tuple[float, int, int]:
    """Pearson statistic for a categorical contingency table."""

    if len(row_labels) != len(column_labels):
        raise ValueError("contingency columns have different lengths")
    row_totals = Counter(row_labels)
    column_totals = Counter(column_labels)
    cells = Counter(zip(row_labels, column_labels))
    sample_size = len(row_labels)
    statistic = 0.0
    for row, row_total in row_totals.items():
        for column, column_total in column_totals.items():
            expected = row_total * column_total / sample_size
            observed = cells[row, column]
            statistic += (observed - expected) ** 2 / expected
    return statistic, len(row_totals), len(column_totals)


def correlation_rows(
    selected_records: Sequence[ZeroRecord], permutations: int, seed: int
) -> list[dict[str, float | int]]:
    """Return small-modulus permutation diagnostics for one record family."""

    rng = random.Random(seed)
    output: list[dict[str, float | int]] = []
    for modulus in SMALL_MODULI:
        selected = [record for record in selected_records if record.prime > modulus]
        if not selected:
            output.append(
                {
                    "q": modulus,
                    "sample_size": 0,
                    "chi_square": 0.0,
                    "degrees_of_freedom": 0,
                    "cramers_v": 0.0,
                    "permutation_p": 1.0,
                }
            )
            continue
        row_labels = [record.prime % modulus for record in selected]
        column_labels = [int(record.gap) % modulus for record in selected]
        observed, row_count, column_count = chi_square_statistic(
            row_labels, column_labels
        )
        exceedances = 0
        shuffled = list(column_labels)
        for _ in range(permutations):
            rng.shuffle(shuffled)
            statistic, _, _ = chi_square_statistic(row_labels, shuffled)
            exceedances += statistic >= observed - 1e-12
        denominator = min(row_count - 1, column_count - 1)
        cramers_v = (
            sqrt(observed / (len(selected) * denominator))
            if denominator > 0
            else 0.0
        )
        output.append(
            {
                "q": modulus,
                "sample_size": len(selected),
                "chi_square": observed,
                "degrees_of_freedom": (row_count - 1) * (column_count - 1),
                "cramers_v": cramers_v,
                "permutation_p": (exceedances + 1) / (permutations + 1),
            }
        )
    return output


def modular_correlations(
    records: Sequence[ZeroRecord], permutations: int, seed: int
) -> dict[str, list[dict[str, float | int]]]:
    """Test outer-gap correlations for all active primes and for doublets.

    The all-active table fulfills the literal computational duty.  The
    doublet table is also reported because only there does the first zero
    determine the entire coefficient zero-set.  Odd q avoid the deterministic
    parity of every reflected outer gap.
    """

    active = [record for record in records if record.zeros]
    doublets = [record for record in active if len(record.zeros) == 2]
    return {
        "all_active_outer_gap": correlation_rows(active, permutations, seed),
        "doublets": correlation_rows(doublets, permutations, seed + 1),
    }


def fixed_anchor_audit(
    records: Sequence[ZeroRecord], limit: int
) -> list[dict[str, object]]:
    """Count finite-range hits and check the Lucas fixed-anchor inclusion."""

    record_by_prime = {record.prime: set(record.zeros) for record in records}
    exact_values = apery_integers(max(FIXED_ANCHORS))
    output: list[dict[str, object]] = []
    for anchor in FIXED_ANCHORS:
        hits = [
            prime
            for prime, zeros in record_by_prime.items()
            if anchor % prime in zeros
        ]
        divisors_in_range = [
            prime
            for prime in record_by_prime
            if exact_values[anchor] % prime == 0
        ]
        if any(prime not in divisors_in_range for prime in hits):
            raise ArithmeticError(f"Lucas fixed-anchor inclusion failed at c={anchor}")
        cumulative = {
            str(cutoff): sum(prime <= min(cutoff, limit) for prime in hits)
            for cutoff in CHECKPOINTS
            if cutoff <= limit
        }
        output.append(
            {
                "anchor": anchor,
                "hits_up_to_limit": hits,
                "hit_count": len(hits),
                "prime_divisors_of_b_anchor_up_to_limit": divisors_in_range,
                "cumulative_counts": cumulative,
                "decimal_digits_of_b_anchor": len(str(exact_values[anchor])),
            }
        )
    return output


def fixed_anchor_growth(
    anchor_rows: Sequence[dict[str, object]], limit: int
) -> list[dict[str, float | int]]:
    """Compare the largest displayed fixed-anchor count with ``pi(X)``."""

    all_primes = primes_up_to(limit)
    output: list[dict[str, float | int]] = []
    for cutoff in CHECKPOINTS:
        if cutoff > limit:
            continue
        prime_count = sum(prime <= cutoff for prime in all_primes)
        maximum = max(
            int(row["cumulative_counts"][str(cutoff)])  # type: ignore[index]
            for row in anchor_rows
        )
        output.append(
            {
                "X": cutoff,
                "pi_X": prime_count,
                "maximum_displayed_anchor_count": maximum,
                "ratio": maximum / prime_count,
            }
        )
    return output


def factor_histogram(factorization: object) -> dict[str, int]:
    """Group Sage factors by irreducible degree."""

    histogram: Counter[int] = Counter()
    for factor, multiplicity in factorization:  # type: ignore[union-attr]
        if int(multiplicity) != 1:
            raise ArithmeticError("the normalized B_p sample is not squarefree")
        histogram[int(factor.degree())] += 1
    return {str(degree): count for degree, count in sorted(histogram.items())}


def sage_hasse_worker(limit: int, mode: str, factor_samples: Sequence[int]) -> None:
    """Run exact Hasse-polynomial checks inside ``sage -python``."""

    try:
        from sage.all import GF, PolynomialRing  # type: ignore[import-not-found]
    except ImportError as error:
        raise SystemExit("--sage-worker must be run by `sage -python`") from error

    primes = [prime for prime in primes_up_to(limit) if prime >= 5]
    if mode == "samples":
        checked_primes = [prime for prime in primes if prime in factor_samples]
    elif mode == "all":
        checked_primes = primes
    else:
        raise SystemExit("the Sage worker accepts only samples/all")

    factor_sample_set = set(factor_samples)
    samples: list[dict[str, object]] = []
    factorization_types: list[dict[str, int]] = []
    reciprocity_signs: Counter[str] = Counter()
    for position, prime in enumerate(checked_primes, start=1):
        inverses = [0] * prime
        inverses[1] = 1
        for value in range(2, prime):
            inverses[value] = (
                prime - (prime // value) * inverses[prime % value] % prime
            )
        coefficients = [0] * prime
        coefficients[0] = 1
        coefficients[1] = 5 % prime
        for index in range(1, prime - 1):
            coefficient = (
                34 * index**3 + 51 * index**2 + 27 * index + 5
            ) % prime
            numerator = (
                coefficient * coefficients[index]
                - (index**3 % prime) * coefficients[index - 1]
            ) % prime
            inverse_cube = pow(inverses[index + 1], 3, prime)
            coefficients[index + 1] = numerator * inverse_cube % prime
        if coefficients[-1] != 1:
            raise ArithmeticError(f"H_p is not monic at p={prime}")

        field = GF(prime)
        polynomial_ring = PolynomialRing(field, "t")
        t = polynomial_ring.gen()
        hasse = polynomial_ring(coefficients)
        delta = t**2 - 34 * t + 1
        residue_class = prime % 24
        if residue_class in EXACT_SQUARE_CLASSES:
            epsilon = 0
        elif residue_class in CORRECTED_SQUARE_CLASSES:
            epsilon = 1
        else:
            raise ArithmeticError(f"unexpected good-prime class p={prime}")
        quotient, remainder = hasse.quo_rem(delta**epsilon)
        if remainder:
            raise ArithmeticError(f"Delta correction does not divide H_{prime}")
        is_square, square_root = quotient.is_square(root=True)
        if not is_square or square_root is None:
            raise ArithmeticError(f"corrected H_{prime} is not a square")
        if square_root[0] != 1:
            square_root = -square_root
        expected_degree = (prime - 1 - 2 * epsilon) // 2
        if int(square_root.degree()) != expected_degree:
            raise ArithmeticError(f"wrong B_p degree at p={prime}")
        if hasse != delta**epsilon * square_root**2:
            raise ArithmeticError(f"square identity failed at p={prime}")
        if square_root.gcd(square_root.derivative()).degree() != 0:
            raise ArithmeticError(f"B_p is not squarefree at p={prime}")
        if square_root.gcd(delta).degree() != 0:
            raise ArithmeticError(f"B_p is not coprime to Delta at p={prime}")

        # A squarefree polynomial divides t^(p^2)-t exactly when every
        # irreducible factor has degree one or two.  This extracts the full
        # degree pattern without performing thousands of costly complete
        # factorizations.  The gcd with t^p-t counts the linear factors.
        frobenius_once = pow(t, prime, square_root)
        linear_factors = int(square_root.gcd(frobenius_once - t).degree())
        frobenius_twice = pow(frobenius_once, prime, square_root)
        if frobenius_twice != t:
            raise ArithmeticError(
                f"B_p has an irreducible factor of degree greater than two at p={prime}"
            )
        if (expected_degree - linear_factors) % 2:
            raise ArithmeticError(f"factor-degree parity failed at p={prime}")
        quadratic_factors = (expected_degree - linear_factors) // 2
        factorization_types.append(
            {
                "p": prime,
                "degree_B": expected_degree,
                "linear_factors": linear_factors,
                "quadratic_factors": quadratic_factors,
            }
        )

        reciprocal = polynomial_ring(list(square_root)[::-1])
        if reciprocal == square_root:
            reciprocity_sign = "+"
        elif reciprocal == -square_root:
            reciprocity_sign = "-"
        else:
            raise ArithmeticError(f"B_p has no reciprocal sign at p={prime}")
        reciprocity_signs[f"{residue_class}:{reciprocity_sign}"] += 1

        if prime in factor_sample_set:
            factorization = square_root.factor()
            histogram = factor_histogram(factorization)
            expected_histogram = {
                str(degree): count
                for degree, count in ((1, linear_factors), (2, quadratic_factors))
                if count
            }
            if histogram != expected_histogram:
                raise ArithmeticError(f"factor-pattern cross-check failed at p={prime}")
            component_count = sum(histogram.values())
            samples.append(
                {
                    "p": prime,
                    "epsilon": epsilon,
                    "degree_B": expected_degree,
                    "reciprocity_sign": reciprocity_sign,
                    "irreducible_degree_histogram": histogram,
                    "irreducible_factor_count": component_count,
                    "geometric_point_count": expected_degree,
                    "hyperelliptic_genus_if_y2_eq_B": (expected_degree - 1) // 2,
                }
            )
        if position % 250 == 0 or position == len(checked_primes):
            print(
                f"Sage Hasse audit {position}/{len(checked_primes)} (p={prime})",
                file=sys.stderr,
                flush=True,
            )

    print(
        json.dumps(
            {
                "mode": mode,
                "checked_prime_count": len(checked_primes),
                "largest_checked_prime": checked_primes[-1] if checked_primes else None,
                "identity_failures": [],
                "all_irreducible_degrees_at_most_two": True,
                "reciprocity_sign_counts_by_p_mod_24": dict(reciprocity_signs),
                "factorization_type_records": factorization_types,
                "factor_samples": samples,
            },
            sort_keys=True,
        )
    )


def run_sage_hasse_audit(
    limit: int, mode: str, factor_samples: Sequence[int]
) -> dict[str, object] | None:
    """Dispatch the exact polynomial phase to the installed Sage runtime."""

    if mode == "skip":
        return None
    sage = shutil.which("sage")
    if sage is None:
        raise RuntimeError("Sage is required unless --hasse-mode skip is used")
    cache = Path(tempfile.gettempdir()) / "oracleB_sage_cache"
    cache.mkdir(parents=True, exist_ok=True)
    environment = dict(__import__("os").environ)
    environment["DOT_SAGE"] = str(cache)
    command = [
        sage,
        "-python",
        str(Path(__file__).resolve()),
        "--sage-worker",
        "--limit",
        str(limit),
        "--hasse-mode",
        mode,
        "--factor-primes",
        ",".join(map(str, factor_samples)),
    ]
    completed = subprocess.run(
        command,
        check=True,
        text=True,
        capture_output=True,
        env=environment,
    )
    if completed.stderr:
        print(completed.stderr, file=sys.stderr, end="")
    return json.loads(completed.stdout)


def zero_set_summary(
    records: Sequence[ZeroRecord], binary_digest: str | None
) -> dict[str, object]:
    """Compute requested coefficient-zero distribution diagnostics."""

    active = [record for record in records if record.zeros]
    doublets = [record for record in active if len(record.zeros) == 2]
    scaled_first = [
        2.0 * int(record.first_zero) / (record.prime - 1) for record in active
    ]
    scaled_gap = [int(record.gap) / (record.prime - 1) for record in active]
    raw_first = [int(record.first_zero) / record.prime for record in active]
    raw_gap = [int(record.gap) / record.prime for record in active]
    doublet_first = [
        2.0 * int(record.first_zero) / (record.prime - 1) for record in doublets
    ]
    doublet_gap = [int(record.gap) / (record.prime - 1) for record in doublets]

    first_mixture_ks = ks_against_cdf(
        scaled_first, lambda value: reflected_min_cdf(active, value)
    )
    gap_mixture_ks = ks_against_cdf(
        scaled_gap,
        lambda value: 1.0
        - reflected_min_cdf(active, nextafter(1.0 - value, float("-inf"))),
    )
    doublet_first_ks, doublet_first_p = uniform_ks(doublet_first)
    doublet_gap_ks, doublet_gap_p = uniform_ks(doublet_gap)

    return {
        "prime_count": len(records),
        "active_prime_count": len(active),
        "zero_count_histogram": {
            str(size): count
            for size, count in sorted(Counter(map(lambda row: len(row.zeros), records)).items())
        },
        "odd_zero_set_primes": [record.prime for record in records if len(record.zeros) % 2],
        "max_zero_count": max(map(lambda row: len(row.zeros), records)),
        "multi_pair_example_Z_181": next(
            (list(record.zeros) for record in records if record.prime == 181), None
        ),
        "r_over_p_quantiles": quantiles(raw_first, (0, 0.1, 0.25, 0.5, 0.75, 0.9, 1)),
        "h_over_p_quantiles": quantiles(raw_gap, (0, 0.1, 0.25, 0.5, 0.75, 0.9, 1)),
        "fixed_cardinality_reflected_null": {
            "scaled_first_zero_KS": first_mixture_ks,
            "scaled_gap_KS": gap_mixture_ks,
            "note": "conditioned on each prime's actual pair count and center status",
        },
        "doublet_uniform_diagnostic": {
            "sample_size": len(doublets),
            "scaled_first_zero_KS": doublet_first_ks,
            "scaled_first_zero_asymptotic_p": doublet_first_p,
            "scaled_gap_KS": doublet_gap_ks,
            "scaled_gap_asymptotic_p": doublet_gap_p,
        },
        "binary_data_sha256": binary_digest,
    }


def render_text(result: dict[str, object]) -> str:
    """Render a compact human-readable version of the JSON certificate."""

    zero = result["zero_sets"]  # type: ignore[index]
    lines = [
        "Oracle B finite exploration",
        f"limit={result['limit']} seed={result['seed']}",
        f"primes={zero['prime_count']} active={zero['active_prime_count']}",  # type: ignore[index]
        f"Z histogram={zero['zero_count_histogram']}",  # type: ignore[index]
        f"odd Z primes={zero['odd_zero_set_primes']}",  # type: ignore[index]
        f"r/p quantiles={zero['r_over_p_quantiles']}",  # type: ignore[index]
        f"h/p quantiles={zero['h_over_p_quantiles']}",  # type: ignore[index]
        f"reflected-null KS={zero['fixed_cardinality_reflected_null']}",  # type: ignore[index]
        f"doublet KS={zero['doublet_uniform_diagnostic']}",  # type: ignore[index]
        "",
        "Fixed anchors:",
    ]
    for row in result["fixed_anchors"]:  # type: ignore[index]
        lines.append(
            f"  c={row['anchor']:4d} hits={row['hits_up_to_limit']} "
            f"digits(b_c)={row['decimal_digits_of_b_anchor']}"
        )
    lines.append(f"fixed-anchor growth={result['fixed_anchor_growth']}")
    lines.append("")
    for family, rows in result["modular_correlations"].items():  # type: ignore[index]
        lines.append(f"Gap correlations ({family}):")
        for row in rows:
            lines.append(
                f"  q={row['q']:2d} n={row['sample_size']:4d} "
                f"V={row['cramers_v']:.4f} perm-p={row['permutation_p']:.4f}"
            )
    hasse = result.get("hasse")
    if hasse:
        lines.extend(
            [
                "",
                f"Hasse audit: mode={hasse['mode']} "
                f"checked={hasse['checked_prime_count']} "
                f"largest={hasse['largest_checked_prime']}",
                f"reciprocity signs={hasse['reciprocity_sign_counts_by_p_mod_24']}",
                "all irreducible factor degrees are at most two: "
                f"{hasse['all_irreducible_degrees_at_most_two']}",
                "factor samples:",
            ]
        )
        for row in hasse["factor_samples"]:
            lines.append(
                f"  p={row['p']:5d} degB={row['degree_B']:5d} "
                f"sign={row['reciprocity_sign']} "
                f"factors={row['irreducible_degree_histogram']}"
            )
    return "\n".join(lines) + "\n"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--limit", type=int, default=DEFAULT_LIMIT)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--permutations", type=int, default=999)
    parser.add_argument(
        "--hasse-mode", choices=("skip", "samples", "all"), default="samples"
    )
    parser.add_argument("--json-output", type=Path)
    parser.add_argument("--sage-worker", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--factor-primes", default=",".join(map(str, FACTOR_SAMPLES)))
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    factor_samples = tuple(
        int(token) for token in args.factor_primes.split(",") if token
    )
    if args.sage_worker:
        sage_hasse_worker(args.limit, args.hasse_mode, factor_samples)
        return
    if args.limit < 5:
        raise SystemExit("--limit must be at least 5")
    if args.permutations < 0:
        raise SystemExit("--permutations must be nonnegative")

    primes = [prime for prime in primes_up_to(args.limit) if prime >= 5]
    records = [ZeroRecord(prime, apery_zero_set(prime)) for prime in primes]
    validate_zero_records(records)

    data_path = Path(__file__).with_name("data_zp_pairs.bin")
    binary_digest: str | None = None
    if data_path.exists():
        binary, binary_digest = load_binary_zero_sets(data_path, args.limit)
        for record in records:
            expected = binary.get(record.prime, ())
            if record.prime == 5 and not expected:
                expected = (1, 3)  # the historical binary starts at p=7
            if record.zeros != expected:
                raise ArithmeticError(
                    f"binary/recomputed zero-set mismatch at p={record.prime}"
                )

    anchor_rows = fixed_anchor_audit(records, args.limit)
    result: dict[str, object] = {
        "limit": args.limit,
        "seed": args.seed,
        "permutations": args.permutations,
        "zero_sets": zero_set_summary(records, binary_digest),
        "fixed_anchors": anchor_rows,
        "fixed_anchor_growth": fixed_anchor_growth(anchor_rows, args.limit),
        "modular_correlations": modular_correlations(
            records, args.permutations, args.seed
        ),
        "hasse": run_sage_hasse_audit(
            args.limit, args.hasse_mode, factor_samples
        ),
    }
    serialized = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.json_output:
        args.json_output.write_text(serialized, encoding="utf-8")
    print(render_text(result), end="")


if __name__ == "__main__":
    main()
