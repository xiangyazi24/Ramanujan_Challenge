#!/usr/bin/env python3
"""E1: zero-position forensics for the zeta(3) Apery numbers modulo p.

The main experiment uses every prime 5 <= p <= 5000, records every index
0 <= j < p for which b_j == 0 (mod p), and studies x = j/p.  The exceptional
high-Z prime 159977 is analysed separately and is not mixed into the aggregate
statistics.

Only the Python standard library is required.  By default the complete report
is written to /tmp/e1_zero_forensics_results.txt.
"""

from __future__ import annotations

import argparse
import bisect
import math
import random
import statistics
import time
from collections import Counter
from dataclasses import dataclass
from datetime import datetime, timezone
from fractions import Fraction
from itertools import combinations
from pathlib import Path
from typing import Iterable, Sequence


DEFAULT_OUTPUT = Path("/tmp/e1_zero_forensics_results.txt")
DEFAULT_SPECIAL_PRIME = 159_977
DEFAULT_SEED = 20260715
DEFAULT_SIMULATIONS = 10_000
HISTOGRAM_BINS = 20
HEIGHT_CUTOFFS = (5, 10, 20, 50, 100, 200, 500, 1_000)


@dataclass(frozen=True)
class PrimeZeros:
    p: int
    zeros: tuple[int, ...]

    @property
    def z(self) -> int:
        return len(self.zeros)

    @property
    def center(self) -> int:
        return (self.p - 1) // 2

    @property
    def has_center_zero(self) -> bool:
        return self.center in self.zeros

    @property
    def pair_count(self) -> int:
        return (self.z - int(self.has_center_zero)) // 2


@dataclass(frozen=True)
class RationalScanResult:
    a: int
    d: int
    observed: int
    expected: float
    variance: float
    z_score: float
    raw_p: float
    bonferroni_p: float


def sieve_primes(limit: int) -> list[int]:
    """Return all primes <= limit."""
    if limit < 2:
        return []
    flags = bytearray(b"\x01") * (limit + 1)
    flags[0:2] = b"\x00\x00"
    for q in range(2, math.isqrt(limit) + 1):
        if flags[q]:
            flags[q * q : limit + 1 : q] = b"\x00" * (
                (limit - q * q) // q + 1
            )
    return [q for q in range(2, limit + 1) if flags[q]]


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    for q in range(3, math.isqrt(n) + 1, 2):
        if n % q == 0:
            return False
    return True


def apery_zero_positions(p: int) -> tuple[int, ...]:
    """Return all j in [0,p) such that b_j == 0 (mod p).

    The recurrence is used only for n=1,...,p-2.  Thus every denominator
    (n+1)^3 is invertible modulo p; no p-adic division or fallback is needed.
    An O(p) inverse table avoids a modular exponentiation at every step.
    """
    if p < 5 or not is_prime(p):
        raise ValueError(f"p must be a prime >= 5, got {p}")

    inverses = [0] * p
    inverses[1] = 1
    for n in range(2, p):
        inverses[n] = (p - (p // n) * inverses[p % n] % p) % p

    b_prev = 1 % p
    b_curr = 5 % p
    zeros: list[int] = []
    if b_prev == 0:
        zeros.append(0)
    if b_curr == 0:
        zeros.append(1)

    for n in range(1, p - 1):
        n2 = n * n
        n3 = n2 * n % p
        coeff = (((34 * n + 51) * n + 27) * n + 5) % p
        inv = inverses[n + 1]
        inv3 = inv * inv % p * inv % p
        b_next = (coeff * b_curr - n3 * b_prev) % p * inv3 % p
        b_prev, b_curr = b_curr, b_next
        if b_curr == 0:
            zeros.append(n + 1)

    return tuple(zeros)


def apery_direct_mod(j: int, p: int) -> int:
    """Slow binomial-sum definition, used only for independent small checks."""
    return sum(
        math.comb(j, k) ** 2 * math.comb(j + k, k) ** 2
        for k in range(j + 1)
    ) % p


def validate_records(records: Sequence[PrimeZeros]) -> dict[str, object]:
    """Check recurrence output against exact small cases and known symmetries."""
    failures: list[str] = []
    direct_checks = 0
    for record in records:
        p = record.p
        zero_set = set(record.zeros)
        if 0 in zero_set or p - 1 in zero_set:
            failures.append(f"endpoint zero at p={p}")
        if any((p - 1 - j) not in zero_set for j in zero_set):
            failures.append(f"reflection failed at p={p}")
        if any(j + 1 in zero_set for j in zero_set):
            failures.append(f"adjacent zeros at p={p}")
        if (record.z - int(record.has_center_zero)) % 2:
            failures.append(f"unpaired zero count at p={p}")
        if p <= 43:
            direct = tuple(j for j in range(p) if apery_direct_mod(j, p) == 0)
            direct_checks += p
            if direct != record.zeros:
                failures.append(
                    f"binomial check failed at p={p}: recurrence={record.zeros}, direct={direct}"
                )
    return {
        "ok": not failures,
        "failures": failures,
        "direct_checks": direct_checks,
        "record_count": len(records),
    }


def continued_fraction_data(numerator: int, denominator: int) -> tuple[list[int], list[tuple[int, int]], int]:
    """Return CF terms, convergents, and H3 for numerator/denominator.

    H3 is the maximum denominator among the first three standard convergents,
    counting the initial 0/1 convergent for a number in (0,1).  If the finite
    continued fraction ends sooner, all available convergents are used.
    """
    n, d = numerator, denominator
    terms: list[int] = []
    while d:
        a = n // d
        terms.append(a)
        n, d = d, n - a * d

    p_nm2, p_nm1 = 0, 1
    q_nm2, q_nm1 = 1, 0
    convergents: list[tuple[int, int]] = []
    for a in terms:
        p_n = a * p_nm1 + p_nm2
        q_n = a * q_nm1 + q_nm2
        convergents.append((p_n, q_n))
        p_nm2, p_nm1 = p_nm1, p_n
        q_nm2, q_nm1 = q_nm1, q_n

    height = max(q for _, q in convergents[:3])
    return terms, convergents, height


def continued_fraction_height3(numerator: int, denominator: int) -> int:
    """Fast H3-only version of continued_fraction_data."""
    n, d = numerator, denominator
    p_nm2, p_nm1 = 0, 1
    q_nm2, q_nm1 = 1, 0
    height = 1
    for _ in range(3):
        if not d:
            break
        a = n // d
        n, d = d, n - a * d
        p_n = a * p_nm1 + p_nm2
        q_n = a * q_nm1 + q_nm2
        height = max(height, q_n)
        p_nm2, p_nm1 = p_nm1, p_n
        q_nm2, q_nm1 = q_nm1, q_n
    return height


def low_denominator_fractions(max_denominator: int = 20) -> list[tuple[int, int]]:
    fractions = {
        (a, d)
        for d in range(2, max_denominator + 1)
        for a in range(1, d)
        if math.gcd(a, d) == 1
    }
    fractions.update({(0, 1), (1, 1)})
    return sorted(fractions, key=lambda item: Fraction(item[0], item[1]))


LOW_FRACTIONS = low_denominator_fractions()
LOW_FRACTION_VALUES = [a / d for a, d in LOW_FRACTIONS]


def nearest_low_fraction(j: int, p: int) -> tuple[int, int, float]:
    """Nearest reduced a/d in [0,1], d<=20, and absolute distance."""
    x = j / p
    insertion = bisect.bisect_left(LOW_FRACTION_VALUES, x)
    candidate_indices = range(
        max(0, insertion - 2), min(len(LOW_FRACTIONS), insertion + 3)
    )
    best_numerator: int | None = None
    best_a = 0
    best_d = 1
    for index in candidate_indices:
        a, d = LOW_FRACTIONS[index]
        numerator = abs(d * j - a * p)
        if (
            best_numerator is None
            or numerator * best_d < best_numerator * d
            or (
                numerator * best_d == best_numerator * d
                and (d, a) < (best_d, best_a)
            )
        ):
            best_numerator = numerator
            best_a = a
            best_d = d
    assert best_numerator is not None
    return best_a, best_d, best_numerator / (p * best_d)


def flatten_points(records: Sequence[PrimeZeros]) -> list[tuple[int, int]]:
    return [(record.p, j) for record in records for j in record.zeros]


def histogram(points: Sequence[tuple[int, int]], bins: int = HISTOGRAM_BINS) -> list[int]:
    counts = [0] * bins
    for p, j in points:
        index = min(bins - 1, (j * bins) // p)
        counts[index] += 1
    return counts


def uniform_histogram(values: Iterable[float], bins: int = HISTOGRAM_BINS) -> list[int]:
    counts = [0] * bins
    for value in values:
        index = min(bins - 1, int(value * bins))
        counts[index] += 1
    return counts


def ks_from_values(values: Iterable[float]) -> float:
    values = sorted(values)
    n = len(values)
    if not n:
        return float("nan")
    d_plus = max((i + 1) / n - x for i, x in enumerate(values))
    d_minus = max(x - i / n for i, x in enumerate(values))
    return max(d_plus, d_minus)


def cramer_von_mises_from_values(values: Iterable[float]) -> float:
    values = sorted(values)
    n = len(values)
    if not n:
        return float("nan")
    return 1 / (12 * n) + sum(
        (x - (2 * i - 1) / (2 * n)) ** 2
        for i, x in enumerate(values, start=1)
    )


def ks_statistic(points: Sequence[tuple[int, int]]) -> float:
    return ks_from_values(j / p for p, j in points)


def cramer_von_mises_statistic(points: Sequence[tuple[int, int]]) -> float:
    return cramer_von_mises_from_values(j / p for p, j in points)


def independent_orbit_values(points: Sequence[tuple[int, int]]) -> list[float]:
    """Map one lower-half representative per reflected pair to [0,1].

    For j<(p-1)/2 use u=2j/(p-1).  Center zeros are excluded and reported
    separately because their orbit has size one.
    """
    return [
        2 * j / (p - 1)
        for p, j in points
        if j < (p - 1) // 2
    ]


def naive_ks_pvalue(d_stat: float, n: int) -> float:
    """Usual continuous-iid KS approximation (reported only as a naive check)."""
    if n <= 0 or not math.isfinite(d_stat):
        return float("nan")
    root_n = math.sqrt(n)
    lam = (root_n + 0.12 + 0.11 / root_n) * d_stat
    total = 0.0
    for k in range(1, 100):
        term = 2 * (-1) ** (k - 1) * math.exp(-2 * k * k * lam * lam)
        total += term
        if abs(term) < 1e-14:
            break
    return min(1.0, max(0.0, total))


def quantile(values: Sequence[float], q: float) -> float:
    if not values:
        return float("nan")
    ordered = sorted(values)
    position = q * (len(ordered) - 1)
    low = math.floor(position)
    high = math.ceil(position)
    if low == high:
        return ordered[low]
    weight = position - low
    return ordered[low] * (1 - weight) + ordered[high] * weight


def bucket_index(value: int, cutoffs: Sequence[int]) -> int:
    return bisect.bisect_left(cutoffs, value)


def height_histogram(points: Sequence[tuple[int, int]]) -> list[int]:
    counts = [0] * (len(HEIGHT_CUTOFFS) + 1)
    for p, j in points:
        height = continued_fraction_height3(j, p)
        counts[bucket_index(height, HEIGHT_CUTOFFS)] += 1
    return counts


def nearest_fraction_metrics(points: Sequence[tuple[int, int]]) -> tuple[list[float], list[float]]:
    absolute: list[float] = []
    scaled: list[float] = []
    for p, j in points:
        _, _, distance = nearest_low_fraction(j, p)
        absolute.append(distance)
        scaled.append(p * distance)
    return absolute, scaled


def sample_reflected_null(records: Sequence[PrimeZeros], rng: random.Random) -> list[tuple[int, int]]:
    """Uniformly sample reflected pairs, conditional on each p's observed Z.

    Endpoints are excluded.  Whether the fixed center is present is also held
    fixed, which handles the p=11 and p=3137 exceptional center zeros without
    treating their deterministic location as evidence of clustering.  Rejection
    sampling additionally enforces the exact no-adjacent-zero constraint.
    """
    sample: list[tuple[int, int]] = []
    for record in records:
        p = record.p
        center = record.center
        fixed = [center] if record.has_center_zero else []
        while True:
            representatives = (
                rng.sample(range(1, center), record.pair_count)
                if record.pair_count
                else []
            )
            candidate = fixed + representatives + [
                p - 1 - j for j in representatives
            ]
            ordered = sorted(candidate)
            if all(b - a > 1 for a, b in zip(ordered, ordered[1:])):
                sample.extend((p, j) for j in candidate)
                break
    return sample


def vector_mean(vectors: Sequence[Sequence[float]]) -> list[float]:
    if not vectors:
        return []
    return [statistics.fmean(column) for column in zip(*vectors)]


def vector_sd(vectors: Sequence[Sequence[float]], means: Sequence[float]) -> list[float]:
    if not vectors:
        return []
    n = len(vectors)
    return [
        math.sqrt(sum((row[i] - means[i]) ** 2 for row in vectors) / n)
        for i in range(len(means))
    ]


def pearson_distance(observed: Sequence[float], expected: Sequence[float]) -> float:
    return sum((o - e) ** 2 / e for o, e in zip(observed, expected) if e > 0)


def empirical_upper_p(observed: float, simulated: Sequence[float]) -> float:
    return (1 + sum(value >= observed for value in simulated)) / (len(simulated) + 1)


def run_null_simulations(
    records: Sequence[PrimeZeros], simulations: int, seed: int
) -> dict[str, object]:
    rng = random.Random(seed)
    histograms: list[list[int]] = []
    heights: list[list[int]] = []
    ks_values: list[float] = []
    cvm_values: list[float] = []
    orbit_histograms: list[list[int]] = []
    orbit_ks_values: list[float] = []
    orbit_cvm_values: list[float] = []
    nearest_counts: dict[str, list[int]] = {
        "scaled<=0.5": [],
        "scaled<=1": [],
        "scaled<=2": [],
        "absolute<=0.001": [],
        "absolute<=0.005": [],
    }

    for _ in range(simulations):
        sample = sample_reflected_null(records, rng)
        histograms.append(histogram(sample))
        heights.append(height_histogram(sample))
        ks_values.append(ks_statistic(sample))
        cvm_values.append(cramer_von_mises_statistic(sample))
        orbit_values = independent_orbit_values(sample)
        orbit_histograms.append(uniform_histogram(orbit_values))
        orbit_ks_values.append(ks_from_values(orbit_values))
        orbit_cvm_values.append(cramer_von_mises_from_values(orbit_values))
        absolute, scaled = nearest_fraction_metrics(sample)
        nearest_counts["scaled<=0.5"].append(sum(x <= 0.5 + 1e-12 for x in scaled))
        nearest_counts["scaled<=1"].append(sum(x <= 1 + 1e-12 for x in scaled))
        nearest_counts["scaled<=2"].append(sum(x <= 2 + 1e-12 for x in scaled))
        nearest_counts["absolute<=0.001"].append(sum(x <= 0.001 + 1e-15 for x in absolute))
        nearest_counts["absolute<=0.005"].append(sum(x <= 0.005 + 1e-15 for x in absolute))

    # Estimate expected bin counts on an independent pilot subset, then use
    # only held-out replicates for Pearson-distance rank calibration.  This
    # avoids letting a simulated row help define its own reference center.
    pilot_count = max(20, simulations // 5)
    pilot_count = min(pilot_count, simulations - 1)
    hist_pilot = histograms[:pilot_count]
    height_pilot = heights[:pilot_count]
    orbit_hist_pilot = orbit_histograms[:pilot_count]
    hist_mean = vector_mean(hist_pilot)
    height_mean = vector_mean(height_pilot)
    orbit_hist_mean = vector_mean(orbit_hist_pilot)
    return {
        "histograms": histograms,
        "hist_mean": hist_mean,
        "hist_sd": vector_sd(hist_pilot, hist_mean),
        "hist_distances": [
            pearson_distance(row, hist_mean) for row in histograms[pilot_count:]
        ],
        "height_histograms": heights,
        "height_mean": height_mean,
        "height_sd": vector_sd(height_pilot, height_mean),
        "height_distances": [
            pearson_distance(row, height_mean) for row in heights[pilot_count:]
        ],
        "ks": ks_values,
        "cvm": cvm_values,
        "orbit_histograms": orbit_histograms,
        "orbit_hist_mean": orbit_hist_mean,
        "orbit_hist_sd": vector_sd(orbit_hist_pilot, orbit_hist_mean),
        "orbit_hist_distances": [
            pearson_distance(row, orbit_hist_mean)
            for row in orbit_histograms[pilot_count:]
        ],
        "orbit_ks": orbit_ks_values,
        "orbit_cvm": orbit_cvm_values,
        "nearest_counts": nearest_counts,
        "pearson_pilot_count": pilot_count,
        "pearson_heldout_count": simulations - pilot_count,
    }


def ceil_div(a: int, b: int) -> int:
    return -((-a) // b)


def rational_window_bounds(p: int, a: int, d: int, mode: str) -> tuple[int, int]:
    """Integer j bounds for |j/p-a/d| <= 1/p or <= 1/200."""
    if mode == "grid":
        low = ceil_div(a * p - d, d)
        high = (a * p + d) // d
    elif mode == "fixed":
        low = ceil_div(200 * a * p - d * p, 200 * d)
        high = (200 * a * p + d * p) // (200 * d)
    else:
        raise ValueError(f"unknown rational-window mode: {mode}")
    return max(1, low), min(p - 2, high)


def in_rational_window(j: int, p: int, a: int, d: int, mode: str) -> bool:
    difference = abs(d * j - a * p)
    if mode == "grid":
        return difference <= d
    if mode == "fixed":
        return 200 * difference <= d * p
    raise ValueError(f"unknown rational-window mode: {mode}")


def scan_rational_clusters(
    records: Sequence[PrimeZeros], mode: str
) -> list[RationalScanResult]:
    """Scan every reduced a/d, d<=20, against a conditional uniform null.

    The expectation and finite-population variance preserve every observed
    pair count and center-zero status.  They do not assume independent members
    inside a reflected pair.
    """
    preliminary: list[tuple[int, int, int, float, float, float, float]] = []
    for a, d in LOW_FRACTIONS:
        observed = 0
        expected = 0.0
        variance = 0.0
        fixed_hits = 0
        total_distribution = [1.0]
        for record in records:
            p = record.p
            observed += sum(
                in_rational_window(j, p, a, d, mode) for j in record.zeros
            )

            center = record.center
            low, high = rational_window_bounds(p, a, d, mode)
            pair_contributions: Counter[int] = Counter()
            fixed_center = 0
            if low <= high:
                for j in range(low, high + 1):
                    if j == center:
                        if record.has_center_zero:
                            fixed_center = 1
                    else:
                        pair_contributions[min(j, p - 1 - j)] += 1

            population_size = center - 1
            draws = record.pair_count
            expected += fixed_center
            fixed_hits += fixed_center
            if population_size <= 0 or draws == 0:
                continue

            first_moment_sum = sum(pair_contributions.values())
            second_moment_sum = sum(c * c for c in pair_contributions.values())
            population_mean = first_moment_sum / population_size
            population_variance = (
                second_moment_sum / population_size - population_mean**2
            )
            expected += draws * population_mean
            if population_size > 1:
                variance += (
                    draws
                    * (population_size - draws)
                    / (population_size - 1)
                    * population_variance
                )

            # Exact finite-population law for this prime.  A reflected pair
            # contributes 0, 1, or 2 hits to a rational window.  Z(p)<=8 in
            # the main range, so enumerating the selected category counts is
            # both cheap and more reliable than a normal tail approximation.
            category_sizes = Counter(pair_contributions.values())
            size_two = category_sizes[2]
            size_one = category_sizes[1]
            size_zero = population_size - size_one - size_two
            denominator = math.comb(population_size, draws)
            local_distribution = [0.0] * (2 * draws + 1)
            for selected_two in range(min(draws, size_two) + 1):
                remaining_after_two = draws - selected_two
                for selected_one in range(min(remaining_after_two, size_one) + 1):
                    selected_zero = remaining_after_two - selected_one
                    if selected_zero > size_zero:
                        continue
                    ways = (
                        math.comb(size_two, selected_two)
                        * math.comb(size_one, selected_one)
                        * math.comb(size_zero, selected_zero)
                    )
                    local_distribution[2 * selected_two + selected_one] += ways / denominator

            convolved = [0.0] * (len(total_distribution) + len(local_distribution) - 1)
            for old_hits, old_probability in enumerate(total_distribution):
                if old_probability == 0:
                    continue
                for new_hits, new_probability in enumerate(local_distribution):
                    if new_probability:
                        convolved[old_hits + new_hits] += old_probability * new_probability
            total_distribution = convolved

        if variance > 0:
            z_score = (observed - expected) / math.sqrt(variance)
        else:
            z_score = 0.0 if math.isclose(observed, expected) else math.copysign(float("inf"), observed - expected)
        random_hits_needed = max(0, observed - fixed_hits)
        raw_p = sum(total_distribution[random_hits_needed:])
        preliminary.append((a, d, observed, expected, variance, z_score, raw_p))

    trials = len(preliminary)
    return [
        RationalScanResult(
            a=a,
            d=d,
            observed=observed,
            expected=expected,
            variance=variance,
            z_score=z_score,
            raw_p=raw_p,
            bonferroni_p=min(1.0, trials * raw_p),
        )
        for a, d, observed, expected, variance, z_score, raw_p in preliminary
    ]


def legendre_symbol(value: int, p: int) -> int:
    value %= p
    if value == 0:
        return 0
    symbol = pow(value, (p - 1) // 2, p)
    return -1 if symbol == p - 1 else symbol


def hypergeometric_two_sided(
    population: int, successes: int, draws: int, observed: int
) -> float:
    low = max(0, draws - (population - successes))
    high = min(draws, successes)
    denominator = math.comb(population, draws)
    probabilities = {
        value: math.comb(successes, value)
        * math.comb(population - successes, draws - value)
        / denominator
        for value in range(low, high + 1)
    }
    observed_probability = probabilities[observed]
    return min(
        1.0,
        sum(
            probability
            for probability in probabilities.values()
            if probability <= observed_probability + 1e-15
        ),
    )


def polynomial_from_roots(roots: Sequence[int], p: int) -> list[int]:
    """Coefficients low-to-high of the monic product of (X-root)."""
    coefficients = [1]
    for root in roots:
        updated = [0] * (len(coefficients) + 1)
        for degree, coefficient in enumerate(coefficients):
            updated[degree] = (updated[degree] - root * coefficient) % p
            updated[degree + 1] = (updated[degree + 1] + coefficient) % p
        coefficients = updated
    return coefficients


def centered_residue(value: int, modulus: int) -> int:
    residue = value % modulus
    return residue if residue <= (modulus - 1) // 2 else residue - modulus


def evaluate_integer_polynomial_mod(
    coefficients: Sequence[int], value: int, modulus: int
) -> int:
    result = 0
    for coefficient in reversed(coefficients):
        result = (result * value + coefficient) % modulus
    return result


def degree_12_height_screen(
    special: PrimeZeros, witness: PrimeZeros
) -> dict[str, object]:
    """Exhaust a concrete class of fixed degree-12 integer polynomials.

    At q=special.p, any polynomial F of degree at most 12 whose reduction is
    nonzero and vanishes at the 12 special zeros must reduce to
    lambda*product(X-j).  If F is a nonzero integer polynomial and every
    coefficient has absolute value <=(q-1)/2, then its reduction is nonzero
    and its coefficients are the unique centered lifts for one of
    lambda=1,...,q-1.  Enumerating all lambda and testing the witness zeros is
    therefore exhaustive for this explicitly bounded-height class.
    """
    q = special.p
    if special.z != 12:
        return {"applicable": False, "reason": f"special Z is {special.z}, not 12"}
    root_polynomial = polynomial_from_roots(special.zeros, q)
    survivors: list[tuple[int, tuple[int, ...]]] = []
    for scalar in range(1, q):
        coefficients = tuple(
            centered_residue(scalar * coefficient, q)
            for coefficient in root_polynomial
        )
        if all(
            evaluate_integer_polynomial_mod(coefficients, j, witness.p) == 0
            for j in witness.zeros
        ):
            survivors.append((scalar, coefficients))
    return {
        "applicable": True,
        "degree": 12,
        "height_bound": (q - 1) // 2,
        "candidates_tested": q - 1,
        "witness_prime": witness.p,
        "witness_zeros": witness.zeros,
        "survivors": survivors,
    }


def degree_6_invariant_height_screen(
    special: PrimeZeros, witness: PrimeZeros
) -> dict[str, object]:
    """Analogous screen in the forced-reflection quotient Y=X(X+1)."""
    q = special.p
    source_representatives = [j for j in special.zeros if j < special.center]
    witness_representatives = [j for j in witness.zeros if j < witness.center]
    source_invariants = tuple(j * (j + 1) % q for j in source_representatives)
    witness_invariants = tuple(
        j * (j + 1) % witness.p for j in witness_representatives
    )
    if len(source_invariants) != 6:
        return {
            "applicable": False,
            "reason": f"special reflected-pair count is {len(source_invariants)}, not 6",
        }

    root_polynomial = polynomial_from_roots(source_invariants, q)
    survivors: list[tuple[int, tuple[int, ...]]] = []
    for scalar in range(1, q):
        coefficients = tuple(
            centered_residue(scalar * coefficient, q)
            for coefficient in root_polynomial
        )
        if all(
            evaluate_integer_polynomial_mod(coefficients, y, witness.p) == 0
            for y in witness_invariants
        ):
            survivors.append((scalar, coefficients))
    return {
        "applicable": True,
        "degree": 6,
        "height_bound": (q - 1) // 2,
        "candidates_tested": q - 1,
        "witness_prime": witness.p,
        "source_invariants": source_invariants,
        "witness_invariants": witness_invariants,
        "survivors": survivors,
    }


def longest_integer_ap(values: Sequence[int]) -> tuple[int, ...]:
    if not values:
        return ()
    if len(values) == 1:
        return (values[0],)
    value_set = set(values)
    best = (values[0], values[1])
    for first, second in combinations(values, 2):
        difference = second - first
        progression = [first, second]
        next_value = second + difference
        while next_value in value_set:
            progression.append(next_value)
            next_value += difference
        if len(progression) > len(best):
            best = tuple(progression)
    return tuple(best)


def three_term_aps(values: Sequence[int], modulus: int | None = None) -> list[tuple[int, int, int]]:
    triples: list[tuple[int, int, int]] = []
    for triple in combinations(values, 3):
        found = False
        for middle in triple:
            endpoints = [x for x in triple if x != middle]
            expression = endpoints[0] + endpoints[1] - 2 * middle
            if (expression == 0) if modulus is None else (expression % modulus == 0):
                found = True
                break
        if found:
            triples.append(triple)
    return triples


def format_polynomial(coefficients: Sequence[int], variable: str = "Y") -> str:
    pieces = []
    for degree, coefficient in enumerate(coefficients):
        if coefficient:
            pieces.append(f"{coefficient}*{variable}^{degree}")
    return " + ".join(pieces) if pieces else "0"


def pattern_lines(record: PrimeZeros) -> list[str]:
    p = record.p
    zeros = list(record.zeros)
    center = record.center
    representatives = [j for j in zeros if j < center]
    mirror_pairs = [(j, p - 1 - j) for j in representatives]
    gaps = [b - a for a, b in zip(zeros, zeros[1:])]
    gap_counts = Counter(gaps)
    repeated_gaps = sorted((gap, count) for gap, count in gap_counts.items() if count > 1)
    integer_aps = three_term_aps(zeros)
    modular_aps = three_term_aps(zeros, p)
    longest_ap = longest_integer_ap(zeros)
    invariants = [j * (j + 1) % p for j in representatives]
    invariant_ap = longest_integer_ap(sorted(invariants))

    position_symbols = [legendre_symbol(j, p) for j in zeros]
    invariant_symbols = [legendre_symbol(value, p) for value in invariants]
    position_qr = position_symbols.count(1)
    invariant_qr = invariant_symbols.count(1)
    invariant_population = center - 1
    invariant_population_qr = sum(
        legendre_symbol(j * (j + 1), p) == 1
        for j in range(1, center)
    )
    invariant_qr_p = hypergeometric_two_sided(
        invariant_population,
        invariant_population_qr,
        record.pair_count,
        invariant_qr,
    )
    gap_gcd = math.gcd(*gaps) if gaps else 0

    lines = [
        f"p={p}, Z(p)={record.z}",
        f"  zero positions: {zeros}",
        "  normalized positions: "
        + ", ".join(f"{j}/{p}={j/p:.12f}" for j in zeros),
        f"  center zero: {record.has_center_zero}; reflected pairs: {mirror_pairs}",
        f"  pair sums: {[a+b for a, b in mirror_pairs]} (forced value p-1={p-1})",
        f"  consecutive gaps: {gaps}",
        f"  gap palindrome: {gaps == list(reversed(gaps))}; gcd(all gaps)={gap_gcd}; repeated gaps={repeated_gaps}",
        f"  full integer AP: {len(set(gaps)) == 1 if gaps else True}; longest integer AP: {list(longest_ap)} (length {len(longest_ap)})",
        f"  nontrivial integer 3-APs: {integer_aps if integer_aps else 'none'}",
        f"  nontrivial mod-p 3-APs: {modular_aps if modular_aps else 'none'}",
        f"  Legendre(j): residues={position_qr}, nonresidues={position_symbols.count(-1)}, zeros={position_symbols.count(0)} (descriptive only; mirrors are dependent)",
        f"  paired invariants y=j(j+1) mod p: {invariants}",
        f"  Legendre(y), one per pair: {invariant_symbols}; residues={invariant_qr}/{len(invariant_symbols)}, conditional hypergeometric p={invariant_qr_p:.4g}",
        f"  longest integer AP among sorted y-values: {list(invariant_ap)} (length {len(invariant_ap)})",
    ]

    if invariants:
        q_coefficients = polynomial_from_roots(invariants, p)
        lines.extend(
            [
                "  tautological paired root polynomial Q_p(Y)=product(Y-j(j+1)) mod p:",
                f"    coefficients low-to-high: {q_coefficients}",
                f"    expanded: {format_polynomial(q_coefficients)} (mod {p})",
                "    This interpolation is recorded for comparison only; every reflected finite set has such a polynomial.",
            ]
        )
    return lines


def format_report(
    records: Sequence[PrimeZeros],
    special: PrimeZeros,
    validation: dict[str, object],
    null_results: dict[str, object],
    special_null_results: dict[str, object],
    grid_scan: Sequence[RationalScanResult],
    fixed_scan: Sequence[RationalScanResult],
    grid_scan_p_gt_20: Sequence[RationalScanResult],
    fixed_scan_p_gt_20: Sequence[RationalScanResult],
    polynomial_screen: dict[str, object],
    invariant_polynomial_screen: dict[str, object],
    pmax: int,
    simulations: int,
    seed: int,
    elapsed: float,
) -> str:
    points = flatten_points(records)
    values = [j / p for p, j in points]
    npoints = len(points)
    observed_hist = histogram(points)
    observed_heights = height_histogram(points)
    heights = [continued_fraction_height3(j, p) for p, j in points]
    absolute_distances, scaled_distances = nearest_fraction_metrics(points)
    observed_ks = ks_statistic(points)
    observed_cvm = cramer_von_mises_statistic(points)
    orbit_values = independent_orbit_values(points)
    observed_orbit_hist = uniform_histogram(orbit_values)
    observed_orbit_ks = ks_from_values(orbit_values)
    observed_orbit_cvm = cramer_von_mises_from_values(orbit_values)

    hist_mean = null_results["hist_mean"]
    hist_sd = null_results["hist_sd"]
    hist_distance = pearson_distance(observed_hist, hist_mean)
    hist_mc_p = empirical_upper_p(hist_distance, null_results["hist_distances"])
    ks_mc_p = empirical_upper_p(observed_ks, null_results["ks"])
    cvm_mc_p = empirical_upper_p(observed_cvm, null_results["cvm"])
    orbit_hist_distance = pearson_distance(
        observed_orbit_hist, null_results["orbit_hist_mean"]
    )
    orbit_hist_mc_p = empirical_upper_p(
        orbit_hist_distance, null_results["orbit_hist_distances"]
    )
    orbit_ks_mc_p = empirical_upper_p(observed_orbit_ks, null_results["orbit_ks"])
    orbit_cvm_mc_p = empirical_upper_p(
        observed_orbit_cvm, null_results["orbit_cvm"]
    )

    height_mean = null_results["height_mean"]
    height_sd = null_results["height_sd"]
    height_distance = pearson_distance(observed_heights, height_mean)
    height_mc_p = empirical_upper_p(height_distance, null_results["height_distances"])

    z_histogram = Counter(record.z for record in records)
    primes_with_zeros = sum(record.z > 0 for record in records)
    center_primes = [record.p for record in records if record.has_center_zero]
    high_records = [record for record in records if record.z >= 8]

    special_points = flatten_points([special])
    special_ks = ks_statistic(special_points)
    special_cvm = cramer_von_mises_statistic(special_points)
    special_orbits = independent_orbit_values(special_points)
    special_orbit_ks = ks_from_values(special_orbits)
    special_orbit_cvm = cramer_von_mises_from_values(special_orbits)
    special_absolute, special_scaled = nearest_fraction_metrics(special_points)
    special_height_hist = height_histogram(special_points)
    special_height_distance = pearson_distance(
        special_height_hist, special_null_results["height_mean"]
    )

    lines: list[str] = [
        "E1: ZERO POSITION FORENSICS FOR APERY NUMBERS MOD p",
        "=" * 78,
        f"Generated (UTC): {datetime.now(timezone.utc).isoformat(timespec='seconds')}",
        f"Range: every prime 5 <= p <= {pmax}",
        f"Special prime (analysed separately): {special.p}",
        f"Monte Carlo: {simulations} reflected-set replicates; seed={seed}",
        f"Runtime: {elapsed:.3f} seconds",
        "",
        "RECURRENCE AND VALIDATION",
        "-" * 78,
        "(j+1)^3 b_(j+1) = (34j^3+51j^2+27j+5)b_j - j^3 b_(j-1),",
        "b_0=1, b_1=5.  The loop ends at j=p-2, so all divided denominators",
        "are nonzero modulo p.",
        f"Validation status: {'PASS' if validation['ok'] else 'FAIL'}",
        f"Prime records checked: {validation['record_count']}; direct binomial values checked: {validation['direct_checks']}",
        f"Failures: {validation['failures'] if validation['failures'] else 'none'}",
        "",
        "NULL MODEL",
        "-" * 78,
        "The primary uniform comparison is conditional on each prime's observed Z(p).",
        "It samples Z(p)/2 representatives uniformly without replacement from",
        "1,...,(p-3)/2, adds their forced mirrors p-1-j, excludes endpoints, and",
        "holds an observed center zero fixed.  This preserves the exact reflection",
        "law, the p=11,3137 center exceptions, and (by rejection sampling) the exact",
        "no-adjacent-zero rule.  The usual iid KS p-value is also",
        "shown but is not used for the conclusion because mirrored observations are",
        "dependent.  Pearson bin expectations use an independent Monte Carlo pilot",
        f"of {null_results['pearson_pilot_count']} replicates and p-values use the remaining",
        f"{null_results['pearson_heldout_count']} held-out replicates.",
        "",
        "BASIC COUNTS",
        "-" * 78,
        f"Number of primes: {len(records)}",
        f"Primes with at least one zero: {primes_with_zeros}",
        f"Total (p,j) zero pairs: {npoints}",
        f"Mean Z(p): {npoints/len(records):.6f}",
        f"Z(p) histogram: {dict(sorted(z_histogram.items()))}",
        f"Center-zero primes: {center_primes}",
        f"Primes with Z(p)>=8: {[(record.p, record.z) for record in high_records]}",
        "",
        "DISTRIBUTION OF x=j/p",
        "-" * 78,
        f"mean={statistics.fmean(values):.9f} (continuous U[0,1] mean=0.5)",
        f"population sd={statistics.pstdev(values):.9f} (continuous U[0,1] sd={math.sqrt(1/12):.9f})",
        f"min={min(values):.9f}; max={max(values):.9f}",
        "quantiles (observed versus continuous-uniform value):",
    ]
    for q in (0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99):
        lines.append(f"  q={q:>4.2f}: observed={quantile(values, q):.9f}, uniform={q:.9f}")
    lines.extend(
        [
            f"KS D={observed_ks:.6f}; naive iid asymptotic p={naive_ks_pvalue(observed_ks, npoints):.6g}; reflected-null MC p={ks_mc_p:.6g}",
            f"Cramer-von Mises W^2={observed_cvm:.6f}; reflected-null MC p={cvm_mc_p:.6g}",
            f"20-bin Pearson distance from reflected-null mean={hist_distance:.6f}; MC p={hist_mc_p:.6g}",
            "",
            "20-bin histogram (bar is observed count; MC mean/sd preserve reflection):",
        ]
    )
    max_count = max(observed_hist) if observed_hist else 1
    for index, count in enumerate(observed_hist):
        low = index / HISTOGRAM_BINS
        high = (index + 1) / HISTOGRAM_BINS
        bar_length = round(40 * count / max_count)
        sd = hist_sd[index]
        z = (count - hist_mean[index]) / sd if sd else 0.0
        lines.append(
            f"  [{low:4.2f},{high:4.2f}) {count:4d} {100*count/npoints:6.2f}% "
            f"MC={hist_mean[index]:6.2f}+/-{sd:5.2f} z={z:6.2f} |{'#'*bar_length}"
        )

    lines.extend(
        [
            "",
            "Independent-orbit sensitivity check:",
            "For each size-two orbit take its lower representative r and map",
            "u=2r/(p-1) to [0,1]; the two fixed center orbits are omitted.",
            f"  pair orbits={len(orbit_values)}; u histogram (20 bins)={observed_orbit_hist}",
            f"  KS D={observed_orbit_ks:.6f}, reflected-null MC p={orbit_ks_mc_p:.6g}",
            f"  Cramer-von Mises W^2={observed_orbit_cvm:.6f}, MC p={orbit_cvm_mc_p:.6g}",
            f"  20-bin Pearson distance={orbit_hist_distance:.6f}, MC p={orbit_hist_mc_p:.6g}",
        ]
    )

    lines.extend(
        [
            "",
            "LOW-DENOMINATOR RATIONAL FORENSICS (all reduced a/d in [0,1], d<=20)",
            "-" * 78,
            f"Number of candidate rationals: {len(LOW_FRACTIONS)}",
            "Two scans are reported: lattice window |x-a/d|<=1/p and fixed window",
            "|x-a/d|<=0.005.  Expected counts and variances condition on each p's",
            "Z(p), reflection pairs, and fixed center status.  These analytic rational",
            "tails do not additionally condition on the no-adjacent rule.  Raw p-values are exact",
            "conditional upper tails; Bonferroni correction covers all candidate",
            "rationals within each scan; a joint two-scan correction is given below.",
            "",
            "Top positive lattice-window deviations:",
            "  rational  observed  expected   z-score  raw upper-p  within-scan Bonf-p",
        ]
    )
    for result in sorted(grid_scan, key=lambda r: (r.z_score, r.observed), reverse=True)[:15]:
        lines.append(
            f"  {result.a:2d}/{result.d:<2d} {result.observed:9d} {result.expected:9.3f} "
            f"{result.z_score:9.3f} {result.raw_p:9.3g} {result.bonferroni_p:14.3g}"
        )
    lines.extend(
        [
            "",
            "Top positive fixed-window deviations:",
            "  rational  observed  expected   z-score  raw upper-p  within-scan Bonf-p",
        ]
    )
    for result in sorted(fixed_scan, key=lambda r: (r.z_score, r.observed), reverse=True)[:15]:
        lines.append(
            f"  {result.a:2d}/{result.d:<2d} {result.observed:9d} {result.expected:9.3f} "
            f"{result.z_score:9.3f} {result.raw_p:9.3g} {result.bonferroni_p:14.3g}"
        )

    small_prime_exact = [
        (p, j)
        for p, j in points
        if p <= 20 and any(d * j == a * p for a, d in LOW_FRACTIONS)
    ]
    best_grid_gt20 = min(
        (result.bonferroni_p for result in grid_scan_p_gt_20), default=1.0
    )
    best_fixed_gt20 = min(
        (result.bonferroni_p for result in fixed_scan_p_gt_20), default=1.0
    )
    best_grid = min((result.bonferroni_p for result in grid_scan), default=1.0)
    best_fixed = min((result.bonferroni_p for result in fixed_scan), default=1.0)
    best_joint = min(1.0, 2 * min(best_grid, best_fixed))
    best_joint_gt20 = min(1.0, 2 * min(best_grid_gt20, best_fixed_gt20))
    lines.extend(
        [
            "",
            "Small-prime sensitivity:",
            f"  Exact low-denominator identities from p<=20: {small_prime_exact}",
            "  These are tautological because the reduced denominator of j/p is p.",
            f"  Repeating both scans with p>20 gives minimum corrected p-values "
            f"{best_grid_gt20:.6g} (lattice) and {best_fixed_gt20:.6g} (fixed).",
            f"  Correcting jointly across all {2*len(LOW_FRACTIONS)} fraction/window tests gives",
            f"  minimum Bonferroni p={best_joint:.6g} (all p) and {best_joint_gt20:.6g} (p>20).",
        ]
    )

    nearest_observed = {
        "scaled<=0.5": sum(x <= 0.5 + 1e-12 for x in scaled_distances),
        "scaled<=1": sum(x <= 1 + 1e-12 for x in scaled_distances),
        "scaled<=2": sum(x <= 2 + 1e-12 for x in scaled_distances),
        "absolute<=0.001": sum(x <= 0.001 + 1e-15 for x in absolute_distances),
        "absolute<=0.005": sum(x <= 0.005 + 1e-15 for x in absolute_distances),
    }
    lines.extend(
        [
            "",
            "Nearest low-denominator rational summaries:",
            f"  Choosing the nearest member already accounts for the {len(LOW_FRACTIONS)}-rational search;",
            "  each threshold below has an MC upper-tail p-value, but the five threshold",
            "  p-values are descriptive and are not jointly familywise-corrected.",
            f"  absolute-distance median={quantile(absolute_distances, 0.5):.8f}, q90={quantile(absolute_distances, 0.9):.8f}",
            f"  p*distance median={quantile(scaled_distances, 0.5):.5f}, q90={quantile(scaled_distances, 0.9):.5f}",
            "  metric                 observed   MC mean     MC sd    upper-tail MC p",
        ]
    )
    for name, observed in nearest_observed.items():
        simulations_for_metric = null_results["nearest_counts"][name]
        lines.append(
            f"  {name:<22s} {observed:8d} {statistics.fmean(simulations_for_metric):10.3f} "
            f"{statistics.pstdev(simulations_for_metric):9.3f} {empirical_upper_p(observed, simulations_for_metric):18.6g}"
        )

    height_labels = [f"H3<={HEIGHT_CUTOFFS[0]}"] + [
        f"{previous + 1}<=H3<={cutoff}"
        for previous, cutoff in zip(HEIGHT_CUTOFFS, HEIGHT_CUTOFFS[1:])
    ] + [f"H3>{HEIGHT_CUTOFFS[-1]}"]
    lines.extend(
        [
            "",
            "CONTINUED-FRACTION HEIGHT",
            "-" * 78,
            "H3=max denominator among the first three standard convergents, counting",
            "0/1 as the first convergent.  The full CF and H3 for every zero are",
            "recorded later in this report.",
            f"H3 median={quantile([float(h) for h in heights], 0.5):.3f}; q90={quantile([float(h) for h in heights], 0.9):.3f}; max={max(heights)}",
            f"Height-bin Pearson distance from reflected-null mean={height_distance:.6f}; MC p={height_mc_p:.6g}",
            "  bin          observed   MC mean     MC sd       z",
        ]
    )
    for index, label in enumerate(height_labels):
        sd = height_sd[index]
        z = (observed_heights[index] - height_mean[index]) / sd if sd else 0.0
        lines.append(
            f"  {label:<12s} {observed_heights[index]:8d} {height_mean[index]:10.3f} {sd:9.3f} {z:8.3f}"
        )

    lines.extend(["", "HIGH-Z PRIME FORENSICS", "-" * 78])
    if not high_records:
        lines.append(f"No p<={pmax} has Z(p)>=8.")
    for record in high_records:
        lines.extend(pattern_lines(record))
        lines.append("")

    lines.extend(
        [
            f"SPECIAL PRIME p={special.p} (excluded from p<={pmax} aggregate statistics)",
            "-" * 78,
        ]
    )
    lines.extend(pattern_lines(special))
    lines.extend(
        [
            "  conditional-uniform diagnostics for this prime alone:",
            f"    all 12 positions: KS D={special_ks:.6f}, MC p={empirical_upper_p(special_ks, special_null_results['ks']):.6g}; "
            f"CvM={special_cvm:.6f}, MC p={empirical_upper_p(special_cvm, special_null_results['cvm']):.6g}",
            f"    6 independent lower-pair coordinates u=2j/(p-1): {[round(x, 9) for x in special_orbits]}",
            f"    orbit KS D={special_orbit_ks:.6f}, MC p={empirical_upper_p(special_orbit_ks, special_null_results['orbit_ks']):.6g}; "
            f"orbit CvM={special_orbit_cvm:.6f}, MC p={empirical_upper_p(special_orbit_cvm, special_null_results['orbit_cvm']):.6g}",
            f"    H3 bins={special_height_hist}; height-distribution MC p={empirical_upper_p(special_height_distance, special_null_results['height_distances']):.6g}",
            "    nearest-rational metric       observed   MC mean    upper-tail MC p",
        ]
    )
    special_nearest_observed = {
        "scaled<=0.5": sum(x <= 0.5 + 1e-12 for x in special_scaled),
        "scaled<=1": sum(x <= 1 + 1e-12 for x in special_scaled),
        "scaled<=2": sum(x <= 2 + 1e-12 for x in special_scaled),
        "absolute<=0.001": sum(x <= 0.001 + 1e-15 for x in special_absolute),
        "absolute<=0.005": sum(x <= 0.005 + 1e-15 for x in special_absolute),
    }
    for name, observed in special_nearest_observed.items():
        simulated = special_null_results["nearest_counts"][name]
        lines.append(
            f"    {name:<24s} {observed:8d} {statistics.fmean(simulated):9.3f} "
            f"{empirical_upper_p(observed, simulated):18.6g}"
        )
    lines.extend(
        [
            "    These single-prime p-values are exploratory: p=159977 was selected",
            "    because it has the largest observed Z, so post-selection inference would",
            "    require calibration over the prime search that found it.",
            "  per-zero continued-fraction data:",
        ]
    )
    for j in special.zeros:
        terms, convergents, height = continued_fraction_data(j, special.p)
        nearest_a, nearest_d, distance = nearest_low_fraction(j, special.p)
        lines.append(
            f"    j={j:6d}, x={j/special.p:.12f}, CF={terms}, first3={convergents[:3]}, "
            f"H3={height}, nearest={nearest_a}/{nearest_d}, |delta|={distance:.9g}, p|delta|={special.p*distance:.6f}"
        )

    lines.extend(["", "FIXED-POLYNOMIAL HEIGHT SCREEN", "-" * 78])
    if polynomial_screen.get("applicable"):
        survivor_count = len(polynomial_screen["survivors"])
        lines.extend(
            [
                "At p=159977, let F be a nonzero integer polynomial of degree <=12",
                "and coefficient height <=79988.  If F vanishes at all 12 special",
                "zeros, its nonzero reduction must be lambda*product(X-j), and each",
                "coefficient is the unique centered lift of that reduction.",
                f"Scalars/candidates exhaustively tested: {polynomial_screen['candidates_tested']}",
                f"Additional witness: p={polynomial_screen['witness_prime']} at zeros={list(polynomial_screen['witness_zeros'])}",
                f"Survivors vanishing at every witness zero: {survivor_count}",
                (
                    "Result: no fixed degree-12 integer polynomial of coefficient height <=79988 fits both high-Z sets."
                    if survivor_count == 0
                    else f"Surviving candidates (scalar, coefficients): {polynomial_screen['survivors'][:10]}"
                ),
                "Together with the root-count argument, this excludes degrees <=11",
                "without a height bound (at good nonzero reductions), and excludes",
                "degree 12 within the stated coefficient-height bound.  It says",
                "nothing conclusive about degree >12, larger coefficients, reduction",
                "identically zero at an exceptional prime, or conditions that vary with p.",
            ]
        )
    else:
        lines.append(f"Screen not applicable: {polynomial_screen.get('reason', 'unknown reason')}")

    if invariant_polynomial_screen.get("applicable"):
        invariant_survivor_count = len(invariant_polynomial_screen["survivors"])
        lines.extend(
            [
                "",
                "Reflection-quotient screen Y=X(X+1):",
                f"  special Y-roots: {list(invariant_polynomial_screen['source_invariants'])}",
                f"  witness Y-roots mod {invariant_polynomial_screen['witness_prime']}: {list(invariant_polynomial_screen['witness_invariants'])}",
                f"  degree-6 centered-lift candidates tested: {invariant_polynomial_screen['candidates_tested']}",
                f"  survivors: {invariant_survivor_count}",
                (
                    "  Result: no fixed degree-6 G(Y) of coefficient height <=79988 fits both quotient sets."
                    if invariant_survivor_count == 0
                    else f"  Surviving candidates: {invariant_polynomial_screen['survivors'][:10]}"
                ),
            ]
        )
    else:
        lines.append(
            f"Reflection-quotient screen not applicable: {invariant_polynomial_screen.get('reason', 'unknown reason')}"
        )

    location_random = ks_mc_p >= 0.05 and hist_mc_p >= 0.05 and cvm_mc_p >= 0.05
    rational_random = best_joint >= 0.05 and best_joint_gt20 >= 0.05
    high_max = max([record.z for record in high_records] + [special.z])
    lines.extend(
        [
            "",
            "INTERPRETATION AND ANSWER TO THE KEY QUESTION",
            "-" * 78,
            (
                "The aggregate locations show no statistically significant departure from the "
                "reflection-preserving uniform null at the 5% level."
                if location_random
                else "At least one aggregate location statistic departs from the reflection-preserving uniform null at the 5% level; inspect the tests above before interpreting it."
            ),
            (
                f"After scanning all {len(LOW_FRACTIONS)} low-denominator rationals in both windows, no excess survives the joint Bonferroni correction."
                if rational_random
                else "At least one low-denominator window survives the nominal Bonferroni 5% threshold; this is a candidate cluster, not by itself an algebraic explanation."
            ),
            "The high-Z sets show the forced reflection j -> p-1-j (and hence palindromic",
            "gaps).  The AP and Legendre diagnostics above test additional elementary",
            "patterns; with only 4 and 6 independent pairs, absence of a 3-AP is expected",
            "for sparse random sets and is not positive evidence.  The tautological Q_p",
            "polynomials are likewise not evidence for a fixed law.",
            f"The {high_max} distinct zeros at p={special.p} rule out a single nonzero fixed",
            f"polynomial of degree <{high_max} whose root set must contain every observed zero",
            "at every good prime (unless its reduction is identically zero at this exceptional",
            "prime).  They do not rule out degree >=12 or a more complicated algebraic condition.",
            (
                "The explicit centered-lift screen also finds no degree-12 fixed integer polynomial of coefficient height <=79988 fitting both p=159977 and p=3727."
                if polynomial_screen.get("applicable") and not polynomial_screen["survivors"]
                else "The explicit degree-12 bounded-height screen is reported above."
            ),
            (
                "In the natural quotient Y=j(j+1), the analogous degree-6 coefficient-height screen also has no survivor."
                if invariant_polynomial_screen.get("applicable") and not invariant_polynomial_screen["survivors"]
                else "The reflection-quotient bounded-height screen is reported above."
            ),
            "Crucially, uniform-looking normalized roots do not distinguish random points from",
            "roots of a fixed polynomial modulo varying primes: algebraic roots can themselves",
            "equidistribute.  Also, these j are zero COEFFICIENTS of the truncated Hasse-Witt",
            "polynomial sum b_j t^j, not roots of that polynomial.  Root-locus ordinarity theory",
            "therefore does not directly apply to this statistic.",
            "Bottom line: within this experiment the positions are consistent with conditional",
            "uniform randomness beyond the exact reflection/center constraints, and no small-",
            "denominator, AP, or quadratic-residue signature is detected.  This is evidence",
            "against fixed real locations, small-coefficient linear factors, and other very simple",
            "position laws, but location statistics alone cannot",
            "exclude a bounded-degree algebraic mechanism.  A defensible verdict is 'random-",
            "looking, not a proof of non-algebraicity.'",
            "",
            "ALL PRIME ZERO SETS AND PER-ZERO RECORDS",
            "-" * 78,
            "Notation: H3 uses the convention defined above; nearest searches reduced a/d in [0,1] with d<=20.",
        ]
    )
    for record in records:
        set_text = "{" + ", ".join(map(str, record.zeros)) + "}"
        lines.append(f"p={record.p:4d} Z={record.z:2d} zeros={set_text}")
        for j in record.zeros:
            terms, convergents, height = continued_fraction_data(j, record.p)
            nearest_a, nearest_d, distance = nearest_low_fraction(j, record.p)
            lines.append(
                f"  j={j:4d} x={j/record.p:.12f} CF={terms} first3={convergents[:3]} "
                f"H3={height} nearest={nearest_a}/{nearest_d} |delta|={distance:.9g} p|delta|={record.p*distance:.6f}"
            )

    return "\n".join(lines) + "\n"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--p-max", type=int, default=5_000, help="largest prime in aggregate experiment")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT, help="text report path")
    parser.add_argument(
        "--simulations",
        type=int,
        default=DEFAULT_SIMULATIONS,
        help="reflection-preserving Monte Carlo replicates",
    )
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED, help="Monte Carlo seed")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.p_max < 5:
        raise SystemExit("--p-max must be at least 5")
    if args.simulations < 100:
        raise SystemExit("--simulations must be at least 100")
    start = time.perf_counter()
    primes = [p for p in sieve_primes(args.p_max) if p >= 5]
    print(f"Computing all zero positions for {len(primes)} primes through {args.p_max}...", flush=True)
    records = [PrimeZeros(p, apery_zero_positions(p)) for p in primes]

    print(f"Computing special prime p={DEFAULT_SPECIAL_PRIME}...", flush=True)
    special = PrimeZeros(
        DEFAULT_SPECIAL_PRIME, apery_zero_positions(DEFAULT_SPECIAL_PRIME)
    )
    validation = validate_records([*records, special])
    if not validation["ok"]:
        raise RuntimeError(f"validation failed: {validation['failures']}")

    print(f"Running {args.simulations} reflected-null simulations...", flush=True)
    null_results = run_null_simulations(records, args.simulations, args.seed)
    special_null_results = run_null_simulations(
        [special], args.simulations, args.seed + 1
    )
    print("Scanning low-denominator rational windows...", flush=True)
    grid_scan = scan_rational_clusters(records, "grid")
    fixed_scan = scan_rational_clusters(records, "fixed")
    records_p_gt_20 = [record for record in records if record.p > 20]
    grid_scan_p_gt_20 = scan_rational_clusters(records_p_gt_20, "grid")
    fixed_scan_p_gt_20 = scan_rational_clusters(records_p_gt_20, "fixed")

    witness = next((record for record in records if record.p == 3727), None)
    if witness is None:
        witness = PrimeZeros(3727, apery_zero_positions(3727))
    print("Screening bounded-height fixed polynomials in X and X(X+1)...", flush=True)
    polynomial_screen = degree_12_height_screen(special, witness)
    invariant_polynomial_screen = degree_6_invariant_height_screen(
        special, witness
    )

    elapsed = time.perf_counter() - start
    report = format_report(
        records=records,
        special=special,
        validation=validation,
        null_results=null_results,
        special_null_results=special_null_results,
        grid_scan=grid_scan,
        fixed_scan=fixed_scan,
        grid_scan_p_gt_20=grid_scan_p_gt_20,
        fixed_scan_p_gt_20=fixed_scan_p_gt_20,
        polynomial_screen=polynomial_screen,
        invariant_polynomial_screen=invariant_polynomial_screen,
        pmax=args.p_max,
        simulations=args.simulations,
        seed=args.seed,
        elapsed=elapsed,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(report, encoding="utf-8")
    print(f"Wrote {args.output} ({len(report):,} bytes) in {elapsed:.2f}s")


if __name__ == "__main__":
    main()
