#!/usr/bin/env python3
"""Search for a closed form of the nonzero tau quarter-point value.

For every prime p == 1 (mod 24), p < 20000, this program computes

    v_p = tau[(p - 1) / 4] (mod p),

together with the quadratic-form, binomial, Jacobsthal-sum, quartic-character,
and midpoint data requested in CODEX_SPEC_CRON_quarterpoint_value.md.  It then
runs three deterministic searches:

* bounded-rational multiples of degree-at-most-two monomials;
* small-coefficient quadratic polynomials for v^2 and v^4;
* linear forms twisted by small quartic characters.

The default run prints the complete search summary.  ``--dump-csv`` emits the
267-row data set instead.  The code uses only Python's standard library,
including a small exact-rational LLL reducer for the dense polynomial probe.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import math
import statistics
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from math import gcd, isqrt
from typing import Callable, Iterable, Sequence


LIMIT = 20_000
EXPECTED_ROW_COUNT = 267
PHI_PARAMETERS = (1, 2, 3, 6, -1, -2, -3, -6)
RATIONAL_BOUND = 64
MONOMIAL_HIT_THRESHOLD = 10
POLYNOMIAL_COEFFICIENT_BOUND = 32
POLYNOMIAL_MULTI_TERM_BOUND = 8
POLYNOMIAL_PROBE_COUNT = 10
QUARTIC_SMALL_BASES = tuple(n for n in range(-12, 13) if n not in (-1, 0, 1))


@dataclass(frozen=True)
class PrimeRow:
    p: int
    v: int
    x: int
    y: int
    a: int
    b: int
    c: int
    d: int
    g1: int
    g2: int
    phi: tuple[int, ...]
    q2: int
    q3: int
    q6: int
    tau_midpoint: int
    quartic_small: tuple[int, ...]


@dataclass(frozen=True)
class Fit:
    family: str
    formula: str
    hits: int
    failures: tuple[int, ...]
    prediction: tuple[int, ...]


def primes_below(limit: int) -> list[int]:
    sieve = bytearray(b"\1") * limit
    sieve[:2] = b"\0\0"
    for divisor in range(2, isqrt(limit - 1) + 1):
        if sieve[divisor]:
            start = divisor * divisor
            count = (limit - 1 - start) // divisor + 1
            sieve[start:limit:divisor] = b"\0" * count
    return [p for p in range(2, limit) if sieve[p]]


def inverses_through(top: int, p: int) -> list[int]:
    inverse = [0] * (top + 1)
    if top >= 1:
        inverse[1] = 1
    for n in range(2, top + 1):
        inverse[n] = -(p // n) * inverse[p % n] % p
    return inverse


def tau_mod(p: int, top: int, inverse: Sequence[int] | None = None) -> list[int]:
    """Return tau_0,...,tau_top from the rank-two recurrence modulo p."""
    if top == 0:
        return [1]
    if inverse is None:
        inverse = inverses_through(top + 1, p)
    values = [1, 5 * pow(2, -1, p) % p]
    inverse_four = pow(4, -1, p)
    for j in range(top - 1):
        middle = 68 * j * j + 170 * j + 107
        numerator = 2 * middle * values[-1] - (2 * j + 1) ** 2 * values[-2]
        denominator_inverse = inverse_four * inverse[j + 2] ** 2
        values.append(numerator * denominator_inverse % p)
    return values


def apery_mod(p: int, top: int, inverse: Sequence[int]) -> list[int]:
    values = [1, 5]
    for n in range(1, top):
        numerator = (2 * n + 1) * (17 * n * n + 17 * n + 5) * values[n]
        numerator -= n**3 * values[n - 1]
        values.append(numerator * inverse[n + 1] ** 3 % p)
    return values[: top + 1]


def direct_sqrt_mod(coefficients: Sequence[int], p: int) -> list[int]:
    root = [1]
    inverse_two = pow(2, -1, p)
    for n in range(1, len(coefficients)):
        convolution = sum(root[j] * root[n - j] for j in range(1, n))
        root.append((coefficients[n] - convolution) * inverse_two % p)
    return root


def binomial_mod(n: int, k: int, p: int, inverse: Sequence[int]) -> int:
    k = min(k, n - k)
    answer = 1
    for j in range(1, k + 1):
        answer = answer * (n - k + j) % p
        answer = answer * inverse[j] % p
    return answer


def positive_representations(p: int, coefficient: int) -> tuple[tuple[int, int], ...]:
    """Positive (r,s) with p = r^2 + coefficient*s^2."""
    answers = []
    for r in range(1, isqrt(p - coefficient) + 1):
        remainder = p - r * r
        if remainder <= 0 or remainder % coefficient:
            continue
        s = isqrt(remainder // coefficient)
        if s > 0 and r * r + coefficient * s * s == p:
            answers.append((r, s))
    return tuple(answers)


def gaussian_representation(p: int) -> tuple[int, int]:
    """p = a^2+b^2 with b>0 and the signed odd a congruent to 1 mod 4."""
    candidates = []
    for odd in range(1, isqrt(p - 1) + 1, 2):
        square = p - odd * odd
        even = isqrt(square)
        if even > 0 and even * even == square and even % 2 == 0:
            signed_odd = odd if odd % 4 == 1 else -odd
            candidates.append((signed_odd, even))
    if len(candidates) != 1:
        raise AssertionError((p, candidates))
    return candidates[0]


def legendre_table(p: int) -> list[int]:
    table = [-1] * p
    table[0] = 0
    for value in range(1, (p + 1) // 2):
        table[value * value % p] = 1
    return table


def jacobsthal_sums(p: int) -> tuple[int, ...]:
    table = legendre_table(p)
    totals = [0] * len(PHI_PARAMETERS)
    parameters_mod_p = tuple(k % p for k in PHI_PARAMETERS)
    for t in range(p):
        square = t * t % p
        for index, parameter in enumerate(parameters_mod_p):
            totals[index] += table[t * (square + parameter) % p]
    return tuple(totals)


def legendre_symbol(value: int, p: int) -> int:
    residue = pow(value % p, (p - 1) // 2, p)
    if residue == p - 1:
        return -1
    return residue


def build_rows(limit: int = LIMIT) -> list[PrimeRow]:
    rows = []
    direct_cross_checks = 0
    for p in primes_below(limit):
        if p % 24 != 1:
            continue
        quarter = (p - 1) // 4
        midpoint = (p - 1) // 2
        inverse = inverses_through(midpoint + 1, p)
        tau = tau_mod(p, midpoint, inverse)

        rep_16 = positive_representations(p, 6)
        rep_12 = positive_representations(p, 2)
        if len(rep_16) != 1 or len(rep_12) != 1:
            raise AssertionError((p, rep_16, rep_12))
        x, y = rep_16[0]
        c, d = rep_12[0]
        a, b = gaussian_representation(p)

        g1 = binomial_mod((p - 1) // 2, quarter, p, inverse)
        g2 = binomial_mod(quarter, (p - 1) // 8, p, inverse)
        q2 = pow(2, quarter, p)
        q3 = pow(3, quarter, p)
        q6 = pow(6, quarter, p)
        quartic_small = tuple(pow(base % p, quarter, p) for base in QUARTIC_SMALL_BASES)

        row = PrimeRow(
            p=p,
            v=tau[quarter],
            x=x,
            y=y,
            a=a,
            b=b,
            c=c,
            d=d,
            g1=g1,
            g2=g2,
            phi=jacobsthal_sums(p),
            q2=q2,
            q3=q3,
            q6=q6,
            tau_midpoint=tau[midpoint],
            quartic_small=quartic_small,
        )

        assert row.v != 0
        assert row.g1 == 2 * row.a % p
        assert row.q2 in (1, p - 1) and row.q3 in (1, p - 1)
        assert row.q6 == row.q2 * row.q3 % p
        assert row.tau_midpoint == legendre_symbol(-2, p) % p == 1

        if p < 250:
            direct = direct_sqrt_mod(apery_mod(p, midpoint, inverse), p)
            assert direct == tau
            direct_cross_checks += 1
        rows.append(row)

    if limit == LIMIT:
        assert len(rows) == EXPECTED_ROW_COUNT
        expected_initial_values = {
            73: 19,
            97: 1,
            193: 187,
            241: 128,
            313: 200,
            337: 175,
            409: 71,
            433: 432,
            457: 338,
            577: 386,
            601: 30,
            673: 432,
        }
        observed = {row.p: row.v for row in rows if row.p in expected_initial_values}
        assert observed == expected_initial_values
        assert direct_cross_checks == 4
    return rows


BASE_FEATURE_NAMES = ("x", "y", "a", "b") + tuple(
    f"phi({parameter})" for parameter in PHI_PARAMETERS
)


def base_features(row: PrimeRow) -> tuple[int, ...]:
    return (row.x, row.y, row.a, row.b) + row.phi


MONOMIALS: tuple[tuple[str, tuple[int, ...]], ...] = (
    (("1", ()),)
    + tuple((name, (index,)) for index, name in enumerate(BASE_FEATURE_NAMES))
    + tuple(
        (
            f"{BASE_FEATURE_NAMES[left]}*{BASE_FEATURE_NAMES[right]}",
            (left, right),
        )
        for left in range(len(BASE_FEATURE_NAMES))
        for right in range(left, len(BASE_FEATURE_NAMES))
    )
)


def monomial_value(row: PrimeRow, indices: tuple[int, ...]) -> int:
    values = base_features(row)
    answer = 1
    for index in indices:
        answer = answer * values[index] % row.p
    return answer


def bounded_fractions(bound: int) -> tuple[Fraction, ...]:
    values = {
        Fraction(numerator, denominator)
        for denominator in range(1, bound + 1)
        for numerator in range(-bound, bound + 1)
        if numerator != 0 and gcd(abs(numerator), denominator) == 1
    }
    return tuple(sorted(values, key=lambda value: (abs(value), value)))


def rational_residue_maps(
    rows: Sequence[PrimeRow], fractions: Sequence[Fraction]
) -> dict[int, dict[int, tuple[Fraction, ...]]]:
    maps = {}
    for row in rows:
        by_residue: dict[int, list[Fraction]] = defaultdict(list)
        inverse_denominator = {
            denominator: pow(denominator, -1, row.p)
            for denominator in {value.denominator for value in fractions}
        }
        for value in fractions:
            residue = value.numerator * inverse_denominator[value.denominator] % row.p
            by_residue[residue].append(value)
        maps[row.p] = {key: tuple(values) for key, values in by_residue.items()}
    return maps


def format_fraction(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"({value.numerator}/{value.denominator})"


def make_fit(
    family: str,
    formula: str,
    rows: Sequence[PrimeRow],
    prediction_function: Callable[[PrimeRow], int],
) -> Fit:
    prediction = tuple(prediction_function(row) % row.p for row in rows)
    failures = tuple(row.p for row, predicted in zip(rows, prediction) if row.v != predicted)
    return Fit(family, formula, len(rows) - len(failures), failures, prediction)


def deduplicate_fits(fits: Iterable[Fit]) -> list[Fit]:
    # The same prediction vector can be compared with different targets
    # (notably v/G1 versus v*G1), so its exact failure set is part of the key.
    best_by_prediction: dict[tuple[tuple[int, ...], tuple[int, ...]], Fit] = {}
    for fit in fits:
        key = (fit.prediction, fit.failures)
        old = best_by_prediction.get(key)
        if old is None or (len(fit.formula), fit.formula) < (len(old.formula), old.formula):
            best_by_prediction[key] = fit
    return sorted(
        best_by_prediction.values(),
        key=lambda fit: (-fit.hits, len(fit.formula), fit.formula),
    )


def deduplicate_supports(fits: Iterable[Fit]) -> list[Fit]:
    """Keep the shortest formula for each exact hit/failure pattern."""
    best_by_failures: dict[tuple[int, ...], Fit] = {}
    for fit in fits:
        old = best_by_failures.get(fit.failures)
        if old is None or (len(fit.formula), fit.formula) < (len(old.formula), old.formula):
            best_by_failures[fit.failures] = fit
    return sorted(
        best_by_failures.values(),
        key=lambda fit: (-fit.hits, len(fit.formula), fit.formula),
    )


def search_scaled_monomials(rows: Sequence[PrimeRow]) -> tuple[list[Fit], list[Fit]]:
    fractions = bounded_fractions(RATIONAL_BOUND)
    residue_maps = rational_residue_maps(rows, fractions)
    threshold_fits = []
    local_best = []

    for exponent in (-1, 0, 1):
        for monomial_name, indices in MONOMIALS:
            counts: Counter[Fraction] = Counter()
            for row in rows:
                if exponent == -1:
                    g_factor = pow(row.g1, -1, row.p)
                elif exponent == 0:
                    g_factor = 1
                else:
                    g_factor = row.g1
                scale = g_factor * monomial_value(row, indices) % row.p
                if scale == 0:
                    continue
                ratio = row.v * pow(scale, -1, row.p) % row.p
                counts.update(residue_maps[row.p].get(ratio, ()))

            if not counts:
                continue
            candidates = counts.most_common(4)
            for constant, count in candidates:
                def predictor(
                    row: PrimeRow,
                    constant: Fraction = constant,
                    exponent: int = exponent,
                    indices: tuple[int, ...] = indices,
                ) -> int:
                    rational = constant.numerator * pow(constant.denominator, -1, row.p)
                    if exponent == -1:
                        g_factor = pow(row.g1, -1, row.p)
                    elif exponent == 0:
                        g_factor = 1
                    else:
                        g_factor = row.g1
                    return rational * g_factor * monomial_value(row, indices) % row.p

                g_text = "G1^-1" if exponent == -1 else ("1" if exponent == 0 else "G1")
                formula = f"v = {format_fraction(constant)} * {g_text} * {monomial_name}"
                fit = make_fit("monomial", formula, rows, predictor)
                assert fit.hits == count
                local_best.append(fit)
            for constant, count in counts.items():
                if count < MONOMIAL_HIT_THRESHOLD or constant in dict(candidates):
                    continue

                def predictor(
                    row: PrimeRow,
                    constant: Fraction = constant,
                    exponent: int = exponent,
                    indices: tuple[int, ...] = indices,
                ) -> int:
                    rational = constant.numerator * pow(constant.denominator, -1, row.p)
                    g_factor = (
                        pow(row.g1, -1, row.p)
                        if exponent == -1
                        else (1 if exponent == 0 else row.g1)
                    )
                    return rational * g_factor * monomial_value(row, indices) % row.p

                g_text = "G1^-1" if exponent == -1 else ("1" if exponent == 0 else "G1")
                formula = f"v = {format_fraction(constant)} * {g_text} * {monomial_name}"
                threshold_fits.append(make_fit("monomial", formula, rows, predictor))

    all_best = deduplicate_fits(local_best + threshold_fits)
    all_threshold = [fit for fit in all_best if fit.hits >= MONOMIAL_HIT_THRESHOLD]
    return all_threshold, all_best


# Modulo p, terms containing p vanish.  The identities x^2+6y^2=p and
# a^2+b^2=p let us eliminate x^2 and a^2.  This is the nonredundant quotient
# basis for the requested degree-at-most-two polynomial search.
POLYNOMIAL_TERMS: tuple[tuple[str, Callable[[PrimeRow], int]], ...] = (
    ("1", lambda row: 1),
    ("x", lambda row: row.x),
    ("y", lambda row: row.y),
    ("a", lambda row: row.a),
    ("b", lambda row: row.b),
    ("x*y", lambda row: row.x * row.y),
    ("x*a", lambda row: row.x * row.a),
    ("x*b", lambda row: row.x * row.b),
    ("y^2", lambda row: row.y * row.y),
    ("y*a", lambda row: row.y * row.a),
    ("y*b", lambda row: row.y * row.b),
    ("a*b", lambda row: row.a * row.b),
    ("b^2", lambda row: row.b * row.b),
)


def polynomial_values(row: PrimeRow) -> tuple[int, ...]:
    return tuple(function(row) % row.p for _, function in POLYNOMIAL_TERMS)


def format_polynomial(coefficients: Sequence[int]) -> str:
    pieces = []
    for coefficient, (name, _) in zip(coefficients, POLYNOMIAL_TERMS):
        if coefficient:
            pieces.append(f"{coefficient:+d}*{name}")
    if not pieces:
        return "0"
    text = " ".join(pieces)
    return text[1:] if text.startswith("+") else text


def polynomial_fit(
    rows: Sequence[PrimeRow], power: int, coefficients: tuple[int, ...], family: str
) -> Fit:
    def predictor(row: PrimeRow) -> int:
        return sum(
            coefficient * value
            for coefficient, value in zip(coefficients, polynomial_values(row))
        ) % row.p

    prediction = tuple(predictor(row) for row in rows)
    failures = tuple(
        row.p
        for row, predicted in zip(rows, prediction)
        if pow(row.v, power, row.p) != predicted
    )
    formula = f"v^{power} = {format_polynomial(coefficients)}"
    return Fit(family, formula, len(rows) - len(failures), failures, prediction)


def target_signature(rows: Sequence[PrimeRow], power: int) -> tuple[int, ...]:
    return tuple(pow(row.v, power, row.p) for row in rows)


def exhaustive_sparse_polynomials(
    rows: Sequence[PrimeRow], power: int
) -> tuple[set[tuple[int, ...]], int]:
    """Find small sparse laws on the probe primes, then let the caller verify all."""
    probes = rows[:POLYNOMIAL_PROBE_COUNT]
    target = target_signature(probes, power)
    term_count = len(POLYNOMIAL_TERMS)
    term_columns = [
        tuple(polynomial_values(row)[index] for row in probes)
        for index in range(term_count)
    ]
    found: set[tuple[int, ...]] = set()
    tested = 0

    # One and two nonzero terms, with the full coefficient bound.
    full_coefficients = tuple(
        coefficient
        for coefficient in range(-POLYNOMIAL_COEFFICIENT_BOUND, POLYNOMIAL_COEFFICIENT_BOUND + 1)
        if coefficient
    )
    for first in range(term_count):
        column_first = term_columns[first]
        for coefficient_first in full_coefficients:
            tested += 1
            signature_first = tuple(
                coefficient_first * value % row.p
                for value, row in zip(column_first, probes)
            )
            if signature_first == target:
                coefficients = [0] * term_count
                coefficients[first] = coefficient_first
                found.add(tuple(coefficients))
        for second in range(first + 1, term_count):
            column_second = term_columns[second]
            for coefficient_first in full_coefficients:
                first_values = tuple(
                    coefficient_first * value % row.p
                    for value, row in zip(column_first, probes)
                )
                for coefficient_second in full_coefficients:
                    tested += 1
                    if all(
                        (left + coefficient_second * right - wanted) % row.p == 0
                        for left, right, wanted, row in zip(
                            first_values, column_second, target, probes
                        )
                    ):
                        coefficients = [0] * term_count
                        coefficients[first] = coefficient_first
                        coefficients[second] = coefficient_second
                        found.add(tuple(coefficients))

    # Three and four nonzero terms, exhaustively at the smaller bound, using
    # meet-in-the-middle signatures across all probe primes.
    multi_coefficients = tuple(
        coefficient
        for coefficient in range(-POLYNOMIAL_MULTI_TERM_BOUND, POLYNOMIAL_MULTI_TERM_BOUND + 1)
        if coefficient
    )
    atoms = []
    for index, column in enumerate(term_columns):
        for coefficient in multi_coefficients:
            signature = tuple(
                coefficient * value % row.p for value, row in zip(column, probes)
            )
            atoms.append((index, coefficient, signature))

    pairs_by_signature: dict[
        tuple[int, ...], list[tuple[int, int, int, int]]
    ] = defaultdict(list)
    pairs = []
    for left in range(term_count):
        for right in range(left + 1, term_count):
            for coefficient_left in multi_coefficients:
                for coefficient_right in multi_coefficients:
                    tested += 1
                    signature = tuple(
                        (
                            coefficient_left * term_columns[left][probe_index]
                            + coefficient_right * term_columns[right][probe_index]
                        )
                        % probes[probe_index].p
                        for probe_index in range(len(probes))
                    )
                    item = (left, coefficient_left, right, coefficient_right)
                    pairs_by_signature[signature].append(item)
                    pairs.append((item, signature))

    for third, coefficient_third, signature_third in atoms:
        complement = tuple(
            (wanted - value) % row.p
            for wanted, value, row in zip(target, signature_third, probes)
        )
        for left, coefficient_left, right, coefficient_right in pairs_by_signature.get(
            complement, ()
        ):
            if right >= third:
                continue
            coefficients = [0] * term_count
            coefficients[left] = coefficient_left
            coefficients[right] = coefficient_right
            coefficients[third] = coefficient_third
            found.add(tuple(coefficients))

    for (third, coefficient_third, fourth, coefficient_fourth), signature_pair in pairs:
        complement = tuple(
            (wanted - value) % row.p
            for wanted, value, row in zip(target, signature_pair, probes)
        )
        for left, coefficient_left, right, coefficient_right in pairs_by_signature.get(
            complement, ()
        ):
            if right >= third:
                continue
            coefficients = [0] * term_count
            coefficients[left] = coefficient_left
            coefficients[right] = coefficient_right
            coefficients[third] = coefficient_third
            coefficients[fourth] = coefficient_fourth
            found.add(tuple(coefficients))

    return found, tested


def crt_combine(residues: Sequence[int], moduli: Sequence[int]) -> tuple[int, int]:
    value = 0
    modulus = 1
    for residue, prime in zip(residues, moduli):
        correction = (residue - value) * pow(modulus, -1, prime) % prime
        value += modulus * correction
        modulus *= prime
    return value % modulus, modulus


def exact_lll(
    input_basis: Sequence[Sequence[int]], delta: Fraction = Fraction(3, 4)
) -> list[list[int]]:
    """Incremental exact-rational LLL for the 15-dimensional CRT lattice.

    This is the standard Lenstra--Lenstra--Lovasz update scheme.  Keeping it
    local avoids a SymPy 1.14 pure-Python postcondition bug on this particular
    large-entry lattice.  Decoded candidates are still verified independently
    modulo all 267 primes.
    """
    y = [list(map(int, row)) for row in input_basis]
    row_count = len(y)
    dimension = len(y[0])
    y_star = [[Fraction(0)] * dimension for _ in range(row_count)]
    mu = [[Fraction(0)] * row_count for _ in range(row_count)]
    g_star = [Fraction(0)] * row_count
    half = Fraction(1, 2)

    def dot(left: Sequence[int | Fraction], right: Sequence[int | Fraction]) -> Fraction:
        return sum((a * b for a, b in zip(left, right)), Fraction(0))

    def closest_integer(value: Fraction) -> int:
        return math.floor(value + half)

    def reduce_row(row: int, previous: int) -> None:
        quotient = closest_integer(mu[row][previous])
        y[row] = [
            value - quotient * old for value, old in zip(y[row], y[previous])
        ]
        for index in range(previous):
            mu[row][index] -= quotient * mu[previous][index]
        mu[row][previous] -= quotient

    for row in range(row_count):
        y_star[row] = [Fraction(value) for value in y[row]]
        for previous in range(row):
            if g_star[previous] == 0:
                raise AssertionError("LLL basis lost rank")
            mu[row][previous] = dot(y[row], y_star[previous]) / g_star[previous]
            y_star[row] = [
                value - mu[row][previous] * old
                for value, old in zip(y_star[row], y_star[previous])
            ]
        g_star[row] = dot(y_star[row], y_star[row])

    row = 1
    iterations = 0
    while row < row_count:
        iterations += 1
        if iterations > 100_000:
            raise RuntimeError("exact LLL did not converge")
        if abs(mu[row][row - 1]) > half:
            reduce_row(row, row - 1)
        if g_star[row] >= (
            delta - mu[row][row - 1] ** 2
        ) * g_star[row - 1]:
            for previous in range(row - 2, -1, -1):
                if abs(mu[row][previous]) > half:
                    reduce_row(row, previous)
            row += 1
            continue

        nu = mu[row][row - 1]
        alpha = g_star[row] + nu**2 * g_star[row - 1]
        if alpha == 0:
            raise AssertionError("LLL basis lost rank")
        beta = g_star[row - 1] / alpha
        mu[row][row - 1] = nu * beta
        g_star[row] *= beta
        g_star[row - 1] = alpha
        y[row], y[row - 1] = y[row - 1], y[row]
        mu[row][: row - 1], mu[row - 1][: row - 1] = (
            mu[row - 1][: row - 1],
            mu[row][: row - 1],
        )
        for following in range(row + 1, row_count):
            old = mu[following][row]
            mu[following][row] = mu[following][row - 1] - nu * old
            mu[following][row - 1] = mu[row][row - 1] * mu[following][row] + old
        row = max(1, row - 1)
    return y


def dense_lll_polynomials(
    rows: Sequence[PrimeRow], power: int
) -> set[tuple[int, ...]]:
    """Use LLL on several CRT embeddings to find dense small laws."""
    term_count = len(POLYNOMIAL_TERMS)
    candidates: set[tuple[int, ...]] = set()
    for probe_count in (10, 12, 16):
        probes = rows[:probe_count]
        moduli = [row.p for row in probes]
        target_crt, modulus = crt_combine(
            [pow(row.v, power, row.p) for row in probes], moduli
        )
        term_crt = []
        for term_index in range(term_count):
            combined, check_modulus = crt_combine(
                [polynomial_values(row)[term_index] for row in probes], moduli
            )
            assert check_modulus == modulus
            term_crt.append(combined)

        weight = 10_000
        dimension = term_count + 2
        basis = [[0] * dimension for _ in range(dimension)]
        basis[0][0] = modulus * weight
        for index, combined in enumerate(term_crt):
            basis[index + 1][0] = combined * weight
            basis[index + 1][index + 1] = 1
        basis[-1][0] = -target_crt * weight
        basis[-1][-1] = 1
        reduced = exact_lll(basis)

        # LLL sometimes puts the affine vector in a short combination rather
        # than a basis row, so inspect coefficients -2..2 of one or two rows.
        vectors = list(reduced)
        for left in range(len(reduced)):
            for right in range(left + 1, len(reduced)):
                for coefficient_left in (-2, -1, 1, 2):
                    for coefficient_right in (-2, -1, 1, 2):
                        vectors.append(
                            [
                                coefficient_left * reduced[left][coordinate]
                                + coefficient_right * reduced[right][coordinate]
                                for coordinate in range(dimension)
                            ]
                        )

        for vector in vectors:
            if vector[0] != 0 or vector[-1] == 0:
                continue
            affine_scale = vector[-1]
            coefficient_coordinates = vector[1:-1]
            if any(coordinate % affine_scale for coordinate in coefficient_coordinates):
                continue
            coefficients = tuple(coordinate // affine_scale for coordinate in coefficient_coordinates)
            if max(map(abs, coefficients), default=0) <= POLYNOMIAL_COEFFICIENT_BOUND:
                candidates.add(coefficients)
    return candidates


def search_polynomials(
    rows: Sequence[PrimeRow], power: int
) -> tuple[list[Fit], int, int]:
    sparse_candidates, tested = exhaustive_sparse_polynomials(rows, power)
    lll_candidates = dense_lll_polynomials(rows, power)
    candidates = sparse_candidates | lll_candidates
    fits = [
        polynomial_fit(rows, power, coefficients, "polynomial")
        for coefficients in candidates
    ]
    qualifying = deduplicate_fits(
        fit for fit in fits if fit.hits >= POLYNOMIAL_PROBE_COUNT
    )
    return qualifying, tested, len(lll_candidates)


def centered(value: int, p: int) -> int:
    value %= p
    return value if value <= p // 2 else value - p


def quartic_factor_families(
    rows: Sequence[PrimeRow], include_fourth_roots: bool = True
) -> list[tuple[str, tuple[int, ...]]]:
    simple_factors: list[tuple[str, tuple[int, ...]]] = []
    simple_generators = (
        ("eps16", tuple(1 if row.p % 16 == 1 else -1 for row in rows)),
        ("q2", tuple(centered(row.q2, row.p) for row in rows)),
        ("q3", tuple(centered(row.q3, row.p) for row in rows)),
    )
    for exponents in itertools.product((0, 1), repeat=3):
        names = [name for exponent, (name, _) in zip(exponents, simple_generators) if exponent]
        name = "*".join(names) if names else "1"
        values = tuple(
            math.prod(
                generator_values[row_index]
                for exponent, (_, generator_values) in zip(exponents, simple_generators)
                if exponent
            )
            % row.p
            for row_index, row in enumerate(rows)
        )
        simple_factors.append((name, values))
    factors = list(simple_factors)
    if include_fourth_roots:
        for base_index, base in enumerate(QUARTIC_SMALL_BASES):
            raw_values = tuple(row.quartic_small[base_index] for row in rows)
            for sign_name, sign_values in simple_factors:
                name = (
                    f"chi4({base})"
                    if sign_name == "1"
                    else f"{sign_name}*chi4({base})"
                )
                factors.append(
                    (
                        name,
                        tuple(
                            sign * raw % row.p
                            for sign, raw, row in zip(sign_values, raw_values, rows)
                        ),
                    )
                )

    # Keep one readable name for each empirically identical factor.
    unique: dict[tuple[int, ...], tuple[str, tuple[int, ...]]] = {}
    for name, values in factors:
        unique.setdefault(values, (name, values))
    return list(unique.values())


def linear_forms() -> tuple[tuple[str, tuple[int, int, int, int]], ...]:
    forms = []
    variables = ("x", "y", "a", "b")
    for coefficients in itertools.product(range(-2, 3), repeat=4):
        if not any(coefficients):
            continue
        nonzero = [(coefficient, name) for coefficient, name in zip(coefficients, variables) if coefficient]
        pieces = [f"{coefficient:+d}*{name}" for coefficient, name in nonzero]
        name = " ".join(pieces)
        if name.startswith("+"):
            name = name[1:]
        forms.append((name, coefficients))
    return tuple(forms)


def linear_value(row: PrimeRow, coefficients: tuple[int, int, int, int]) -> int:
    return sum(
        coefficient * value
        for coefficient, value in zip(coefficients, (row.x, row.y, row.a, row.b))
    ) % row.p


def search_quartic_twists(rows: Sequence[PrimeRow]) -> tuple[list[Fit], list[Fit]]:
    factors = quartic_factor_families(rows)
    forms = linear_forms()
    all_fits = []
    threshold_fits = []
    targets = {
        "v/G1": tuple(row.v * pow(row.g1, -1, row.p) % row.p for row in rows),
        "v*G1": tuple(row.v * row.g1 % row.p for row in rows),
    }
    for target_name, target_values in targets.items():
        for form_name, coefficients in forms:
            form_values = tuple(linear_value(row, coefficients) for row in rows)
            for factor_name, factor_values in factors:
                prediction = tuple(
                    form * factor % row.p
                    for form, factor, row in zip(form_values, factor_values, rows)
                )
                failures = tuple(
                    row.p
                    for row, target, predicted in zip(rows, target_values, prediction)
                    if target != predicted
                )
                fit = Fit(
                    "quartic twist",
                    f"{target_name} = ({form_name}) * {factor_name}",
                    len(rows) - len(failures),
                    failures,
                    prediction,
                )
                all_fits.append(fit)
                if fit.hits >= MONOMIAL_HIT_THRESHOLD:
                    threshold_fits.append(fit)
    return deduplicate_fits(threshold_fits), deduplicate_fits(all_fits)


def distribution_after_quartic_sign(rows: Sequence[PrimeRow]) -> dict[str, object]:
    sign_factors = quartic_factor_families(rows, include_fourth_roots=False)
    scored = []
    for name, values in sign_factors:
        normalized = []
        for row, factor in zip(rows, values):
            residue = row.v * pow(row.g1 * factor % row.p, -1, row.p) % row.p
            normalized.append(centered(residue, row.p) / row.p)
        mean = statistics.fmean(normalized)
        scored.append((abs(mean), mean, name, values, normalized))
    _, mean, name, values, normalized = max(scored, key=lambda item: (item[0], item[2]))
    orientation = 1 if mean >= 0 else -1
    if orientation == -1:
        normalized = [-value for value in normalized]
        name = f"-({name})"

    bins = [0] * 10
    for value in normalized:
        index = min(9, max(0, int((value + 0.5) * 10)))
        bins[index] += 1
    ordered = sorted(normalized)

    def empirical_quantile(fraction: Fraction) -> float:
        index = (len(ordered) - 1) * fraction.numerator // fraction.denominator
        return ordered[index]

    return {
        "guess": name,
        "selection": "maximal absolute empirical mean among the 8 Walsh sign factors",
        "mean": statistics.fmean(normalized),
        "stdev": statistics.pstdev(normalized),
        "minimum": min(normalized),
        "q1": empirical_quantile(Fraction(1, 4)),
        "median": empirical_quantile(Fraction(1, 2)),
        "q3": empirical_quantile(Fraction(3, 4)),
        "maximum": max(normalized),
        "bins": tuple(bins),
    }


def dataset_digest(rows: Sequence[PrimeRow]) -> str:
    fields = []
    for row in rows:
        fields.append(
            ",".join(
                map(
                    str,
                    (
                        row.p,
                        row.v,
                        row.x,
                        row.y,
                        row.a,
                        row.b,
                        row.c,
                        row.d,
                        row.g1,
                        row.g2,
                        *row.phi,
                        row.q2,
                        row.q3,
                        row.q6,
                        row.tau_midpoint,
                    ),
                )
            )
        )
    return hashlib.sha256(("\n".join(fields) + "\n").encode()).hexdigest()


def dump_csv(rows: Sequence[PrimeRow]) -> None:
    writer = csv.writer(sys.stdout, lineterminator="\n")
    writer.writerow(
        (
            "p",
            "v",
            "x",
            "y",
            "a",
            "b",
            "c",
            "d",
            "G1",
            "G2",
            *(f"phi_{parameter}" for parameter in PHI_PARAMETERS),
            "quartic_2",
            "quartic_3",
            "quartic_6",
            "tau_midpoint",
        )
    )
    for row in rows:
        writer.writerow(
            (
                row.p,
                row.v,
                row.x,
                row.y,
                row.a,
                row.b,
                row.c,
                row.d,
                row.g1,
                row.g2,
                *row.phi,
                centered(row.q2, row.p),
                centered(row.q3, row.p),
                centered(row.q6, row.p),
                row.tau_midpoint,
            )
        )


def print_fit(label: str, fit: Fit) -> None:
    print(f"{label}: {fit.formula}")
    print(f"  hits/misses = {fit.hits}/{len(fit.failures)}")
    print(f"  failure primes = {list(fit.failures)}")


def run_analysis(rows: Sequence[PrimeRow]) -> None:
    print("DATA")
    print(f"  primes p == 1 (mod 24), p < {LIMIT}: {len(rows)}")
    print(f"  first/last prime: {rows[0].p}/{rows[-1].p}")
    print(f"  canonical data SHA-256: {dataset_digest(rows)}")
    print("  direct sqrt(Apery) cross-check: 4/4 primes below 250")
    print(f"  midpoint tau_((p-1)/2)=(-2|p): {len(rows)}/{len(rows)}")
    print(f"  Gauss normalization G1=2a: {len(rows)}/{len(rows)}")
    for name in ("q2", "q3", "q6"):
        counts = Counter(centered(getattr(row, name), row.p) for row in rows)
        print(f"  quartic {name[1:]} distribution: {dict(sorted(counts.items()))}")

    print("\nBOUNDED-RATIONAL MONOMIAL SEARCH")
    monomial_threshold, monomial_best = search_scaled_monomials(rows)
    print(
        f"  searched {len(MONOMIALS)} monomials * 3 G1 exponents; "
        f"|num|,den <= {RATIONAL_BOUND}"
    )
    print(f"  distinct candidates with >= {MONOMIAL_HIT_THRESHOLD} hits: {len(monomial_threshold)}")
    if monomial_best:
        print(f"  best support: {monomial_best[0].hits}/{len(rows)} for {monomial_best[0].formula}")
    for index, fit in enumerate(monomial_threshold, 1):
        print_fit(f"  qualifying {index}", fit)

    print("\nSMALL-COEFFICIENT POLYNOMIAL SEARCH")
    polynomial_fits = []
    for power in (2, 4):
        qualifying, tested, lll_count = search_polynomials(rows, power)
        polynomial_fits.extend(qualifying)
        print(
            f"  v^{power}: {tested} sparse coefficient tuples screened on "
            f"{POLYNOMIAL_PROBE_COUNT} primes; {lll_count} bounded LLL candidates; "
            f"{len(qualifying)} candidates retained with >= {POLYNOMIAL_PROBE_COUNT} hits"
        )
        for index, fit in enumerate(qualifying, 1):
            print_fit(f"    v^{power} qualifying {index}", fit)

    print("\nQUARTIC-TWISTED LINEAR SEARCH")
    quartic_threshold, quartic_best = search_quartic_twists(rows)
    print(
        f"  searched {len(linear_forms())} linear forms for each of v/G1 and v*G1 "
        f"against {len(quartic_factor_families(rows))} distinct quartic factors"
    )
    print(f"  distinct candidates with >= {MONOMIAL_HIT_THRESHOLD} hits: {len(quartic_threshold)}")
    if quartic_best:
        print(f"  best support: {quartic_best[0].hits}/{len(rows)} for {quartic_best[0].formula}")
    for index, fit in enumerate(quartic_threshold, 1):
        print_fit(f"  qualifying {index}", fit)

    complete = deduplicate_fits(
        fit
        for fit in monomial_threshold + polynomial_fits + quartic_threshold
        if fit.hits == len(rows)
    )
    print("\nVERDICT")
    print(f"  full identities found: {len(complete)}")
    for index, fit in enumerate(complete, 1):
        print_fit(f"  identity {index}", fit)

    if not complete:
        # Dense LLL candidates are deliberately fitted on ten probes and are
        # already reported above.  For theory guidance, rank the a-priori
        # monomial/twist families and require distinct exact support sets.
        near_misses = deduplicate_supports(monomial_best + quartic_best)[:3]
        print("  three best a-priori near-misses with distinct support sets:")
        for index, fit in enumerate(near_misses, 1):
            print_fit(f"    near-miss {index}", fit)

    distribution = distribution_after_quartic_sign(rows)
    print("\nEMPIRICAL DISTRIBUTION OF v/(2a*quartic-sign-guess)")
    print(f"  sign guess: {distribution['guess']}")
    print(f"  selection rule: {distribution['selection']}")
    print(
        "  min/q1/median/q3/max = "
        f"{distribution['minimum']:.6f}, {distribution['q1']:.6f}, "
        f"{distribution['median']:.6f}, {distribution['q3']:.6f}, "
        f"{distribution['maximum']:.6f}"
    )
    print(f"  mean/stdev = {distribution['mean']:.6f}/{distribution['stdev']:.6f}")
    print(f"  10 bins on [-1/2,1/2): {list(distribution['bins'])}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--dump-csv",
        action="store_true",
        help="emit all computed prime rows as CSV instead of running the fits",
    )
    arguments = parser.parse_args()
    rows = build_rows()
    if arguments.dump_csv:
        dump_csv(rows)
    else:
        run_analysis(rows)


if __name__ == "__main__":
    main()
