#!/usr/bin/env python3
"""Round 4, Task A: exhaustive corrected-boundary recurrence test.

For p in {31, 73, 97}, this script computes the zeta(3) Apery numbers

    b_n = sum_{k=0}^n binom(n,k)^2 binom(n+k,k)^2

modulo p, labels the conifold roots using the least integer square root of 2,
and exhausts every (c1, d1) in F_p^2 in the proposed correction

    B_j = (c1 + d1*j)t1^(-j) + (c1 - d1 - d1*j)t2^(-j),
    R_j = b_j + B_j.

The score is the number of 1 <= j <= p-2 at which R satisfies the Apery
recurrence.  Results are written to /tmp/round4_task_a.txt by default.

The exhaustive search includes (0, 0), as the task requests.  Since b itself
satisfies the recurrence, (0, 0) has the full score p-2 a priori.  To make the
search diagnostic, the report also gives the best nonzero pair and the best
pair after fixing c1=H_p(t1).
"""

from __future__ import annotations

import argparse
from collections import Counter
from math import comb, factorial, isqrt
from pathlib import Path
from typing import Iterable, Sequence


PRIMES = (31, 73, 97)
DEFAULT_OUTPUT = Path("/tmp/round4_task_a.txt")


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    return all(n % divisor for divisor in range(2, isqrt(n) + 1))


def apery_binomial_mod_p(n: int, p: int) -> int:
    """Compute b_n modulo p directly from its defining binomial sum."""
    return sum(
        comb(n, k) ** 2 * comb(n + k, k) ** 2 for k in range(n + 1)
    ) % p


def recurrence_polynomial(n: int) -> int:
    """P(n)=34*n^3+51*n^2+27*n+5 over the integers."""
    return ((34 * n + 51) * n + 27) * n + 5


def apery_recurrence_mod_p(p: int) -> list[int]:
    """Return b_0,...,b_{p-1} modulo p from the Apery recurrence."""
    assert is_prime(p) and p > 2
    b = [0] * p
    b[0] = 1
    b[1] = 5 % p
    for n in range(1, p - 1):
        numerator = (
            recurrence_polynomial(n) * b[n] - n**3 * b[n - 1]
        ) % p
        denominator = pow(n + 1, 3, p)
        b[n + 1] = numerator * pow(denominator, -1, p) % p
    return b


def evaluate_polynomial(coefficients: Sequence[int], t: int, p: int) -> int:
    value = 0
    for coefficient in reversed(coefficients):
        value = (value * t + coefficient) % p
    return value


def conifold_roots(p: int) -> tuple[tuple[int, int], int, int]:
    """Return ((sqrt2_small, sqrt2_large), t1, t2).

    The canonical labeling uses the least integer representative r of a square
    root of 2 and sets t1=17-12*r, t2=17+12*r.  The other square root swaps t1
    and t2.  There is no intrinsic ordering of elements of F_p.
    """
    square_roots = tuple(r for r in range(p) if r * r % p == 2)
    assert len(square_roots) == 2
    small, large = square_roots
    assert small < large and (small + large) % p == 0
    t1 = (17 - 12 * small) % p
    t2 = (17 + 12 * small) % p
    assert t1 != t2
    assert t1 * t2 % p == 1
    assert (t1 * t1 - 34 * t1 + 1) % p == 0
    assert (t2 * t2 - 34 * t2 + 1) % p == 0
    return (small, large), t1, t2


def correction_vector(
    p: int, t1: int, t2: int, c1: int, d1: int
) -> list[int]:
    return [
        (
            (c1 + d1 * j) * pow(t1, -j, p)
            + (c1 - d1 - d1 * j) * pow(t2, -j, p)
        )
        % p
        for j in range(p)
    ]


def corrected_residual(
    b: Sequence[int], p: int, t1: int, t2: int, c1: int, d1: int
) -> list[int]:
    correction = correction_vector(p, t1, t2, c1, d1)
    return [(b[j] + correction[j]) % p for j in range(p)]


def recurrence_residuals(values: Sequence[int], p: int) -> list[int]:
    """Return recurrence defects at j=1,...,p-2."""
    assert len(values) == p
    return [
        (
            (j + 1) ** 3 * values[j + 1]
            - recurrence_polynomial(j) * values[j]
            + j**3 * values[j - 1]
        )
        % p
        for j in range(1, p - 1)
    ]


def score(values: Sequence[int], p: int) -> int:
    return sum(defect == 0 for defect in recurrence_residuals(values, p))


def chunks(items: Sequence[str], size: int) -> Iterable[str]:
    for start in range(0, len(items), size):
        yield " ".join(items[start : start + size])


def vector_lines(label: str, values: Sequence[int], width: int = 10) -> list[str]:
    entries = [f"{j}:{value}" for j, value in enumerate(values)]
    return [label, *(f"  {line}" for line in chunks(entries, width))]


def pair_lines(label: str, pairs: Sequence[tuple[int, int]], width: int = 8) -> list[str]:
    entries = [f"({c},{d})" for c, d in pairs]
    return [f"{label} ({len(pairs)} pair(s)):", *(f"  {line}" for line in chunks(entries, width))]


def swap_label_pair(pair: tuple[int, int], p: int) -> tuple[int, int]:
    """Transform parameters when (t1,t2) is relabeled as (t2,t1)."""
    c1, d1 = pair
    return ((c1 - d1) % p, (-d1) % p)


def projective_parameter_direction(
    coefficient_c: int, coefficient_d: int, p: int
) -> tuple[int, int]:
    """Normalize the line coefficient_c*c+coefficient_d*d=0 in F_p^2."""
    assert coefficient_c or coefficient_d
    # A spanning vector for the solution line is (coefficient_d,-coefficient_c).
    vector_c = coefficient_d % p
    vector_d = (-coefficient_c) % p
    if vector_c:
        inverse = pow(vector_c, -1, p)
        return (1, vector_d * inverse % p)
    assert vector_d
    return (0, 1)


def occupancy_maximum_tail(number_bins: int, number_balls: int, observed: int) -> float:
    """P(max occupancy >= observed) for iid uniform balls in bins.

    Uses
      P(max < observed) = n!/m^n [x^n](sum_{r=0}^{observed-1}x^r/r!)^m.
    The largest case here is only (m,n)=(98,47).
    """
    polynomial = [1.0] + [0.0] * number_balls
    factor = [1.0 / factorial(r) for r in range(observed)]
    for _ in range(number_bins):
        product = [0.0] * (number_balls + 1)
        for degree, coefficient in enumerate(polynomial):
            if coefficient == 0.0:
                continue
            for shift, factor_coefficient in enumerate(factor):
                if degree + shift <= number_balls:
                    product[degree + shift] += coefficient * factor_coefficient
        polynomial = product
    probability_below = (
        factorial(number_balls)
        * polynomial[number_balls]
        / number_bins**number_balls
    )
    return 1.0 - probability_below


def analyze_prime(p: int) -> dict:
    b = apery_recurrence_mod_p(p)
    b_direct = [apery_binomial_mod_p(j, p) for j in range(p)]
    assert b == b_direct, f"binomial/recurrence mismatch for p={p}"
    assert recurrence_residuals(b, p) == [0] * (p - 2)
    assert all(b[j] == b[p - 1 - j] for j in range(p))

    sqrt2_roots, t1, t2 = conifold_roots(p)
    h1 = evaluate_polynomial(b, t1, p)
    h2 = evaluate_polynomial(b, t2, p)
    assert h1 == h2  # follows here also from palindromy and t1*t2=1
    assert h1 != 0

    scores: dict[tuple[int, int], int] = {}
    score_histogram: Counter[int] = Counter()
    for c1 in range(p):
        for d1 in range(p):
            residual = corrected_residual(b, p, t1, t2, c1, d1)
            pair_score = score(residual, p)
            scores[(c1, d1)] = pair_score
            score_histogram[pair_score] += 1

    maximum = max(scores.values())
    maximizers = sorted(pair for pair, value in scores.items() if value == maximum)
    # This is an important sanity check: the uncorrected Apery sequence is an
    # exact solution, so the zero boundary term must attain the full score.
    assert scores[(0, 0)] == p - 2
    assert maximum == p - 2

    nonzero_maximum = max(
        value for pair, value in scores.items() if pair != (0, 0)
    )
    nonzero_maximizers = sorted(
        pair
        for pair, value in scores.items()
        if pair != (0, 0) and value == nonzero_maximum
    )

    fixed_c_maximum = max(scores[(h1, d1)] for d1 in range(p))
    fixed_c_maximizers = [
        (h1, d1) for d1 in range(p) if scores[(h1, d1)] == fixed_c_maximum
    ]
    fixed_c_r_zeros = {
        pair: [
            j
            for j, value in enumerate(
                corrected_residual(b, p, t1, t2, *pair)
            )
            if value == 0
        ]
        for pair in fixed_c_maximizers
    }
    naive_pair = (h1, 0)

    # The recurrence defect is linear and homogeneous in (c1,d1), because b
    # already has zero defect.  Locate positions at which it vanishes for every
    # pair, and use this to check the exact mean score over F_p^2.
    defect_c = recurrence_residuals(
        correction_vector(p, t1, t2, 1, 0), p
    )
    defect_d = recurrence_residuals(
        correction_vector(p, t1, t2, 0, 1), p
    )
    for (c1, d1), direct_score in scores.items():
        linear_score = sum(
            (c1 * coefficient_c + d1 * coefficient_d) % p == 0
            for coefficient_c, coefficient_d in zip(defect_c, defect_d)
        )
        assert direct_score == linear_score
    universal_positions = [
        j
        for j, (coefficient_c, coefficient_d) in enumerate(
            zip(defect_c, defect_d), start=1
        )
        if coefficient_c == coefficient_d == 0
    ]
    universal_count = len(universal_positions)
    score_sum = sum(value * count for value, count in score_histogram.items())
    expected_numerator = universal_count * p + (p - 2 - universal_count)
    # expected_numerator/p is the exact all-pair mean.
    assert score_sum * p == expected_numerator * p * p

    # Reflection pairs all noncentral recurrence positions.  Each pair votes
    # for one of the p+1 projective parameter directions at which its defect
    # vanishes.  This supplies a multiple-search-aware random comparison.
    direction_counts: Counter[tuple[int, int]] = Counter()
    position_directions: dict[int, tuple[int, int]] = {}
    for j, (coefficient_c, coefficient_d) in enumerate(
        zip(defect_c, defect_d), start=1
    ):
        if coefficient_c == coefficient_d == 0:
            continue
        direction = projective_parameter_direction(
            coefficient_c, coefficient_d, p
        )
        direction_counts[direction] += 1
        position_directions[j] = direction
    assert universal_positions == [(p - 1) // 2]
    assert all(
        position_directions[j] == position_directions[p - 1 - j]
        for j in position_directions
    )
    assert all(count % 2 == 0 for count in direction_counts.values())
    maximum_occupancy = max(count // 2 for count in direction_counts.values())
    assert nonzero_maximum == 1 + 2 * maximum_occupancy
    projective_tail_probability = occupancy_maximum_tail(
        number_bins=p + 1,
        number_balls=(p - 3) // 2,
        observed=maximum_occupancy,
    )

    # Scaling a nonzero parameter pair does not change which homogeneous
    # defects vanish.  Check this projective invariance throughout F_p^2.
    for pair, pair_score in scores.items():
        if pair == (0, 0):
            continue
        c1, d1 = pair
        for scalar in (2, p - 1):
            scaled = (scalar * c1 % p, scalar * d1 % p)
            assert scores[scaled] == pair_score

    selected_best = maximizers[0]
    best_residual = corrected_residual(b, p, t1, t2, *selected_best)
    best_zero_positions = [j for j, value in enumerate(best_residual) if value == 0]

    # Relabeling roots must preserve B and R after (c,d)->(c-d,-d).
    sample_pairs = {(0, 0), naive_pair, nonzero_maximizers[0]}
    for pair in sample_pairs:
        swapped_pair = swap_label_pair(pair, p)
        original = corrected_residual(b, p, t1, t2, *pair)
        swapped = corrected_residual(b, p, t2, t1, *swapped_pair)
        assert original == swapped

    return {
        "p": p,
        "b": b,
        "sqrt2_roots": sqrt2_roots,
        "t1": t1,
        "t2": t2,
        "h1": h1,
        "h2": h2,
        "b_zeros": [j for j, value in enumerate(b) if value == 0],
        "scores": scores,
        "score_histogram": score_histogram,
        "maximum": maximum,
        "maximizers": maximizers,
        "nonzero_maximum": nonzero_maximum,
        "nonzero_maximizers": nonzero_maximizers,
        "fixed_c_maximum": fixed_c_maximum,
        "fixed_c_maximizers": fixed_c_maximizers,
        "fixed_c_r_zeros": fixed_c_r_zeros,
        "naive_pair": naive_pair,
        "naive_score": scores[naive_pair],
        "universal_positions": universal_positions,
        "expected_numerator": expected_numerator,
        "maximum_occupancy": maximum_occupancy,
        "projective_tail_probability": projective_tail_probability,
        "selected_best": selected_best,
        "best_residual": best_residual,
        "best_zero_positions": best_zero_positions,
    }


def render_report(rows: Sequence[dict]) -> str:
    lines = [
        "Round 4, Task A: corrected-boundary exhaustive recurrence test",
        "=" * 72,
        "",
        "Verified definitions and conventions",
        "------------------------------------",
        "b_n = sum_{k=0}^n binom(n,k)^2 binom(n+k,k)^2.",
        "b_0=1, b_1=5, and",
        "  (n+1)^3 b_(n+1) - P(n)b_n + n^3 b_(n-1) = 0,",
        "  P(n)=34n^3+51n^2+27n+5.",
        "Every b_j below was computed by this recurrence and independently",
        "checked against the defining binomial sum.",
        "",
        "The phrase 'smaller sqrt(2)' is interpreted using least nonnegative",
        "integer representatives: choose the smaller r in {0,...,p-1} with",
        "r^2=2, then t1=17-12r and t2=17+12r.  Finite fields have no",
        "intrinsic order.  Choosing the other square root swaps t1,t2 and sends",
        "  (c1,d1) -> (c1-d1,-d1),",
        "without changing B(j), R(j), or any score.",
        "",
        "Logical warning",
        "---------------",
        "The global search includes (c1,d1)=(0,0).  At that pair R=b, and b",
        "satisfies the scored recurrence by definition.  Therefore the global",
        "maximum is forced to be p-2 before any boundary hypothesis is tested.",
        "The nonzero-pair and fixed-c1 results below are the informative parts.",
    ]

    for row in rows:
        p = row["p"]
        small, large = row["sqrt2_roots"]
        nominal = (p - 2) / p
        exact_mean = row["expected_numerator"] / p
        nonzero_mean = (
            len(row["universal_positions"])
            + (p - 2 - len(row["universal_positions"])) / (p + 1)
        )
        alternate_maximizers = sorted(
            swap_label_pair(pair, p) for pair in row["maximizers"]
        )
        alternate_nonzero_maximizers = sorted(
            swap_label_pair(pair, p) for pair in row["nonzero_maximizers"]
        )

        lines.extend(
            [
                "",
                f"p={p}",
                "-" * len(f"p={p}"),
                f"sqrt(2) representatives: {small}, {large}; chosen r={small}",
                f"canonical roots: t1={row['t1']}, t2={row['t2']}",
                f"alternate labeling: t1={row['t2']}, t2={row['t1']}",
                f"checks: t1*t2 mod p={(row['t1'] * row['t2']) % p}; "
                f"H_p(t1)={row['h1']}, H_p(t2)={row['h2']}",
                f"zeros of b: {row['b_zeros']}",
                "",
                f"Global maximum score: {row['maximum']}/{p-2}",
                f"Nominal single-random-pair baseline: (p-2)/p={nominal:.9f}",
                f"Exact mean over all p^2 pairs: {exact_mean:.9f}",
                f"Here the exact mean is 1+(p-3)/p, not (p-2)/p, because",
                "the center position is forced by reflection symmetry.",
                f"Universal recurrence positions for every pair: "
                f"{row['universal_positions'] or 'none'}",
            ]
        )
        lines.extend(pair_lines("All global maximizers, canonical labeling", row["maximizers"]))
        lines.extend(pair_lines("All global maximizers, alternate labeling", alternate_maximizers))
        lines.extend(
            [
                "The global maximum is numerically far above the single-pair",
                "baseline, but this is tautological: it contains the zero correction.",
                "",
                f"Best score with (c1,d1)!=(0,0): "
                f"{row['nonzero_maximum']}/{p-2}",
                f"Exact mean over nonzero pairs: {nonzero_mean:.9f}",
                f"After reflection pairing this is max occupancy "
                f"{row['maximum_occupancy']} among p+1 projective directions.",
                f"Independent-uniform occupancy null: "
                f"Pr(max >= observed)={row['projective_tail_probability']:.9f}.",
                "Thus the nonzero peak is not significant under this",
                "multiple-search-aware benchmark.",
            ]
        )
        lines.extend(
            pair_lines(
                "All nonzero maximizers, canonical labeling",
                row["nonzero_maximizers"],
            )
        )
        lines.extend(
            pair_lines(
                "All nonzero maximizers, alternate labeling",
                alternate_nonzero_maximizers,
            )
        )
        lines.extend(
            [
                "",
                f"With c1 fixed to H_p(t1)={row['h1']}: best score "
                f"{row['fixed_c_maximum']}/{p-2}",
            ]
        )
        lines.extend(
            pair_lines(
                "All fixed-c1 maximizers", row["fixed_c_maximizers"]
            )
        )
        lines.append("R-zero audit for the fixed-c1 maximizers:")
        lines.extend(
            f"  {pair}: {zeros if zeros else 'none (R is nonzero everywhere)'}"
            for pair, zeros in row["fixed_c_r_zeros"].items()
        )
        lines.extend(
            [
                f"Naive pair (H_p(t1),0)={row['naive_pair']}: "
                f"score {row['naive_score']}/{p-2}",
                "",
                "Score histogram (score:number of parameter pairs):",
                "  "
                + " ".join(
                    f"{value}:{row['score_histogram'][value]}"
                    for value in sorted(row["score_histogram"])
                ),
                "",
                f"Selected best pair for the requested R table: "
                f"{row['selected_best']}",
            ]
        )
        if row["selected_best"] == (0, 0):
            lines.append("At this pair R(j)=b_j, so one combined table suffices.")
        lines.extend(vector_lines("b(j)=R_best(j) (j:value):", row["best_residual"]))
        lines.extend(
            [
                f"R_best nonzero everywhere: "
                f"{'YES' if not row['best_zero_positions'] else 'NO'}",
                f"zero positions of R_best: {row['best_zero_positions']}",
            ]
        )

    lines.extend(
        [
            "",
            "Overall verdict",
            "---------------",
            "For all three primes, the requested unrestricted sweep is dominated",
            "by the forced zero-boundary solution (0,0), whose R=b has full score",
            "and retains the known paired zeros.  Thus its R is not nonzero",
            "everywhere.  A full recurrence score in this sweep does not identify",
            "a nontrivial corrected trace term; the test needs an independent",
            "normalization excluding (0,0), for example a prescribed nonzero c1.",
        ]
    )
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    rows = [analyze_prime(p) for p in PRIMES]
    report = render_report(rows)
    args.output.write_text(report, encoding="utf-8")
    print(f"Wrote {args.output} ({args.output.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
