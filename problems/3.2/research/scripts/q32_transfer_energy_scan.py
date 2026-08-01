#!/usr/bin/env python3
"""Exact mod-p collision-energy scan for the Apéry zeta(3) numbers.

The scanned multiset is ``{b_r mod p : 1 <= r <= p-2}``.  Its off-diagonal
energy counts ordered pairs.  Gap counts use one orientation, ``r < s``, so
the reflection contribution is one pair per two-cycle of ``r -> p-1-r``.

Only the requested normalized ratios use floating point.  The recurrence,
histogram, energies, reflection separation, and gap counts are integer-exact.
"""

from __future__ import annotations

import argparse
from array import array
from collections import Counter
from dataclasses import dataclass
from functools import cmp_to_key
import math
import os
from pathlib import Path
import random
import sys
import time
from typing import Iterable, Sequence


SCRIPT_PATH = Path(__file__).resolve()
PROJECT_ROOT = SCRIPT_PATH.parents[2]
DEFAULT_PROGRESS = SCRIPT_PATH.with_name("te_scan_progress.txt")
DEFAULT_REPORT = PROJECT_ROOT / "CODEX_TE_SCAN.md"


@dataclass
class PrimeRecord:
    p: int
    energy: int
    reflection_energy: int
    nonforced_energy: int
    zero_count: int
    max_multiplicity: int
    popular_value: int
    energy_over_p: float
    energy_over_p_log2: float
    energy_over_p_125: float
    nonforced_over_p: float
    nonforced_over_p_log2: float
    nonforced_over_p_125: float


@dataclass
class Inspection:
    p: int
    top_nonreflection_gaps: tuple[tuple[int, int], ...]
    max_nonreflection_gap: int
    max_nonreflection_gap_count: int
    nonforced_unordered_pairs: int
    danger: bool
    popular_positions: tuple[int, ...]
    popular_reflection_pairs: int
    popular_central_fixed_point: bool
    longest_ap_length: int
    longest_ap_step: int
    adjacent_gap_mode: int
    adjacent_gap_mode_count: int


def primes_in_interval(lower: int, upper: int) -> list[int]:
    """Return all primes in the inclusive interval [lower, upper]."""
    if upper < 2 or lower > upper:
        return []
    sieve = bytearray(b"\x01") * (upper + 1)
    sieve[0:2] = b"\x00\x00"
    for q in range(2, math.isqrt(upper) + 1):
        if sieve[q]:
            start = q * q
            sieve[start : upper + 1 : q] = b"\x00" * (
                (upper - start) // q + 1
            )
    return [p for p in range(max(2, lower), upper + 1) if sieve[p]]


def normalized(value: int, p: int) -> tuple[float, float, float]:
    log_p = math.log(p)
    return value / p, value / (p * log_p * log_p), value / (p ** 1.25)


def apery_energy(p: int) -> tuple[list[int], list[int], int, int, int]:
    """Compute b_0,...,b_{p-2}, its histogram, and exact energy data."""
    inverses = [0] * p
    inverses[1] = 1
    for i in range(2, p):
        inverses[i] = p - (p // i) * inverses[p % i] % p

    values = [0] * (p - 1)
    values[0] = 1
    values[1] = 5 % p
    counts = [0] * p
    counts[values[1]] = 1
    energy = 0
    previous = values[0]
    current = values[1]

    # n=1,...,p-3 produces b_2,...,b_{p-2}.
    for n in range(1, p - 2):
        n2 = n * n % p
        n3 = n2 * n % p
        coefficient = (34 * n3 + 51 * n2 + 27 * n + 5) % p
        numerator = (coefficient * current - n3 * previous) % p
        inverse = inverses[n + 1]
        inverse_cube = inverse * inverse % p * inverse % p
        next_value = numerator * inverse_cube % p
        values[n + 1] = next_value

        old_count = counts[next_value]
        energy += 2 * old_count
        counts[next_value] = old_count + 1
        previous, current = current, next_value

    max_multiplicity = max(counts)
    popular_value = counts.index(max_multiplicity)
    return values, counts, energy, max_multiplicity, popular_value


def sanity_gate(
    p: int, values: Sequence[int], counts: Sequence[int], sample_number: int
) -> str:
    """Run the required histogram-zero and three-point reflection checks."""
    histogram_zeros = counts[0]
    direct_zeros = sum(value == 0 for value in values[1:])
    if histogram_zeros != direct_zeros:
        raise AssertionError(
            f"p={p}: histogram zero count {histogram_zeros} != {direct_zeros}"
        )

    rng = random.Random((p << 20) ^ sample_number ^ 0x5A17E2)
    sampled_positions = rng.sample(range(1, p - 1), 3)
    for r in sampled_positions:
        reflected = p - 1 - r
        if values[r] != values[reflected]:
            raise AssertionError(
                f"p={p}: reflection failed at r={r}, reflected={reflected}"
            )
    return (
        f"SANITY prime_index={sample_number} p={p} "
        f"zeros={histogram_zeros} direct={direct_zeros} reflection=3/3 PASS"
    )


def longest_arithmetic_progression(positions: Sequence[int]) -> tuple[int, int]:
    """Return the length and step of a longest positive-step AP subsequence."""
    size = len(positions)
    if size < 2:
        return size, 0
    tables: list[dict[int, int]] = [{} for _ in positions]
    best_length, best_step = 2, positions[1] - positions[0]
    for right in range(size):
        table = tables[right]
        for left in range(right):
            step = positions[right] - positions[left]
            length = tables[left].get(step, 1) + 1
            if length > table.get(step, 0):
                table[step] = length
            if length > best_length or (
                length == best_length and step < best_step
            ):
                best_length, best_step = length, step
    return best_length, best_step


def inspect_prime(
    record: PrimeRecord, values: Sequence[int], counts: Sequence[int]
) -> Inspection:
    """Compute the separated gap spectrum and popular-fiber structure."""
    positions_by_value: dict[int, list[int]] = {
        value: [] for value, count in enumerate(counts) if count >= 2
    }
    for position in range(1, record.p - 1):
        value = values[position]
        if value in positions_by_value:
            positions_by_value[value].append(position)

    nonreflection_gaps: Counter[int] = Counter()
    forced_unordered = 0
    for positions in positions_by_value.values():
        for right_index in range(1, len(positions)):
            right = positions[right_index]
            for left in positions[:right_index]:
                if left + right == record.p - 1:
                    forced_unordered += 1
                else:
                    gap = right - left
                    # C_p(h) is directed and modular: the two orientations of
                    # an unordered pair contribute h and p-h.
                    nonreflection_gaps[gap] += 1
                    nonreflection_gaps[record.p - gap] += 1

    expected_forced = (record.p - 3) // 2
    if forced_unordered != expected_forced:
        raise AssertionError(
            f"p={record.p}: forced gap pairs {forced_unordered} != {expected_forced}"
        )
    nonforced_ordered = sum(nonreflection_gaps.values())
    if nonforced_ordered != record.nonforced_energy:
        raise AssertionError(
            f"p={record.p}: nonreflection gap energy mismatch"
        )
    nonforced_unordered = nonforced_ordered // 2

    top = tuple(
        sorted(nonreflection_gaps.items(), key=lambda item: (-item[1], item[0]))[
            :5
        ]
    )
    if top:
        max_gap, max_gap_count = top[0]
    else:
        max_gap, max_gap_count = 0, 0
    # count > p^(1/4), tested without a floating-point threshold.
    danger = max_gap_count**4 > record.p

    popular_positions = tuple(
        position
        for position in range(1, record.p - 1)
        if values[position] == record.popular_value
    )
    popular_set = set(popular_positions)
    popular_reflection_pairs = sum(
        1
        for position in popular_positions
        if position < record.p - 1 - position
        and record.p - 1 - position in popular_set
    )
    central = (record.p - 1) // 2 in popular_set
    ap_length, ap_step = longest_arithmetic_progression(popular_positions)
    adjacent_gaps = Counter(
        popular_positions[index + 1] - popular_positions[index]
        for index in range(len(popular_positions) - 1)
    )
    if adjacent_gaps:
        adjacent_gap, adjacent_count = min(
            adjacent_gaps.items(), key=lambda item: (-item[1], item[0])
        )
    else:
        adjacent_gap, adjacent_count = 0, 0

    return Inspection(
        p=record.p,
        top_nonreflection_gaps=top,
        max_nonreflection_gap=max_gap,
        max_nonreflection_gap_count=max_gap_count,
        nonforced_unordered_pairs=nonforced_unordered,
        danger=danger,
        popular_positions=popular_positions,
        popular_reflection_pairs=popular_reflection_pairs,
        popular_central_fixed_point=central,
        longest_ap_length=ap_length,
        longest_ap_step=ap_step,
        adjacent_gap_mode=adjacent_gap,
        adjacent_gap_mode_count=adjacent_count,
    )


def compare_energy_over_p(left: PrimeRecord, right: PrimeRecord) -> int:
    lhs = left.energy * right.p
    rhs = right.energy * left.p
    return (lhs > rhs) - (lhs < rhs)


def dyadic_lower(p: int) -> int:
    return 1 << (p.bit_length() - 1)


def ratio_max(
    records: Sequence[PrimeRecord], numerator: str, denominator: str
) -> PrimeRecord:
    def compare(left: PrimeRecord, right: PrimeRecord) -> int:
        lv = getattr(left, numerator)
        rv = getattr(right, numerator)
        if denominator == "p":
            lhs, rhs = lv * right.p, rv * left.p
        elif denominator == "p_log2":
            # This ratio is report-only; energy remains exact.
            lhs = lv / (left.p * math.log(left.p) ** 2)
            rhs = rv / (right.p * math.log(right.p) ** 2)
            return (lhs > rhs) - (lhs < rhs)
        elif denominator == "p_125":
            lhs = lv / (left.p ** 1.25)
            rhs = rv / (right.p ** 1.25)
            return (lhs > rhs) - (lhs < rhs)
        else:
            raise ValueError(denominator)
        return (lhs > rhs) - (lhs < rhs)

    return max(records, key=cmp_to_key(compare))


def regression_slope(blocks: Sequence[Sequence[PrimeRecord]]) -> tuple[float, int]:
    points: list[tuple[float, float]] = []
    for block in blocks:
        raw_max = max(block, key=lambda record: (record.energy, record.p))
        if raw_max.energy > 0:
            points.append((math.log(raw_max.p), math.log(raw_max.energy)))
    if len(points) < 2:
        return float("nan"), len(points)
    mean_x = sum(point[0] for point in points) / len(points)
    mean_y = sum(point[1] for point in points) / len(points)
    denominator = sum((point[0] - mean_x) ** 2 for point in points)
    numerator = sum(
        (point[0] - mean_x) * (point[1] - mean_y) for point in points
    )
    return numerator / denominator, len(points)


def complete_blocks(
    blocks: Sequence[Sequence[PrimeRecord]], minimum: int, maximum: int
) -> list[Sequence[PrimeRecord]]:
    """Select dyadic blocks not truncated by either requested endpoint."""
    result: list[Sequence[PrimeRecord]] = []
    for block in blocks:
        lower = dyadic_lower(block[0].p)
        if minimum <= lower + 1 and maximum >= 2 * lower:
            result.append(block)
    return result


def fmt_ratio(value: float) -> str:
    return f"{value:.8g}"


def fmt_gap_top(top: Sequence[tuple[int, int]]) -> str:
    if not top:
        return "none"
    return ", ".join(f"({gap}, {count})" for gap, count in top)


def fmt_positions(positions: Sequence[int]) -> str:
    if len(positions) <= 20:
        return ", ".join(map(str, positions))
    head = ", ".join(map(str, positions[:10]))
    tail = ", ".join(map(str, positions[-10:]))
    return f"{head}, ..., {tail}"


def build_report(
    records: Sequence[PrimeRecord],
    inspections: dict[int, Inspection],
    minimum: int,
    maximum: int,
    started_at: str,
    elapsed: float,
) -> str:
    blocks_by_lower: dict[int, list[PrimeRecord]] = {}
    for record in records:
        blocks_by_lower.setdefault(dyadic_lower(record.p), []).append(record)
    ordered_blocks = [blocks_by_lower[lower] for lower in sorted(blocks_by_lower)]
    fitted_blocks = complete_blocks(ordered_blocks, minimum, maximum)
    slope, regression_points = regression_slope(fitted_blocks)
    worst = sorted(records, key=cmp_to_key(compare_energy_over_p), reverse=True)[:5]
    gap_primes = sorted(inspections)
    largest_inspected_primes = {record.p for record in records[-50:]}
    envelope_violations = [
        record
        for record in records
        if record.energy > 100 * record.p * math.log(record.p) ** 2
    ]
    danger_primes = [p for p in gap_primes if inspections[p].danger]
    large_range_danger_primes = [
        p for p in danger_primes if p in largest_inspected_primes
    ]
    concentrated_danger_primes = [
        p
        for p in danger_primes
        if inspections[p].nonforced_unordered_pairs
        and inspections[p].max_nonreflection_gap_count
        / inspections[p].nonforced_unordered_pairs
        >= 0.25
    ]

    lines = [
        "# TE_{5/4} Apéry collision-energy scan",
        "",
        "## Coverage and conventions",
        "",
        f"- Exact range tested: every prime `p` in `[{minimum}, {maximum}]` "
        f"({len(records)} primes; first `{records[0].p}`, last `{records[-1].p}`).",
        f"- Run started: `{started_at}`; elapsed: `{elapsed:.1f}` seconds.",
        "- The histogram uses exactly `1 <= r <= p-2`. `E_off` counts ordered "
        "pairs `(r,s)` with `r != s`; all energy computations are integer-exact.",
        "- The involution `r -> p-1-r` has one fixed point and `(p-3)/2` "
        "two-cycles, so its ordered forced contribution is `E_refl=p-3`.",
        "- Gap counts are directed and modular: an unordered pair `r<s` adds "
        "one count at `h=s-r` and one at `p-h`. Both are removed as forced "
        "exactly when `r+s=p-1`.",
        "",
        "The incremental per-prime ledger, including all six requested "
        "normalizations, is `research/scripts/te_scan_progress.txt`.",
        "",
        "## Dyadic block maxima: total energy",
        "",
        "| block | coverage | primes | max E/p^(5/4) (p) | max E/p (p) | "
        "max E/(p log^2 p) (p) |",
        "|---|---:|---:|---:|---:|---:|",
    ]

    for block in ordered_blocks:
        lower = dyadic_lower(block[0].p)
        upper = 2 * lower
        max_125 = ratio_max(block, "energy", "p_125")
        max_p = ratio_max(block, "energy", "p")
        max_log = ratio_max(block, "energy", "p_log2")
        coverage = f"{block[0].p}-{block[-1].p}"
        lines.append(
            f"| ({lower}, {upper}] | {coverage} | {len(block)} | "
            f"{fmt_ratio(max_125.energy_over_p_125)} ({max_125.p}) | "
            f"{fmt_ratio(max_p.energy_over_p)} ({max_p.p}) | "
            f"{fmt_ratio(max_log.energy_over_p_log2)} ({max_log.p}) |"
        )

    lines.extend(
        [
            "",
            "The first and last rows are partial dyadic blocks when the stated "
            "scan endpoints cut them. The regression excludes those partial "
            "rows. Each complete block contributes "
            "the prime attaining the largest raw `E_off` in that block; ordinary "
            "least squares on `(log p, log E_off)` gives",
            "",
            f"- fitted slope: **{slope:.6f}** from {regression_points} complete-block maxima;",
            f"- danger threshold `slope > 1.25`: **{'TRIGGERED' if slope > 1.25 else 'not triggered'}**.",
            "",
            "## Reflection-separated block maxima",
            "",
            "| block | max (E-E_refl)/p^(5/4) (p) | max (E-E_refl)/p (p) | "
            "max (E-E_refl)/(p log^2 p) (p) |",
            "|---|---:|---:|---:|",
        ]
    )
    for block in ordered_blocks:
        lower = dyadic_lower(block[0].p)
        upper = 2 * lower
        max_125 = ratio_max(block, "nonforced_energy", "p_125")
        max_p = ratio_max(block, "nonforced_energy", "p")
        max_log = ratio_max(block, "nonforced_energy", "p_log2")
        lines.append(
            f"| ({lower}, {upper}] | "
            f"{fmt_ratio(max_125.nonforced_over_p_125)} ({max_125.p}) | "
            f"{fmt_ratio(max_p.nonforced_over_p)} ({max_p.p}) | "
            f"{fmt_ratio(max_log.nonforced_over_p_log2)} ({max_log.p}) |"
        )

    lines.extend(
        [
            "",
            "## Gap spectrum",
            "",
            "Inspected primes are the largest 50 tested primes together with the "
            "final five spike primes by `E_off/p`. The danger test is exact: "
            "`C_p(h)^4 > p`, equivalent to `C_p(h) > p^(1/4)`.",
            "",
            "| p | top five non-reflection (h, C_p(h)) | max > p^(1/4)? |",
            "|---:|---|:---:|",
        ]
    )
    for p in gap_primes:
        inspection = inspections[p]
        lines.append(
            f"| {p} | {fmt_gap_top(inspection.top_nonreflection_gaps)} | "
            f"{'YES' if inspection.danger else 'no'} |"
        )
    lines.extend(
        [
            "",
            f"Danger-signal primes: "
            f"{'none' if not danger_primes else ', '.join(map(str, danger_primes))}.",
            f"Among the largest 50 tested primes: "
            f"{'none' if not large_range_danger_primes else ', '.join(map(str, large_range_danger_primes))}.",
            "",
            "## Spike autopsies",
            "",
            "A spike is labeled *concentrated* when one non-reflection gap "
            "accounts for at least 25% of all nonforced unordered pairs, and "
            "*diffuse* otherwise. An AP is reported as structural only at length "
            "at least 3; length 2 is tautological.",
            "",
            "| rank | p | E/p | Z_p | max N_p(a), a | gap diagnosis | "
            "popular-fiber structure |",
            "|---:|---:|---:|---:|---|---|---|",
        ]
    )
    for rank, record in enumerate(worst, 1):
        inspection = inspections[record.p]
        if inspection.nonforced_unordered_pairs:
            share = (
                inspection.max_nonreflection_gap_count
                / inspection.nonforced_unordered_pairs
            )
            diagnosis = "concentrated" if share >= 0.25 else "diffuse"
            gap_text = (
                f"{diagnosis}; h={inspection.max_nonreflection_gap}, "
                f"count={inspection.max_nonreflection_gap_count}, "
                f"share={share:.3%}"
            )
        else:
            gap_text = "forced-only; no non-reflection gap"
        ap_text = (
            f"positions [{fmt_positions(inspection.popular_positions)}]; "
            f"reflection pairs={inspection.popular_reflection_pairs}, "
            f"center={'yes' if inspection.popular_central_fixed_point else 'no'}; "
            f"longest AP={inspection.longest_ap_length} "
            f"(step {inspection.longest_ap_step}); adjacent-gap mode="
            f"{inspection.adjacent_gap_mode} x{inspection.adjacent_gap_mode_count}"
        )
        lines.append(
            f"| {rank} | {record.p} | {fmt_ratio(record.energy_over_p)} | "
            f"{record.zero_count} | {record.max_multiplicity}, "
            f"a={record.popular_value} | {gap_text} | {ap_text} |"
        )

    lines.extend(
        [
            "",
            "## Envelope test and verdict",
            "",
            f"Violations of `E_off <= 100 p log^2 p`: {len(envelope_violations)}.",
        ]
    )
    if envelope_violations:
        lines.append(
            "Violating primes: "
            + ", ".join(str(record.p) for record in envelope_violations)
            + "."
        )
    else:
        lines.append("Violating primes: none.")

    max_125 = max(records, key=lambda record: record.energy_over_p_125)
    max_p = max(records, key=cmp_to_key(compare_energy_over_p))
    max_nonforced = max(records, key=lambda record: record.nonforced_over_p_125)
    genuine_resonance = (
        slope > 1.25
        or bool(envelope_violations)
        or bool(large_range_danger_primes)
        or bool(concentrated_danger_primes)
    )
    lines.extend(
        [
            "",
            f"Globally, `max E_off/p^(5/4) = "
            f"{fmt_ratio(max_125.energy_over_p_125)}` at `p={max_125.p}`, "
            f"while `max E_off/p = {fmt_ratio(max_p.energy_over_p)}` at "
            f"`p={max_p.p}`. After subtracting reflection, "
            f"`max (E_off-E_refl)/p^(5/4) = "
            f"{fmt_ratio(max_nonforced.nonforced_over_p_125)}` at "
            f"`p={max_nonforced.p}`.",
            "",
        ]
    )
    if genuine_resonance:
        verdict = (
            "**Verdict.** At least one prespecified resonance alarm fired, so "
            "this range does not justify prioritizing a p^(1+eps) theorem; the "
            "5/4 target should be treated cautiously pending analysis of the "
            "listed signal."
        )
    else:
        if danger_primes:
            if concentrated_danger_primes:
                concentration_text = (
                    "include concentrated cases at p="
                    + ", ".join(map(str, concentrated_danger_primes))
                )
            else:
                concentration_text = "are all diffuse under the stated 25% test"
            danger_qualification = (
                f" The {len(danger_primes)} gap-threshold crossings are confined "
                f"to low-p spike primes (largest p={max(danger_primes)}), "
                f"{concentration_text}, and do not persist among the largest 50 "
                "primes; they are not a growing resonance."
            )
        else:
            danger_qualification = " No inspected gap crosses p^(1/4)."
        verdict = (
            "**Verdict.** No genuine resonance was found: the fitted block-max "
            "slope stays below 1.25, the 100 p log^2 p envelope has no violation, "
            "and no p^(1/4) gap signal survives at the top of the tested range."
            + danger_qualification
            + " The data therefore support attacking the stronger p^(1+eps) "
            "collision-energy theorem; the 5/4 target is comfortably supported "
            "as the fallback."
        )
    lines.extend([verdict, ""])
    return "\n".join(lines)


def progress_header(minimum: int, maximum: int, started_at: str) -> str:
    columns = (
        "p\tE_off\tE_refl\tE_nonforced\tZ_p\tmax_N\tpopular_a\t"
        "E_over_p\tE_over_p_log2\tE_over_p_1.25\t"
        "Enf_over_p\tEnf_over_p_log2\tEnf_over_p_1.25"
    )
    return (
        "# q32 transfer-energy scan; integer-exact recurrence and energies\n"
        f"# requested_range=[{minimum},{maximum}] started={started_at}\n"
        f"{columns}\n"
    )


def progress_row(record: PrimeRecord) -> str:
    return (
        f"{record.p}\t{record.energy}\t{record.reflection_energy}\t"
        f"{record.nonforced_energy}\t{record.zero_count}\t"
        f"{record.max_multiplicity}\t{record.popular_value}\t"
        f"{record.energy_over_p:.17g}\t"
        f"{record.energy_over_p_log2:.17g}\t"
        f"{record.energy_over_p_125:.17g}\t"
        f"{record.nonforced_over_p:.17g}\t"
        f"{record.nonforced_over_p_log2:.17g}\t"
        f"{record.nonforced_over_p_125:.17g}\n"
    )


def parse_args(argv: Sequence[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--min-prime", type=int, default=1_000)
    parser.add_argument("--max-prime", type=int, default=100_000)
    parser.add_argument("--progress", type=Path, default=DEFAULT_PROGRESS)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--progress-interval", type=float, default=30.0)
    parser.add_argument("--sanity-every", type=int, default=100)
    parser.add_argument("--largest-gap-primes", type=int, default=50)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(sys.argv[1:] if argv is None else argv)
    if args.min_prime < 7:
        raise SystemExit("--min-prime must be at least 7")
    if args.max_prime < args.min_prime:
        raise SystemExit("--max-prime must be at least --min-prime")
    if args.sanity_every <= 0 or args.largest_gap_primes < 0:
        raise SystemExit("sanity/gap-prime counts must be nonnegative (sanity > 0)")

    primes = primes_in_interval(args.min_prime, args.max_prime)
    if not primes:
        raise SystemExit("no primes in the requested interval")
    largest_gap_set = set(primes[-args.largest_gap_primes :])
    started_wall = time.strftime("%Y-%m-%d %H:%M:%S %Z")
    started = time.monotonic()
    next_progress = started + args.progress_interval
    records: list[PrimeRecord] = []
    inspections: dict[int, Inspection] = {}
    running_spikes: list[PrimeRecord] = []
    spike_snapshots: dict[int, tuple[array, array]] = {}

    args.progress.parent.mkdir(parents=True, exist_ok=True)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    print(
        f"START range=[{args.min_prime},{args.max_prime}] primes={len(primes)} "
        f"first={primes[0]} last={primes[-1]}",
        flush=True,
    )
    with args.progress.open("w", encoding="utf-8", buffering=1) as progress:
        progress.write(
            progress_header(args.min_prime, args.max_prime, started_wall)
        )
        for index, p in enumerate(primes, 1):
            values, counts, energy, max_multiplicity, popular_value = apery_energy(p)
            reflection_energy = p - 3
            if energy < reflection_energy:
                raise AssertionError(
                    f"p={p}: E_off={energy} smaller than reflection {reflection_energy}"
                )
            nonforced_energy = energy - reflection_energy
            e_norm = normalized(energy, p)
            nf_norm = normalized(nonforced_energy, p)
            record = PrimeRecord(
                p=p,
                energy=energy,
                reflection_energy=reflection_energy,
                nonforced_energy=nonforced_energy,
                zero_count=counts[0],
                max_multiplicity=max_multiplicity,
                popular_value=popular_value,
                energy_over_p=e_norm[0],
                energy_over_p_log2=e_norm[1],
                energy_over_p_125=e_norm[2],
                nonforced_over_p=nf_norm[0],
                nonforced_over_p_log2=nf_norm[1],
                nonforced_over_p_125=nf_norm[2],
            )
            records.append(record)
            progress.write(progress_row(record))

            if index % args.sanity_every == 0:
                message = sanity_gate(p, values, counts, index)
                print(message, flush=True)
                progress.write(f"# {message}\n")
                progress.flush()
                os.fsync(progress.fileno())

            enters_spike_list = len(running_spikes) < 5
            if not enters_spike_list:
                worst_running = min(
                    running_spikes, key=cmp_to_key(compare_energy_over_p)
                )
                enters_spike_list = energy * worst_running.p > worst_running.energy * p
            if enters_spike_list:
                running_spikes.append(record)
                running_spikes.sort(
                    key=cmp_to_key(compare_energy_over_p), reverse=True
                )
                del running_spikes[5:]

                surviving_spikes = {item.p for item in running_spikes}
                if p in surviving_spikes:
                    spike_snapshots[p] = (array("I", values), array("I", counts))
                for old_p in list(spike_snapshots):
                    if old_p not in surviving_spikes:
                        del spike_snapshots[old_p]

            if p in largest_gap_set:
                inspections[p] = inspect_prime(record, values, counts)

            now = time.monotonic()
            if now >= next_progress:
                elapsed = now - started
                rate = index / elapsed
                remaining = (len(primes) - index) / rate if rate else float("inf")
                line = (
                    f"PROGRESS {index}/{len(primes)} p={p} "
                    f"elapsed={elapsed:.1f}s eta={remaining:.1f}s "
                    f"E/p={record.energy_over_p:.6g}"
                )
                print(line, flush=True)
                progress.write(f"# {line}\n")
                next_progress = now + args.progress_interval

    elapsed = time.monotonic() - started
    final_worst = sorted(
        records, key=cmp_to_key(compare_energy_over_p), reverse=True
    )[:5]
    for record in final_worst:
        if record.p not in inspections:
            values, counts = spike_snapshots[record.p]
            inspections[record.p] = inspect_prime(record, values, counts)
    missing = [record.p for record in final_worst if record.p not in inspections]
    if missing:
        raise AssertionError(f"missing spike inspections: {missing}")

    report = build_report(
        records,
        inspections,
        args.min_prime,
        args.max_prime,
        started_wall,
        elapsed,
    )
    args.report.write_text(report, encoding="utf-8")
    summary_blocks = [
            [record for record in records if dyadic_lower(record.p) == lower]
            for lower in sorted({dyadic_lower(record.p) for record in records})
        ]
    slope, points = regression_slope(
        complete_blocks(summary_blocks, args.min_prime, args.max_prime)
    )
    max_125 = max(records, key=lambda record: record.energy_over_p_125)
    max_p = max(records, key=cmp_to_key(compare_energy_over_p))
    violations = sum(
        record.energy > 100 * record.p * math.log(record.p) ** 2
        for record in records
    )
    danger_count = sum(inspection.danger for inspection in inspections.values())
    print(
        f"DONE primes={len(records)} elapsed={elapsed:.1f}s "
        f"slope={slope:.6f} blocks={points} "
        f"max_E/p^1.25={max_125.energy_over_p_125:.8g}@{max_125.p} "
        f"max_E/p={max_p.energy_over_p:.8g}@{max_p.p} "
        f"envelope_violations={violations} gap_dangers={danger_count}",
        flush=True,
    )
    print(f"REPORT {args.report}", flush=True)
    print(f"PROGRESS_FILE {args.progress}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
