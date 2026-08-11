#!/usr/bin/env python3
"""Q7310: exact fixed-q cutoff-intersection scan.

This is a mechanical audit built on the existing modules

    q32_fixed_q_content.py
    q32_multicutoff_content_intersection.py
    q32_newton.py
    q32_strehl_gcd.py

For n=m and q in a requested set, put

    H = floor((m-q)/(2q+1)),
    I_w = gcd(Gamma_(m,H), ..., Gamma_(m,H+w-1)),
    E   = binom(m,H+1) binom(m+H+1,H+1),
    C_w = gcd(I_w,E),                  w in {1,2,4,8}.

The identity used by q32_multicutoff_content_intersection.py gives the exact
fast update

    I_w = gcd(Gamma_(m,H), kappa_(m,H+1), ..., kappa_(m,H+w-1)).

The default scan writes exact decimal I_w and C_w for every row carrying a
lower-digit hit or a high-digit singular prime in the requested dyadic
shells.  Use --all-candidate-rows to scan every m having at least one prime
p in (X,2X] with floor(m/p)=q.

No asymptotic or proof claim is made by this program.  It asserts the local
support identities at every scanned candidate prime and reports raw data.
"""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from math import comb, gcd, prod
from pathlib import Path
from typing import Iterable

from q32_fixed_q_content import truncation_content
from q32_multicutoff_content_intersection import kappa_against
from q32_newton import apery_numbers
from q32_strehl_gcd import franel_numbers, primes_up_to


DEFAULT_XS = (32, 64, 128)
DEFAULT_QS = (1, 2, 3)
DEFAULT_WIDTHS = (1, 2, 4, 8)


def parse_int_tuple(text: str) -> tuple[int, ...]:
    values = tuple(int(piece) for piece in text.split(",") if piece.strip())
    if not values:
        raise argparse.ArgumentTypeError("expected a nonempty comma-separated list")
    if any(value <= 0 for value in values):
        raise argparse.ArgumentTypeError("all values must be positive")
    return values


def log_int(value: int) -> float:
    """Natural logarithm without converting a large integer to float."""
    if value <= 1:
        return 0.0
    shift = max(0, value.bit_length() - 53)
    return math.log(value >> shift) + shift * math.log(2.0)


def factor_bounded(value: int, primes: list[int], bound: int) -> dict[int, int]:
    """Factor value, knowing that every prime factor is at most bound."""
    if value < 1:
        raise ValueError("factor_bounded expects a positive integer")
    remaining = value
    factors: dict[int, int] = {}
    for prime in primes:
        if prime > bound or prime * prime > remaining:
            break
        exponent = 0
        while remaining % prime == 0:
            remaining //= prime
            exponent += 1
        if exponent:
            factors[prime] = exponent
    if remaining > 1:
        if remaining > bound:
            raise AssertionError(("unbounded residual factor", remaining, bound))
        factors[remaining] = factors.get(remaining, 0) + 1
    check = prod(prime**exponent for prime, exponent in factors.items())
    assert check == value, (value, factors, check)
    return factors


def factor_string(factors: dict[int, int]) -> str:
    if not factors:
        return "1"
    return "*".join(
        str(prime) if exponent == 1 else f"{prime}^{exponent}"
        for prime, exponent in sorted(factors.items())
    )


def support_string(values: Iterable[int]) -> str:
    items = sorted(set(values))
    return ";".join(str(value) for value in items) if items else "-"


def product_from_factors(
    factors: dict[int, int], selected_primes: set[int]
) -> int:
    return prod(
        prime**exponent
        for prime, exponent in factors.items()
        if prime in selected_primes
    )


def exact_intersections(
    n: int,
    base_cutoff: int,
    gamma: int,
    widths: tuple[int, ...],
) -> dict[int, int]:
    requested = set(widths)
    result: dict[int, int] = {}
    current = gamma
    if 1 in requested:
        result[1] = current
    for offset in range(1, max(widths)):
        current = kappa_against(n, base_cutoff + offset, current)
        width = offset + 1
        if width in requested:
            result[width] = current
    assert set(result) == requested
    return result


def write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def scan(
    xs: tuple[int, ...],
    qs: tuple[int, ...],
    widths: tuple[int, ...],
    output_dir: Path,
    all_candidate_rows: bool,
) -> None:
    if 1 not in widths:
        raise ValueError("width list must contain 1")
    if tuple(sorted(set(widths))) != widths:
        raise ValueError("widths must be strictly increasing")
    if max(widths) > 8:
        raise ValueError("this audit is intentionally capped at eight cutoffs")
    if min(xs) < 32:
        raise ValueError("use X>=32 so H+7<p is uniform in the requested slices")

    output_dir.mkdir(parents=True, exist_ok=True)
    max_x = max(xs)
    max_q = max(qs)
    max_n = 2 * (max_q + 1) * max_x - 1
    max_endpoint_index = max_n + (max_n // 3) + 2

    print(f"precompute max_n={max_n} max_endpoint_index={max_endpoint_index}")
    franel = franel_numbers(max_n)
    apery = apery_numbers(2 * max_x)
    all_primes = primes_up_to(max_endpoint_index)

    gamma_cache: dict[tuple[int, int], int] = {}
    raw_rows: list[dict[str, object]] = []
    internal_rows: list[dict[str, object]] = []
    meta: dict[tuple[int, int], dict[str, object]] = {}

    for x in xs:
        shell_primes = [prime for prime in all_primes if x < prime <= 2 * x]
        zero_residues = {
            prime: [residue for residue in range(prime) if apery[residue] % prime == 0]
            for prime in shell_primes
        }

        for quotient in qs:
            low_targets: dict[int, list[int]] = defaultdict(list)
            singular_targets: dict[int, list[int]] = defaultdict(list)
            candidate_rows: set[int] = set()

            for prime in shell_primes:
                for residue in range(prime):
                    n = quotient * prime + residue
                    assert n < x * x
                    candidate_rows.add(n)
                    if apery[quotient] % prime == 0:
                        singular_targets[n].append(prime)
                for residue in zero_residues[prime]:
                    n = quotient * prime + residue
                    low_targets[n].append(prime)

            event_rows = set(low_targets) | set(singular_targets)
            scanned_rows = candidate_rows if all_candidate_rows else event_rows
            scan_mode = "all-candidate" if all_candidate_rows else "event-supported"
            meta[(x, quotient)] = {
                "scan_mode": scan_mode,
                "shell_prime_count": len(shell_primes),
                "candidate_row_count": len(candidate_rows),
                "event_row_count": len(event_rows),
                "low_active_row_count": len(low_targets),
                "singular_active_row_count": len(singular_targets),
                "low_target_event_count": sum(len(values) for values in low_targets.values()),
                "singular_event_count": sum(len(values) for values in singular_targets.values()),
                "distinct_low_target_primes": sorted(
                    {prime for values in low_targets.values() for prime in values}
                ),
                "distinct_singular_primes": sorted(
                    {prime for values in singular_targets.values() for prime in values}
                ),
            }

            print(
                f"X={x} q={quotient} mode={scan_mode} "
                f"candidate_rows={len(candidate_rows)} event_rows={len(event_rows)} "
                f"scanned_rows={len(scanned_rows)}"
            )

            for row_number, n in enumerate(sorted(scanned_rows), start=1):
                candidates = [
                    prime for prime in shell_primes if n // prime == quotient
                ]
                assert candidates
                low = sorted(low_targets.get(n, []))
                singular = sorted(singular_targets.get(n, []))
                full_hits = sorted(set(low) | set(singular))
                cutoff = (n - quotient) // (2 * quotient + 1)
                last_cutoff = cutoff + max(widths) - 1
                assert last_cutoff < min(candidates), (
                    x,
                    quotient,
                    n,
                    cutoff,
                    min(candidates),
                )

                cache_key = (n, cutoff)
                gamma = gamma_cache.get(cache_key)
                if gamma is None:
                    gamma = truncation_content(n, quotient, franel)
                    gamma_cache[cache_key] = gamma
                intersections = exact_intersections(
                    n, cutoff, gamma, widths
                )

                endpoint_index = cutoff + 1
                endpoint = comb(n, endpoint_index) * comb(
                    n + endpoint_index, endpoint_index
                )
                factor_bound = n + endpoint_index
                endpoint_factors = factor_bounded(
                    endpoint, all_primes, factor_bound
                )

                low_set = set(low)
                candidate_set = set(candidates)
                for prime in candidates:
                    residue = n % prime
                    folded = min(residue, prime - 1 - residue)
                    is_low_hit = apery[folded] % prime == 0
                    assert is_low_hit == (prime in low_set)
                    assert (gamma % prime == 0) == is_low_hit, (
                        "Gamma support mismatch",
                        x,
                        quotient,
                        n,
                        prime,
                        folded,
                    )
                    assert endpoint % prime == 0, (
                        "endpoint misses candidate",
                        x,
                        quotient,
                        n,
                        prime,
                    )
                    for width, intersection in intersections.items():
                        assert (intersection % prime == 0) == is_low_hit, (
                            "intersection support mismatch",
                            x,
                            quotient,
                            n,
                            prime,
                            width,
                        )

                width_one_carrier = gcd(intersections[1], endpoint)
                width_one_log = log_int(width_one_carrier)

                for width in widths:
                    intersection = intersections[width]
                    carrier = gcd(intersection, endpoint)
                    carrier_factors = factor_bounded(
                        carrier, all_primes, factor_bound
                    )
                    carrier_support = set(carrier_factors)
                    assert carrier_support & candidate_set == low_set, (
                        "capped candidate support mismatch",
                        x,
                        quotient,
                        n,
                        width,
                        carrier_factors,
                        candidates,
                        low,
                    )

                    target_part = product_from_factors(carrier_factors, low_set)
                    target_radical = prod(low) if low else 1
                    assert target_part % target_radical == 0
                    nuisance_part = carrier // target_part
                    nuisance_factors = {
                        prime: exponent
                        for prime, exponent in carrier_factors.items()
                        if prime not in low_set
                    }
                    nuisance_support = sorted(nuisance_factors)
                    candidate_phantoms = sorted(
                        (carrier_support & candidate_set) - low_set
                    )
                    assert not candidate_phantoms

                    internal = {
                        "X": x,
                        "q": quotient,
                        "m": n,
                        "H": cutoff,
                        "width": width,
                        "cutoff_last": cutoff + width - 1,
                        "candidate_primes_list": candidates,
                        "low_target_primes_list": low,
                        "singular_high_primes_list": singular,
                        "full_hit_primes_list": full_hits,
                        "carrier_factors_dict": carrier_factors,
                        "nuisance_factors_dict": nuisance_factors,
                        "intersection": intersection,
                        "endpoint": endpoint,
                        "carrier": carrier,
                        "target_part": target_part,
                        "target_radical": target_radical,
                        "nuisance_part": nuisance_part,
                        "log_intersection_over_X": log_int(intersection) / x,
                        "log_endpoint_over_X": log_int(endpoint) / x,
                        "log_carrier_over_X": log_int(carrier) / x,
                        "log_target_radical_over_X": log_int(target_radical) / x,
                        "log_target_part_over_X": log_int(target_part) / x,
                        "log_nuisance_part_over_X": log_int(nuisance_part) / x,
                        "log_removed_from_width1_over_X": (
                            width_one_log - log_int(carrier)
                        )
                        / x,
                        "target_support_preserved": carrier_support & candidate_set == low_set,
                    }
                    internal_rows.append(internal)

                    raw_rows.append(
                        {
                            "X": x,
                            "q": quotient,
                            "m": n,
                            "H": cutoff,
                            "width": width,
                            "cutoff_last": cutoff + width - 1,
                            "candidate_primes": support_string(candidates),
                            "low_target_primes": support_string(low),
                            "singular_high_primes": support_string(singular),
                            "full_hit_primes": support_string(full_hits),
                            "Gamma_intersection_exact": str(intersection),
                            "Gamma_intersection_digits": len(str(intersection)),
                            "endpoint_exact": str(endpoint),
                            "endpoint_factorization": factor_string(endpoint_factors),
                            "capped_carrier_exact": str(carrier),
                            "capped_carrier_factorization": factor_string(carrier_factors),
                            "target_part_exact": str(target_part),
                            "target_part_factorization": factor_string(
                                {
                                    prime: exponent
                                    for prime, exponent in carrier_factors.items()
                                    if prime in low_set
                                }
                            ),
                            "target_radical_exact": str(target_radical),
                            "nuisance_part_exact": str(nuisance_part),
                            "nuisance_factorization": factor_string(nuisance_factors),
                            "nuisance_support": support_string(nuisance_support),
                            "candidate_phantoms": support_string(candidate_phantoms),
                            "log_intersection_over_X": f"{internal['log_intersection_over_X']:.12f}",
                            "log_endpoint_over_X": f"{internal['log_endpoint_over_X']:.12f}",
                            "log_carrier_over_X": f"{internal['log_carrier_over_X']:.12f}",
                            "log_target_radical_over_X": f"{internal['log_target_radical_over_X']:.12f}",
                            "log_target_part_over_X": f"{internal['log_target_part_over_X']:.12f}",
                            "log_nuisance_part_over_X": f"{internal['log_nuisance_part_over_X']:.12f}",
                            "log_removed_from_width1_over_X": f"{internal['log_removed_from_width1_over_X']:.12f}",
                            "target_support_preserved": "YES",
                        }
                    )

                if row_number % 10 == 0 or row_number == len(scanned_rows):
                    print(
                        f"  X={x} q={quotient}: {row_number}/{len(scanned_rows)} rows"
                    )

    summary_rows: list[dict[str, object]] = []
    for x in xs:
        for quotient in qs:
            group_meta = meta[(x, quotient)]
            q_rows = [
                row
                for row in internal_rows
                if row["X"] == x and row["q"] == quotient
            ]
            by_width = {
                width: [row for row in q_rows if row["width"] == width]
                for width in widths
            }
            width_one_by_m = {
                int(row["m"]): row for row in by_width[1]
            }

            for width in widths:
                rows = by_width[width]
                nuisance_primes = sorted(
                    {
                        prime
                        for row in rows
                        for prime in row["nuisance_factors_dict"]
                    }
                )
                carrier_changed = sum(
                    int(row["carrier"] != width_one_by_m[int(row["m"])]["carrier"])
                    for row in rows
                )
                intersection_changed = sum(
                    int(
                        row["intersection"]
                        != width_one_by_m[int(row["m"])]["intersection"]
                    )
                    for row in rows
                )
                low_active_rows = [
                    row for row in rows if row["low_target_primes_list"]
                ]
                max_hits = max(
                    (len(row["low_target_primes_list"]) for row in rows),
                    default=0,
                )
                multi_hit_rows = sum(
                    len(row["low_target_primes_list"]) >= 2 for row in rows
                )
                target_preserved = sum(
                    bool(row["target_support_preserved"]) for row in rows
                )
                summary_rows.append(
                    {
                        "X": x,
                        "q": quotient,
                        "width": width,
                        "scan_mode": group_meta["scan_mode"],
                        "shell_prime_count": group_meta["shell_prime_count"],
                        "candidate_row_count": group_meta["candidate_row_count"],
                        "event_row_count": group_meta["event_row_count"],
                        "rows_scanned": len(rows),
                        "low_active_rows": group_meta["low_active_row_count"],
                        "singular_active_rows": group_meta["singular_active_row_count"],
                        "low_target_events": group_meta["low_target_event_count"],
                        "singular_events": group_meta["singular_event_count"],
                        "distinct_low_target_primes": support_string(
                            group_meta["distinct_low_target_primes"]
                        ),
                        "distinct_singular_primes": support_string(
                            group_meta["distinct_singular_primes"]
                        ),
                        "multi_low_hit_rows": multi_hit_rows,
                        "max_low_hits_in_row": max_hits,
                        "target_support_preserved_rows": target_preserved,
                        "intersection_changed_from_width1_rows": intersection_changed,
                        "carrier_changed_from_width1_rows": carrier_changed,
                        "nuisance_free_rows": sum(
                            int(row["nuisance_part"] == 1) for row in rows
                        ),
                        "distinct_nuisance_primes": support_string(nuisance_primes),
                        "max_log_Gamma_intersection_over_X": f"{max((row['log_intersection_over_X'] for row in rows), default=0.0):.12f}",
                        "mean_log_Gamma_intersection_over_X": f"{(sum(row['log_intersection_over_X'] for row in rows) / len(rows) if rows else 0.0):.12f}",
                        "max_log_capped_carrier_over_X": f"{max((row['log_carrier_over_X'] for row in rows), default=0.0):.12f}",
                        "mean_log_capped_carrier_over_X": f"{(sum(row['log_carrier_over_X'] for row in rows) / len(rows) if rows else 0.0):.12f}",
                        "max_log_target_radical_over_X": f"{max((row['log_target_radical_over_X'] for row in low_active_rows), default=0.0):.12f}",
                        "mean_log_target_radical_over_X_on_low_active_rows": f"{(sum(row['log_target_radical_over_X'] for row in low_active_rows) / len(low_active_rows) if low_active_rows else 0.0):.12f}",
                        "max_log_nuisance_part_over_X": f"{max((row['log_nuisance_part_over_X'] for row in rows), default=0.0):.12f}",
                        "mean_log_nuisance_part_over_X": f"{(sum(row['log_nuisance_part_over_X'] for row in rows) / len(rows) if rows else 0.0):.12f}",
                        "mean_log_removed_from_width1_over_X": f"{(sum(row['log_removed_from_width1_over_X'] for row in rows) / len(rows) if rows else 0.0):.12f}",
                        "max_log_removed_from_width1_over_X": f"{max((row['log_removed_from_width1_over_X'] for row in rows), default=0.0):.12f}",
                    }
                )

    summary_fields = list(summary_rows[0].keys()) if summary_rows else []
    raw_fields = list(raw_rows[0].keys()) if raw_rows else []
    write_csv(output_dir / "q7310_summary.csv", summary_fields, summary_rows)
    write_csv(output_dir / "q7310_raw_rows.csv", raw_fields, raw_rows)

    for x in xs:
        for quotient in qs:
            subset = [
                row
                for row in raw_rows
                if row["X"] == x and row["q"] == quotient
            ]
            write_csv(
                output_dir / f"q7310_raw_X{x}_q{quotient}.csv",
                raw_fields,
                subset,
            )

    # A concise Markdown table for direct inclusion in the report.
    markdown_columns = [
        "X",
        "q",
        "width",
        "rows_scanned",
        "low_target_events",
        "multi_low_hit_rows",
        "max_low_hits_in_row",
        "nuisance_free_rows",
        "carrier_changed_from_width1_rows",
        "max_log_capped_carrier_over_X",
        "max_log_target_radical_over_X",
        "max_log_nuisance_part_over_X",
    ]
    with (output_dir / "q7310_summary.md").open("w", encoding="utf-8") as handle:
        handle.write("|" + "|".join(markdown_columns) + "|\n")
        handle.write("|" + "|".join(["---"] * len(markdown_columns)) + "|\n")
        for row in summary_rows:
            handle.write(
                "|"
                + "|".join(str(row[column]) for column in markdown_columns)
                + "|\n"
            )

    with (output_dir / "q7310_manifest.txt").open("w", encoding="utf-8") as handle:
        handle.write(f"xs={','.join(map(str, xs))}\n")
        handle.write(f"qs={','.join(map(str, qs))}\n")
        handle.write(f"widths={','.join(map(str, widths))}\n")
        handle.write(f"all_candidate_rows={all_candidate_rows}\n")
        handle.write(f"raw_record_count={len(raw_rows)}\n")
        handle.write(f"summary_record_count={len(summary_rows)}\n")
        handle.write(f"gamma_cache_entries={len(gamma_cache)}\n")

    print(f"wrote {len(raw_rows)} raw records and {len(summary_rows)} summary records")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--xs",
        type=parse_int_tuple,
        default=DEFAULT_XS,
        help="comma-separated dyadic X values",
    )
    parser.add_argument(
        "--qs",
        type=parse_int_tuple,
        default=DEFAULT_QS,
        help="comma-separated quotient slices",
    )
    parser.add_argument(
        "--widths",
        type=parse_int_tuple,
        default=DEFAULT_WIDTHS,
        help="comma-separated cutoff-intersection widths",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
    )
    parser.add_argument(
        "--all-candidate-rows",
        action="store_true",
        help="scan all candidate rows rather than rows carrying an event",
    )
    args = parser.parse_args()
    scan(
        tuple(args.xs),
        tuple(args.qs),
        tuple(args.widths),
        args.output_dir,
        args.all_candidate_rows,
    )


if __name__ == "__main__":
    main()
