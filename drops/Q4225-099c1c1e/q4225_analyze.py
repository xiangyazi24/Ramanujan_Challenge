#!/usr/bin/env python3
"""Analyze Q4225 canonical popular-edge ledger outputs.

The script never substitutes an unselected mass for M_sigma.  Fits are emitted
only when at least three rows have both M_sigma>0 and E_pop>0.  An event-free
or zero-denominator range is reported as censored/undefined, not as exponent
infinity and not as a theorem.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
from statistics import mean
from typing import Iterable


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def linear_fit(points: list[tuple[float, float]]) -> tuple[float, float, float]:
    """Return intercept, slope, R^2 for y=intercept+slope*x."""
    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    xb, yb = mean(xs), mean(ys)
    sxx = sum((x - xb) ** 2 for x in xs)
    if sxx == 0:
        raise ValueError("degenerate x values")
    slope = sum((x - xb) * (y - yb) for x, y in points) / sxx
    intercept = yb - slope * xb
    sst = sum((y - yb) ** 2 for y in ys)
    ssr = sum((y - (intercept + slope * x)) ** 2 for x, y in points)
    r2 = 1.0 - ssr / sst if sst else 1.0
    return intercept, slope, r2


def fmt(value: float | None, digits: int = 6) -> str:
    return "NA" if value is None or not math.isfinite(value) else f"{value:.{digits}g}"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    src = Path(args.input)
    dst = Path(args.output)
    dst.mkdir(parents=True, exist_ok=True)

    manifest = json.loads((src / "manifest.json").read_text())
    ledger = read_csv(src / "ledger_summary.csv")
    dprm = read_csv(src / "dprm_cells.csv")
    leaves = read_csv(src / "raw_leaves.csv")

    fit_rows: list[dict[str, object]] = []
    groups: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in ledger:
        groups[(row["sign"], f"{float(row['beta']):.3f}")].append(row)

    for (sign, beta_text), rows in sorted(groups.items()):
        beta = float(beta_text)
        positive: list[tuple[float, float]] = []
        for row in rows:
            M = int(row["M"])
            e = float(row["Epop"])
            if M > 0 and e > 0:
                positive.append((math.log(int(row["T"])), math.log(e / M)))
        record: dict[str, object] = {
            "sign": sign,
            "beta": beta,
            "positive_points": len(positive),
            "fit_status": "undefined: fewer than three positive actual-mass ledger points",
        }
        if len(positive) >= 3:
            intercept, slope, r2 = linear_fit(positive)
            total_exponent = -slope
            record.update(
                fit_status="fit",
                intercept=intercept,
                slope=slope,
                r2=r2,
                total_exponent=total_exponent,
                extra_kappa_over_beta3=total_exponent - beta / 3.0,
            )
        fit_rows.append(record)

    with (dst / "fits.csv").open("w", newline="") as handle:
        fields = [
            "sign", "beta", "positive_points", "fit_status", "intercept",
            "slope", "r2", "total_exponent", "extra_kappa_over_beta3",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in fit_rows:
            writer.writerow(row)

    gap_counter = Counter(row["delta"] for row in dprm)
    type_counter = Counter(row["type"] for row in dprm)
    R_counter = Counter(row["R"] for row in dprm)
    content_counter = Counter(row["primitive_content"] for row in dprm)

    standard_rows = [row for row in ledger if float(row["beta"]) < 0.1]
    max_scales: dict[str, int] = {}
    for field in ("Y", "Lambda", "H", "Jstar"):
        max_scales[field] = max((int(row[field]) for row in standard_rows), default=0)

    actual_M_rows = [row for row in ledger if int(row["M"]) > 0]
    positive_E_rows = [row for row in ledger if int(row["M"]) > 0 and float(row["Epop"]) > 0]

    lines: list[str] = []
    lines.append("# Q4225 computed-law analysis")
    lines.append("")
    lines.append("## Execution census")
    lines.append("")
    lines.append(f"- Prime cutoff: `{manifest['pmax']}`")
    lines.append(f"- Prime count: `{manifest['prime_count']}`")
    lines.append(f"- Apéry zero records: `{manifest['total_zero_records']}`")
    lines.append(f"- Selected upper states `(q,t)`: `{manifest['selected_states']}`")
    lines.append(f"- Raw plus leaves before minus-first: `{manifest['raw_plus_before_minus_first']}`")
    lines.append(f"- Raw minus leaves: `{manifest['raw_minus']}`")
    lines.append(f"- Actual plus leaves after minus-first: `{manifest['actual_plus_after_minus_first']}`")
    lines.append(f"- Actual minus leaves: `{manifest['actual_minus']}`")
    lines.append(f"- Runtime: `{manifest['elapsed_seconds']:.3f}` seconds on `{manifest['threads']}` threads")
    lines.append("")

    checks = manifest["checks"]
    lines.append("## Exact correctness gates")
    lines.append("")
    for name in (
        "barrett_failures", "recurrence_failures", "reflection_failures",
        "consecutive_failures",
    ):
        lines.append(f"- `{name}`: `{checks[name]}`")
    lines.append(f"- q<=5000 selected-state checkpoint: `{checks['q5000_checkpoint_ok']}`")
    lines.append("")

    lines.append("## Measured popular-edge law")
    lines.append("")
    if not actual_M_rows:
        lines.append(
            "The denominator `M_sigma(T)` is zero in every computed shell and both signs. "
            "Consequently `E_pop/M_sigma` is **undefined**, not zero, and no decay exponent "
            "can be fitted honestly.  The data are pre-first-event evidence only."
        )
    elif not positive_E_rows:
        lines.append(
            "Some shells have positive actual marked mass, but no canonical DPRM edge was "
            "found.  The ratios are right-censored at zero; a power-law fit is not reported."
        )
    else:
        lines.append("Positive ledger points exist; fitted exponents are in `fits.csv`.")
    lines.append("")
    lines.append(
        f"Across the standard beta<0.1 rows, the largest finite scales reached were "
        f"`Y={max_scales['Y']}`, `Lambda={max_scales['Lambda']}`, "
        f"`H={max_scales['H']}`, `J*={max_scales['Jstar']}`."
    )
    lines.append("")

    lines.append("## Canonical DPRM distributions")
    lines.append("")
    lines.append(f"- DPRM cell records: `{len(dprm)}`")
    lines.append(f"- Canonical gap histogram: `{dict(gap_counter)}`")
    lines.append(f"- Matching-size histogram: `{dict(R_counter)}`")
    lines.append(f"- Fold/order type histogram: `{dict(type_counter)}`")
    lines.append(f"- Common primitive-content statuses: `{dict(content_counter)}`")
    lines.append("")

    lines.append("## Interpretation guard")
    lines.append("")
    lines.append(
        "An empty actual ledger is consistent with every proposed upper bound, including "
        "no extra saving, exactly `T^{-beta/3}`, and `T^{-beta/3-kappa}`.  It distinguishes "
        "none of them.  Likewise the tiny integer values of `Lambda` and `H` in the reachable "
        "standard-beta range are a severe finite-size effect.  Neither emptiness nor a fitted "
        "finite exponent is promoted to a theorem."
    )
    lines.append("")

    lines.append("## Files")
    lines.append("")
    lines.append("- `zero_sets.csv`: every computed prime and its full Apéry zero set.")
    lines.append("- `selected_states.csv`: every actual upper mark `(q,t)`.")
    lines.append("- `raw_leaves.csv`: exact raw two-row leaves after minus-first.")
    lines.append("- `ledger_summary.csv`: all `(T,beta,sign)` normalizations and ledger masses.")
    lines.append("- `dprm_cells.csv`: every canonical DPRM cell and matching.")
    lines.append("- `gap_hist.csv`, `R_hist.csv`, `type_hist.csv`, `residue_hist.csv`: raw structural histograms.")
    lines.append("- `manifest.json`: cutoffs, conventions, runtime, and correctness gates.")
    lines.append("")

    (dst / "analysis.md").write_text("\n".join(lines) + "\n")
    (dst / "analysis.json").write_text(
        json.dumps(
            {
                "manifest": manifest,
                "actual_mass_rows": len(actual_M_rows),
                "positive_popular_edge_rows": len(positive_E_rows),
                "dprm_cells": len(dprm),
                "raw_leaves": len(leaves),
                "max_standard_scales": max_scales,
                "fits": fit_rows,
                "gap_hist": dict(gap_counter),
                "R_hist": dict(R_counter),
                "type_hist": dict(type_counter),
                "primitive_content_hist": dict(content_counter),
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
