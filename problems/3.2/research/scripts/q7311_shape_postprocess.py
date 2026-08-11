#!/usr/bin/env python3
"""Mechanical classification table for Q7311 raw output.

No mathematical claim is attached to the labels.  They are deterministic
threshold flags recorded in the CSV itself:

* prime-triple diffuse: at least 20 active triples, top-10 absolute share
  <= 0.60, and participation ratio >= 20;
* frequency diffuse: for every examined leading actual triple, the absolute
  residual after its ten largest conjugate frequency pairs is >= 0.50 |C|;
* few coherent modes: for every examined leading actual triple, that residual
  is <= 0.25 |C|.
"""

from __future__ import annotations

import csv
import statistics
from pathlib import Path
from typing import Dict, List


OUT = Path("q7311-output")
REPORT = Path("Q7311_RESULTS.md")
MARKER = "\n## excess_shape.csv\n"


def read_csv(name: str) -> List[dict]:
    with (OUT / name).open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def write_csv(name: str, rows: List[dict], fields: List[str]) -> None:
    with (OUT / name).open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def csv_block(rows: List[dict], fields: List[str]) -> str:
    import io

    buf = io.StringIO()
    writer = csv.DictWriter(buf, fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return "```csv\n" + buf.getvalue() + "```\n"


def main() -> None:
    concentration = {int(r["X"]): r for r in read_csv("triple_concentration.csv")}
    modes: Dict[int, List[dict]] = {}
    for row in read_csv("mode_capture_summary.csv"):
        if row["set_label"] == "actual":
            modes.setdefault(int(row["X"]), []).append(row)

    random_rows = read_csv("random_comparison.csv")
    random_lookup = {
        (int(r["X"]), r["metric"]): r
        for r in random_rows
    }

    fields = [
        "X",
        "active_prime_triples",
        "prime_top10_abs_share",
        "prime_participation_ratio",
        "prime_diffuse_min_triples_threshold",
        "prime_diffuse_top10_share_threshold",
        "prime_diffuse_participation_threshold",
        "prime_triple_diffuse_flag",
        "examined_leading_triples_for_modes",
        "mode_abs_residual_after_top10_min",
        "mode_abs_residual_after_top10_median",
        "mode_abs_residual_after_top10_max",
        "frequency_diffuse_residual_threshold",
        "few_coherent_residual_threshold",
        "frequency_diffuse_flag",
        "few_coherent_modes_flag",
        "F3_actual",
        "F3_random_mean",
        "F3_random_zscore",
        "F3_empirical_ge_fraction",
        "ordered_primitive_abs_mass_actual",
        "ordered_primitive_abs_mass_random_mean",
        "ordered_primitive_abs_mass_random_zscore",
        "mechanical_shape_label",
    ]

    rows: List[dict] = []
    for X in sorted(concentration):
        c = concentration[X]
        active = int(c["active_triples"])
        top10 = float(c["top10_abs_share"])
        participation = float(c["participation_ratio"])
        prime_diffuse = active >= 20 and top10 <= 0.60 and participation >= 20.0

        residuals = [abs(float(r["residual_after_top10_over_C"])) for r in modes.get(X, [])]
        if residuals:
            residual_min = min(residuals)
            residual_med = statistics.median(residuals)
            residual_max = max(residuals)
            frequency_diffuse = all(v >= 0.50 for v in residuals)
            few_coherent = all(v <= 0.25 for v in residuals)
        else:
            residual_min = residual_med = residual_max = float("nan")
            frequency_diffuse = False
            few_coherent = False

        if active == 0:
            label = "no_primitive_prime_triples"
        elif prime_diffuse and frequency_diffuse:
            label = "diffuse_over_prime_triples_and_examined_modes"
        elif active < 20 and frequency_diffuse:
            label = "few_available_prime_triples_but_examined_modes_diffuse"
        elif few_coherent:
            label = "examined_excess_carried_by_few_coherent_modes"
        elif prime_diffuse:
            label = "diffuse_over_prime_triples_mode_result_mixed"
        else:
            label = "mixed_under_recorded_thresholds"

        f3 = random_lookup[(X, "F3")]
        mass = random_lookup[(X, "primitive_abs_mass")]
        rows.append(
            {
                "X": X,
                "active_prime_triples": active,
                "prime_top10_abs_share": c["top10_abs_share"],
                "prime_participation_ratio": c["participation_ratio"],
                "prime_diffuse_min_triples_threshold": 20,
                "prime_diffuse_top10_share_threshold": 0.60,
                "prime_diffuse_participation_threshold": 20.0,
                "prime_triple_diffuse_flag": int(prime_diffuse),
                "examined_leading_triples_for_modes": len(residuals),
                "mode_abs_residual_after_top10_min": "NA" if not residuals else f"{residual_min:.17g}",
                "mode_abs_residual_after_top10_median": "NA" if not residuals else f"{residual_med:.17g}",
                "mode_abs_residual_after_top10_max": "NA" if not residuals else f"{residual_max:.17g}",
                "frequency_diffuse_residual_threshold": 0.50,
                "few_coherent_residual_threshold": 0.25,
                "frequency_diffuse_flag": int(frequency_diffuse),
                "few_coherent_modes_flag": int(few_coherent),
                "F3_actual": f3["actual"],
                "F3_random_mean": f3["random_mean"],
                "F3_random_zscore": f3["zscore"],
                "F3_empirical_ge_fraction": f3["empirical_ge_fraction"],
                "ordered_primitive_abs_mass_actual": mass["actual"],
                "ordered_primitive_abs_mass_random_mean": mass["random_mean"],
                "ordered_primitive_abs_mass_random_zscore": mass["zscore"],
                "mechanical_shape_label": label,
            }
        )

    write_csv("excess_shape.csv", rows, fields)

    text = REPORT.read_text(encoding="utf-8")
    if MARKER in text:
        text = text.split(MARKER, 1)[0]
    addition = [
        "## excess_shape.csv",
        csv_block(rows, fields),
        "## shape_postprocessor.py",
        "```python\n" + Path(__file__).read_text(encoding="utf-8") + "\n```\n",
    ]
    REPORT.write_text(text.rstrip() + "\n\n" + "\n".join(addition), encoding="utf-8")


if __name__ == "__main__":
    main()
