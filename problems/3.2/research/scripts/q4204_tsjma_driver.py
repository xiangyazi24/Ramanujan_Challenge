#!/usr/bin/env python3
"""Build and run the exact Q4204 TSJMA/SRT computation.

This is the reproducible Python entry point for the campaign.  The numerical
kernel is C++20 because the large-range computation is dominated by billions
of exact modular-recurrence steps.  No floating-point value is used in Apéry
zero generation, row construction, fibre assignment, cell occupancy, prefix
support, or valuation tests; floating point is used only for displayed logs,
correlations, and regression slopes.

The launcher performs a checkpoint run before the requested full run.  The
checkpoint must reproduce the previously banked exact counts:

* 163 selected states through q=1000;
* 68 lower columns with at least two zeros through p=1000;
* 605 selected states for 17 <= q <= 5000;
* zero clean actual q6 rows through q=5000.

Example
-------
python3 problems/3.2/research/scripts/q4204_tsjma_driver.py \
  --qmax 300000 \
  --betas 0.10,0.15,0.30,0.50 \
  --out-dir q4204_results
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
from typing import Sequence


SCRIPT = Path(__file__).resolve()
SOURCE = SCRIPT.with_name("q4204_tsjma_scan.cpp")
PROJECT_ROOT = SCRIPT.parents[4]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--qmax", type=int, default=300_000)
    parser.add_argument("--sanity-qmax", type=int, default=5_000)
    parser.add_argument("--betas", default="0.10,0.15,0.30,0.50")
    parser.add_argument("--threads", type=int, default=0)
    parser.add_argument("--top-cells", type=int, default=24)
    parser.add_argument("--out-dir", type=Path, default=Path("q4204_results"))
    parser.add_argument("--build-dir", type=Path, default=Path(".q4204_build"))
    parser.add_argument("--compiler", default=os.environ.get("CXX", "g++"))
    parser.add_argument("--skip-sanity", action="store_true")
    return parser.parse_args()


def run_checked(command: Sequence[str], *, cwd: Path) -> None:
    print("+", " ".join(command), flush=True)
    subprocess.run(list(command), cwd=cwd, check=True)


def compile_backend(compiler: str, build_dir: Path) -> Path:
    if not SOURCE.is_file():
        raise FileNotFoundError(SOURCE)
    resolved = shutil.which(compiler)
    if resolved is None:
        raise RuntimeError(f"C++ compiler not found: {compiler}")
    build_dir.mkdir(parents=True, exist_ok=True)
    binary = build_dir / "q4204_tsjma_scan"
    # -include array repairs the intentionally minimal backend include list and
    # makes the build independent of transitive standard-library includes.
    command = [
        resolved,
        "-O3",
        "-DNDEBUG",
        "-std=c++20",
        "-pthread",
        "-march=native",
        "-include",
        "array",
        str(SOURCE),
        "-o",
        str(binary),
    ]
    run_checked(command, cwd=PROJECT_ROOT)
    return binary.resolve()


def invoke_backend(
    binary: Path,
    *,
    qmax: int,
    betas: str,
    threads: int,
    top_cells: int,
    out_dir: Path,
) -> None:
    if out_dir.exists():
        shutil.rmtree(out_dir)
    command = [
        str(binary),
        "--qmax",
        str(qmax),
        "--betas",
        betas,
        "--threads",
        str(threads),
        "--top-cells",
        str(top_cells),
        "--out-dir",
        str(out_dir),
    ]
    run_checked(command, cwd=PROJECT_ROOT)


def read_summary(out_dir: Path) -> dict[str, object]:
    summary_path = out_dir / "summary.json"
    report_path = out_dir / "report.md"
    if not summary_path.is_file() or not report_path.is_file():
        raise RuntimeError(f"backend did not create complete output in {out_dir}")
    with summary_path.open("r", encoding="utf-8") as source:
        return json.load(source)


def verify_sanity(out_dir: Path) -> None:
    summary = read_summary(out_dir)
    report = (out_dir / "report.md").read_text(encoding="utf-8")
    required_text = {
        "states q<=1000             = 163": "selected-state q<=1000 checkpoint",
        "lower z_p>=2, p<=1000      = 68": "lower zero-column checkpoint",
        "states 17<=q<=5000         = 605": "selected-state q<=5000 checkpoint",
        "clean actual rows q<=5000      = 0": "actual q6 emptiness checkpoint",
    }
    missing = [label for text, label in required_text.items() if text not in report]
    if missing:
        raise AssertionError(f"sanity checkpoint failure: {missing}")
    if int(summary["qmax"]) != 5_000:
        raise AssertionError(summary)
    print("Q4204_SANITY_OK", json.dumps(summary, sort_keys=True), flush=True)


def verify_full(out_dir: Path, qmax: int) -> dict[str, object]:
    summary = read_summary(out_dir)
    if int(summary["qmax"]) != qmax:
        raise AssertionError((summary.get("qmax"), qmax))
    for name in ("actual_rows.csv", "unselected_overcarrier_rows.csv", "translated_metrics.csv"):
        if not (out_dir / name).is_file():
            raise RuntimeError(f"missing output: {name}")
    return summary


def main() -> None:
    args = parse_args()
    if args.qmax < 1_000:
        raise ValueError("qmax must be at least 1000")
    if args.sanity_qmax != 5_000 and not args.skip_sanity:
        raise ValueError("the exact checkpoint run is fixed at qmax=5000")

    started = time.perf_counter()
    build_dir = args.build_dir.resolve()
    full_dir = args.out_dir.resolve()
    binary = compile_backend(args.compiler, build_dir)

    if not args.skip_sanity:
        sanity_dir = full_dir.parent / f"{full_dir.name}_sanity"
        invoke_backend(
            binary,
            qmax=args.sanity_qmax,
            betas="0.10,0.30,0.50",
            threads=args.threads,
            top_cells=min(args.top_cells, 8),
            out_dir=sanity_dir,
        )
        verify_sanity(sanity_dir)

    invoke_backend(
        binary,
        qmax=args.qmax,
        betas=args.betas,
        threads=args.threads,
        top_cells=args.top_cells,
        out_dir=full_dir,
    )
    summary = verify_full(full_dir, args.qmax)
    elapsed = time.perf_counter() - started
    print("Q4204_FULL_OK", json.dumps(summary, sort_keys=True), flush=True)
    print(f"Q4204_DRIVER_SECONDS={elapsed:.6f}", flush=True)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:  # fail loudly in Actions and local reproductions
        print(f"Q4204_DRIVER_FATAL: {exc}", file=sys.stderr, flush=True)
        raise
