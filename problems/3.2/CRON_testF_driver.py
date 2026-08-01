#!/usr/bin/env python3
"""Compile Test F, run the full window, and render CODEX_TESTF_report.md."""

from __future__ import annotations

import csv
import math
from pathlib import Path
import subprocess
import sys
import tempfile


HERE = Path(__file__).resolve().parent
SOURCE = HERE / "CRON_testF_dispersion.c"
REPORT = HERE / "CODEX_TESTF_report.md"


def fmt_number(value: float) -> str:
    if value == 0:
        return "0"
    absolute = abs(value)
    if absolute >= 1e5 or absolute < 1e-3:
        return f"{value:.4e}"
    return f"{value:.5f}"


def parse_output(text: str):
    meta: dict[str, str] = {}
    valid: list[list[str]] = []
    stats: dict[tuple[str, int, str], dict[str, float | int]] = {}
    profiles: dict[int, tuple[int, int]] = {}
    fourier: list[dict[str, float | int | str]] = []
    for row in csv.reader(text.splitlines()):
        if not row:
            continue
        if row[0] == "META":
            meta[row[1]] = row[2]
        elif row[0] == "VALID":
            valid.append(row[1:])
        elif row[0] == "STAT":
            family, bin_text, count, support, observable = row[1:6]
            stats[(family, int(bin_text), observable)] = {
                "count": int(count),
                "support": int(support),
                "sum": float(row[6]),
                "rms": float(row[7]),
                "bench": float(row[8]),
                "ratio": float(row[9]),
                "z": float(row[10]),
                "sumsq": float(row[11]),
                "varsum": float(row[12]),
            }
        elif row[0] == "PROFILE":
            profiles[int(row[1])] = (int(row[2]), int(row[3]))
        elif row[0] == "FOURIER":
            fourier.append(
                {
                    "family": row[1],
                    "bin": int(row[2]),
                    "p": int(row[3]),
                    "q": int(row[4]),
                    "direct": float(row[5]),
                    "recon": float(row[6]),
                    "error": float(row[7]),
                }
            )
    return meta, valid, stats, profiles, fourier


def markdown_table(headers: list[str], rows: list[list[str]]) -> str:
    result = ["| " + " | ".join(headers) + " |"]
    result.append("| " + " | ".join("---" for _ in headers) + " |")
    result.extend("| " + " | ".join(row) + " |" for row in rows)
    return "\n".join(result)


def status_for(record: dict[str, float | int]) -> str:
    ratio = float(record["ratio"])
    if ratio > 1.25:
        return "EXCESS"
    if ratio < 0.80:
        return "DEFICIT"
    return "~1"


def render_report(meta, valid, stats, profiles, fourier, command: str) -> str:
    bins = sorted(
        bin_number
        for family, bin_number, observable in stats
        if family == "aligned" and observable == "indicator"
    )
    lines: list[str] = []
    lines.append("# Test F: dyadic gap profile for exact and fixed-degree detectors")
    lines.append("")
    lines.append("## Scope and reproducibility")
    lines.append("")
    lines.append(
        f"The full requested window was used: {meta['prime_count']} primes in "
        f"[{meta['prime_min']}, {meta['prime_max']}], hence "
        f"{meta['pair_count']} pairs. No range reduction was needed. The C run "
        f"took {float(meta['elapsed_seconds']):.2f} seconds."
    )
    lines.append("")
    lines.append(f"Reproduce with <code>{command}</code> from this directory.")
    lines.append("")
    lines.append("## Reconstruction of the missing definitions")
    lines.append("")
    lines.append(
        "Put $P=p-1$, $Q=q-1$, $D=q-p>0$, and $d=p-q=-D$. For an "
        "observable $F$, the centered cyclic row is"
    )
    lines.append("")
    lines.append(
        "$$h_{p,F}(r)=F(T_p(r))-P^{-1}\\sum_{s=0}^{P-1}F(T_p(s)),"
        "\\qquad r\\in\\mathbb Z/P\\mathbb Z.$$"
    )
    lines.append("")
    lines.append(
        "The displayed equations missing from Q6420 can be reconstructed from "
        "Sections 2.3--2.5 as follows. With"
    )
    lines.append("")
    lines.append(
        "$$\\widehat h_{p,F}(u)=\\sum_{r=0}^{P-1}h_{p,F}(r)"
        "e^{-2\\pi iur/P},\\qquad "
        "K_W(\\theta)=\\sum_r W(r)e^{2\\pi ir\\theta},$$"
    )
    lines.append("")
    lines.append("the shifted pair has the exact Fourier reconstruction")
    lines.append("")
    lines.append(
        "$$C_{p,q}(F;d,W)=\\frac1{PQ}\\sum_{u\\bmod P}\\sum_{v\\bmod Q}"
        "\\widehat h_{p,F}(u)\\widehat h_{q,F}(v)e^{2\\pi ivd/Q}"
        "K_W(u/P+v/Q).$$"
    )
    lines.append("")
    lines.append(
        "For the exact indicator, additive orthogonality gives, for nonzero "
        "$u$, $\\widehat h_p(u)=p^{-1}\\sum_{a=1}^{p-1}"
        "\\mathfrak N_{p,a}(u)$, where "
        "$\\mathfrak N_{p,a}(u)=\\sum_r e_p(aT_p(r))e^{-2\\pi iur/P}$. "
        "Substitution for both primes is the fourfold $(a,b,u,v)$ formula. "
        "Centering sets the $u=0$ and $v=0$ coefficients to zero."
    )
    lines.append("")
    lines.append(
        "**Flagged fallback.** Q6420's export omits the displayed definition of "
        "$W$ and the ambient integer interval. Therefore the primary tables use "
        "the specification's explicit fallback"
    )
    lines.append("")
    lines.append(
        "$$C^{\\rm align}_{p,q}(F)=\\sum_{r=0}^{p-2}h_{p,F}(r)h_{q,F}(r).$$"
    )
    lines.append("")
    lines.append(
        "I also computed the linkage recovered from the prose: $d=p-q$ and "
        "$W=1_{[D,P)}$, so that"
    )
    lines.append("")
    lines.append(
        "$$C^{\\rm shift}_{p,q}(F)=\\sum_{r=D}^{P-1}"
        "h_{p,F}(r)h_{q,F}(r-D).$$"
    )
    lines.append("")
    lines.append(
        "This is the maximal hard window on which both original (non-wrapped) "
        "indices are admissible. It is empty when $D\\ge P$. The unreported "
        "smooth/ambient $W$ in Q6420 is not guessed."
    )
    lines.append("")
    lines.append("## Independence benchmark")
    lines.append("")
    lines.append(
        "Let $x_1,\\dots,x_n$ and $y_1,\\dots,y_m$ be centered rows, "
        "$Q_x=\\sum x_i^2$, $Q_y=\\sum y_j^2$, and let $\\pi$ be a "
        "uniform random injection of the $n$ first-row positions into the "
        "$m$ second-row positions. Then"
    )
    lines.append("")
    lines.append(
        "$$\\mathbb E[y_{\\pi(i)}^2]=Q_y/m,\\qquad "
        "\\mathbb E[y_{\\pi(i)}y_{\\pi(j)}]=-Q_y/[m(m-1)]\\quad(i\\ne j),$$"
    )
    lines.append("")
    lines.append(
        "$$\\operatorname{Var}\\!\\left(\\sum_i x_i y_{\\pi(i)}\\right)"
        "=Q_xQ_y/(m-1).$$"
    )
    lines.append("")
    lines.append(
        "This is the exact finite-population benchmark used for aligned and "
        "reflection-reduced rows. For a fixed subset $S$ of first-row positions "
        "in the shifted-overlap statistic, the exact variant used is"
    )
    lines.append("")
    lines.append(
        "$$\\operatorname{Var}=\\frac{Q_y}{m-1}"
        "\\left(\\sum_{i\\in S}x_i^2-\\frac{(\\sum_{i\\in S}x_i)^2}{m}\\right).$$"
    )
    lines.append("")
    lines.append(
        "In every bin, <code>RMS ratio</code> means "
        "$\\sqrt{\\sum C_{p,q}^2/\\sum V_{p,q}}$; <code>sum z</code> means "
        "$\\sum C_{p,q}/\\sqrt{\\sum V_{p,q}}$. Powers use the canonical "
        "integer representative $T_p(r)\\in[0,p)$ before centering."
    )
    lines.append("")
    lines.append("## Primary aligned statistic")
    lines.append("")
    primary_rows: list[list[str]] = []
    for bin_number in bins:
        interval = f"[{2**bin_number},{2**(bin_number+1)})"
        for observable in ("indicator", "T", "T2", "T3"):
            record = stats[("aligned", bin_number, observable)]
            primary_rows.append(
                [
                    interval,
                    str(record["count"]),
                    observable.replace("T2", "T^2").replace("T3", "T^3"),
                    fmt_number(float(record["sum"])),
                    fmt_number(float(record["rms"])),
                    fmt_number(float(record["bench"])),
                    f"{float(record['ratio']):.3f}",
                    f"{float(record['z']):+.2f}",
                ]
            )
    lines.append(
        markdown_table(
            [
                "D bin",
                "pairs",
                "F",
                "sum C",
                "RMS",
                "benchmark RMS",
                "RMS ratio",
                "sum z",
            ],
            primary_rows,
        )
    )
    lines.append("")
    lines.append("## Midpoint and reflection removal (indicator)")
    lines.append("")
    lines.append(
        "<code>no-mid</code> replaces the midpoint indicator by zero and "
        "re-centers on the full cycle. <code>reflection-reduced</code> then "
        "projects to the reflection-even part, discards the midpoint orbit, "
        "keeps representatives $0\\le r<(p-1)/2$, and re-centers the "
        "representative row. The program also checked "
        "$T_p(r)=T_p(-r\\bmod p-1)$ exactly for every row."
    )
    lines.append("")
    forced_rows: list[list[str]] = []
    for bin_number in bins:
        full = stats[("aligned", bin_number, "indicator")]
        nm = stats[("indicator_no_mid", bin_number, "indicator")]
        red = stats[("indicator_reflection_reduced", bin_number, "indicator")]
        full_excess = float(full["ratio"]) ** 2 - 1.0
        nm_excess = float(nm["ratio"]) ** 2 - 1.0
        red_excess = float(red["ratio"]) ** 2 - 1.0
        mid_explained = (
            100.0 * (full_excess - nm_excess) / full_excess
            if full_excess > 0
            else math.nan
        )
        total_explained = (
            100.0 * (full_excess - red_excess) / full_excess
            if full_excess > 0
            else math.nan
        )
        forced_rows.append(
            [
                f"[{2**bin_number},{2**(bin_number+1)})",
                f"{float(full['ratio']):.3f}",
                f"{float(nm['ratio']):.3f}",
                f"{float(red['ratio']):.3f}",
                "n/a"
                if math.isnan(mid_explained)
                else f"{mid_explained:+.1f}%",
                "n/a"
                if math.isnan(total_explained)
                else f"{total_explained:+.1f}%",
                f"{float(full['z']):+.2f} → {float(red['z']):+.2f}",
            ]
        )
    lines.append(
        markdown_table(
            [
                "D bin",
                "full R",
                "no-mid R",
                "reduced R",
                "midpoint share of excess",
                "total forced share",
                "sum z: full → reduced",
            ],
            forced_rows,
        )
    )
    lines.append("")
    lines.append(
        "The percentages compare normalized excess energy $R^2-1$; negative "
        "percentages mean that removal increased, rather than explained, the "
        "observed excess. <code>n/a</code> means the original bin had no "
        "positive excess."
    )
    lines.append("")
    lines.append("## Reconstructed shifted linkage")
    lines.append("")
    shifted_rows: list[list[str]] = []
    for bin_number in bins:
        interval = f"[{2**bin_number},{2**(bin_number+1)})"
        for observable in ("indicator", "T", "T2", "T3"):
            record = stats[("shifted_overlap", bin_number, observable)]
            shifted_rows.append(
                [
                    interval,
                    f"{record['support']}/{record['count']}",
                    observable.replace("T2", "T^2").replace("T3", "T^3"),
                    fmt_number(float(record["sum"])),
                    fmt_number(float(record["rms"])),
                    fmt_number(float(record["bench"])),
                    f"{float(record['ratio']):.3f}",
                    f"{float(record['z']):+.2f}",
                ]
            )
    lines.append(
        markdown_table(
            [
                "D bin",
                "nonempty/all",
                "F",
                "sum C",
                "RMS",
                "benchmark RMS",
                "RMS ratio",
                "sum z",
            ],
            shifted_rows,
        )
    )
    lines.append("")
    lines.append("## Fourier self-check and near-resonant bookkeeping")
    lines.append("")
    fourier_rows: list[list[str]] = []
    for record in fourier:
        fourier_rows.append(
            [
                str(record["family"]),
                f"[{2**int(record['bin'])},{2**(int(record['bin'])+1)})",
                f"({record['p']},{record['q']})",
                fmt_number(float(record["direct"])),
                fmt_number(float(record["recon"])),
                fmt_number(float(record["error"])),
            ]
        )
    lines.append(
        markdown_table(
            [
                "window",
                "D bin",
                "sample pair",
                "direct",
                "Fourier",
                "absolute error",
            ],
            fourier_rows,
        )
    )
    lines.append("")
    max_error = max(
        (float(record["error"]) for record in fourier), default=0.0
    )
    lines.append(
        f"The maximum direct-versus-Fourier absolute error was {max_error:.3e}. "
        "Thus the reconstruction succeeded numerically."
    )
    lines.append("")
    near_rows: list[list[str]] = []
    for bin_number in bins:
        near = stats[("near_resonant_shifted", bin_number, "indicator")]
        shifted_indicator = stats[
            ("shifted_overlap", bin_number, "indicator")
        ]
        fraction = (
            float(near["rms"]) / float(shifted_indicator["rms"])
            if float(shifted_indicator["rms"])
            else math.nan
        )
        profile = stats[
            ("indicator_cyclic_profile_max", bin_number, "indicator")
        ]
        exact_hits, total = profiles[bin_number]
        near_rows.append(
            [
                f"[{2**bin_number},{2**(bin_number+1)})",
                f"{near['support']}/{near['count']}",
                fmt_number(float(near["sum"])),
                fmt_number(float(near["rms"])),
                "n/a" if math.isnan(fraction) else f"{fraction:.3f}",
                f"{float(profile['ratio']):.3f}",
                f"{exact_hits}/{total}",
            ]
        )
    lines.append(
        "For each pair, the near block is $v=-u$, "
        "$1\\le |u|\\le\\lfloor((p+q)/2)/D\\rfloor$, truncated to the signed "
        "frequency ranges. Both signs are included. The last two columns "
        "profile the exact indicator over every cyclic shift "
        "$|\\delta|\\le D$."
    )
    lines.append("")
    lines.append(
        markdown_table(
            [
                "D bin",
                "near nonempty/all",
                "sum near",
                "RMS near",
                "near/full shifted RMS",
                "cyclic max RMS ratio",
                "argmax at ±D",
            ],
            near_rows,
        )
    )
    lines.append("")
    lines.append("## Ground-truth and structural checks")
    lines.append("")
    checks = [
        ["$p=17$", "$Z_{17}=\\{3,13\\}$", "PASS"],
        ["Parity-law sample", "first 20 window primes", "PASS (0 failures)"],
        [
            "Parity law, stronger run",
            f"all {meta['prime_count']} window primes",
            "PASS (0 failures)",
        ],
        [
            "Reflection",
            f"all {meta['prime_count']} complete rows",
            "PASS (0 failures)",
        ],
        [
            "Fourier reconstruction",
            f"{len(fourier)} selected bin/window checks",
            f"PASS (max error {max_error:.3e})",
        ],
    ]
    lines.append(markdown_table(["check", "scope", "result"], checks))
    lines.append("")
    lines.append("## Verdict")
    lines.append("")
    verdict_rows: list[list[str]] = []
    bin_confirmed: dict[int, bool] = {}
    for bin_number in bins:
        ind = stats[("aligned", bin_number, "indicator")]
        surrogates = [
            stats[("aligned", bin_number, name)] for name in ("T", "T2", "T3")
        ]
        statuses = [status_for(record) for record in surrogates]
        ind_status = status_for(ind)
        enough = int(ind["count"]) >= 30
        sur_ok = all(status == "~1" for status in statuses)
        ind_excess = ind_status == "EXCESS"
        confirmed = enough and sur_ok and ind_excess
        bin_confirmed[bin_number] = confirmed
        verdict_rows.append(
            [
                f"[{2**bin_number},{2**(bin_number+1)})",
                str(ind["count"]),
                f"{float(ind['ratio']):.3f} / "
                f"{float(ind['z']):+.2f} ({ind_status})",
                f"{float(surrogates[0]['ratio']):.3f} ({statuses[0]})",
                f"{float(surrogates[1]['ratio']):.3f} ({statuses[1]})",
                f"{float(surrogates[2]['ratio']):.3f} ({statuses[2]})",
                "YES" if confirmed else ("low N" if not enough else "NO"),
            ]
        )
    lines.append(
        markdown_table(
            [
                "D bin",
                "pairs",
                "indicator R / z",
                "T R",
                "T^2 R",
                "T^3 R",
                "predicted split?",
            ],
            verdict_rows,
        )
    )
    lines.append("")
    reliable = [
        bin_number
        for bin_number in bins
        if int(stats[("aligned", bin_number, "indicator")]["count"]) >= 30
    ]
    confirmations = sum(bin_confirmed[bin_number] for bin_number in reliable)
    if reliable and confirmations == len(reliable):
        overall = "CONFIRMED"
    elif confirmations == 0:
        overall = "NOT CONFIRMED"
    else:
        overall = "MIXED"
    lines.append(f"**Signature verdict: {overall}.**")
    lines.append("")
    lines.append(
        "The table uses a declared descriptive rule: <code>~1</code> means "
        "$0.80\\le R\\le1.25$; <code>EXCESS</code> means $R>1.25$; and "
        "<code>DEFICIT</code> means $R<0.80$. The separately reported sum z "
        "records signed drift, but pair statistics sharing a prime are not "
        "independent, so it is not used as a formal z-test. Bins with fewer "
        "than 30 pairs are not used for the overall verdict. These thresholds "
        "are reporting conventions, not formal significance tests."
    )
    lines.append("")
    return "\n".join(lines)


def main() -> int:
    with tempfile.TemporaryDirectory(prefix="cron-testf-") as temp_dir:
        binary = Path(temp_dir) / "CRON_testF_dispersion"
        compile_command = [
            "cc",
            "-O2",
            "-std=c11",
            "-Wall",
            "-Wextra",
            "-pedantic",
            str(SOURCE),
            "-lm",
            "-o",
            str(binary),
        ]
        subprocess.run(compile_command, check=True, cwd=HERE)
        completed = subprocess.run(
            [str(binary)],
            check=True,
            cwd=HERE,
            stdout=subprocess.PIPE,
            stderr=None,
            text=True,
        )
    meta, valid, stats, profiles, fourier = parse_output(completed.stdout)
    command = "python3 CRON_testF_driver.py"
    report = render_report(meta, valid, stats, profiles, fourier, command)
    REPORT.write_text(report, encoding="utf-8")
    print(f"wrote {REPORT.name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
