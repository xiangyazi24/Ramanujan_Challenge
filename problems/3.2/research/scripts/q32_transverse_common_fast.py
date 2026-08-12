#!/usr/bin/env python3
"""Sage driver for the fast exact Apéry transverse common-pair scan.

Run under Sage, for example

    sage -python problems/3.2/research/scripts/q32_transverse_common_fast.py all \
        --pmax 100000 --workers 16

The C++/FLINT worker makes every support decision modulo the prime being
scanned.  This driver compiles the worker, balances primes across processes,
resumes from per-prime atomic checkpoints, verifies the known r<=10000 high
pairs, and computes p-adic multiplicities of the hits by independent Sage
power-series arithmetic modulo p^e.

No floating-point number is used to decide a common pair or a multiplicity.
Floating point occurs only in workload balancing, elapsed-time reporting, and
human-readable p/r output.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import heapq
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
from typing import Iterable, Sequence

try:
    from sage.all import Integers, PowerSeriesRing, ZZ, prime_range
except Exception as exc:  # pragma: no cover - diagnostic outside Sage
    raise SystemExit(
        "This driver needs Sage.  Run it as `sage -python "
        "problems/3.2/research/scripts/q32_transverse_common_fast.py ...`.\n"
        f"Import error: {exc}"
    )


SCRIPT = Path(__file__).resolve()
SCRIPT_DIR = SCRIPT.parent
CPP_SOURCE = SCRIPT.with_suffix(".cpp")
DEFAULT_BINARY = SCRIPT_DIR / ".build" / "q32_transverse_common_fast"
DEFAULT_CHECKPOINTS = SCRIPT_DIR / "q32_transverse_fast_checkpoints"
ALGO_VERSION = "q32-transverse-common-fast-v1"
KNOWN_R10000_HIGH = {(17, 13), (2237, 492)}


def apery_mod(modulus: int, N: int) -> list[int]:
    """b_0,...,b_{N-1} modulo modulus, assuming every 1,...,N is a unit."""
    if N <= 0:
        return []
    b = [0] * N
    b[0] = 1 % modulus
    if N == 1:
        return b
    b[1] = 5 % modulus
    for n in range(1, N - 1):
        n2 = n * n % modulus
        n3 = n2 * n % modulus
        coeff = (34 * n3 + 51 * n2 + 27 * n + 5) % modulus
        numerator = (coeff * b[n] - n3 * b[n - 1]) % modulus
        inv = pow(n + 1, -1, modulus)
        b[n + 1] = numerator * pow(inv, 3, modulus) % modulus
    return b


def inv_sqrt_D_mod(modulus: int, N: int) -> list[int]:
    """Coefficients of (1-34t+t^2)^(-1/2) modulo modulus."""
    if N <= 0:
        return []
    q = [0] * N
    q[0] = 1 % modulus
    if N == 1:
        return q
    q[1] = 17 % modulus
    for n in range(1, N - 1):
        numerator = ((34 * n + 17) * q[n] - n * q[n - 1]) % modulus
        q[n + 1] = numerator * pow(n + 1, -1, modulus) % modulus
    return q


def xi_residue_mod_prime_power(p: int, r: int, exponent: int) -> tuple[int, int]:
    """Return (b_r, Xi_r) modulo p^exponent by an independent Sage path.

    Since r<p for every scanner hit, all recurrence denominators are units
    modulo p^e.  F(0)=1 makes the power-series inverse valid over Z/p^e Z.
    """
    if not (0 <= r < p):
        raise ValueError("valuation verifier requires r < p")
    modulus = int(p) ** int(exponent)
    N = r + 1
    b = apery_mod(modulus, N)
    q = inv_sqrt_D_mod(modulus, N)
    ring = Integers(modulus)
    series = PowerSeriesRing(ring, "t", default_prec=N)
    F = series([ring(x) for x in b]).add_bigoh(N)
    Q = series([ring(x) for x in q]).add_bigoh(N)
    Finv = (~F).add_bigoh(N)
    G = (Finv * Finv * Q).add_bigoh(N)
    xi = ring(-1)
    for m in range(1, N):
        xi -= ring(5) * G[m] * ring(b[m - 1])
    return int(ring(b[r])), int(xi)


def residue_valuation(residue: int, p: int) -> int:
    if residue == 0:
        raise ValueError("zero residue does not determine a finite valuation")
    v = 0
    while residue % p == 0:
        residue //= p
        v += 1
    return v


def exact_hit_valuation(p: int, r: int) -> dict[str, int]:
    """Determine v_p(b_r), v_p(Xi_r), and their minimum exactly.

    The exponent doubles until both residues are nonzero.  b_r is a positive
    integer and Xi_r=-1-5*sum_{m<=r} g_m b_{m-1}<0 over Z, so both valuations
    are finite; the loop therefore terminates.  Once x mod p^e is nonzero,
    with v_p(x)<e, its residue has exactly the true p-adic valuation.
    """
    exponent = 2
    vb = None
    vx = None
    while vb is None or vx is None:
        br, xir = xi_residue_mod_prime_power(p, r, exponent)
        if vb is None and br != 0:
            vb = residue_valuation(br, p)
        if vx is None and xir != 0:
            vx = residue_valuation(xir, p)
        exponent *= 2
    return {"vp_b": int(vb), "vp_xi": int(vx), "vp_gcd": int(min(vb, vx))}


def self_test() -> None:
    """Small exact characteristic-zero test of the formulas used by both paths."""
    b = [ZZ(1), ZZ(5)]
    for n in range(1, 2):
        numerator = (34 * n**3 + 51 * n**2 + 27 * n + 5) * b[n] - n**3 * b[n - 1]
        denominator = (n + 1) ** 3
        assert numerator % denominator == 0
        b.append(numerator // denominator)
    assert b == [1, 5, 73]

    q = [ZZ(1), ZZ(17)]
    for n in range(1, 2):
        numerator = (34 * n + 17) * q[n] - n * q[n - 1]
        assert numerator % (n + 1) == 0
        q.append(numerator // (n + 1))
    assert q == [1, 17, 433]

    S = PowerSeriesRing(ZZ, "t", default_prec=3)
    F = S(b).add_bigoh(3)
    Q = S(q).add_bigoh(3)
    G = ((~F) * (~F) * Q).add_bigoh(3)
    assert [ZZ(G[i]) for i in range(3)] == [1, 7, 192]
    xi = ZZ(-1)
    xis = [xi]
    for m in range(1, 3):
        xi -= 5 * ZZ(G[m]) * b[m - 1]
        xis.append(xi)
    assert xis == [-1, -36, -4836]


def compile_worker(binary: Path, *, force: bool = False) -> None:
    binary = binary.resolve()
    binary.parent.mkdir(parents=True, exist_ok=True)
    if binary.exists() and not force and binary.stat().st_mtime >= CPP_SOURCE.stat().st_mtime:
        return

    cxx = os.environ.get("CXX") or shutil.which("g++") or shutil.which("c++")
    if not cxx:
        raise SystemExit("No C++ compiler found (set CXX or install g++).")
    cmd = [cxx, "-O3", "-DNDEBUG", "-std=c++17", str(CPP_SOURCE), "-o", str(binary)]

    sage_local = os.environ.get("SAGE_LOCAL")
    if sage_local:
        include = Path(sage_local) / "include"
        lib = Path(sage_local) / "lib"
        cmd += [f"-I{include}", f"-L{lib}", f"-Wl,-rpath,{lib}"]
    else:
        try:
            flags = subprocess.check_output(
                ["pkg-config", "--cflags", "--libs", "flint"], text=True
            ).split()
            cmd += flags
        except Exception:
            pass
    cmd += ["-lflint", "-lgmp", "-lmpfr"]
    print("BUILD", " ".join(map(str, cmd)), flush=True)
    subprocess.run(cmd, check=True)


def read_checkpoint(path: Path) -> dict:
    data = json.loads(path.read_text())
    if data.get("version") != ALGO_VERSION:
        raise RuntimeError(f"checkpoint {path} has wrong version: {data.get('version')}")
    return data


def completed_primes(checkpoint_dir: Path) -> set[int]:
    done = set()
    if not checkpoint_dir.exists():
        return done
    for path in checkpoint_dir.glob("p_*.json"):
        try:
            data = read_checkpoint(path)
        except Exception:
            continue
        done.add(int(data["p"]))
    return done


def balanced_bins(primes: Sequence[int], workers: int) -> list[list[int]]:
    """Greedy LPT partition using p*log(p) as a conservative FLINT cost proxy."""
    workers = max(1, min(workers, len(primes) or 1))
    heap: list[tuple[float, int, list[int]]] = [(0.0, i, []) for i in range(workers)]
    heapq.heapify(heap)
    for p in sorted(primes, reverse=True):
        load, i, values = heapq.heappop(heap)
        values.append(p)
        load += p * max(1.0, math.log2(p))
        heapq.heappush(heap, (load, i, values))
    bins = [[] for _ in range(workers)]
    while heap:
        _load, i, values = heapq.heappop(heap)
        bins[i] = sorted(values)
    return [values for values in bins if values]


def run_scan(binary: Path, checkpoint_dir: Path, pmax: int, workers: int) -> None:
    if pmax < 2:
        return
    all_primes = [int(p) for p in prime_range(2, pmax + 1)]
    done = completed_primes(checkpoint_dir)
    pending = [p for p in all_primes if p not in done]
    print(
        f"SCAN p<= {pmax}: total={len(all_primes)} done={len(done & set(all_primes))} "
        f"pending={len(pending)} workers={workers}",
        flush=True,
    )
    if not pending:
        return

    jobs_dir = checkpoint_dir / ".jobs"
    jobs_dir.mkdir(parents=True, exist_ok=True)
    processes: list[tuple[subprocess.Popen, Path]] = []
    for index, values in enumerate(balanced_bins(pending, workers)):
        job = jobs_dir / f"job_{index:03d}.txt"
        job.write_text("".join(f"{p}\n" for p in values))
        cmd = [
            str(binary),
            "--prime-file",
            str(job),
            "--pmax",
            str(pmax),
            "--checkpoint-dir",
            str(checkpoint_dir),
        ]
        processes.append((subprocess.Popen(cmd), job))

    failures = []
    for proc, job in processes:
        code = proc.wait()
        if code:
            failures.append((job, code))
    if failures:
        raise RuntimeError(f"worker failures: {failures}")


def collect_summary(checkpoint_dir: Path, pmax: int, *, allow_partial: bool = False) -> dict:
    expected_primes = [int(p) for p in prime_range(2, pmax + 1)]
    records: dict[int, dict] = {}
    for path in checkpoint_dir.glob("p_*.json"):
        try:
            row = read_checkpoint(path)
        except Exception:
            continue
        p = int(row["p"])
        if p <= pmax:
            records[p] = row

    missing = [p for p in expected_primes if p not in records]
    if missing and not allow_partial:
        raise RuntimeError(
            f"scan incomplete: {len(missing)} missing primes; first={missing[:10]}"
        )

    pairs = sorted(
        (p, int(r))
        for p, row in records.items()
        for r in row.get("pairs", [])
    )
    detailed = []
    for p, r in pairs:
        valuation = exact_hit_valuation(p, r)
        detailed.append(
            {
                "p": p,
                "r": r,
                "ratio": str(Fraction(p, r)),
                "ratio_float": p / r,
                "barrier_violation": p > 5 * r,
                **valuation,
            }
        )

    max_pair = None
    if pairs:
        p, r = max(pairs, key=lambda pr: Fraction(pr[0], pr[1]))
        max_pair = {
            "p": p,
            "r": r,
            "ratio": str(Fraction(p, r)),
            "ratio_float": p / r,
        }

    baseline_observed = {(p, r) for p, r in pairs if r <= 10000}
    baseline_checked = pmax >= 100000 and not missing
    baseline_ok = baseline_checked and baseline_observed == KNOWN_R10000_HIGH
    if baseline_checked and not baseline_ok:
        raise RuntimeError(
            "r<=10000 high-pair regression: "
            f"observed={sorted(baseline_observed)} expected={sorted(KNOWN_R10000_HIGH)}"
        )

    summary = {
        "version": ALGO_VERSION,
        "pmax": pmax,
        "prime_count_expected": len(expected_primes),
        "prime_count_complete": len(records),
        "missing_primes": missing,
        "pair_count": len(pairs),
        "pairs": detailed,
        "max_ratio_pair": max_pair,
        "barrier_violations": [row for row in detailed if row["barrier_violation"]],
        "known_r10000_reference": sorted([list(x) for x in KNOWN_R10000_HIGH]),
        "known_r10000_reference_checked": baseline_checked,
        "known_r10000_reference_ok": baseline_ok,
        "total_b_zero_count": sum(int(row.get("b_zero_count", 0)) for row in records.values()),
        "total_worker_seconds": sum(float(row.get("seconds", 0.0)) for row in records.values()),
    }
    return summary


def print_summary(summary: dict) -> None:
    print(json.dumps(summary, indent=2, sort_keys=True))
    if summary["max_ratio_pair"]:
        row = summary["max_ratio_pair"]
        print(
            f"MAX p/r: p={row['p']} r={row['r']} ratio={row['ratio']} "
            f"(~{row['ratio_float']:.12f})"
        )
    print("COMMON p>r PAIRS:")
    for row in summary["pairs"]:
        print(
            f"  p={row['p']} r={row['r']} p/r={row['ratio']} "
            f"vp(b)={row['vp_b']} vp(Xi)={row['vp_xi']} "
            f"vp(gcd)={row['vp_gcd']}"
        )
    if summary["barrier_violations"]:
        print("P>5R COUNTEREXAMPLES:")
        for row in summary["barrier_violations"]:
            print(f"  p={row['p']} r={row['r']} ratio={row['ratio']}")
    else:
        print("P>5R COUNTEREXAMPLES: none in the completed scan range")


def parse_args(argv: Sequence[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "command",
        nargs="?",
        choices=("build", "scan", "summarize", "all", "self-test"),
        default="all",
    )
    parser.add_argument("--pmax", type=int, default=100000)
    parser.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) // 2))
    parser.add_argument("--binary", type=Path, default=DEFAULT_BINARY)
    parser.add_argument("--checkpoint-dir", type=Path, default=DEFAULT_CHECKPOINTS)
    parser.add_argument("--force-build", action="store_true")
    parser.add_argument("--allow-partial", action="store_true")
    parser.add_argument("--summary-json", type=Path)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(sys.argv[1:] if argv is None else argv)
    self_test()
    if args.command == "self-test":
        print("self-test: OK")
        return
    if args.command in ("build", "scan", "all"):
        compile_worker(args.binary, force=args.force_build)
        if args.command == "build":
            return
    if args.command in ("scan", "all"):
        run_scan(args.binary.resolve(), args.checkpoint_dir.resolve(), args.pmax, args.workers)
    if args.command in ("summarize", "all"):
        summary = collect_summary(
            args.checkpoint_dir.resolve(), args.pmax, allow_partial=args.allow_partial
        )
        print_summary(summary)
        output = args.summary_json or (args.checkpoint_dir / "summary.json")
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
        print(f"summary written to {output}")


if __name__ == "__main__":
    main()
