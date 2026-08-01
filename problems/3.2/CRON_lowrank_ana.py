#!/usr/bin/env -S sage -python
"""Exact clock-Fourier low-rank tests from Q6523 section 3.

The script computes the mirror-subtracted collision channel ``Y`` and the
mandatory linear determinant control ``Z``.  It then scans all 24 affine
twists from Q6523, four dispersed seeds, h0 in {1,3,7}, s in {1,2,4,8}, and
every dyadic L allowed by h0+2L-2 <= floor(sqrt(p)).  Ranks are exact over
two independently chosen finite fields F_q containing a primitive p-th root
of unity.

Default full run::

    sage -python CRON_lowrank_ana.py

Output is JSON Lines so every raw rank is machine-readable.
"""

from __future__ import annotations

import argparse
import json
import math
import time
from collections import defaultdict
from dataclasses import dataclass
from itertools import count
from typing import Callable, Iterable

import numpy as np
from sympy import isprime, primitive_root

try:
    from sage.all import GF, matrix
except ImportError as exc:  # pragma: no cover - exercised only outside Sage
    raise SystemExit(
        "CRON_lowrank_ana.py requires SageMath; run "
        "`sage -python CRON_lowrank_ana.py`."
    ) from exc


PRIMES = (1009, 3001, 10007)
H0S = (1, 3, 7)
SVALS = (1, 2, 4, 8)
TWIST_AS = (1, -1, 2, -2, 3, -3, 6, -6)
TWIST_BS = (0, 1, -1)
TWISTS = tuple((a, b) for a in TWIST_AS for b in TWIST_BS)
OBSERVABLES = ("Y", "Z")


@dataclass(frozen=True)
class OrbitData:
    b: np.ndarray
    c: np.ndarray
    points: tuple[tuple[int, int], ...]


def projective_key(x: int, y: int, p: int) -> tuple[int, int]:
    if x % p:
        return (1, y * pow(x, -1, p) % p)
    return (0, 1)


def orbit_projective(p: int) -> OrbitData:
    """Apéry companion pair and its projective orbit, exactly as in Q6523."""
    b = np.zeros(p - 1, dtype=np.int64)
    c = np.zeros(p - 1, dtype=np.int64)
    b[0], b[1] = 1, 5 % p
    c[0], c[1] = 0, 1
    for n in range(2, p - 1):
        A = (34 * n**3 - 51 * n**2 + 27 * n - 5) % p
        B = (n - 1) ** 3 % p
        inv = pow(n**3 % p, -1, p)
        b[n] = (A * int(b[n - 1]) - B * int(b[n - 2])) * inv % p
        c[n] = (A * int(c[n - 1]) - B * int(c[n - 2])) * inv % p
    points = tuple(projective_key(int(b[n]), int(c[n]), p) for n in range(p - 1))
    return OrbitData(b=b, c=c, points=points)


def collision_and_delta_arrays(p: int, orbit: OrbitData, Hmax: int):
    """Return dense clock arrays for C^circ and the linear Delta control.

    ``Delta`` is represented by its canonical Python residue in 0,...,p-1,
    exactly as obtained by the ``% p`` arithmetic in the Q6523 skeleton,
    before reducing the resulting integer Fourier sums modulo q.
    """
    collision = np.zeros((Hmax + 1, p), dtype=np.int64)
    delta = np.zeros((Hmax + 1, p), dtype=np.int64)
    counts = {}
    for h in range(1, Hmax + 1):
        valid = p - 1 - h
        br = orbit.b[:valid]
        cr = orbit.c[:valid]
        bh = orbit.b[h : h + valid]
        ch = orbit.c[h : h + valid]
        # Products are <p^2 (p<=10007), safely inside int64.
        dh = (br * ch - bh * cr) % p
        equal_from_delta = dh == 0
        equal_from_points = np.fromiter(
            (orbit.points[r] == orbit.points[r + h] for r in range(valid)),
            dtype=np.bool_,
            count=valid,
        )
        if not np.array_equal(equal_from_delta, equal_from_points):
            raise AssertionError(("projective/determinant collision mismatch", p, h))
        # Store by the centered clock z=2r+h+1, not by r.  This is the
        # Fourier coordinate fixed in Q6523 section 3.1.
        zcoords = (2 * np.arange(valid, dtype=np.int64) + h + 1) % p
        delta[h, zcoords] = dh
        collision[h, zcoords] = equal_from_delta.astype(np.int64)

        removed = 0
        if h % 2 == 0:
            mirror_r = (p - 1 - h) // 2
            mirror_z = (2 * mirror_r + h + 1) % p
            if not (0 <= mirror_r < valid and mirror_z == 0 and collision[h, mirror_z] == 1):
                raise AssertionError(("missing forced mirror collision", p, h, mirror_r))
            collision[h, mirror_z] = 0
            removed = 1
        counts[h] = {
            "primitive_collisions": int(collision[h].sum()),
            "forced_mirror_removed": removed,
        }
    return collision, delta, counts


def two_aux_primes(p: int, start: int = 2) -> tuple[int, int]:
    out = []
    for m in count(start):
        q = m * p + 1
        if isprime(q):
            out.append(q)
            if len(out) == 2:
                return out[0], out[1]
    raise AssertionError("unreachable")


def phi_iter(k: int, a: int, b: int, p: int, n: int) -> int:
    for _ in range(n):
        k = (a * k + b) % p
    return k


def orbit_period(k: int, a: int, b: int, p: int, cap: int):
    x = k
    for n in range(1, cap + 1):
        x = (a * x + b) % p
        if x == k:
            return n
    return None


def effective_s(seeds: Iterable[int], a: int, b: int, p: int, s: int) -> int:
    # This is the conservative multi-seed rule in the Q6523 skeleton.
    periods = [orbit_period(k, a, b, p, 2 * s + 2) for k in seeds]
    finite = [period for period in periods if period is not None]
    return min([s] + finite) if finite else s


def dyadic_windows(Hmax: int, h0: int) -> tuple[int, ...]:
    values = []
    L = 2
    while h0 + 2 * L - 2 <= Hmax:
        values.append(L)
        L *= 2
    return tuple(values)


def needed_frequencies(p: int, seeds: tuple[int, ...], windows: dict[int, tuple[int, ...]]):
    """Collect every frequency used by any requested diagonal/bispectral matrix."""
    needed = {0}
    max_L = max((max(values) for values in windows.values() if values), default=0)
    for a, b in TWISTS:
        for k0 in seeds:
            for n in range(max(0, 2 * max_L - 1)):
                needed.add(phi_iter(k0, a, b, p, n))
        for s in SVALS:
            seff = effective_s(seeds, a, b, p, s)
            for k0 in seeds:
                for n in range(2 * seff - 1):
                    needed.add(phi_iter(k0, a, b, p, n))
    return tuple(sorted(needed))


def exact_dot_mod(row: np.ndarray, omega_powers: np.ndarray, k: int, p: int, q: int) -> int:
    total = 0
    for z, coefficient in enumerate(row):
        if coefficient:
            total += int(coefficient) * int(omega_powers[(k * z) % p])
    return total % q


def fourier_tables(
    p: int,
    q: int,
    collision: np.ndarray,
    delta: np.ndarray,
    frequencies: tuple[int, ...],
    chunk_size: int = 128,
):
    """Compute the requested DFT columns exactly, using checked float GEMM.

    Every dot product is an integer of absolute value below 2^53 for the
    mandated ranges, so binary64 represents every product and partial sum
    exactly.  Direct Python-integer spot checks below guard this optimization.
    """
    generator = int(primitive_root(q))
    omega = pow(generator, (q - 1) // p, q)
    if pow(omega, p, q) != 1 or omega == 1:
        raise AssertionError(("not a primitive p-th root", p, q, omega))

    omega_powers = np.empty(p, dtype=np.int64)
    omega_powers[0] = 1
    for z in range(1, p):
        omega_powers[z] = int(omega_powers[z - 1]) * omega % q
    assert int(omega_powers[-1]) * omega % q == 1

    # Worst possible unsigned Delta dot product.  This also bounds every
    # BLAS partial sum because all entries are nonnegative.
    worst = (p - 1) * (p - 1) * (q - 1)
    if worst >= 2**53:
        raise AssertionError(("binary64 exactness bound exceeded", p, q, worst))

    stacked = np.vstack((collision, delta))
    stacked_float = stacked.astype(np.float64)
    nrows = stacked.shape[0]
    result = np.empty((nrows, len(frequencies)), dtype=np.int64)
    zcoords = np.arange(p, dtype=np.int64)[:, None]
    for start in range(0, len(frequencies), chunk_size):
        stop = min(start + chunk_size, len(frequencies))
        ks = np.asarray(frequencies[start:stop], dtype=np.int64)[None, :]
        exponents = (zcoords * ks) % p
        weights = omega_powers[exponents]
        raw = stacked_float @ weights.astype(np.float64)
        rounded = np.rint(raw)
        if not np.array_equal(raw, rounded):
            raise AssertionError(("nonintegral binary64 Fourier dot", p, q, start))
        result[:, start:stop] = rounded.astype(np.int64) % q

    Hrows = collision.shape[0]
    Y = result[:Hrows]
    Z = result[Hrows:]
    index = {k: i for i, k in enumerate(frequencies)}

    # Independent exact-integer checks at boundary and interior samples.
    sample_h = sorted({1, min(3, Hrows - 1), Hrows - 1})
    sample_k = sorted({frequencies[0], frequencies[len(frequencies) // 2], frequencies[-1]})
    for h in sample_h:
        for k in sample_k:
            column = index[k]
            assert int(Y[h, column]) == exact_dot_mod(collision[h], omega_powers, k, p, q)
            assert int(Z[h, column]) == exact_dot_mod(delta[h], omega_powers, k, p, q)
    return {"Y": Y, "Z": Z}, index, omega


def rank_mod(A, q: int) -> int:
    """Exact rank over F_q using Sage's finite-field matrix backend."""
    A = np.asarray(A, dtype=np.int64)
    if A.ndim != 2:
        raise ValueError("rank_mod expects a matrix")
    if A.size == 0:
        return 0
    return int(matrix(GF(q), A.tolist()).rank())


def diag_hankel(values: np.ndarray, index: dict[int, int], p: int, k0: int, h0: int, L: int, a: int, b: int):
    sequence = np.asarray(
        [values[h0 + n, index[phi_iter(k0, a, b, p, n)]] for n in range(2 * L - 1)],
        dtype=np.int64,
    )
    indices = np.add.outer(np.arange(L), np.arange(L))
    return sequence[indices]


def bispectral_hankel(
    values: np.ndarray,
    index: dict[int, int],
    p: int,
    seeds: tuple[int, ...],
    h0: int,
    L: int,
    s: int,
    a: int,
    b: int,
):
    seff = effective_s(seeds, a, b, p, s)
    rows = np.empty((len(seeds) * seff * L, seff * L), dtype=np.int64)
    uv = np.add.outer(np.arange(L), np.arange(L))
    for seed_index, k0 in enumerate(seeds):
        for aa in range(seff):
            row_slice = slice((seed_index * seff + aa) * L, (seed_index * seff + aa + 1) * L)
            for bb in range(seff):
                kk = phi_iter(k0, a, b, p, aa + bb)
                block = values[h0 + uv, index[kk]]
                col_slice = slice(bb * L, (bb + 1) * L)
                rows[row_slice, col_slice] = block
    return rows, seff


def self_test():
    # Exact rank sanity checks, including a dependent rectangular matrix.
    assert rank_mod([[1, 2], [2, 4]], 101) == 1
    assert rank_mod([[1, 0, 3], [0, 1, 4]], 101) == 2

    p = 101
    Hmax = math.isqrt(p)
    orbit = orbit_projective(p)
    collision, delta, counts = collision_and_delta_arrays(p, orbit, Hmax)
    assert counts[2]["primitive_collisions"] == 0
    q = two_aux_primes(p)[0]
    frequencies = (0, 1, 2, p - 1)
    tables, index, omega = fourier_tables(p, q, collision, delta, frequencies, chunk_size=4)
    assert pow(omega, p, q) == 1 and omega != 1
    for h in range(1, Hmax + 1):
        assert int(tables["Y"][h, index[0]]) == counts[h]["primitive_collisions"] % q
        # The mirror-subtracted event array is even in centered clock z.
        assert int(tables["Y"][h, index[1]]) == int(tables["Y"][h, index[p - 1]])
    print(json.dumps({"record": "self_test", "test": "analysis", "status": "OK"}))


def rank_records_for_field(
    p: int,
    q: int,
    tables: dict[str, np.ndarray],
    index: dict[int, int],
    seeds: tuple[int, ...],
    windows: dict[int, tuple[int, ...]],
):
    records = []
    for observable in OBSERVABLES:
        values = tables[observable]
        for h0 in H0S:
            for L in windows[h0]:
                for a, b in TWISTS:
                    diag_ranks = []
                    for k0 in seeds:
                        Hd = diag_hankel(values, index, p, k0, h0, L, a, b)
                        diag_ranks.append(rank_mod(Hd, q))
                    for s in SVALS:
                        Hb, seff = bispectral_hankel(
                            values, index, p, seeds, h0, L, s, a, b
                        )
                        rank_bis = rank_mod(Hb, q)
                        rmax = seff * L
                        record = {
                            "record": "analysis_rank",
                            "observable": observable,
                            "p": p,
                            "q": q,
                            "h0": h0,
                            "L": L,
                            "a": a,
                            "b": b,
                            "s": s,
                            "seff": seff,
                            "rank_diag": diag_ranks[0],
                            "diag_ranks": diag_ranks,
                            "rank_bis": rank_bis,
                            "rmax": rmax,
                            "rho": rank_bis / rmax,
                            "rho2": ((rank_bis - 2) / (rmax - 2)) if rmax > 2 else 0.0,
                            "full_bispectral_rank": rank_bis == rmax,
                            "full_diagonal_rank_all_seeds": all(rank == L for rank in diag_ranks),
                        }
                        records.append(record)
                        print(json.dumps(record, sort_keys=True), flush=True)
    return records


def paired_q_records(records: list[dict], qs: tuple[int, int]):
    by_key = defaultdict(dict)
    for record in records:
        key = tuple(
            record[name]
            for name in ("observable", "p", "h0", "L", "a", "b", "s")
        )
        by_key[key][record["q"]] = record
    paired = []
    for key, qmap in sorted(by_key.items()):
        if set(qmap) != set(qs):
            raise AssertionError(("missing q specialization", key, sorted(qmap)))
        left, right = qmap[qs[0]], qmap[qs[1]]
        stable = (
            left["rank_bis"] == right["rank_bis"]
            and left["diag_ranks"] == right["diag_ranks"]
            and left["seff"] == right["seff"]
        )
        item = {
            "record": "analysis_q_pair",
            "observable": key[0],
            "p": key[1],
            "h0": key[2],
            "L": key[3],
            "a": key[4],
            "b": key[5],
            "s": key[6],
            "q1": qs[0],
            "q2": qs[1],
            "seff": left["seff"],
            "rank_bis": [left["rank_bis"], right["rank_bis"]],
            "diag_ranks": [left["diag_ranks"], right["diag_ranks"]],
            "rmax": left["rmax"],
            "specialization_stable": stable,
        }
        paired.append(item)
        print(json.dumps(item, sort_keys=True), flush=True)
    return paired


def summarize_prime(p: int, paired: list[dict], windows: dict[int, tuple[int, ...]]):
    summaries = []
    for observable in OBSERVABLES:
        subset = [r for r in paired if r["observable"] == observable]
        unstable = [r for r in subset if not r["specialization_stable"]]
        full_bis = [r for r in subset if r["rank_bis"][0] == r["rmax"]]
        full_diag = [
            r
            for r in subset
            if all(rank == r["L"] for rank in r["diag_ranks"][0])
        ]

        # A threshold-independent sufficient form of Q6523's strong-negative
        # test: maximal rank on all of the final three scales for every actual
        # long orbit (seff=s), excluding the two stated mirror sanity twists.
        curves = defaultdict(list)
        for r in subset:
            curves[(r["h0"], r["a"], r["b"], r["s"])].append(r)
        eligible_curves = []
        maximal_last3 = []
        for key, curve in curves.items():
            curve.sort(key=lambda r: r["L"])
            h0, a, b, s = key
            if len(curve) < 3 or s != 8 or (a, b) in ((1, 0), (-1, 0)):
                continue
            if any(r["seff"] != s for r in curve):
                continue
            eligible_curves.append(key)
            last3 = curve[-3:]
            maximal_last3.append(
                all(
                    r["specialization_stable"]
                    and r["rank_bis"][0] == r["rmax"]
                    for r in last3
                )
            )

        summary = {
            "record": "analysis_prime_summary",
            "observable": observable,
            "p": p,
            "windows": {str(h0): list(windows[h0]) for h0 in H0S},
            "rank_records": len(subset),
            "specialization_unstable": len(unstable),
            "full_bispectral_records": len(full_bis),
            "full_diagonal_records": len(full_diag),
            "eligible_long_orbit_last3_curves": len(eligible_curves),
            "maximal_rank_on_every_eligible_last3_curve": (
                bool(eligible_curves) and all(maximal_last3)
            ),
        }
        summaries.append(summary)
        print(json.dumps(summary, sort_keys=True), flush=True)
    return summaries


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--prime", type=int, action="append", dest="primes", help="run only this prime (repeatable)"
    )
    parser.add_argument("--self-test-only", action="store_true")
    parser.add_argument(
        "--chunk-size", type=int, default=128, help="Fourier column chunk size (default: 128)"
    )
    return parser.parse_args()


def main():
    args = parse_args()
    if args.chunk_size < 1:
        raise SystemExit("--chunk-size must be positive")
    self_test()
    if args.self_test_only:
        return
    primes = tuple(args.primes or PRIMES)
    unknown = set(primes) - set(PRIMES)
    if unknown:
        raise SystemExit(f"full Q6523 implementation is fixed to primes {PRIMES}; unknown: {sorted(unknown)}")

    started = time.time()
    global_summaries = []
    for p in primes:
        p_started = time.time()
        Hmax = math.isqrt(p)
        windows = {h0: dyadic_windows(Hmax, h0) for h0 in H0S}
        if any(len(windows[h0]) < 3 for h0 in H0S):
            raise AssertionError(("fewer than three admissible dyadic scales", p, windows))
        seeds = (1, p // 7, 2 * p // 7, 3 * p // 7)
        if len(set(seeds)) != 4 or 0 in seeds:
            raise AssertionError(("bad dispersed seeds", p, seeds))
        frequencies = needed_frequencies(p, seeds, windows)
        qs = two_aux_primes(p)
        print(
            json.dumps(
                {
                    "record": "analysis_prime_start",
                    "p": p,
                    "Hmax": Hmax,
                    "windows": {str(h0): list(windows[h0]) for h0 in H0S},
                    "seeds": seeds,
                    "twists": len(TWISTS),
                    "frequencies": len(frequencies),
                    "aux_primes": qs,
                },
                sort_keys=True,
            ),
            flush=True,
        )

        orbit = orbit_projective(p)
        collision, delta, counts = collision_and_delta_arrays(p, orbit, Hmax)
        print(
            json.dumps(
                {
                    "record": "collision_counts",
                    "p": p,
                    "Hmax": Hmax,
                    "total_primitive_collisions": int(collision.sum()),
                    "per_h": counts,
                },
                sort_keys=True,
            ),
            flush=True,
        )

        all_records = []
        for q in qs:
            field_started = time.time()
            tables, index, omega = fourier_tables(
                p, q, collision, delta, frequencies, chunk_size=args.chunk_size
            )
            print(
                json.dumps(
                    {
                        "record": "fourier_table",
                        "p": p,
                        "q": q,
                        "omega": omega,
                        "frequencies": len(frequencies),
                        "status": "exact_spot_checks_OK",
                        "seconds": round(time.time() - field_started, 3),
                    },
                    sort_keys=True,
                ),
                flush=True,
            )
            all_records.extend(rank_records_for_field(p, q, tables, index, seeds, windows))

        paired = paired_q_records(all_records, qs)
        summaries = summarize_prime(p, paired, windows)
        global_summaries.extend(summaries)
        print(
            json.dumps(
                {
                    "record": "analysis_prime_done",
                    "p": p,
                    "seconds": round(time.time() - p_started, 3),
                },
                sort_keys=True,
            ),
            flush=True,
        )

    final = {
        "record": "analysis_final",
        "primes": list(primes),
        "all_specializations_stable": all(
            summary["specialization_unstable"] == 0 for summary in global_summaries
        ),
        "Y_maximal_last3_all_primes": all(
            summary["maximal_rank_on_every_eligible_last3_curve"]
            for summary in global_summaries
            if summary["observable"] == "Y"
        ),
        "Z_maximal_last3_all_primes": all(
            summary["maximal_rank_on_every_eligible_last3_curve"]
            for summary in global_summaries
            if summary["observable"] == "Z"
        ),
        # Q6523 demands at least four primary primes for a strong positive;
        # the campaign spec fixes three, so a positive cannot be certified in
        # this mandated run even before inspecting ranks.
        "strong_positive_possible_with_prime_count": len(primes) >= 4,
        "elapsed_seconds": round(time.time() - started, 3),
    }
    print(json.dumps(final, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
