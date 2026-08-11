#!/usr/bin/env python3
"""Q7311: exact centered CRT/Fourier verifier for Apéry zero sets.

The exact arithmetic uses the integer-centered residues

    u_p(a) = p * 1_{Z_p}(a) - |Z_p|.

For distinct p,q,r and N=X^2, put H_S for the number of n<N hit by every
prime in S and A_S = sum_{n<N} prod_{p in S} u_p(n).  Then

    p q r H_{pqr}
      = N z_p z_q z_r
      + A_p z_q z_r + A_q z_p z_r + A_r z_p z_q
      + A_{pq} z_r + A_{pr} z_q + A_{qr} z_p
      + A_{pqr}.

All reconstruction assertions are integer equalities.  The primitive mixed
correlation is C_{pqr}=A_{pqr}/(p q r).  Its Fourier expansion uses only
frequencies nonzero in every local coordinate.  Conjugate frequency pairs are
ranked with a rigorous unseen-mode envelope and 80/140-digit recomputation.

The script is dependency-free except for mpmath, used only for the displayed
high-precision frequency values.  All moment and centered-correlation tables
are exact Python integer/Fraction computations.
"""

from __future__ import annotations

import argparse
import cmath
import csv
import heapq
import io
import itertools
import math
import random
import statistics
import time
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, MutableMapping, Sequence, Tuple


ANSWER_LINE = "ANSWER Q7311 93679a56"
DEFAULT_XS = (32, 64, 128, 256)


@dataclass(frozen=True)
class PrimeData:
    p: int
    zeros: Tuple[int, ...]
    bits: int
    hits: int
    center_hit: bool

    @property
    def z(self) -> int:
        return len(self.zeros)


@dataclass
class InstanceAnalysis:
    X: int
    N: int
    label: str
    pdata: List[PrimeData]
    active: List[PrimeData]
    K: List[int]
    F1: int
    F2: int
    F3: int
    M3: int
    pair_hits: Dict[Tuple[int, int], int]
    triple_rows: List[dict]
    centered_degree: Tuple[Fraction, Fraction, Fraction, Fraction]
    pair_centered_degree: Tuple[Fraction, Fraction, Fraction]
    one_centered_degree: Tuple[Fraction, Fraction]
    metrics: dict


def sieve_primes(n: int) -> List[int]:
    mark = bytearray(b"\x01") * (n + 1)
    if n >= 0:
        mark[0] = 0
    if n >= 1:
        mark[1] = 0
    for d in range(2, math.isqrt(n) + 1):
        if mark[d]:
            start = d * d
            mark[start : n + 1 : d] = b"\x00" * (((n - start) // d) + 1)
    return [i for i in range(2, n + 1) if mark[i]]


def primes_in_shell(X: int) -> List[int]:
    return [p for p in sieve_primes(2 * X) if X < p <= 2 * X]


def apery_mod_recurrence(p: int) -> List[int]:
    """b_0,...,b_{p-1} mod p from the three-term recurrence."""
    b = [0] * p
    b[0] = 1
    if p == 1:
        return b
    b[1] = 5 % p
    for n in range(1, p - 1):
        Pn = (2 * n + 1) * (17 * n * n + 17 * n + 5)
        num = (Pn * b[n] - n**3 * b[n - 1]) % p
        den = pow(n + 1, 3, p)
        b[n + 1] = num * pow(den, -1, p) % p
    return b


def apery_mod_binomial(p: int) -> List[int]:
    """Independent direct binomial sum for every n<p."""
    out: List[int] = []
    for n in range(p):
        term = 1
        total = 1
        for j in range(n):
            inv = pow(j + 1, -1, p)
            base = (n - j) * (n + j + 1) % p
            base = base * inv % p
            base = base * inv % p
            term = term * base % p
            term = term * base % p
            total = (total + term) % p
        out.append(total)
    return out


def actual_zero_set(p: int) -> Tuple[int, ...]:
    rec = apery_mod_recurrence(p)
    direct = apery_mod_binomial(p)
    if rec != direct:
        bad = next(i for i, (a, b) in enumerate(zip(rec, direct)) if a != b)
        raise AssertionError(f"p={p}: recurrence/binomial mismatch at n={bad}: {rec[bad]} != {direct[bad]}")
    zeros = tuple(i for i, value in enumerate(rec) if value == 0)
    reflected = tuple(sorted((p - 1 - r) % p for r in zeros))
    if zeros != reflected:
        raise AssertionError(f"p={p}: zero set is not reflection invariant: {zeros}")
    if 0 in zeros or p - 1 in zeros:
        raise AssertionError(f"p={p}: unexpected endpoint zero: {zeros}")
    return zeros


def hit_bits(zeros: Sequence[int], p: int, N: int) -> int:
    bits = 0
    for r in zeros:
        for n in range(r, N, p):
            bits |= 1 << n
    return bits


def make_prime_data(p: int, zeros: Sequence[int], N: int) -> PrimeData:
    zt = tuple(sorted(zeros))
    bits = hit_bits(zt, p, N)
    return PrimeData(
        p=p,
        zeros=zt,
        bits=bits,
        hits=bits.bit_count(),
        center_hit=((p - 1) // 2 in zt),
    )


def load_vector(pdata: Sequence[PrimeData], N: int) -> List[int]:
    K = [0] * N
    for item in pdata:
        for r in item.zeros:
            for n in range(r, N, item.p):
                K[n] += 1
    return K


def fraction_text(x: Fraction) -> str:
    if x.denominator == 1:
        return str(x.numerator)
    return f"{x.numerator}/{x.denominator}"


def decimal_text(x: Fraction | float | int, digits: int = 15) -> str:
    return f"{float(x):.{digits}g}"


def pair_centered_data(a: PrimeData, b: PrimeData, N: int, pair_hit: int) -> Tuple[int, int, int]:
    A_a = a.p * a.hits - N * a.z
    A_b = b.p * b.hits - N * b.z
    A_ab = (
        a.p * b.p * pair_hit
        - a.p * b.z * a.hits
        - b.p * a.z * b.hits
        + N * a.z * b.z
    )
    return A_a, A_b, A_ab


def triple_decomposition(
    a: PrimeData,
    b: PrimeData,
    c: PrimeData,
    N: int,
    pair_hits: Mapping[Tuple[int, int], int],
) -> dict:
    p, q, r = a.p, b.p, c.p
    hp, hq, hr = a.hits, b.hits, c.hits
    hpq = pair_hits[(p, q)]
    hpr = pair_hits[(p, r)]
    hqr = pair_hits[(q, r)]
    hpqr = (a.bits & b.bits & c.bits).bit_count()

    Ap = p * hp - N * a.z
    Aq = q * hq - N * b.z
    Ar = r * hr - N * c.z
    Apq = p * q * hpq - p * b.z * hp - q * a.z * hq + N * a.z * b.z
    Apr = p * r * hpr - p * c.z * hp - r * a.z * hr + N * a.z * c.z
    Aqr = q * r * hqr - q * c.z * hq - r * b.z * hr + N * b.z * c.z
    Apqr = (
        p * q * r * hpqr
        - p * q * c.z * hpq
        - p * r * b.z * hpr
        - q * r * a.z * hqr
        + p * b.z * c.z * hp
        + q * a.z * c.z * hq
        + r * a.z * b.z * hr
        - N * a.z * b.z * c.z
    )

    base_num = N * a.z * b.z * c.z
    single_p_num = Ap * b.z * c.z
    single_q_num = Aq * a.z * c.z
    single_r_num = Ar * a.z * b.z
    pair_pq_num = Apq * c.z
    pair_pr_num = Apr * b.z
    pair_qr_num = Aqr * a.z
    primitive_num = Apqr
    den = p * q * r
    reconstructed_num = (
        base_num
        + single_p_num
        + single_q_num
        + single_r_num
        + pair_pq_num
        + pair_pr_num
        + pair_qr_num
        + primitive_num
    )
    direct_num = den * hpqr
    if reconstructed_num != direct_num:
        raise AssertionError(
            f"triple reconstruction failed for {(p,q,r)}: {reconstructed_num} != {direct_num}"
        )

    C = Fraction(primitive_num, den)
    lower = Fraction(reconstructed_num - primitive_num, den)
    return {
        "p": p,
        "q": q,
        "r": r,
        "zp": a.z,
        "zq": b.z,
        "zr": c.z,
        "Hp": hp,
        "Hq": hq,
        "Hr": hr,
        "Hpq": hpq,
        "Hpr": hpr,
        "Hqr": hqr,
        "Hpqr": hpqr,
        "den": den,
        "base_num": base_num,
        "single_p_num": single_p_num,
        "single_q_num": single_q_num,
        "single_r_num": single_r_num,
        "pair_pq_num": pair_pq_num,
        "pair_pr_num": pair_pr_num,
        "pair_qr_num": pair_qr_num,
        "primitive_num": primitive_num,
        "direct_num": direct_num,
        "reconstructed_num": reconstructed_num,
        "C_frac": C,
        "lower_frac": lower,
        "abs_C_frac": abs(C),
        "total_excess_frac": Fraction(direct_num - base_num, den),
    }


def analyze_instance(X: int, zero_sets: Mapping[int, Sequence[int]], label: str) -> InstanceAnalysis:
    N = X * X
    pdata = [make_prime_data(p, zero_sets[p], N) for p in sorted(zero_sets)]
    active = [d for d in pdata if d.z]
    K = load_vector(pdata, N)
    F1 = sum(K)
    F2 = sum(k * (k - 1) for k in K)
    F3 = sum(k * (k - 1) * (k - 2) for k in K)
    M3 = sum(k**3 for k in K)
    if M3 != F3 + 3 * F2 + F1:
        raise AssertionError("ordinary/factorial third-moment diagonal identity failed")

    pair_hits: Dict[Tuple[int, int], int] = {}
    for a, b in itertools.combinations(active, 2):
        pair_hits[(a.p, b.p)] = (a.bits & b.bits).bit_count()

    triple_rows = [
        triple_decomposition(a, b, c, N, pair_hits)
        for a, b, c in itertools.combinations(active, 3)
    ]

    deg0 = sum((Fraction(row["base_num"], row["den"]) for row in triple_rows), Fraction())
    deg1 = sum(
        (
            Fraction(
                row["single_p_num"] + row["single_q_num"] + row["single_r_num"],
                row["den"],
            )
            for row in triple_rows
        ),
        Fraction(),
    )
    deg2 = sum(
        (
            Fraction(
                row["pair_pq_num"] + row["pair_pr_num"] + row["pair_qr_num"],
                row["den"],
            )
            for row in triple_rows
        ),
        Fraction(),
    )
    deg3 = sum((row["C_frac"] for row in triple_rows), Fraction())
    centered_degree = (6 * deg0, 6 * deg1, 6 * deg2, 6 * deg3)
    if sum(centered_degree, Fraction()) != F3:
        raise AssertionError(
            f"X={X} {label}: centered triple aggregate does not reconstruct F3: "
            f"{sum(centered_degree, Fraction())} != {F3}"
        )

    # Ordered F2 = 2 * sum_{p<q} H_{pq}, split by centered degree 0,1,2.
    f2d0 = Fraction()
    f2d1 = Fraction()
    f2d2 = Fraction()
    for a, b in itertools.combinations(active, 2):
        hpq = pair_hits[(a.p, b.p)]
        Aa, Ab, Aab = pair_centered_data(a, b, N, hpq)
        den = a.p * b.p
        f2d0 += Fraction(N * a.z * b.z, den)
        f2d1 += Fraction(Aa * b.z + Ab * a.z, den)
        f2d2 += Fraction(Aab, den)
    pair_centered_degree = (2 * f2d0, 2 * f2d1, 2 * f2d2)
    if sum(pair_centered_degree, Fraction()) != F2:
        raise AssertionError("centered pair aggregate does not reconstruct F2")

    # F1 split by centered degree 0,1.
    f1d0 = sum((Fraction(N * a.z, a.p) for a in active), Fraction())
    f1d1 = sum((Fraction(a.p * a.hits - N * a.z, a.p) for a in active), Fraction())
    one_centered_degree = (f1d0, f1d1)
    if sum(one_centered_degree, Fraction()) != F1:
        raise AssertionError("centered one-prime aggregate does not reconstruct F1")

    abs_values = sorted((abs(row["C_frac"]) for row in triple_rows), reverse=True)
    abs_mass = sum(abs_values, Fraction())
    ordered_abs_mass = 6 * abs_mass
    signed_primitive = centered_degree[3]
    if ordered_abs_mass < abs(signed_primitive):
        raise AssertionError("ordered primitive absolute mass is below absolute signed mass")
    top1 = abs_values[0] if abs_values else Fraction()
    top5 = sum(abs_values[:5], Fraction())
    top10 = sum(abs_values[:10], Fraction())
    sum_sq = sum((float(v) ** 2 for v in abs_values), 0.0)
    participation = (float(abs_mass) ** 2 / sum_sq) if sum_sq else 0.0
    hist: Dict[int, int] = {}
    for k in K:
        hist[k] = hist.get(k, 0) + 1
    metrics = {
        "F1": F1,
        "F2": F2,
        "F3": F3,
        "M3": M3,
        "lower_diagonal_3F2_plus_F1": 3 * F2 + F1,
        "max_load": max(K) if K else 0,
        "rows_K_ge_3": sum(v for k, v in hist.items() if k >= 3),
        "active_triples": len(triple_rows),
        "ordered_primitive_signed_total": signed_primitive,
        "ordered_primitive_abs_mass": ordered_abs_mass,
        "triple_top1_abs_share": (top1 / abs_mass if abs_mass else Fraction()),
        "triple_top5_abs_share": (top5 / abs_mass if abs_mass else Fraction()),
        "triple_top10_abs_share": (top10 / abs_mass if abs_mass else Fraction()),
        "triple_participation_ratio": participation,
        "load_histogram": hist,
    }

    return InstanceAnalysis(
        X=X,
        N=N,
        label=label,
        pdata=pdata,
        active=active,
        K=K,
        F1=F1,
        F2=F2,
        F3=F3,
        M3=M3,
        pair_hits=pair_hits,
        triple_rows=triple_rows,
        centered_degree=centered_degree,
        pair_centered_degree=pair_centered_degree,
        one_centered_degree=one_centered_degree,
        metrics=metrics,
    )


def reflection_matched_random_set(p: int, actual: Sequence[int], rng: random.Random) -> Tuple[int, ...]:
    center = (p - 1) // 2
    actual_set = set(actual)
    center_hit = center in actual_set
    remaining = len(actual) - int(center_hit)
    if remaining % 2:
        raise AssertionError(f"p={p}: reflection signature has odd nonfixed cardinality")
    pair_count = remaining // 2
    # Match actual endpoint exclusion as well as reflection and fixed-point status.
    representatives = list(range(1, center))
    if pair_count > len(representatives):
        raise AssertionError("random reflection pool too small")
    chosen = rng.sample(representatives, pair_count)
    out = set()
    for r in chosen:
        out.add(r)
        out.add(p - 1 - r)
    if center_hit:
        out.add(center)
    result = tuple(sorted(out))
    if len(result) != len(actual):
        raise AssertionError("random reflection cardinality mismatch")
    if tuple(sorted(p - 1 - r for r in result)) != result:
        raise AssertionError("random reflection symmetry mismatch")
    if (center in result) != center_hit:
        raise AssertionError("random reflection fixed-point mismatch")
    if 0 in result or p - 1 in result:
        raise AssertionError("random reflection endpoint mismatch")
    return result


def random_zero_sets(actual_sets: Mapping[int, Sequence[int]], seed: int) -> Dict[int, Tuple[int, ...]]:
    rng = random.Random(seed)
    return {
        p: reflection_matched_random_set(p, zeros, rng)
        for p, zeros in sorted(actual_sets.items())
    }


def metric_float(value: object) -> float:
    if isinstance(value, Fraction):
        return float(value)
    return float(value)  # type: ignore[arg-type]


def random_ensemble(
    X: int,
    actual_sets: Mapping[int, Sequence[int]],
    actual: InstanceAnalysis,
    reps: int,
) -> Tuple[List[dict], List[dict], Dict[int, Dict[int, Tuple[int, ...]]]]:
    metric_names = [
        "F3",
        "ordered_primitive_signed_total",
        "ordered_primitive_abs_mass",
        "triple_top1_abs_share",
        "triple_top5_abs_share",
        "triple_top10_abs_share",
        "triple_participation_ratio",
        "max_load",
        "rows_K_ge_3",
    ]
    rep_rows: List[dict] = []
    rep_sets: Dict[int, Dict[int, Tuple[int, ...]]] = {}
    for rep in range(reps):
        seed = 731100000 + X * 10000 + rep
        zsets = random_zero_sets(actual_sets, seed)
        if rep < 8:
            rep_sets[rep] = zsets
        trial = analyze_instance(X, zsets, f"random_{rep}")
        row = {"X": X, "rep": rep, "seed": seed}
        for name in metric_names:
            row[name] = metric_float(trial.metrics[name])
        rep_rows.append(row)

    summary_rows: List[dict] = []
    for name in metric_names:
        values = [row[name] for row in rep_rows]
        mean = statistics.fmean(values) if values else 0.0
        sd = statistics.stdev(values) if len(values) >= 2 else 0.0
        observed = metric_float(actual.metrics[name])
        zscore = (observed - mean) / sd if sd else math.nan
        ge = sum(v >= observed for v in values)
        le = sum(v <= observed for v in values)
        summary_rows.append(
            {
                "X": X,
                "metric": name,
                "actual": f"{observed:.17g}",
                "random_reps": reps,
                "random_mean": f"{mean:.17g}",
                "random_sd": f"{sd:.17g}",
                "random_min": f"{min(values):.17g}" if values else "",
                "random_max": f"{max(values):.17g}" if values else "",
                "zscore": f"{zscore:.17g}" if math.isfinite(zscore) else "NA",
                "empirical_ge_fraction": f"{ge / reps:.17g}" if reps else "NA",
                "empirical_le_fraction": f"{le / reps:.17g}" if reps else "NA",
            }
        )
    return rep_rows, summary_rows, rep_sets


def fourier_array_double(p: int, zeros: Sequence[int]) -> List[complex]:
    arr = [0j] * p
    two_pi_over_p = 2.0 * math.pi / p
    for a in range(1, p):
        arr[a] = sum(cmath.exp(-1j * two_pi_over_p * a * z) for z in zeros)
    return arr


def local_frequencies(k: int, p: int, q: int, r: int) -> Tuple[int, int, int]:
    a = k * pow((q * r) % p, -1, p) % p
    b = k * pow((p * r) % q, -1, q) % q
    c = k * pow((p * q) % r, -1, r) % r
    return a, b, c


def dirichlet_kernel_double(N: int, k: int, M: int) -> complex:
    theta = math.pi * k / M
    den = math.sin(theta)
    amp = math.sin(N * theta) / den
    phase = (N - 1) * theta
    return complex(math.cos(phase), math.sin(phase)) * amp


def conjugate_pair_contribution_double(
    k: int,
    p: int,
    q: int,
    r: int,
    N: int,
    Fp: Sequence[complex],
    Fq: Sequence[complex],
    Fr: Sequence[complex],
) -> Tuple[float, int, int, int]:
    M = p * q * r
    a, b, c = local_frequencies(k, p, q, r)
    D = dirichlet_kernel_double(N, k, M)
    term = Fp[a] * Fq[b] * Fr[c] * D / M
    return 2.0 * term.real, a, b, c


def unseen_pair_bound(zprod: int, N: int, M: int, first_unseen: int) -> float:
    if first_unseen > M // 2:
        return 0.0
    theta = math.pi * first_unseen / M
    return min(2.0 * zprod * N / M, 2.0 * zprod / (M * math.sin(theta)))


def high_precision_mode_value(
    k: int,
    p: int,
    q: int,
    r: int,
    N: int,
    zp: Sequence[int],
    zq: Sequence[int],
    zr: Sequence[int],
    dps: int,
):
    import mpmath as mp

    with mp.workdps(dps):
        M = p * q * r
        a, b, c = local_frequencies(k, p, q, r)
        j = mp.j
        twopi = 2 * mp.pi
        fp = mp.fsum([mp.e ** (-j * twopi * a * z / p) for z in zp])
        fq = mp.fsum([mp.e ** (-j * twopi * b * z / q) for z in zq])
        fr = mp.fsum([mp.e ** (-j * twopi * c * z / r) for z in zr])
        theta = mp.pi * k / M
        D = mp.e ** (j * (N - 1) * theta) * mp.sin(N * theta) / mp.sin(theta)
        value = 2 * mp.re(fp * fq * fr * D / M)
        return +value


def top_primitive_frequency_pairs(
    X: int,
    label: str,
    triple_row: Mapping[str, object],
    zero_sets: Mapping[int, Sequence[int]],
    top_modes: int,
    max_scan: int,
) -> Tuple[List[dict], dict]:
    p, q, r = int(triple_row["p"]), int(triple_row["q"]), int(triple_row["r"])
    N = X * X
    M = p * q * r
    zp, zq, zr = zero_sets[p], zero_sets[q], zero_sets[r]
    zprod = len(zp) * len(zq) * len(zr)
    Fp = fourier_array_double(p, zp)
    Fq = fourier_array_double(q, zq)
    Fr = fourier_array_double(r, zr)

    # Retain every scanned primitive pair; high-precision reranking then
    # covers the complete scanned prefix rather than a floating heap.
    scanned_candidates: List[Tuple[float, int, float, int, int, int]] = []
    certified_by_envelope = False
    final_bound = math.inf
    scanned_through = 0
    for k in range(1, min(M // 2, max_scan) + 1):
        scanned_through = k
        if math.gcd(k, M) == 1:
            value, a, b, c = conjugate_pair_contribution_double(k, p, q, r, N, Fp, Fq, Fr)
            item = (abs(value), k, value, a, b, c)
            scanned_candidates.append(item)
        if k % 64 == 0 and len(scanned_candidates) >= top_modes:
            ranked = sorted(scanned_candidates, reverse=True)
            kth = ranked[top_modes - 1][0]
            final_bound = unseen_pair_bound(zprod, N, M, k + 1)
            if final_bound + 1e-12 < kth:
                certified_by_envelope = True
                break
    else:
        final_bound = unseen_pair_bound(zprod, N, M, scanned_through + 1)
        if scanned_through >= M // 2:
            certified_by_envelope = True

    candidates = sorted(scanned_candidates, reverse=True)
    hp_candidates = []
    for _, k, _, a, b, c in candidates:
        v80 = high_precision_mode_value(k, p, q, r, N, zp, zq, zr, 80)
        v140 = high_precision_mode_value(k, p, q, r, N, zp, zq, zr, 140)
        hp_candidates.append((abs(v140), k, v140, abs(v140 - v80), a, b, c))
    hp_candidates.sort(reverse=True)
    selected = hp_candidates[:top_modes]

    import mpmath as mp

    exact_C = triple_row["C_frac"]
    assert isinstance(exact_C, Fraction)
    C_float = float(exact_C)
    cumulative = mp.mpf("0")
    rows: List[dict] = []
    max_hp_delta = max((entry[3] for entry in selected), default=mp.mpf("0"))
    kth_abs = selected[-1][0] if selected else mp.mpf("0")
    inclusion_margin = kth_abs - final_bound
    certified_tolerance = mp.mpf("1e-11") + max_hp_delta
    ranking_certified = bool(certified_by_envelope and inclusion_margin > certified_tolerance)
    for rank, (absval, k, value, hp_delta, a, b, c) in enumerate(selected, 1):
        cumulative += value
        rows.append(
            {
                "X": X,
                "set_label": label,
                "p": p,
                "q": q,
                "r": r,
                "M": M,
                "rank": rank,
                "k": k,
                "conjugate_k": M - k,
                "a_mod_p": a,
                "b_mod_q": b,
                "c_mod_r": c,
                "pair_contribution": mp.nstr(value, 45),
                "abs_pair_contribution": mp.nstr(absval, 45),
                "hp_80_140_delta": mp.nstr(hp_delta, 8),
                "cumulative_signed": mp.nstr(cumulative, 45),
                "cumulative_over_C": (
                    mp.nstr(cumulative / C_float, 25) if C_float != 0.0 else "NA"
                ),
                "unseen_single_mode_bound": f"{final_bound:.17g}",
                "ranking_inclusion_margin": mp.nstr(inclusion_margin, 20),
                "ranking_certified_tolerance": int(ranking_certified),
            }
        )

    top1_abs = selected[0][0] if selected else mp.mpf("0")
    top5_signed = mp.fsum([v for _, _, v, _, _, _, _ in selected[:5]])
    top10_signed = mp.fsum([v for _, _, v, _, _, _, _ in selected[:10]])
    top5_abs = mp.fsum([av for av, *_ in selected[:5]])
    top10_abs = mp.fsum([av for av, *_ in selected[:10]])
    mode_summary = {
        "X": X,
        "set_label": label,
        "p": p,
        "q": q,
        "r": r,
        "C_exact": fraction_text(exact_C),
        "C_decimal": f"{C_float:.17g}",
        "scanned_through_k": scanned_through,
        "M_half": M // 2,
        "unseen_single_mode_bound": f"{final_bound:.17g}",
        "ranking_certified_tolerance": int(ranking_certified),
        "top1_abs_over_abs_C": (
            mp.nstr(top1_abs / abs(C_float), 25) if C_float else "NA"
        ),
        "top5_abs_over_abs_C": (
            mp.nstr(top5_abs / abs(C_float), 25) if C_float else "NA"
        ),
        "top10_abs_over_abs_C": (
            mp.nstr(top10_abs / abs(C_float), 25) if C_float else "NA"
        ),
        "top5_signed_over_C": (
            mp.nstr(top5_signed / C_float, 25) if C_float else "NA"
        ),
        "top10_signed_over_C": (
            mp.nstr(top10_signed / C_float, 25) if C_float else "NA"
        ),
        "residual_after_top10_over_C": (
            mp.nstr((C_float - top10_signed) / C_float, 25) if C_float else "NA"
        ),
        "max_hp_80_140_delta": mp.nstr(max_hp_delta, 8),
    }
    return rows, mode_summary


def full_fourier_sum_double(
    X: int,
    triple_row: Mapping[str, object],
    zero_sets: Mapping[int, Sequence[int]],
) -> Tuple[float, float, int]:
    p, q, r = int(triple_row["p"]), int(triple_row["q"]), int(triple_row["r"])
    N = X * X
    M = p * q * r
    Fp = fourier_array_double(p, zero_sets[p])
    Fq = fourier_array_double(q, zero_sets[q])
    Fr = fourier_array_double(r, zero_sets[r])
    contributions = []
    for k in range(1, M // 2 + 1):
        if math.gcd(k, M) == 1:
            value, _, _, _ = conjugate_pair_contribution_double(k, p, q, r, N, Fp, Fq, Fr)
            contributions.append(value)
    value = math.fsum(contributions)
    abs_mass = math.fsum(abs(v) for v in contributions)
    return value, abs_mass, len(contributions)


def write_csv(path: Path, rows: Sequence[Mapping[str, object]], fields: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(fields), extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def csv_block(rows: Sequence[Mapping[str, object]], fields: Sequence[str]) -> str:
    out = io.StringIO()
    writer = csv.DictWriter(out, fieldnames=list(fields), extrasaction="ignore", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return "```csv\n" + out.getvalue() + "```\n"


def public_triple_row(X: int, row: Mapping[str, object], rank: int) -> dict:
    C = row["C_frac"]
    lower = row["lower_frac"]
    excess = row["total_excess_frac"]
    assert isinstance(C, Fraction) and isinstance(lower, Fraction) and isinstance(excess, Fraction)
    return {
        "X": X,
        "rank_abs_C": rank,
        "p": row["p"],
        "q": row["q"],
        "r": row["r"],
        "zp": row["zp"],
        "zq": row["zq"],
        "zr": row["zr"],
        "Hpqr": row["Hpqr"],
        "den": row["den"],
        "base_num": row["base_num"],
        "single_total_num": int(row["single_p_num"]) + int(row["single_q_num"]) + int(row["single_r_num"]),
        "pair_total_num": int(row["pair_pq_num"]) + int(row["pair_pr_num"]) + int(row["pair_qr_num"]),
        "primitive_num": row["primitive_num"],
        "lower_exact": fraction_text(lower),
        "C_exact": fraction_text(C),
        "C_decimal": decimal_text(C, 17),
        "abs_C_decimal": decimal_text(abs(C), 17),
        "total_excess_exact": fraction_text(excess),
        "scaled_reconstruction_error": int(row["reconstructed_num"]) - int(row["direct_num"]),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--xs", nargs="+", type=int, default=list(DEFAULT_XS))
    parser.add_argument("--output", type=Path, default=Path("q7311-output"))
    parser.add_argument("--random-reps", type=int, default=64)
    parser.add_argument("--mode-random-reps", type=int, default=4)
    parser.add_argument("--mode-triples", type=int, default=3)
    parser.add_argument("--top-triples", type=int, default=25)
    parser.add_argument("--top-modes", type=int, default=12)
    parser.add_argument("--mode-max-scan", type=int, default=2_000_000)
    parser.add_argument("--full-fourier-max-M", type=int, default=350_000)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outdir: Path = args.output
    outdir.mkdir(parents=True, exist_ok=True)

    zero_rows: List[dict] = []
    moment_rows: List[dict] = []
    centered_rows: List[dict] = []
    reconstruction_rows: List[dict] = []
    concentration_rows: List[dict] = []
    all_top_triples: List[dict] = []
    all_triples_public: List[dict] = []
    random_rep_rows: List[dict] = []
    random_summary_rows: List[dict] = []
    mode_rows: List[dict] = []
    mode_summary_rows: List[dict] = []
    full_fourier_rows: List[dict] = []
    random_signature_rows: List[dict] = []
    runtime_rows: List[dict] = []

    actual_analyses: Dict[int, InstanceAnalysis] = {}
    actual_sets_by_X: Dict[int, Dict[int, Tuple[int, ...]]] = {}
    random_sets_for_modes: Dict[int, Dict[int, Dict[int, Tuple[int, ...]]]] = {}

    total_start = time.perf_counter()
    for X in args.xs:
        start = time.perf_counter()
        primes = primes_in_shell(X)
        actual_sets = {p: actual_zero_set(p) for p in primes}
        actual_sets_by_X[X] = actual_sets
        analysis = analyze_instance(X, actual_sets, "actual")
        actual_analyses[X] = analysis

        for item in analysis.pdata:
            zero_rows.append(
                {
                    "X": X,
                    "p": item.p,
                    "Z": item.z,
                    "center_hit": int(item.center_hit),
                    "reflection_ok": int(tuple(sorted(item.p - 1 - z for z in item.zeros)) == item.zeros),
                    "endpoint_free": int(0 not in item.zeros and item.p - 1 not in item.zeros),
                    "zeros": " ".join(map(str, item.zeros)),
                }
            )

        hist = analysis.metrics["load_histogram"]
        moment_rows.append(
            {
                "X": X,
                "N": X * X,
                "prime_count": len(analysis.pdata),
                "active_prime_count": len(analysis.active),
                "sum_Z": sum(d.z for d in analysis.pdata),
                "F1": analysis.F1,
                "F2": analysis.F2,
                "F3": analysis.F3,
                "M3": analysis.M3,
                "all_equal_diagonal_F1": analysis.F1,
                "exactly_two_equal_diagonal_3F2": 3 * analysis.F2,
                "lower_diagonal_total": 3 * analysis.F2 + analysis.F1,
                "M3_minus_F3_minus_3F2_minus_F1": analysis.M3 - analysis.F3 - 3 * analysis.F2 - analysis.F1,
                "max_K": max(analysis.K),
                "rows_K_ge_3": analysis.metrics["rows_K_ge_3"],
                "load_histogram": " ".join(f"{k}:{hist[k]}" for k in sorted(hist)),
            }
        )

        for order, degrees, direct in [
            (1, analysis.one_centered_degree, analysis.F1),
            (2, analysis.pair_centered_degree, analysis.F2),
            (3, analysis.centered_degree, analysis.F3),
        ]:
            row = {"X": X, "factorial_order": order, "direct": direct}
            for degree in range(4):
                value = degrees[degree] if degree < len(degrees) else Fraction()
                row[f"centered_degree_{degree}_exact"] = fraction_text(value)
                row[f"centered_degree_{degree}_decimal"] = decimal_text(value, 17)
            reconstructed = sum(degrees, Fraction())
            row["reconstructed_exact"] = fraction_text(reconstructed)
            row["reconstruction_error_exact"] = fraction_text(reconstructed - direct)
            centered_rows.append(row)

        all_scaled_errors = [int(row["reconstructed_num"]) - int(row["direct_num"]) for row in analysis.triple_rows]
        reconstruction_rows.append(
            {
                "X": X,
                "triple_count": len(analysis.triple_rows),
                "max_abs_scaled_triple_error": max(map(abs, all_scaled_errors), default=0),
                "F3_centered_error_exact": fraction_text(sum(analysis.centered_degree, Fraction()) - analysis.F3),
                "F2_centered_error_exact": fraction_text(sum(analysis.pair_centered_degree, Fraction()) - analysis.F2),
                "F1_centered_error_exact": fraction_text(sum(analysis.one_centered_degree, Fraction()) - analysis.F1),
                "ordinary_diagonal_error": analysis.M3 - analysis.F3 - 3 * analysis.F2 - analysis.F1,
            }
        )

        metrics = analysis.metrics
        concentration_rows.append(
            {
                "X": X,
                "active_triples": metrics["active_triples"],
                "ordered_primitive_signed_total_exact": fraction_text(metrics["ordered_primitive_signed_total"]),
                "ordered_primitive_signed_total_decimal": decimal_text(metrics["ordered_primitive_signed_total"], 17),
                "ordered_primitive_abs_mass_exact": fraction_text(metrics["ordered_primitive_abs_mass"]),
                "ordered_primitive_abs_mass_decimal": decimal_text(metrics["ordered_primitive_abs_mass"], 17),
                "top1_abs_share": decimal_text(metrics["triple_top1_abs_share"], 17),
                "top5_abs_share": decimal_text(metrics["triple_top5_abs_share"], 17),
                "top10_abs_share": decimal_text(metrics["triple_top10_abs_share"], 17),
                "participation_ratio": f"{metrics['triple_participation_ratio']:.17g}",
            }
        )

        ranked = sorted(
            analysis.triple_rows,
            key=lambda row: (row["abs_C_frac"], row["p"], row["q"], row["r"]),
            reverse=True,
        )
        for rank, row in enumerate(ranked, 1):
            pub = public_triple_row(X, row, rank)
            all_triples_public.append(pub)
            if rank <= args.top_triples:
                all_top_triples.append(pub)

        reps, summaries, retained_sets = random_ensemble(
            X, actual_sets, analysis, args.random_reps
        )
        random_rep_rows.extend(reps)
        random_summary_rows.extend(summaries)
        random_sets_for_modes[X] = retained_sets

        if retained_sets:
            first = retained_sets[min(retained_sets)]
            for p in primes:
                actual_z = actual_sets[p]
                random_z = first[p]
                random_signature_rows.append(
                    {
                        "X": X,
                        "rep": min(retained_sets),
                        "p": p,
                        "actual_Z": len(actual_z),
                        "random_Z": len(random_z),
                        "actual_center_hit": int((p - 1) // 2 in actual_z),
                        "random_center_hit": int((p - 1) // 2 in random_z),
                        "random_reflection_ok": int(tuple(sorted(p - 1 - z for z in random_z)) == random_z),
                        "random_endpoint_free": int(0 not in random_z and p - 1 not in random_z),
                        "random_zeros": " ".join(map(str, random_z)),
                    }
                )

        runtime_rows.append(
            {
                "X": X,
                "status": "completed",
                "seconds_through_exact_and_random": f"{time.perf_counter() - start:.6f}",
            }
        )

    # Frequency rankings for the largest |C| actual triples and the largest |C|
    # triple in a few random reflection-matched realizations.
    for X in args.xs:
        analysis = actual_analyses[X]
        actual_sets = actual_sets_by_X[X]
        ranked_actual = sorted(
            analysis.triple_rows,
            key=lambda row: row["abs_C_frac"],
            reverse=True,
        )
        for row in ranked_actual[: args.mode_triples]:
            rows, summary = top_primitive_frequency_pairs(
                X,
                "actual",
                row,
                actual_sets,
                args.top_modes,
                args.mode_max_scan,
            )
            mode_rows.extend(rows)
            mode_summary_rows.append(summary)
            M = int(row["den"])
            if M <= args.full_fourier_max_M:
                fsum, fabs, count = full_fourier_sum_double(X, row, actual_sets)
                C = row["C_frac"]
                assert isinstance(C, Fraction)
                err = abs(fsum - float(C))
                tolerance = 5e-10 * (1.0 + fabs)
                if err > tolerance:
                    raise AssertionError(
                        f"full Fourier sum outside tolerance for X={X}, triple={(row['p'],row['q'],row['r'])}: "
                        f"error={err}, tolerance={tolerance}"
                    )
                full_fourier_rows.append(
                    {
                        "X": X,
                        "set_label": "actual",
                        "p": row["p"],
                        "q": row["q"],
                        "r": row["r"],
                        "M": M,
                        "conjugate_pair_count": count,
                        "C_exact": fraction_text(C),
                        "C_decimal": f"{float(C):.17g}",
                        "full_fourier_sum_double": f"{fsum:.17g}",
                        "absolute_error": f"{err:.17g}",
                        "asserted_tolerance": f"{tolerance:.17g}",
                        "within_tolerance": int(err <= tolerance),
                        "frequency_abs_mass_double": f"{fabs:.17g}",
                    }
                )

        retained = random_sets_for_modes[X]
        for rep in sorted(retained)[: args.mode_random_reps]:
            zsets = retained[rep]
            trial = analyze_instance(X, zsets, f"random_{rep}")
            if not trial.triple_rows:
                continue
            row = max(trial.triple_rows, key=lambda t: t["abs_C_frac"])
            rows, summary = top_primitive_frequency_pairs(
                X,
                f"random_{rep}",
                row,
                zsets,
                args.top_modes,
                args.mode_max_scan,
            )
            mode_rows.extend(rows)
            mode_summary_rows.append(summary)

    runtime_rows.append(
        {
            "X": "ALL",
            "status": "completed",
            "seconds_through_exact_and_random": f"{time.perf_counter() - total_start:.6f}",
        }
    )

    zero_fields = ["X", "p", "Z", "center_hit", "reflection_ok", "endpoint_free", "zeros"]
    moment_fields = [
        "X", "N", "prime_count", "active_prime_count", "sum_Z", "F1", "F2", "F3", "M3",
        "all_equal_diagonal_F1", "exactly_two_equal_diagonal_3F2", "lower_diagonal_total",
        "M3_minus_F3_minus_3F2_minus_F1", "max_K", "rows_K_ge_3", "load_histogram",
    ]
    centered_fields = [
        "X", "factorial_order", "direct",
        "centered_degree_0_exact", "centered_degree_1_exact", "centered_degree_2_exact", "centered_degree_3_exact",
        "centered_degree_0_decimal", "centered_degree_1_decimal", "centered_degree_2_decimal", "centered_degree_3_decimal",
        "reconstructed_exact", "reconstruction_error_exact",
    ]
    reconstruction_fields = [
        "X", "triple_count", "max_abs_scaled_triple_error", "F3_centered_error_exact",
        "F2_centered_error_exact", "F1_centered_error_exact", "ordinary_diagonal_error",
    ]
    concentration_fields = [
        "X", "active_triples", "ordered_primitive_signed_total_exact", "ordered_primitive_signed_total_decimal",
        "ordered_primitive_abs_mass_exact", "ordered_primitive_abs_mass_decimal", "top1_abs_share", "top5_abs_share",
        "top10_abs_share", "participation_ratio",
    ]
    triple_fields = [
        "X", "rank_abs_C", "p", "q", "r", "zp", "zq", "zr", "Hpqr", "den", "base_num",
        "single_total_num", "pair_total_num", "primitive_num", "lower_exact", "C_exact", "C_decimal",
        "abs_C_decimal", "total_excess_exact", "scaled_reconstruction_error",
    ]
    random_rep_fields = [
        "X", "rep", "seed", "F3", "ordered_primitive_signed_total", "ordered_primitive_abs_mass",
        "triple_top1_abs_share", "triple_top5_abs_share", "triple_top10_abs_share",
        "triple_participation_ratio", "max_load", "rows_K_ge_3",
    ]
    random_summary_fields = [
        "X", "metric", "actual", "random_reps", "random_mean", "random_sd", "random_min", "random_max",
        "zscore", "empirical_ge_fraction", "empirical_le_fraction",
    ]
    mode_fields = [
        "X", "set_label", "p", "q", "r", "M", "rank", "k", "conjugate_k", "a_mod_p", "b_mod_q",
        "c_mod_r", "pair_contribution", "abs_pair_contribution", "hp_80_140_delta", "cumulative_signed",
        "cumulative_over_C", "unseen_single_mode_bound", "ranking_inclusion_margin",
        "ranking_certified_tolerance",
    ]
    mode_summary_fields = [
        "X", "set_label", "p", "q", "r", "C_exact", "C_decimal", "scanned_through_k", "M_half",
        "unseen_single_mode_bound", "ranking_certified_tolerance", "top1_abs_over_abs_C",
        "top5_abs_over_abs_C", "top10_abs_over_abs_C", "top5_signed_over_C", "top10_signed_over_C",
        "residual_after_top10_over_C", "max_hp_80_140_delta",
    ]
    full_fourier_fields = [
        "X", "set_label", "p", "q", "r", "M", "conjugate_pair_count", "C_exact", "C_decimal",
        "full_fourier_sum_double", "absolute_error", "asserted_tolerance", "within_tolerance",
        "frequency_abs_mass_double",
    ]
    random_signature_fields = [
        "X", "rep", "p", "actual_Z", "random_Z", "actual_center_hit", "random_center_hit",
        "random_reflection_ok", "random_endpoint_free", "random_zeros",
    ]
    runtime_fields = ["X", "status", "seconds_through_exact_and_random"]

    write_csv(outdir / "zero_sets.csv", zero_rows, zero_fields)
    write_csv(outdir / "moment_summary.csv", moment_rows, moment_fields)
    write_csv(outdir / "centered_decomposition.csv", centered_rows, centered_fields)
    write_csv(outdir / "reconstruction_audit.csv", reconstruction_rows, reconstruction_fields)
    write_csv(outdir / "triple_concentration.csv", concentration_rows, concentration_fields)
    write_csv(outdir / "all_prime_triples.csv", all_triples_public, triple_fields)
    write_csv(outdir / "top_prime_triples.csv", all_top_triples, triple_fields)
    write_csv(outdir / "random_replicates.csv", random_rep_rows, random_rep_fields)
    write_csv(outdir / "random_comparison.csv", random_summary_rows, random_summary_fields)
    write_csv(outdir / "primitive_frequency_modes.csv", mode_rows, mode_fields)
    write_csv(outdir / "mode_capture_summary.csv", mode_summary_rows, mode_summary_fields)
    write_csv(outdir / "full_fourier_checks.csv", full_fourier_rows, full_fourier_fields)
    write_csv(outdir / "random_signature_audit.csv", random_signature_rows, random_signature_fields)
    write_csv(outdir / "runtime.csv", runtime_rows, runtime_fields)

    report: List[str] = [ANSWER_LINE, "", "# Q7311 exact centered CRT Fourier verifier", ""]
    report.extend([
        "## run_config.csv",
        csv_block(
            [{
                "X_values": " ".join(map(str, args.xs)),
                "random_reps": args.random_reps,
                "mode_random_reps": args.mode_random_reps,
                "mode_triples_per_X": args.mode_triples,
                "top_prime_triples_per_X": args.top_triples,
                "top_frequency_pairs_per_triple": args.top_modes,
                "mode_max_scan": args.mode_max_scan,
                "full_fourier_max_M": args.full_fourier_max_M,
                "exact_arithmetic": "integers_and_Fraction",
                "frequency_precision": "80_and_140_decimal_digits",
            }],
            [
                "X_values", "random_reps", "mode_random_reps", "mode_triples_per_X",
                "top_prime_triples_per_X", "top_frequency_pairs_per_triple", "mode_max_scan",
                "full_fourier_max_M", "exact_arithmetic", "frequency_precision",
            ],
        ),
        "## moment_summary.csv",
        csv_block(moment_rows, moment_fields),
        "## centered_decomposition.csv",
        csv_block(centered_rows, centered_fields),
        "## reconstruction_audit.csv",
        csv_block(reconstruction_rows, reconstruction_fields),
        "## triple_concentration.csv",
        csv_block(concentration_rows, concentration_fields),
        "## random_comparison.csv",
        csv_block(random_summary_rows, random_summary_fields),
        "## full_fourier_checks.csv",
        csv_block(full_fourier_rows, full_fourier_fields),
        "## mode_capture_summary.csv",
        csv_block(mode_summary_rows, mode_summary_fields),
        "## top_prime_triples.csv",
        csv_block(all_top_triples, triple_fields),
        "## primitive_frequency_modes.csv",
        csv_block(mode_rows, mode_fields),
        "## zero_sets.csv",
        csv_block(zero_rows, zero_fields),
        "## random_signature_audit.csv",
        csv_block(random_signature_rows, random_signature_fields),
        "## runtime.csv",
        csv_block(runtime_rows, runtime_fields),
        "## reproducible_code.py",
        "```python\n" + Path(__file__).read_text(encoding="utf-8") + "\n```\n",
        "## generated_raw_files.csv",
        csv_block(
            [{"path": str(path), "bytes": path.stat().st_size} for path in sorted(outdir.glob("*.csv"))],
            ["path", "bytes"],
        ),
    ])
    Path("Q7311_RESULTS.md").write_text("\n".join(report), encoding="utf-8")


if __name__ == "__main__":
    main()
