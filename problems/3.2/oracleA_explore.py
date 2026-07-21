#!/usr/bin/env python3
"""Reproducible Fourier diagnostics for Oracle A in Problem 3.2.

The experiment has two deliberately separate parts.

1.  For doublet primes p in (N/2,N], use the *complete* reflected pair

        m_1=p+r_p,  m_2=2p-1-r_p

    and sample S(theta)=sum_p(e(theta*m_1)+e(theta*m_2)) on a fine FFT
    grid.  This is the sum in Remark ``rem:palindromic-fourier``.

2.  For the true centered covariance, use the fixed shell I_N=(N,2N]
    and the actual shell columns

        Omega_p={m in (N,2N] : m mod p is in Z_p}.

    All lifts r+kp in the shell are enumerated.  Since p>N/2, a residue has
    one or two such lifts, but the relevant k need not equal one.

    The direct covariance is independently reconstructed from its finite
    Fourier bilinear form on Z/NZ.

The default invocation regenerates ``oracleA_exploration.md``.  NumPy is
used only for FFTs and vectorized evaluation; all zero-set input is exact.
No numerical observation is presented as a proof of AMTD.
"""

from __future__ import annotations

import argparse
import hashlib
import math
import os
import struct
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import DefaultDict, Dict, List, Mapping, Sequence, Tuple

try:
    import numpy as np
except ImportError as error:  # pragma: no cover - exercised only off the run host
    raise SystemExit("oracleA_explore.py requires NumPy for its fine-grid FFTs") from error


DEFAULT_NS = (8192, 16384, 32768, 65536)
HERE = Path(__file__).resolve().parent
DEFAULT_DATA = HERE / "data_zp_pairs.bin"
DEFAULT_OUTPUT = HERE / "oracleA_exploration.md"


@dataclass(frozen=True)
class Doublet:
    p: int
    r: int
    h: int
    m1: int
    m2: int


@dataclass(frozen=True)
class Peak:
    k: int
    theta: float
    magnitude: float
    root_t: float
    arc: str
    rational: str
    scaled_offset: float


@dataclass(frozen=True)
class CorrelationRow:
    a: int
    q: int
    modulation_corr: float
    gap_corr: float
    relevant_amplitude: float
    phase_variance: float
    cosine_variance: float


@dataclass(frozen=True)
class FourierResult:
    n: int
    doublets: Tuple[Doublet, ...]
    grid_size: int
    t: int
    collision: int
    l2_exact: int
    l2_numeric: float
    peaks: Tuple[Peak, ...]
    minor_peak: Peak
    correlations: Tuple[CorrelationRow, ...]
    uniform_gap_collision: float
    iid_collision: float


@dataclass(frozen=True)
class CovarianceRow:
    weight_name: str
    columns: int
    hits: int
    diagonal: float
    e_direct: float
    e_fourier: float
    variance: float
    random_sd: float


@dataclass(frozen=True)
class ShellResult:
    n: int
    rows: Tuple[CovarianceRow, ...]


def load_zero_sets(path: Path, maximum_prime: int) -> Tuple[Dict[int, Tuple[int, ...]], str, int]:
    """Read the documented headerless little-endian uint32 pair file."""

    raw = path.read_bytes()
    if len(raw) % 8:
        raise ValueError(f"{path} has {len(raw)} bytes, not a multiple of 8")
    digest = hashlib.sha256(raw).hexdigest()
    grouped: DefaultDict[int, List[int]] = defaultdict(list)
    previous = (-1, -1)
    records = 0
    for p, r in struct.iter_unpack("<II", raw):
        if (p, r) <= previous:
            raise ValueError(f"pair file is not strictly sorted at {(p, r)}")
        if r >= p:
            raise ValueError(f"invalid zero record {(p, r)}")
        previous = (p, r)
        records += 1
        if p <= maximum_prime:
            grouped[p].append(r)
    return {p: tuple(rs) for p, rs in grouped.items()}, digest, records


def selected_zero_sets(
    zeros: Mapping[int, Tuple[int, ...]], n: int
) -> List[Tuple[int, Tuple[int, ...]]]:
    return [(p, rs) for p, rs in zeros.items() if n // 2 < p <= n]


def validate_zero_sets(zeros: Mapping[int, Tuple[int, ...]]) -> int:
    """Check the structural invariants relevant to this experiment."""

    for p, rs in zeros.items():
        if tuple(sorted(rs)) != rs or len(set(rs)) != len(rs):
            raise AssertionError(f"unsorted or repeated residues at p={p}: {rs}")
        residue_set = set(rs)
        for r in rs:
            if p - 1 - r not in residue_set:
                raise AssertionError(f"reflection failure at p={p}, r={r}")
            if (r + 1) % p in residue_set:
                raise AssertionError(f"consecutive zeros at p={p}, r={r}")
    return len(zeros)


def make_doublets(
    zeros: Mapping[int, Tuple[int, ...]], n: int
) -> Tuple[Doublet, ...]:
    answer: List[Doublet] = []
    for p, rs in selected_zero_sets(zeros, n):
        if len(rs) != 2:
            continue
        r, reflected = rs
        if r + reflected != p - 1 or not (0 < r < reflected < p):
            raise AssertionError(f"purported doublet is not reflected: p={p}, Z={rs}")
        h = p - 1 - 2 * r
        m1 = p + r
        m2 = 2 * p - 1 - r
        # These are exact integer checks of the phase factorization used below.
        assert m2 == p + reflected
        assert m1 + m2 == 3 * p - 1
        assert m2 - m1 == h
        answer.append(Doublet(p, r, h, m1, m2))
    return tuple(answer)


def next_power_of_two(value: int) -> int:
    return 1 << (value - 1).bit_length()


def nearest_small_rational(
    k: int, grid: int, n: int, q_max: int
) -> Tuple[str, str, float]:
    """Classify k/grid using |theta-a/q| <= Q/(qN), q<=Q.

    The returned offset is N times the signed circular displacement from the
    nearest (in ordinary distance on [0,1/2]) small-denominator rational.
    Peaks are reported only on this half-circle, so no wrap ambiguity occurs.
    """

    best = None
    major = False
    for q in range(1, q_max + 1):
        a = (2 * k * q + grid) // (2 * grid)
        difference_numerator = k * q - a * grid
        distance = abs(difference_numerator) / (grid * q)
        if best is None or distance < best[0] or (
            distance == best[0] and q < best[2]
        ):
            best = (distance, a, q, difference_numerator / (grid * q))
        if abs(difference_numerator) * n <= q_max * grid:
            major = True
    assert best is not None
    _, a, q, signed = best
    divisor = math.gcd(a, q)
    a //= divisor
    q //= divisor
    return ("major" if major else "minor", f"{a}/{q}", n * signed)


def peak_from_index(
    k: int,
    magnitudes: "np.ndarray",
    grid: int,
    n: int,
    t: int,
    q_major: int,
) -> Peak:
    arc, rational, offset = nearest_small_rational(k, grid, n, q_major)
    magnitude = float(magnitudes[k])
    return Peak(
        k=k,
        theta=k / grid,
        magnitude=magnitude,
        root_t=magnitude / math.sqrt(t),
        arc=arc,
        rational=rational,
        scaled_offset=offset,
    )


def find_peaks(
    magnitudes: "np.ndarray",
    grid: int,
    n: int,
    t: int,
    q_major: int,
    peak_count: int,
) -> Tuple[Tuple[Peak, ...], Peak]:
    """Find sampled local maxima on 0<=theta<=1/2 and the largest minor one."""

    half = magnitudes[: grid // 2 + 1]
    interior = np.nonzero(
        (half[1:-1] >= half[:-2]) & (half[1:-1] >= half[2:])
    )[0] + 1
    candidates = np.concatenate((np.array([0, grid // 2]), interior))
    candidates = np.unique(candidates)

    # Vectorized exact major-arc test on the local maxima.  For
    # theta=k/G and nearest a/q, q*|theta-a/q|=|kq-aG|/G.
    major = np.zeros(candidates.shape, dtype=bool)
    for q in range(1, q_major + 1):
        a = np.floor((candidates * q + grid / 2) / grid).astype(np.int64)
        difference = np.abs(candidates * q - a * grid)
        major |= difference * n <= q_major * grid

    order = np.argsort(half[candidates])[::-1]
    chosen = candidates[order[:peak_count]]
    peaks = tuple(
        peak_from_index(int(k), magnitudes, grid, n, t, q_major)
        for k in chosen
    )

    minor_candidates = candidates[~major]
    if len(minor_candidates) == 0:
        raise AssertionError("major arcs covered every sampled local maximum")
    minor_k = int(minor_candidates[np.argmax(half[minor_candidates])])
    minor_peak = peak_from_index(minor_k, magnitudes, grid, n, t, q_major)
    assert minor_peak.arc == "minor"
    return peaks, minor_peak


def correlation_rows(
    doublets: Sequence[Doublet], q_max: int
) -> Tuple[CorrelationRow, ...]:
    """Connected correlations at reduced theta=a/q, 0<a<q<=q_max.

    The relevant connected statistic is Corr(e(3p theta/2),
    cos(pi theta h_p)); the raw-gap diagnostic replaces the cosine by h_p/p.
    ``relevant_amplitude`` is |S(theta)|/sqrt(T).
    """

    primes = np.array([item.p for item in doublets], dtype=np.float64)
    gaps = np.array([item.h for item in doublets], dtype=np.float64)
    normalized_gaps = gaps / primes
    centered_gaps = normalized_gaps - normalized_gaps.mean()
    gap_variance = float(np.mean(centered_gaps * centered_gaps))
    t = 2 * len(doublets)
    rows: List[CorrelationRow] = []
    for q in range(2, q_max + 1):
        for a in range(1, q):
            if math.gcd(a, q) != 1:
                continue
            theta = a / q
            phase = np.exp(1j * math.pi * 3.0 * primes * theta)
            cosine = np.cos(math.pi * gaps * theta)
            centered_phase = phase - phase.mean()
            centered_cosine = cosine - cosine.mean()
            phase_variance = float(np.mean(np.abs(centered_phase) ** 2))
            cosine_variance = float(np.mean(centered_cosine * centered_cosine))
            if phase_variance > 1e-24 and cosine_variance > 1e-24:
                modulation_corr = abs(
                    np.mean(centered_phase * centered_cosine)
                ) / math.sqrt(phase_variance * cosine_variance)
            else:
                modulation_corr = float("nan")
            if phase_variance > 1e-24 and gap_variance > 1e-24:
                gap_corr = abs(np.mean(centered_phase * centered_gaps)) / math.sqrt(
                    phase_variance * gap_variance
                )
            else:
                gap_corr = float("nan")
            relevant_amplitude = 2.0 * abs(np.sum(phase * cosine)) / math.sqrt(t)
            rows.append(
                CorrelationRow(
                    a,
                    q,
                    float(modulation_corr),
                    float(gap_corr),
                    float(relevant_amplitude),
                    phase_variance,
                    cosine_variance,
                )
            )
    return tuple(rows)


def count_shift_pairs(left: int, right: int, shift: int) -> int:
    """Count 1<=r<=left, 1<=s<=right with s=r+shift."""

    lower = max(1, 1 - shift)
    upper = min(left, right - shift)
    return max(0, upper - lower + 1)


def count_sum_pairs(left: int, right: int, target: int) -> int:
    """Count 1<=r<=left, 1<=s<=right with r+s=target."""

    lower = max(1, target - right)
    upper = min(left, target - 1)
    return max(0, upper - lower + 1)


def check_counting_formulas() -> None:
    """Exhaustively check the O(1) interval formulas on small boxes."""

    for left in range(1, 9):
        for right in range(1, 9):
            for shift in range(-10, 11):
                brute = sum(
                    s == r + shift
                    for r in range(1, left + 1)
                    for s in range(1, right + 1)
                )
                assert count_shift_pairs(left, right, shift) == brute
            for target in range(-2, 20):
                brute = sum(
                    r + s == target
                    for r in range(1, left + 1)
                    for s in range(1, right + 1)
                )
                assert count_sum_pairs(left, right, target) == brute


def residue_lifts_in_shell(n: int, p: int, residue: int) -> Tuple[int, ...]:
    """Return every residue+kp in the integer shell (n,2n]."""

    first_k = (n - residue) // p + 1
    last_k = (2 * n - residue) // p
    return tuple(residue + k * p for k in range(first_k, last_k + 1))


def check_shell_lift_enumerator() -> None:
    """Compare the quotient-bound enumerator with brute force on small boxes."""

    for n in range(8, 33):
        for p in range(n // 2 + 1, n + 1):
            for residue in range(p):
                expected = tuple(
                    m for m in range(n + 1, 2 * n + 1) if m % p == residue
                )
                assert residue_lifts_in_shell(n, p, residue) == expected


def uniform_gap_collision_prediction(doublets: Sequence[Doublet]) -> float:
    """Independent-uniform first-zero prediction for ordered collisions.

    For each observed pair of doublet prime labels p<q, independently replace
    r_p and r_q by uniforms on 1..(p-3)/2 and 1..(q-3)/2.  The four equality
    types m_i(p)=m_j(q) are counted exactly.  Multiplication by two converts
    the unordered-prime sum to the ordered collision convention in F.
    """

    expected_unordered = 0.0
    primes = [item.p for item in doublets]
    for left_index, p in enumerate(primes):
        hp = (p - 3) // 2
        for q in primes[left_index + 1 :]:
            hq = (q - 3) // 2
            count_a = count_shift_pairs(hp, hq, p - q)
            count_b = count_sum_pairs(hp, hq, 2 * q - 1 - p)
            count_c = count_sum_pairs(hp, hq, 2 * p - 1 - q)
            count_d = count_shift_pairs(hp, hq, 2 * (q - p))
            expected_unordered += (count_a + count_b + count_c + count_d) / (
                hp * hq
            )
    return 2.0 * expected_unordered


def analyze_fourier(
    zeros: Mapping[int, Tuple[int, ...]],
    n: int,
    grid_factor: int,
    q_major: int,
    q_corr: int,
    peak_count: int,
) -> FourierResult:
    doublets = make_doublets(zeros, n)
    if not doublets:
        raise ValueError(f"no doublets in ({n // 2},{n}]")
    hits = [position for item in doublets for position in (item.m1, item.m2)]
    multiplicities = Counter(hits)
    t = len(hits)
    collision = sum(value * (value - 1) for value in multiplicities.values())
    l2_exact = sum(value * value for value in multiplicities.values())
    assert l2_exact == t + collision

    grid = next_power_of_two(grid_factor * n)
    if grid <= max(hits) - min(hits):
        raise ValueError("FFT grid is too short for alias-free discrete Parseval")
    incidence = np.zeros(grid, dtype=np.float64)
    for position, multiplicity in multiplicities.items():
        incidence[position] = multiplicity
    transform = np.fft.fft(incidence)
    magnitudes = np.abs(transform)
    l2_numeric = float(np.vdot(transform, transform).real / grid)
    peaks, minor_peak = find_peaks(
        magnitudes, grid, n, t, q_major, peak_count
    )
    correlations = correlation_rows(doublets, q_corr)
    uniform_prediction = uniform_gap_collision_prediction(doublets)
    k = len(doublets)
    iid_collision = 4.0 * k * (k - 1) / n
    return FourierResult(
        n,
        doublets,
        grid,
        t,
        collision,
        l2_exact,
        l2_numeric,
        peaks,
        minor_peak,
        correlations,
        uniform_prediction,
        iid_collision,
    )


def shell_columns(
    zeros: Mapping[int, Tuple[int, ...]], n: int
) -> List[Tuple[int, Tuple[int, ...]]]:
    columns: List[Tuple[int, Tuple[int, ...]]] = []
    for p, rs in selected_zero_sets(zeros, n):
        hits = tuple(
            sorted(
                m
                for residue in rs
                for m in residue_lifts_in_shell(n, p, residue)
            )
        )
        if hits:
            assert len(hits) == len(set(hits))
            assert all(n < m <= 2 * n and m % p in rs for m in hits)
            assert len(hits) <= 2 * len(rs), (n, p, rs, hits)
            columns.append((p, hits))
    return columns


def covariance_for_weights(
    n: int,
    columns: Sequence[Tuple[int, Tuple[int, ...]]],
    weight_name: str,
) -> CovarianceRow:
    if weight_name == "1":
        weight = lambda p: 1.0
    elif weight_name == "log p":
        weight = lambda p: math.log(p)
    else:  # pragma: no cover - internal programming error
        raise ValueError(weight_name)

    load = np.zeros(n, dtype=np.float64)
    sum_weighted_hits = 0.0
    sum_weight2_hits = 0.0
    sum_weight2_hits2 = 0.0
    diagonal = 0.0
    fourth_term = 0.0
    hit_count = 0
    for p, hits in columns:
        w = weight(p)
        a_p = len(hits)
        hit_count += a_p
        indices = np.array([m - n - 1 for m in hits], dtype=np.int64)
        load[indices] += w
        sum_weighted_hits += w * a_p
        sum_weight2_hits += w * w * a_p
        sum_weight2_hits2 += w * w * a_p * a_p
        d_p = a_p * (1.0 - a_p / n)
        diagonal += w * w * d_p
        fourth_term += w**4 * d_p * d_p

    raw_cross = float(np.dot(load, load)) - sum_weight2_hits
    centered_main = (sum_weighted_hits**2 - sum_weight2_hits2) / n
    e_direct = raw_cross - centered_main
    mean = sum_weighted_hits / n
    variance = float(np.dot(load - mean, load - mean))
    if not math.isclose(variance, diagonal + e_direct, rel_tol=2e-12, abs_tol=2e-8):
        raise AssertionError((n, weight_name, variance, diagonal + e_direct))

    transform = np.fft.fft(load)
    aggregate_nonzero = float(np.sum(np.abs(transform[1:]) ** 2))
    column_nonzero = n * sum_weight2_hits - sum_weight2_hits2
    e_fourier = (aggregate_nonzero - column_nonzero) / n
    if not math.isclose(e_direct, e_fourier, rel_tol=2e-9, abs_tol=2e-7):
        raise AssertionError((n, weight_name, e_direct, e_fourier))

    radicand = 2.0 * max(0.0, diagonal * diagonal - fourth_term) / (n - 1)
    random_sd = math.sqrt(radicand)
    return CovarianceRow(
        weight_name,
        len(columns),
        hit_count,
        diagonal,
        e_direct,
        e_fourier,
        variance,
        random_sd,
    )


def analyze_shell(
    zeros: Mapping[int, Tuple[int, ...]], n: int
) -> ShellResult:
    columns = shell_columns(zeros, n)
    if not columns:
        raise ValueError(f"no active shell columns in ({n // 2},{n}]")
    return ShellResult(
        n,
        tuple(covariance_for_weights(n, columns, name) for name in ("1", "log p")),
    )


def finite(value: float) -> bool:
    return not (math.isnan(value) or math.isinf(value))


def fraction_label(row: CorrelationRow) -> str:
    return f"{row.a}/{row.q}"


def fmt(value: float, digits: int = 6) -> str:
    if math.isnan(value):
        return "degenerate"
    return f"{value:.{digits}f}"


def render_report(
    results: Sequence[FourierResult],
    shells: Sequence[ShellResult],
    data_path: Path,
    digest: str,
    record_count: int,
    validated_sets: int,
    grid_factor: int,
    q_major: int,
    q_corr: int,
    peak_count: int,
    command: str,
) -> str:
    lines: List[str] = []
    try:
        displayed_data_path = data_path.resolve().relative_to(Path.cwd().resolve()).as_posix()
    except ValueError:
        displayed_data_path = data_path.as_posix()
    lines.extend(
        [
            "# Oracle A: top-block Fourier exploration",
            "",
            "This is a reproducible diagnostic for the analytic oracle. It is evidence, not a proof of `eq:short-arc` or AMTD.",
            "",
            "## Definitions and reproducibility",
            "",
            f"The data source is `{displayed_data_path}` ({record_count:,} records, SHA-256 `{digest}`). Its {validated_sets:,} stored nonempty zero sets passed sortedness, reflection, and no-consecutive-zero checks. The analyzed prime block is always `N/2 < p <= N`; the default scales are `{', '.join(str(item.n) for item in results)}`.",
            "",
            "For a doublet `Z_p={r_p,p-1-r_p}`, with `r_p<(p-1)/2`, set",
            "",
            "```text",
            "h_p  = p-1-2r_p",
            "m_1  = p+r_p              = (3p-1-h_p)/2",
            "m_2  = 2p-1-r_p           = (3p-1+h_p)/2",
            "S(theta) = sum_p [e(theta m_1)+e(theta m_2)]",
            "         = 2 sum_p e(theta(3p-1)/2) cos(pi theta h_p).",
            "```",
            "",
            f"The FFT grid has size the next power of two at least `{grid_factor}N`. Because it is longer than the diameter of the integer hit set, sampled Parseval is alias-free and equals the exact integral up to floating-point roundoff. A sampled point is called a major-arc point when `|theta-a/q| <= {q_major}/(qN)` for some reduced `a/q` with `q <= {q_major}`; this cutoff is only a declared diagnostic convention.",
            "",
            "The true covariance experiment is separate: on `I_N=(N,2N]`, it uses `Omega_p={m in I_N: m mod p in Z_p}` and both unit and `log p` weights. Every lift `r+kp` in the shell is included; because `p>N/2`, each residue has one or two lifts. Its exact finite Fourier comparison is",
            "",
            "```text",
            "E^o = sum_{p!=q} w_p w_q (|Omega_p cap Omega_q|-A_p A_q/N)",
            "    = (1/N) sum_{k=1}^{N-1} (|sum_p w_p F_p(k)|^2",
            "                                  -sum_p w_p^2 |F_p(k)|^2),",
            "F_p(k)=sum_{m in Omega_p} e(-k(m-N-1)/N).",
            "```",
            "",
            "Run:",
            "",
            "```sh",
            command,
            "```",
            "",
            f"Runtime dependencies: Python `{sys.version.split()[0]}`, NumPy `{np.__version__}`. The report used `{os.cpu_count()}` visible CPUs, although NumPy decides its own FFT threading.",
            "",
            "## Doublet Fourier norm and minor-arc maxima",
            "",
            "Here `T=2K`, `F=sum_m lambda(m)(lambda(m)-1)`, and the exact identity is `integral |S|^2 = T+F`. The minor maximum is over sampled local maxima outside the declared major arcs.",
            "",
            "| N | K | complete T | complete positions in I_N | all doublet shell lifts | F | (T+F)/T | FFT grid | Parseval rel. error | max minor theta | max minor / sqrt(T) | nearest a/q | N offset |",
            "|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|:---:|---:|",
        ]
    )
    for item in results:
        error = abs(item.l2_numeric - item.l2_exact) / item.l2_exact
        peak = item.minor_peak
        retained = sum(
            item.n < position <= 2 * item.n
            for doublet in item.doublets
            for position in (doublet.m1, doublet.m2)
        )
        shell_lifts = sum(
            len(residue_lifts_in_shell(item.n, doublet.p, residue))
            for doublet in item.doublets
            for residue in (doublet.r, doublet.p - 1 - doublet.r)
        )
        lines.append(
            f"| {item.n} | {len(item.doublets)} | {item.t} | {retained} | {shell_lifts} | {item.collision} | "
            f"{item.l2_exact / item.t:.6f} | {item.grid_size:,} | {error:.3e} | "
            f"{peak.theta:.9f} | {peak.root_t:.4f} | {peak.rational} | {peak.scaled_offset:+.3f} |"
        )

    lines.extend(["", f"Top `{peak_count}` sampled local maxima on `0 <= theta <= 1/2` at each scale:", ""])
    for item in results:
        lines.extend(
            [
                f"### Peaks for N={item.n}",
                "",
                "| theta | |S| | |S|/sqrt(T) | arc | nearest a/q | N offset |",
                "|---:|---:|---:|:---:|:---:|---:|",
            ]
        )
        for peak in item.peaks:
            lines.append(
                f"| {peak.theta:.9f} | {peak.magnitude:.4f} | {peak.root_t:.4f} | "
                f"{peak.arc} | {peak.rational} | {peak.scaled_offset:+.3f} |"
            )
        lines.append("")

    lines.extend(
        [
            "## Small-denominator gap/phase correlations",
            "",
            "At each reduced `theta=a/q`, `0<a<q<=%d`, the relevant statistic is the connected complex correlation `Corr(e(3p theta/2), cos(pi theta h_p))`. The raw-gap column instead uses `h_p/p`. Cases where the prime phase or cosine is constant are marked degenerate and excluded from maxima. `|S|/sqrt(T)` is also recorded because a large value can come from a nonzero marginal even when connected correlation is small." % q_corr,
            "",
            "| N | max modulation corr. | theta | median modulation corr. | max raw-gap corr. | theta | max rational |S|/sqrt(T) | theta |",
            "|---:|---:|:---:|---:|---:|:---:|---:|:---:|",
        ]
    )
    for item in results:
        modulation = [row for row in item.correlations if finite(row.modulation_corr)]
        gaps = [row for row in item.correlations if finite(row.gap_corr)]
        max_mod = max(modulation, key=lambda row: row.modulation_corr)
        max_gap = max(gaps, key=lambda row: row.gap_corr)
        max_amp = max(item.correlations, key=lambda row: row.relevant_amplitude)
        median_mod = float(np.median([row.modulation_corr for row in modulation]))
        lines.append(
            f"| {item.n} | {max_mod.modulation_corr:.4f} | {fraction_label(max_mod)} | "
            f"{median_mod:.4f} | {max_gap.gap_corr:.4f} | {fraction_label(max_gap)} | "
            f"{max_amp.relevant_amplitude:.4f} | {fraction_label(max_amp)} |"
        )

    largest = results[-1]
    lines.extend(
        [
            "",
            f"For the largest scale `N={largest.n}`, maxima over reduced numerators at each denominator are:",
            "",
            "| q | max modulation corr. (a/q) | max raw-gap corr. (a/q) | max |S|/sqrt(T) (a/q) |",
            "|---:|:---|:---|:---|",
        ]
    )
    for q in range(2, q_corr + 1):
        q_rows = [row for row in largest.correlations if row.q == q]
        mod_rows = [row for row in q_rows if finite(row.modulation_corr)]
        gap_rows = [row for row in q_rows if finite(row.gap_corr)]
        mod_text = "degenerate"
        gap_text = "degenerate"
        if mod_rows:
            row = max(mod_rows, key=lambda value: value.modulation_corr)
            mod_text = f"{row.modulation_corr:.4f} ({fraction_label(row)})"
        if gap_rows:
            row = max(gap_rows, key=lambda value: value.gap_corr)
            gap_text = f"{row.gap_corr:.4f} ({fraction_label(row)})"
        row = max(q_rows, key=lambda value: value.relevant_amplitude)
        lines.append(
            f"| {q} | {mod_text} | {gap_text} | "
            f"{row.relevant_amplitude:.4f} ({fraction_label(row)}) |"
        )

    lines.extend(
        [
            "",
            "## Collision bilinear model for the complete doublets",
            "",
            "For each observed pair of prime labels `p<q`, the bilinear model independently replaces the first zeros by uniforms `r_p in [1,(p-3)/2]`, `r_q in [1,(q-3)/2]`. It exactly counts the four linear collision equations for `m_i(p)=m_j(q)` and doubles the unordered sum. Thus",
            "",
            "```text",
            "F_bil = 2 sum_{p<q} (C_A(p,q)+C_B(p,q)+C_C(p,q)+C_D(p,q))",
            "                         / (((p-3)/2)((q-3)/2)).",
            "```",
            "",
            "The iid fixed-margin comparator is `F_iid=4K(K-1)/N`. The centered entries below subtract that same comparator; `F_bil` is a geometric independent-gap prediction, while `F` is the observed ordered collision count.",
            "",
            "| N | observed F | F_bil | F_iid | F-F_iid | F_bil-F_iid | observed/bilinear |",
            "|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for item in results:
        ratio = item.collision / item.uniform_gap_collision if item.uniform_gap_collision else float("nan")
        lines.append(
            f"| {item.n} | {item.collision} | {item.uniform_gap_collision:.4f} | "
            f"{item.iid_collision:.4f} | {item.collision - item.iid_collision:+.4f} | "
            f"{item.uniform_gap_collision - item.iid_collision:+.4f} | {fmt(ratio, 4)} |"
        )

    lines.extend(
        [
            "",
            "## True shell off-diagonal and Fourier reconstruction",
            "",
            "The random-CRT prediction has mean zero. Its displayed standard deviation is the exact fixed-margin formula from `prop:random-crt-moment`, with the observed `A_p` and weights. Agreement of `E_direct` and `E_Fourier` is an identity check, not statistical evidence by itself.",
            "",
            "| N | weight | active columns | hits | diagonal S_c | E_direct | E_Fourier | difference | E/S_c | random sd | E/sd | V/S_c |",
            "|---:|:---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for shell in shells:
        for row in shell.rows:
            difference = row.e_direct - row.e_fourier
            lines.append(
                f"| {shell.n} | {row.weight_name} | {row.columns} | {row.hits} | "
                f"{row.diagonal:.4f} | {row.e_direct:+.4f} | {row.e_fourier:+.4f} | "
                f"{difference:+.2e} | {row.e_direct / row.diagonal:+.5f} | "
                f"{row.random_sd:.4f} | {row.e_direct / row.random_sd:+.3f} | "
                f"{row.variance / row.diagonal:.5f} |"
            )

    l2_ratios = [item.l2_exact / item.t for item in results]
    minor_ratios = [item.minor_peak.root_t for item in results]
    weighted_rows = [shell.rows[1] for shell in shells]
    weighted_ratios = [abs(row.e_direct / row.diagonal) for row in weighted_rows]
    weighted_zscores = [abs(row.e_direct / row.random_sd) for row in weighted_rows]
    observed_prediction_ratios = [
        item.collision / item.uniform_gap_collision for item in results
    ]
    lines.extend(
        [
            "",
            "## What the computation says (and does not say)",
            "",
            f"- Complete-doublet Parseval ratios `(T+F)/T` lie in `{min(l2_ratios):.4f}..{max(l2_ratios):.4f}`. Thus the integrated norm is cleanly at the diagonal scale on every tested block.",
            f"- With the explicitly declared `q<={q_major}` major arcs removed, the largest sampled minor-arc values are `{min(minor_ratios):.2f}..{max(minor_ratios):.2f}` times `sqrt(T)`. This is compatible with square-root-scale cancellation plus an extreme-value factor; it is not a uniform minor-arc bound.",
            f"- The independent-uniform-gap bilinear model predicts the observed collision energy only to constant-factor/fluctuation accuracy: observed/predicted ranges from `{min(observed_prediction_ratios):.3f}` to `{max(observed_prediction_ratios):.3f}`. It is a useful null geometry, not an arithmetic theorem about `r_p`.",
            f"- In the true shell with logarithmic weights, `|E^o|/S_c` ranges from `{min(weighted_ratios):.5f}` to `{max(weighted_ratios):.5f}`, the ratios are nonmonotone, and the largest fixed-margin null-model magnitude is `{max(weighted_zscores):.2f}` standard deviations. The direct/Fourier computations agree to the displayed precision.",
            "- Small-denominator correlations must be read denominator by denominator: exact or near degeneracies of the prime phase are major-arc structure, while the connected nondegenerate correlations are finite-sample diagnostics. None supplies the horizontal cross-prime estimate required by `eq:short-arc`.",
            "- The two shell-count columns quantify a boundary/lifting issue: the complete reflected-pair sum contains only the designated `k=1` lifts `p+r`, while a true shell column can omit these below `N` and include higher lifts (`2p+r`, and sometimes `3p+r`). Thus the complete sum is not the exact Fourier transform of the shell columns. Any analytic reduction from this `S(theta)` to shell AMTD must insert the interval-dependent lift sums (or average/translate the interval); the complete-pair Parseval identity alone does not make that reduction.",
            "- Consequently these computations support the doublet-cancellation heuristic but do not remove the generic `P^2` large-sieve barrier and do not prove `hyp:amtd`.",
            "",
        ]
    )
    return "\n".join(lines)


def parse_args(argv: Sequence[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=DEFAULT_DATA)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--n", type=int, nargs="+", default=list(DEFAULT_NS))
    parser.add_argument("--grid-factor", type=int, default=32)
    parser.add_argument("--q-major", type=int, default=16)
    parser.add_argument("--q-corr", type=int, default=16)
    parser.add_argument("--peak-count", type=int, default=8)
    return parser.parse_args(argv)


def main(argv: Sequence[str] = ()) -> int:
    args = parse_args(argv or sys.argv[1:])
    ns = sorted(set(args.n))
    if not ns or min(ns) < 8192:
        raise SystemExit("all requested N must be at least 8192")
    if args.grid_factor < 4:
        raise SystemExit("--grid-factor must be at least 4")
    if args.q_major < 1 or args.q_corr < 2 or args.peak_count < 1:
        raise SystemExit("invalid q/peak parameter")

    check_counting_formulas()
    check_shell_lift_enumerator()
    zeros, digest, record_count = load_zero_sets(args.data, max(ns))
    validated_sets = validate_zero_sets(zeros)
    print(
        f"loaded {record_count:,} records; analyzing N={','.join(map(str, ns))}",
        file=sys.stderr,
        flush=True,
    )
    results: List[FourierResult] = []
    shells: List[ShellResult] = []
    for n in ns:
        print(f"N={n}: doublet FFT and correlations", file=sys.stderr, flush=True)
        results.append(
            analyze_fourier(
                zeros,
                n,
                args.grid_factor,
                args.q_major,
                args.q_corr,
                args.peak_count,
            )
        )
        print(f"N={n}: shell covariance", file=sys.stderr, flush=True)
        shells.append(analyze_shell(zeros, n))

    command = "python3 problems/3.2/oracleA_explore.py"
    report = render_report(
        results,
        shells,
        args.data,
        digest,
        record_count,
        validated_sets,
        args.grid_factor,
        args.q_major,
        args.q_corr,
        args.peak_count,
        command,
    )
    args.output.write_text(report, encoding="utf-8")
    print(f"wrote {args.output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
