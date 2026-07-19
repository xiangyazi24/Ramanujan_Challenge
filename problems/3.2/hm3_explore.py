#!/usr/bin/env python3
"""Direct exploration of the third factorial moment in Problem 3.2.

For each requested X, this program computes the actual Apéry zero sets for
primes in (X, 2X] and enumerates every *canonical* zero triple

    p < q < s,  r_p in Z_p, r_q in Z_q, r_s in Z_s

whose CRT representative lies in [0, X^2).  Each canonical triple represents
exactly 3! ordered triples.  The enumeration is over prime/residue tuples, not
over an array of m-values.  A separate m-scatter is used only as an independent
cross-check and for conditional statistics.

The default invocation reproduces ``hm3_exploration.md``:

    python3 problems/3.2/hm3_explore.py

The implementation uses only the Python standard library.
"""

from __future__ import annotations

import argparse
import bisect
import math
import statistics
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from decimal import Decimal, getcontext
from pathlib import Path
from typing import Dict, FrozenSet, Iterable, List, Mapping, Sequence, Tuple


getcontext().prec = 50

DEFAULT_XS = (256, 512, 1024, 2048, 4096)
# Three-decimal values printed by scripts/p32_hm_check.py.  The TeX table
# subsequently displays these at two decimals (including 0.455 as 0.46).
CITED_R3_3DP = {256: "0.455", 512: "0.717", 1024: "0.833", 2048: "0.933"}


@dataclass(frozen=True)
class PrimeZeroSet:
    p: int
    zeros: Tuple[int, ...]
    zero_set: FrozenSet[int]


@dataclass(frozen=True)
class Triple:
    """A canonical (in increasing-prime order) contributing zero triple."""

    m: int
    p: int
    rp: int
    q: int
    rq: int
    s: int
    rs: int


@dataclass(frozen=True)
class PairData:
    p: int
    q: int
    zp: int
    zq: int
    reps: Tuple[int, ...]
    main: float
    error: float


@dataclass
class XResult:
    x: int
    prime_count: int
    active_count: int
    zero_count: int
    lam: Decimal
    expected3: Decimal
    s1: int
    s2: int
    s3: int
    r2: Decimal
    r3: Decimal
    max_k: int
    canonical_triples: int
    diagnostics: Mapping[str, float]
    section: str


def primes_upto(n: int) -> List[int]:
    """Return all primes <= n by an ordinary bytearray sieve."""

    if n < 2:
        return []
    sieve = bytearray(b"\x01") * (n + 1)
    sieve[0:2] = b"\x00\x00"
    for p in range(2, math.isqrt(n) + 1):
        if sieve[p]:
            start = p * p
            sieve[start : n + 1 : p] = b"\x00" * (((n - start) // p) + 1)
    return [p for p in range(2, n + 1) if sieve[p]]


def apery_zeros(p: int) -> Tuple[int, ...]:
    """Compute {0 <= r < p : b_r = 0 (mod p)} in O(p) operations.

    We iterate Y_n=(n!)^3 b_n, for which

        Y_{n+1} = P(n)Y_n - n^6 Y_{n-1}.

    Since n! is a unit for n < p, Y_n and b_n have the same zero set.
    """

    if p < 2:
        raise ValueError("p must be prime")
    zeros: List[int] = []
    ym1, y = 1 % p, 5 % p
    if ym1 == 0:
        zeros.append(0)
    if p > 1 and y == 0:
        zeros.append(1)
    for n in range(1, p - 1):
        poly = (34 * n * n * n + 51 * n * n + 27 * n + 5) % p
        ynext = (poly * y - pow(n, 6, p) * ym1) % p
        ym1, y = y, ynext
        if y == 0:
            zeros.append(n + 1)
    return tuple(zeros)


def apery_zeros_division(p: int) -> Tuple[int, ...]:
    """Slow, independent check using the original divided recurrence."""

    bm1, b = 1 % p, 5 % p
    zeros: List[int] = []
    if bm1 == 0:
        zeros.append(0)
    if b == 0:
        zeros.append(1)
    for n in range(1, p - 1):
        poly = (34 * n * n * n + 51 * n * n + 27 * n + 5) % p
        numerator = (poly * b - pow(n, 3, p) * bm1) % p
        inv = pow(n + 1, -1, p)
        bnext = numerator * pow(inv, 3, p) % p
        bm1, b = b, bnext
        if b == 0:
            zeros.append(n + 1)
    return tuple(zeros)


def build_zero_sets(x: int, deep_check: bool) -> Tuple[List[PrimeZeroSet], int]:
    primes = [p for p in primes_upto(2 * x) if p > x and p >= 7]
    data: List[PrimeZeroSet] = []
    for p in primes:
        zeros = apery_zeros(p)
        data.append(PrimeZeroSet(p, zeros, frozenset(zeros)))

    if deep_check:
        check_indices = range(len(data))
    else:
        candidates = {
            0,
            len(data) // 4,
            len(data) // 2,
            (3 * len(data)) // 4,
            len(data) - 1,
        }
        active_indices = [i for i, item in enumerate(data) if item.zeros]
        if active_indices:
            candidates.update((active_indices[0], active_indices[-1]))
        check_indices = sorted(i for i in candidates if 0 <= i < len(data))

    checked = 0
    for i in check_indices:
        item = data[i]
        reference = apery_zeros_division(item.p)
        if reference != item.zeros:
            raise AssertionError(
                f"zero-set recurrence mismatch at p={item.p}: "
                f"cleared={item.zeros}, divided={reference}"
            )
        checked += 1
    return data, checked


def enumerate_pairs_and_triples(
    active: Sequence[PrimeZeroSet], m_limit: int
) -> Tuple[List[PairData], List[Triple]]:
    """Directly enumerate pair representatives and canonical zero triples.

    For p<q, pq>X^2.  Thus a residue pair has at most one representative
    ``a`` below X^2.  A triple p<q<s contributes exactly when that same a has
    a mod s in Z_s.  This is precisely a CRT-tuple enumeration, with no scan
    over m.
    """

    pairs: List[PairData] = []
    triples: List[Triple] = []
    for i, left in enumerate(active):
        p = left.p
        for j in range(i + 1, len(active)):
            right = active[j]
            q = right.p
            inv_p_mod_q = pow(p, -1, q)
            reps: List[int] = []
            for rp in left.zeros:
                for rq in right.zeros:
                    a = rp + p * (((rq - rp) * inv_p_mod_q) % q)
                    if a < m_limit:
                        reps.append(a)

            zprod = len(left.zeros) * len(right.zeros)
            main = m_limit * zprod / (p * q)
            pairs.append(
                PairData(p, q, len(left.zeros), len(right.zeros), tuple(reps), main, len(reps) - main)
            )

            # The outer indices explicitly impose p<q<s, so every canonical
            # prime triple and its unique residue triple are visited once.
            for a in reps:
                for k in range(j + 1, len(active)):
                    third = active[k]
                    rs = a % third.p
                    if rs in third.zero_set:
                        triples.append(
                            Triple(a, p, a % p, q, a % q, third.p, rs)
                        )
    return pairs, triples


def scatter_k(active: Sequence[PrimeZeroSet], m_limit: int) -> bytearray:
    """Independent m-scatter, used only after direct tuple enumeration."""

    k_values = bytearray(m_limit)
    for item in active:
        for r in item.zeros:
            for m in range(r, m_limit, item.p):
                if k_values[m] == 255:
                    raise OverflowError("K_X(m) does not fit in a byte")
                k_values[m] += 1
    return k_values


def decimal_sum_ratio(terms: Iterable[Tuple[int, int]]) -> Decimal:
    return sum((Decimal(a) / Decimal(b) for a, b in terms), Decimal(0))


def side(p: int, r: int) -> str:
    twice = 2 * r
    centre = p - 1
    if twice < centre:
        return "L"
    if twice > centre:
        return "R"
    return "C"


def reflected_side_pattern(pattern: str) -> str:
    return pattern.translate(str.maketrans({"L": "R", "R": "L", "C": "C"}))


def pearson(xs: Sequence[float], ys: Sequence[float]) -> float:
    if len(xs) < 2 or len(xs) != len(ys):
        return float("nan")
    mx = statistics.fmean(xs)
    my = statistics.fmean(ys)
    vx = math.fsum((x - mx) ** 2 for x in xs)
    vy = math.fsum((y - my) ** 2 for y in ys)
    if vx == 0.0 or vy == 0.0:
        return float("nan")
    cov = math.fsum((x - mx) * (y - my) for x, y in zip(xs, ys))
    return cov / math.sqrt(vx * vy)


def bin_label(value: float, edges: Sequence[float]) -> str:
    index = bisect.bisect_right(edges, value) - 1
    index = max(0, min(index, len(edges) - 2))

    def fmt(v: float) -> str:
        if math.isinf(v):
            return "-inf" if v < 0 else "+inf"
        return f"{v:g}"

    return f"[{fmt(edges[index])}, {fmt(edges[index + 1])})"


def dyadic_span_label(span: int) -> str:
    if span == 0:
        return "0"
    if span == 1:
        return "1"
    low = 1 << (span.bit_length() - 1)
    return f"{low}-{2 * low - 1}"


def md_table(headers: Sequence[str], rows: Iterable[Sequence[object]]) -> str:
    rows_list = [[str(value) for value in row] for row in rows]
    output = ["| " + " | ".join(headers) + " |"]
    output.append("|" + "|".join("---" for _ in headers) + "|")
    output.extend("| " + " | ".join(row) + " |" for row in rows_list)
    return "\n".join(output)


def fmt_decimal(value: Decimal, digits: int = 9) -> str:
    return f"{value:.{digits}f}"


def hit_data(m: int, active: Sequence[PrimeZeroSet]) -> List[Tuple[int, int, int]]:
    return [
        (item.p, m % item.p, m // item.p)
        for item in active
        if m % item.p in item.zero_set
    ]


def analyze_x(x: int, deep_check: bool, verbose: bool) -> XResult:
    if x < 8:
        raise ValueError(f"X must be at least 8, got {x}")
    if verbose:
        print(f"X={x}: computing zero sets", file=sys.stderr, flush=True)
    all_data, recurrence_checks = build_zero_sets(x, deep_check)
    active = [item for item in all_data if item.zeros]
    m_limit = x * x

    # Reflection is a theorem available in proof.tex, but checking it here is
    # a useful guard against recurrence/indexing mistakes.
    reflection_failures = [
        (item.p, r)
        for item in all_data
        for r in item.zeros
        if item.p - 1 - r not in item.zero_set
    ]
    if reflection_failures:
        raise AssertionError(f"zero-set reflection failures: {reflection_failures[:5]}")

    if verbose:
        print(
            f"X={x}: direct CRT enumeration over {len(active)} active primes",
            file=sys.stderr,
            flush=True,
        )
    pairs, triples = enumerate_pairs_and_triples(active, m_limit)

    if verbose:
        print(f"X={x}: independent m-scatter cross-check", file=sys.stderr, flush=True)
    k_values = scatter_k(active, m_limit)
    k_hist = Counter(k_values)
    s1 = sum(k * count for k, count in k_hist.items())
    s2 = sum(k * (k - 1) * count for k, count in k_hist.items())
    s3 = sum(k * (k - 1) * (k - 2) * count for k, count in k_hist.items())
    max_k = max(k_hist)

    canonical_count = len(triples)
    pair_rep_count = sum(len(pair.reps) for pair in pairs)
    triple_m_hist = Counter(triple.m for triple in triples)
    combination_count = sum(math.comb(k, 3) * count for k, count in k_hist.items() if k >= 3)

    # Four independent combinatorial decompositions must agree exactly.
    if 6 * canonical_count != s3:
        raise AssertionError(
            f"X={x}: direct triple count {6 * canonical_count} != scatter S3 {s3}"
        )
    if 2 * pair_rep_count != s2:
        raise AssertionError(
            f"X={x}: direct pair count {2 * pair_rep_count} != scatter S2 {s2}"
        )
    if combination_count != canonical_count:
        raise AssertionError(
            f"X={x}: sum_m C(K,3)={combination_count} != {canonical_count}"
        )
    for m, count in triple_m_hist.items():
        expected = math.comb(k_values[m], 3)
        if count != expected:
            raise AssertionError(
                f"X={x}, m={m}: direct multiplicity {count} != C(K,3)={expected}"
            )

    lam = decimal_sum_ratio((len(item.zeros), item.p) for item in all_data)
    expected2 = Decimal(m_limit) * lam * lam
    expected3 = expected2 * lam
    r2 = Decimal(s2) / expected2 if expected2 else Decimal("NaN")
    r3 = Decimal(s3) / expected3 if expected3 else Decimal("NaN")
    zero_count = sum(len(item.zeros) for item in all_data)

    cited_note = "not tabulated"
    if x in CITED_R3_3DP:
        cited = CITED_R3_3DP[x]
        rounded = f"{r3:.3f}"
        if rounded != cited:
            raise AssertionError(
                f"X={x}: R3={r3} rounds to {rounded}, not cited {cited}"
            )
        cited_note = f"PASS (p32_hm_check.py prints {cited})"

    lines: List[str] = []
    lines.append(f"## X = {x}\n")
    lines.append(
        md_table(
            ("quantity", "value"),
            (
                ("primes in (X,2X]", len(all_data)),
                ("active primes Z(p)>0", len(active)),
                ("sum Z(p)", zero_count),
                ("lambda_X", fmt_decimal(lam, 12)),
                ("direct canonical triples p<q<s", canonical_count),
                ("exact ordered S3", s3),
                ("X^2 lambda_X^3", fmt_decimal(expected3, 6)),
                ("R3", fmt_decimal(r3, 9)),
                ("R3 table cross-check", cited_note),
                ("exact ordered S2", s2),
                ("R2", fmt_decimal(r2, 9)),
                ("max K_X(m)", max_k),
            ),
        )
    )
    lines.append(
        "\nThe direct enumerator found one canonical tuple for each `p<q<s`; "
        "multiplication by 6 gives the ordered factorial moment. The separate "
        f"m-scatter agrees exactly: `S3 = 6*{canonical_count} = {s3}` and "
        f"`S2 = 2*{pair_rep_count} = {s2}`. The original divided Apéry "
        f"recurrence was also checked at {recurrence_checks} deterministic sample "
        "primes (at every prime under `--deep-check`), and every zero set passed "
        "reflection closure.\n"
    )

    # Quotient-pattern analysis.
    quotient_hist: Counter[Tuple[int, int, int]] = Counter()
    quotient_class = Counter()
    quotient_spans = Counter()
    for triple in triples:
        qt = (triple.m // triple.p, triple.m // triple.q, triple.m // triple.s)
        if not (qt[0] >= qt[1] >= qt[2]):
            raise AssertionError(f"nonmonotone quotient tuple: {qt}")
        quotient_hist[qt] += 1
        if qt[0] == qt[2]:
            quotient_class["all equal"] += 1
        elif qt[0] == qt[1]:
            quotient_class["first two equal"] += 1
        elif qt[1] == qt[2]:
            quotient_class["last two equal"] += 1
        else:
            quotient_class["all distinct"] += 1
        quotient_spans[dyadic_span_label(qt[0] - qt[2])] += 1

    lines.append("\n### Quotient patterns\n")
    if canonical_count:
        qclass_order = ("all equal", "first two equal", "last two equal", "all distinct")
        lines.append(
            md_table(
                ("pattern for q_p >= q_q >= q_s", "canonical triples", "fraction"),
                (
                    (label, quotient_class[label], f"{quotient_class[label] / canonical_count:.4f}")
                    for label in qclass_order
                ),
            )
        )
        lines.append("\nMost frequent exact quotient triples:\n")
        top_q = sorted(quotient_hist.items(), key=lambda item: (-item[1], item[0]))[:12]
        lines.append(
            md_table(
                ("(q_p,q_q,q_s)", "canonical triples"),
                ((str(qt), count) for qt, count in top_q),
            )
        )
        lines.append("\nDyadic histogram of quotient span `q_p-q_s`:\n")

        def span_key(label: str) -> int:
            return int(label.split("-")[0])

        lines.append(
            md_table(
                ("span", "canonical triples", "fraction"),
                (
                    (label, quotient_spans[label], f"{quotient_spans[label] / canonical_count:.4f}")
                    for label in sorted(quotient_spans, key=span_key)
                ),
            )
        )
    else:
        lines.append("No contributing triples.\n")

    # Residue and reflection-orbit analysis.
    side_hist = Counter()
    distance_edges = (0.0, 0.02, 0.05, 0.10, 0.20, 0.30, 0.40, 0.5000001)
    distance_hist = Counter()
    orbit_patterns: Dict[Tuple[Tuple[int, int], ...], set] = defaultdict(set)
    near_thresholds = (0.01, 0.02, 0.05, 0.10)
    near_counts = Counter()
    aligned_signed = Counter()
    aligned_base = Counter()
    signed_coords = ([], [], [])
    base_coords = ([], [], [])

    for triple in triples:
        coords = ((triple.p, triple.rp), (triple.q, triple.rq), (triple.s, triple.rs))
        pattern = "".join(side(p, r) for p, r in coords)
        side_hist[pattern] += 1
        orbit_key = tuple((p, min(r, p - 1 - r)) for p, r in coords)
        orbit_patterns[orbit_key].add(pattern)
        nearest = [min(r, p - r) / p for p, r in coords]
        bases = [min(r, p - 1 - r) / (p - 1) for p, r in coords]
        signed = [r / p for p, r in coords]
        for value in bases:
            distance_hist[bin_label(value, distance_edges)] += 1
        for threshold in near_thresholds:
            if min(nearest) < threshold:
                near_counts[threshold] += 1
        for threshold in (0.02, 0.05, 0.10):
            if max(signed) - min(signed) < threshold:
                aligned_signed[threshold] += 1
            if max(bases) - min(bases) < threshold:
                aligned_base[threshold] += 1
        for i in range(3):
            signed_coords[i].append(signed[i])
            base_coords[i].append(bases[i])

    occupancy = Counter(len(patterns) for patterns in orbit_patterns.values())
    complement_pairs = 0
    for patterns in orbit_patterns.values():
        for pattern in patterns:
            complement = reflected_side_pattern(pattern)
            if pattern < complement and complement in patterns:
                complement_pairs += 1

    lines.append("\n### Residues and reflection\n")
    if canonical_count:
        lines.append("Reflection-side patterns (`L/R` relative to `(p-1)/2`):\n")
        lines.append(
            md_table(
                ("pattern", "canonical triples", "fraction"),
                (
                    (pattern, count, f"{count / canonical_count:.4f}")
                    for pattern, count in sorted(side_hist.items())
                ),
            )
        )
        lines.append("\nMarginal reflected endpoint distance `min(r,p-1-r)/(p-1)`:\n")
        lines.append(
            md_table(
                ("distance bin", "residue coordinates", "fraction"),
                (
                    (label, count, f"{count / (3 * canonical_count):.4f}")
                    for label, count in sorted(
                        distance_hist.items(), key=lambda item: float(item[0].split(",")[0][1:])
                    )
                ),
            )
        )
        lines.append("\nNear-multiple and normalized-alignment diagnostics:\n")
        diagnostics = []
        for threshold in near_thresholds:
            diagnostics.append(
                (
                    f"some min(r,p-r)/p < {threshold:.2f}",
                    near_counts[threshold],
                    f"{near_counts[threshold] / canonical_count:.4f}",
                )
            )
        for threshold in (0.02, 0.05, 0.10):
            diagnostics.append(
                (
                    f"range(r_i/p_i) < {threshold:.2f}",
                    aligned_signed[threshold],
                    f"{aligned_signed[threshold] / canonical_count:.4f}",
                )
            )
            diagnostics.append(
                (
                    f"range(reflected distances) < {threshold:.2f}",
                    aligned_base[threshold],
                    f"{aligned_base[threshold] / canonical_count:.4f}",
                )
            )
        lines.append(md_table(("event", "canonical triples", "fraction"), diagnostics))
        lines.append(
            "\nPairwise Pearson correlations of normalized coordinates:\n\n"
            + md_table(
                ("coordinates", "(p,q)", "(p,s)", "(q,s)"),
                (
                    (
                        "signed r/p",
                        f"{pearson(signed_coords[0], signed_coords[1]):+.4f}",
                        f"{pearson(signed_coords[0], signed_coords[2]):+.4f}",
                        f"{pearson(signed_coords[1], signed_coords[2]):+.4f}",
                    ),
                    (
                        "reflected distance",
                        f"{pearson(base_coords[0], base_coords[1]):+.4f}",
                        f"{pearson(base_coords[0], base_coords[2]):+.4f}",
                        f"{pearson(base_coords[1], base_coords[2]):+.4f}",
                    ),
                ),
            )
        )
        lines.append(
            "\nReflection-orbit occupancy (an orbit fixes each pair "
            "`{r,p-1-r}`):\n\n"
            + md_table(
                ("contributing sign variants", "orbits", "canonical triples"),
                (
                    (size, count, size * count)
                    for size, count in sorted(occupancy.items())
                ),
            )
        )
        lines.append(
            f"\nThere are **{complement_pairs}** nontrivial fully reflected pairs "
            "among contributing variants. (For a tuple with modulus product `P`, "
            "full reflection sends `m` to `P-1-m`, which is normally far above "
            "`X^2`.)\n"
        )
    else:
        lines.append("No contributing triples.\n")

    # Special-m clustering and the window radical c_m.
    special_ms = [m for m, k in enumerate(k_values) if k >= 3]
    special_details: Dict[int, List[Tuple[int, int, int]]] = {
        m: hit_data(m, active) for m in special_ms
    }
    decile_ms = Counter(int(10 * m / m_limit) for m in special_ms)
    decile_triples = Counter(int(10 * triple.m / m_limit) for triple in triples)

    lines.append("\n### Clustering on special m\n")
    lines.append(
        md_table(
            ("K", "number of m", "canonical C(K,3) contribution", "ordered contribution"),
            (
                (
                    k,
                    k_hist[k],
                    math.comb(k, 3) * k_hist[k] if k >= 3 else 0,
                    k * (k - 1) * (k - 2) * k_hist[k] if k >= 3 else 0,
                )
                for k in sorted(k_hist)
                if k >= 1
            ),
        )
    )
    lines.append("\nLocation in `[0,X^2)`:\n")
    lines.append(
        md_table(
            ("m/X^2 bin", "distinct m with K>=3", "canonical triples"),
            (
                (f"[{d/10:.1f},{(d+1)/10:.1f})", decile_ms[d], decile_triples[d])
                for d in range(10)
            ),
        )
    )

    def nearest_fraction(m: int) -> float:
        return min(min(r, p - r) / p for p, r, _ in special_details[m])

    top_ms = sorted(special_ms, key=lambda m: (-k_values[m], m))[:12]
    near_ms = sorted(special_ms, key=lambda m: (nearest_fraction(m), -k_values[m], m))[:12]

    def special_rows(ms: Sequence[int]) -> Iterable[Sequence[object]]:
        for m in ms:
            hits = special_details[m]
            radical = math.prod(p for p, _, _ in hits)
            hit_string = "; ".join(f"{p}:{r}/{q}" for p, r, q in hits)
            yield (
                m,
                f"{m / m_limit:.6f}",
                k_values[m],
                f"{nearest_fraction(m):.6g}",
                radical,
                hit_string,
            )

    if top_ms:
        lines.append(
            "\nLargest-K representatives (`p:r/q` means prime, residue, quotient):\n\n"
            + md_table(
                ("m", "m/X^2", "K", "nearest multiple", "c_m", "p:r/q hits"),
                special_rows(top_ms),
            )
        )
        lines.append(
            "\nRepresentatives closest to a prime multiple:\n\n"
            + md_table(
                ("m", "m/X^2", "K", "nearest multiple", "c_m", "p:r/q hits"),
                special_rows(near_ms),
            )
        )
    lines.append(
        "\nHere `c_m` is the exactly known squarefree radical contributed by "
        "bad primes in `(X,2X]`. The computation does not construct or factor "
        "the enormous integer `b_m`, so it cannot assert smoothness of the full "
        "Apéry number; the K-histogram is the reproducible window-smoothness "
        "diagnostic.\n"
    )

    # Pair errors and conditional third moment.
    raw_edges = (-math.inf, -8.0, -4.0, -2.0, -1.0, 0.0, 1.0, 2.0, 4.0, 8.0, math.inf)
    norm_edges = (-1.0000001, -0.75, -0.50, -0.25, 0.0, 0.25, 0.50, 0.75, 1.0000001)
    raw_error_hist = Counter(bin_label(pair.error, raw_edges) for pair in pairs)
    norm_error_hist = Counter(
        bin_label(pair.error / (pair.zp * pair.zq), norm_edges) for pair in pairs
    )
    a_hist = Counter(len(pair.reps) for pair in pairs)
    errors = [pair.error for pair in pairs]
    normalized_errors = [pair.error / (pair.zp * pair.zq) for pair in pairs]
    total_error = math.fsum(errors)

    conditional: Dict[int, Counter] = defaultdict(Counter)
    extension_by_pair: List[int] = []
    lam_by_prime = {
        item.p: Decimal(len(item.zeros)) / Decimal(item.p) for item in active
    }
    active_counts = {
        item.p: sum((m_limit - 1 - r) // item.p + 1 for r in item.zeros)
        for item in active
    }
    rho_by_prime = {
        p: Decimal(count) / Decimal(m_limit) for p, count in active_counts.items()
    }
    mu = Decimal(s1) / Decimal(m_limit)
    palm_periodic_unordered = Decimal(0)
    palm_finite_unordered = Decimal(0)
    baseline_periodic_unordered = Decimal(0)
    baseline_finite_unordered = Decimal(0)
    for pair in pairs:
        a = len(pair.reps)
        extensions = sum(k_values[m] - 2 for m in pair.reps)
        extension_by_pair.append(extensions)
        conditional[a]["pairs"] += 1
        conditional[a]["pair_hits"] += a
        conditional[a]["extensions"] += extensions
        conditional[a]["positive_pairs"] += int(extensions > 0)
        lambda_pq = max(
            lam - lam_by_prime[pair.p] - lam_by_prime[pair.q], Decimal(0)
        )
        remaining_count = s1 - active_counts[pair.p] - active_counts[pair.q]
        if remaining_count < 0:
            raise AssertionError(f"X={x}: negative finite pair marginal")
        mu_pq = Decimal(remaining_count) / Decimal(m_limit)
        periodic_baseline = Decimal(a) * lambda_pq
        finite_baseline = Decimal(a) * mu_pq
        baseline_periodic_unordered += periodic_baseline
        baseline_finite_unordered += finite_baseline
        palm_periodic_unordered += max(
            Decimal(extensions) - periodic_baseline, Decimal(0)
        )
        palm_finite_unordered += max(
            Decimal(extensions) - finite_baseline, Decimal(0)
        )

    total_extensions = sum(extension_by_pair)
    if 2 * total_extensions != s3 or total_extensions != 3 * canonical_count:
        raise AssertionError(
            f"X={x}: conditional extensions {total_extensions} disagree with S3={s3}"
        )

    # The theorem in hm3_result.tex uses ordered pairs.  The direct pair
    # enumerator uses p<q, so every aggregate below is doubled exactly.
    palm_periodic = 2 * palm_periodic_unordered
    palm_finite = 2 * palm_finite_unordered
    baseline_periodic = 2 * baseline_periodic_unordered
    baseline_finite = 2 * baseline_finite_unordered
    s3_decimal = Decimal(s3)
    if not (
        palm_periodic <= s3_decimal <= palm_periodic + baseline_periodic
    ):
        raise AssertionError(f"X={x}: periodic pair-Palm sandwich failed")
    if not (palm_finite <= s3_decimal <= palm_finite + baseline_finite):
        raise AssertionError(f"X={x}: finite pair-Palm sandwich failed")
    if baseline_periodic > 5 * expected3:
        raise AssertionError(f"X={x}: periodic Palm baseline exceeds HM2 bound")
    finite_hm2_bound = 5 * (Decimal(1) + Decimal(2) / Decimal(x)) * expected3
    if baseline_finite > finite_hm2_bound:
        raise AssertionError(f"X={x}: finite Palm baseline exceeds HM2 bound")

    pair_main_decimal = Decimal(m_limit) * (
        lam * lam
        - sum(
            (Decimal(len(item.zeros)) / Decimal(item.p)) ** 2
            for item in active
        )
    )
    error_from_identity = (Decimal(s2) - pair_main_decimal) / Decimal(2)
    if abs(float(error_from_identity) - total_error) > 1e-7 * max(1.0, abs(total_error)):
        raise AssertionError(
            f"X={x}: pair-error sum mismatch {total_error} vs {error_from_identity}"
        )

    total_prime_pairs = math.comb(len(all_data), 2)
    inactive_pairs = total_prime_pairs - len(pairs)
    lines.append("\n### Pair errors and conditional third moment\n")
    lines.append(
        md_table(
            ("pair statistic", "value"),
            (
                ("all unordered prime pairs", total_prime_pairs),
                ("active pairs Z(p)Z(q)>0", len(pairs)),
                ("inactive pairs (E=0)", inactive_pairs),
                ("sum A(p,q)", pair_rep_count),
                ("sum E(p,q), active pairs", f"{total_error:+.9f}"),
                ("mean E", f"{statistics.fmean(errors):+.9f}" if errors else "nan"),
                (
                    "population sd(E)",
                    f"{statistics.pstdev(errors):.9f}" if len(errors) >= 2 else "nan",
                ),
                (
                    "mean E/(ZpZq)",
                    f"{statistics.fmean(normalized_errors):+.9f}" if normalized_errors else "nan",
                ),
                ("sum conditional extensions", total_extensions),
                ("corr(E, extensions per pair)", f"{pearson(errors, extension_by_pair):+.5f}"),
                ("ordered periodic Palm excess", fmt_decimal(palm_periodic, 9)),
                (
                    "periodic Palm excess / S3",
                    fmt_decimal(palm_periodic / s3_decimal, 9) if s3 else "--",
                ),
                ("ordered finite Palm excess", fmt_decimal(palm_finite, 9)),
                (
                    "finite Palm excess / S3",
                    fmt_decimal(palm_finite / s3_decimal, 9) if s3 else "--",
                ),
            ),
        )
    )
    lines.append(
        "\nHere `A(p,q)` is the exact pair intersection count and "
        "`E(p,q)=A(p,q)-X^2 Z(p)Z(q)/(pq)`. Pairs involving an empty zero "
        "set have exactly zero error and are listed separately. Raw active-pair "
        "error histogram:\n\n"
        + md_table(
            ("E bin", "active pairs", "fraction"),
            (
                (label, raw_error_hist[label], f"{raw_error_hist[label] / len(pairs):.4f}")
                for label in [bin_label((raw_edges[i] + raw_edges[i + 1]) / 2 if not math.isinf(raw_edges[i]) and not math.isinf(raw_edges[i + 1]) else (-9 if i == 0 else 9), raw_edges) for i in range(len(raw_edges) - 1)]
            )
            if pairs
            else (),
        )
    )
    lines.append(
        "\nNormalized error histogram:\n\n"
        + md_table(
            ("E/(ZpZq) bin", "active pairs", "fraction"),
            (
                (
                    f"[{norm_edges[i]:g}, {norm_edges[i + 1]:g})",
                    norm_error_hist[f"[{norm_edges[i]:g}, {norm_edges[i + 1]:g})"],
                    f"{norm_error_hist[f'[{norm_edges[i]:g}, {norm_edges[i + 1]:g})'] / len(pairs):.4f}",
                )
                for i in range(len(norm_edges) - 1)
            )
            if pairs
            else (),
        )
    )
    lines.append("\nExact histogram of pair counts `A(p,q)`:\n")
    lines.append(
        md_table(
            ("A", "active pairs", "fraction"),
            (
                (a, count, f"{count / len(pairs):.4f}")
                for a, count in sorted(a_hist.items())
            )
            if pairs
            else (),
        )
    )
    lines.append(
        "\nConditional third moment: for an unordered pair `{p,q}`, its "
        "extension count is `T_pq=sum_{m in Omega_p cap Omega_q}(K_X(m)-2)`. "
        "Thus `2 sum T_pq=S3`. Grouping by the exact pair count A gives:\n"
    )
    cond_rows = []
    for a in sorted(conditional):
        row = conditional[a]
        hits = row["pair_hits"]
        extensions = row["extensions"]
        cond_rows.append(
            (
                a,
                row["pairs"],
                row["positive_pairs"],
                hits,
                extensions,
                f"{extensions / hits:.6f}" if hits else "--",
                f"{(2 * extensions / s3):.4f}" if s3 else "--",
            )
        )
    lines.append(
        md_table(
            (
                "A",
                "pairs",
                "pairs with T>0",
                "pair hits",
                "T extensions",
                "T/pair hit",
                "fraction of S3",
            ),
            cond_rows,
        )
    )

    # Release the X^2-sized bytearray before the next scale.
    del k_values
    section = "\n".join(lines).rstrip() + "\n"
    return XResult(
        x=x,
        prime_count=len(all_data),
        active_count=len(active),
        zero_count=zero_count,
        lam=lam,
        expected3=expected3,
        s1=s1,
        s2=s2,
        s3=s3,
        r2=r2,
        r3=r3,
        max_k=max_k,
        canonical_triples=canonical_count,
        diagnostics={
            "repeated_quotient_triples": float(
                canonical_count - quotient_class["all distinct"]
            ),
            "all_equal_quotient_triples": float(quotient_class["all equal"]),
            "maximum_reflection_orbit_occupancy": float(max(occupancy, default=0)),
            "near_one_percent_triples": float(near_counts[0.01]),
            "maximum_abs_signed_residue_correlation": max(
                (
                    abs(pearson(signed_coords[0], signed_coords[1])),
                    abs(pearson(signed_coords[0], signed_coords[2])),
                    abs(pearson(signed_coords[1], signed_coords[2])),
                ),
                default=float("nan"),
            ),
            "special_m_count": float(len(special_ms)),
            "k_at_least_four_m_count": float(
                sum(count for k, count in k_hist.items() if k >= 4)
            ),
            "special_decile_min": float(min(decile_ms.values(), default=0)),
            "special_decile_max": float(max(decile_ms.values(), default=0)),
            "pair_error_sum": total_error,
            "pair_error_mean": statistics.fmean(errors) if errors else float("nan"),
            "pair_error_sd": statistics.pstdev(errors) if len(errors) >= 2 else float("nan"),
            "conditional_extensions_per_pair_hit": (
                total_extensions / pair_rep_count if pair_rep_count else float("nan")
            ),
            "palm_periodic": float(palm_periodic),
            "palm_periodic_ratio": (
                float(palm_periodic / s3_decimal) if s3 else float("nan")
            ),
            "palm_finite": float(palm_finite),
            "palm_finite_ratio": (
                float(palm_finite / s3_decimal) if s3 else float("nan")
            ),
        },
        section=section,
    )


def render_report(results: Sequence[XResult], command: str) -> str:
    first = results[0]
    last = results[-1]
    d = last.diagnostics
    lines = [
        "# Direct CRT exploration of (HM)_3",
        "",
        "This report was generated by `hm3_explore.py`. It enumerates all zero",
        "triples directly in prime/residue space. For `p<q`, the pair CRT",
        "representative `a` lies in `[0,pq)` and `pq>X^2`; only `a<X^2` can",
        "extend to a triple below `X^2`, and it extends to `p<q<s` exactly when",
        "`a mod s` belongs to `Z_s`. Hence no contributing residue tuple is",
        "sampled or omitted. The m-scatter is computed afterward only for exact",
        "cross-checks and conditional statistics.",
        "",
        f"Reproduction command: `{command}`",
        "",
        "## Summary",
        "",
        md_table(
            (
                "X",
                "#p",
                "active",
                "sum Z(p)",
                "lambda",
                "canonical triples",
                "ordered S3",
                "X^2 lambda^3",
                "R3",
                "R2",
                "max K",
            ),
            (
                (
                    result.x,
                    result.prime_count,
                    result.active_count,
                    result.zero_count,
                    fmt_decimal(result.lam, 9),
                    result.canonical_triples,
                    result.s3,
                    fmt_decimal(result.expected3, 6),
                    fmt_decimal(result.r3, 9),
                    fmt_decimal(result.r2, 9),
                    result.max_k,
                )
                for result in results
            ),
        ),
        "",
        "The exact ordered counts agree with the independent identity "
        "`sum_m (K_X(m))_3=6*(number of canonical p<q<s tuples)` at every scale. "
        "For X=256,512,1024,2048, the exact ratios reproduce the three-decimal "
        "output of `scripts/p32_hm_check.py` (from which the two-decimal TeX "
        "table was prepared). X=4096 is included here as the additional "
        "feasible scale.",
        "",
        "## Cross-scale findings",
        "",
        f"- The ordered count at the largest requested scale is exactly "
        f"{last.s3}, versus `X^2 lambda_X^3 = {fmt_decimal(last.expected3, 6)}`; "
        f"thus `R3 = {fmt_decimal(last.r3, 9)}`. Across the requested scales, "
        f"R3 rises from {fmt_decimal(first.r3, 4)} to "
        f"{fmt_decimal(last.r3, 4)}, while R2 reaches "
        f"{fmt_decimal(last.r2, 4)}.",
        "",
        f"- At X={last.x}, only "
        f"{int(d['repeated_quotient_triples'])} of "
        f"{last.canonical_triples} canonical triples have a repeated quotient, "
        f"including {int(d['all_equal_quotient_triples'])} all-equal triples. "
        "The eight L/R reflection-side "
        "patterns are close in size; the largest absolute pairwise correlation "
        f"among signed normalized residues is {d['maximum_abs_signed_residue_correlation']:.4f}. "
        f"The maximum contributing reflection-orbit occupancy is "
        f"{int(d['maximum_reflection_orbit_occupancy'])}.",
        "",
        f"- The {last.canonical_triples} canonical contributions at X={last.x} "
        f"are supported on {int(d['special_m_count'])} integers: "
        f"{int(d['k_at_least_four_m_count'])} have K>=4 and all others have "
        "K=3. Their counts across the ten equal m-intervals range from "
        f"{int(d['special_decile_min'])} to {int(d['special_decile_max'])}; "
        f"{int(d['near_one_percent_triples'])} canonical triples have some "
        "residue within 1% of a prime multiple.",
        "",
        f"- On active pairs at X={last.x}, the signed discrepancy sum is "
        f"{d['pair_error_sum']:+.6f}, the mean is "
        f"{d['pair_error_mean']:+.6f}, and its population standard deviation is "
        f"{d['pair_error_sd']:.6f}. Conditional third extensions occur at rate "
        f"{d['conditional_extensions_per_pair_hit']:.6f} per pair hit, close to "
        f"`lambda_X = {fmt_decimal(last.lam, 6)}`; the detailed A(p,q)-conditional "
        "rows show no monotone growth with the pair count. The ordered "
        f"periodic pair-Palm positive excess is {d['palm_periodic']:.6f} "
        f"({d['palm_periodic_ratio']:.4f} of S3), and its exact finite-marginal "
        f"version is {d['palm_finite']:.6f} "
        f"({d['palm_finite_ratio']:.4f} of S3); both theorem sandwiches pass.",
        "",
        "These are diagnostics, not an arithmetic estimate uniform in X.",
        "",
    ]
    for result in results:
        lines.append(result.section)
    lines.extend(
        [
            "## Interpretation limited to the computation",
            "",
            "The tables test clustering and signed pair discrepancies, but they "
            "are not used as proof. In particular, bounded empirical R3, small "
            "coordinate correlations, or cancellation in the sum of E(p,q) do "
            "not provide the uniform cross-prime estimate required by (HM)_3.",
            "",
        ]
    )
    return "\n".join(lines)


def parse_args(argv: Sequence[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--xs",
        metavar="X",
        type=int,
        nargs="+",
        default=list(DEFAULT_XS),
        help="dyadic lower endpoints to analyze (default: 256 512 1024 2048 4096)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).with_name("hm3_exploration.md"),
        help="Markdown output path; use '-' for stdout",
    )
    parser.add_argument(
        "--deep-check",
        action="store_true",
        help="cross-check cleared and divided Apéry recurrences at every prime",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="suppress progress messages on stderr",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(sys.argv[1:] if argv is None else argv)
    xs = sorted(set(args.xs))
    if not xs:
        raise ValueError("at least one X is required")
    results = [analyze_x(x, args.deep_check, not args.quiet) for x in xs]
    command = "python3 problems/3.2/hm3_explore.py"
    if xs != list(DEFAULT_XS):
        command += " --xs " + " ".join(str(x) for x in xs)
    if args.deep_check:
        command += " --deep-check"
    report = render_report(results, command)
    if str(args.output) == "-":
        sys.stdout.write(report)
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(report, encoding="utf-8")
        if not args.quiet:
            print(f"wrote {args.output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (AssertionError, ValueError, OverflowError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
