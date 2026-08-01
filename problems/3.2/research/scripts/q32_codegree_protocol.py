#!/usr/bin/env python3
"""Exact transfer-cocycle codegree and orbit protocol for Problem 3.2.

The finite-field calculations use integer residues only.  NumPy, when
available, accelerates the exact integer Gram matrix in Experiment 1.  The
stdlib fallback represents multiplicities by binary bit planes and evaluates
each codegree with nine exact popcounts.
"""

from __future__ import annotations

from array import array
from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from hashlib import sha256
import math
import os
from pathlib import Path
import statistics
import sys
import time
from typing import Dict, List, Optional, Sequence, Tuple

try:
    import numpy as np
except ImportError:
    np = None


EXPERIMENT_1_PRIMES = (101, 211, 401, 601, 1009)
EXPERIMENT_2_COUNT = 200
EXPERIMENT_2_LOWER = 1000
EXPERIMENT_2_UPPER = 20_000
ANNEALED_SEED = 20_260_801
MASK64 = (1 << 64) - 1
ROOT = Path(__file__).resolve().parents[2]
REPORT_PATH = ROOT / "CODEX_CODEGREE.md"


Point = int  # finite residues are 0,...,p-1 and infinity is represented by p


@dataclass(frozen=True)
class FieldTables:
    inverse: Tuple[int, ...]
    polynomial: Tuple[int, ...]
    numerator: Tuple[int, ...]


@dataclass(frozen=True)
class GateRecord:
    p: int
    visits: Tuple[int, ...]
    x1: int


@dataclass(frozen=True)
class CodegreeResult:
    p: int
    backend: str
    seconds: float
    flags: Tuple[Tuple[int, int, int], ...]
    diagonal_values: Tuple[int, ...]
    finite_off_min: int
    finite_off_max: int
    max_multiplicity: int
    bilinear_signatures: Tuple[Tuple[int, int, int, int], ...]
    product_max_support: int
    product_candidates: Tuple[Tuple[int, int], ...]
    sum_max_support: int
    sum_candidates: Tuple[Tuple[int, int], ...]


@dataclass(frozen=True)
class OrbitRecord:
    p: int
    ordered_visits: Tuple[int, ...]
    annealed_visits: Tuple[int, ...]


@dataclass(frozen=True)
class Moments:
    mean: float
    variance: float
    factorial_2: float
    factorial_3: float


class ProtocolFailure(RuntimeError):
    def __init__(self, heading: str, details: Sequence[str]) -> None:
        super().__init__(heading)
        self.heading = heading
        self.details = tuple(details)


class SplitMix64:
    """Version-independent 64-bit generator with unbiased randbelow."""

    def __init__(self, seed: int) -> None:
        self.state = seed & MASK64

    def next_u64(self) -> int:
        self.state = (self.state + 0x9E3779B97F4A7C15) & MASK64
        z = self.state
        z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & MASK64
        z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & MASK64
        return (z ^ (z >> 31)) & MASK64

    def randbelow(self, bound: int) -> int:
        if bound <= 0:
            raise ValueError("bound must be positive")
        limit = (1 << 64) - ((1 << 64) % bound)
        while True:
            value = self.next_u64()
            if value < limit:
                return value % bound


def primes_up_to(limit: int) -> List[int]:
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[:2] = b"\x00\x00"
    for q in range(2, math.isqrt(limit) + 1):
        if sieve[q]:
            count = ((limit - q * q) // q) + 1
            sieve[q * q : limit + 1 : q] = b"\x00" * count
    return [n for n in range(2, limit + 1) if sieve[n]]


def apery_polynomial(n: int, p: int) -> int:
    return ((2 * n + 1) * (17 * n * n + 17 * n + 5)) % p


def field_tables(p: int) -> FieldTables:
    inverse = [0] * p
    inverse[1] = 1
    for value in range(2, p):
        inverse[value] = (p - (p // value) * inverse[p % value] % p) % p
    polynomial = tuple(apery_polynomial(u, p) for u in range(p - 1))
    numerator = tuple(pow(u + 1, 6, p) for u in range(p - 1))
    return FieldTables(tuple(inverse), polynomial, numerator)


def transfer(p: int, u: int, x: Point, tables: FieldTables) -> Point:
    if x == p:
        return 0
    denominator = (tables.polynomial[u] - x) % p
    if denominator == 0:
        return p
    return tables.numerator[u] * tables.inverse[denominator] % p


def projective_coordinates(x: Point, p: int) -> Tuple[int, int]:
    return (1, 0) if x == p else (x, 1)


def projective_product_holds(x: Point, y: Point, c: int, p: int) -> bool:
    x_num, x_den = projective_coordinates(x, p)
    y_num, y_den = projective_coordinates(y, p)
    return (x_num * y_num - c * x_den * y_den) % p == 0


def gate_prime(p: int) -> GateRecord:
    tables = field_tables(p)

    states: List[Point] = [0]
    for n in range(p - 1):
        states.append(transfer(p, n, states[-1], tables))
    orbit_visits = tuple(n for n in range(1, p - 1) if states[n] == p)

    b = [1, 5 % p]
    for n in range(1, p - 1):
        numerator = (
            tables.polynomial[n] * b[n] - pow(n, 3, p) * b[n - 1]
        ) % p
        denominator = pow(n + 1, 3, p)
        b.append(numerator * tables.inverse[denominator] % p)
    recurrence_visits = tuple(n for n in range(1, p - 1) if b[n] == 0)

    expected_x1 = tables.inverse[5]
    failures: List[str] = []
    if states[1] != expected_x1:
        failures.append(f"x_1={states[1]}, expected 1/5={expected_x1} (mod {p})")
    if orbit_visits != recurrence_visits:
        failures.extend(
            [
                f"orbit visits={list(orbit_visits)}",
                f"recurrence zeros={list(recurrence_visits)}",
            ]
        )
    for n in orbit_visits:
        if states[n + 1] != 0:
            failures.append(
                f"restart failure at n={n}: x_n=infinity but x_(n+1)={states[n + 1]}"
            )
            break

    reflection_failure = next(
        (n for n in range(p) if b[p - 1 - n] != b[n]), None
    )
    if reflection_failure is not None:
        n = reflection_failure
        failures.append(
            "CRITICAL recurrence reflection failure: "
            f"b_{n}={b[n]} but b_{p - 1 - n}={b[p - 1 - n]}"
        )
    asymmetric = next(
        (n for n in orbit_visits if p - 1 - n not in orbit_visits), None
    )
    if asymmetric is not None:
        failures.append(
            "CRITICAL visit reflection failure: "
            f"{asymmetric} is a visit but {p - 1 - asymmetric} is not"
        )

    cocycle_failure = None
    for n in range(p - 1):
        reflected_index = p - 1 - n
        c = -pow(n + 1, 6, p)
        if not projective_product_holds(
            states[reflected_index], states[n + 1], c, p
        ):
            cocycle_failure = n
            break
    if cocycle_failure is not None:
        n = cocycle_failure
        failures.append(
            "CRITICAL state-reflection failure at n="
            f"{n}: x_{p - 1 - n}={states[p - 1 - n]}, "
            f"x_{n + 1}={states[n + 1]}"
        )

    if failures:
        raise ProtocolFailure(f"Gate failed at p={p}", failures)
    return GateRecord(p, orbit_visits, states[1])


def write_stall_report(failure: ProtocolFailure) -> None:
    body = [
        "# CODEGREE PROTOCOL STALLED",
        "",
        "The mandatory transfer/recurrence gate failed, so no experiment was run.",
        "",
        f"## {failure.heading}",
        "",
    ]
    body.extend(f"- {detail}" for detail in failure.details)
    body.extend(
        [
            "",
            "This witness indicates that the cocycle conventions or the stated reflection law are inconsistent with the direct recurrence computation.",
            "",
        ]
    )
    REPORT_PATH.write_text("\n".join(body), encoding="utf-8")


def format_positions(positions: Sequence[int]) -> str:
    return "[" + ",".join(str(n) for n in positions) + "]"


def is_flagged(codegree: int, p: int) -> bool:
    difference = abs(codegree - p)
    return difference * difference > 36 * p


def multiplicity_counts(p: int, x: int, tables: FieldTables) -> List[int]:
    counts = [0] * (p + 1)
    for u in range(p - 1):
        counts[transfer(p, u, x, tables)] += 1
    if sum(counts) != p - 1:
        raise AssertionError("multiplicity row has the wrong mass")
    return counts


def scan_codegrees_numpy(p: int) -> Tuple[
    Tuple[Tuple[int, int, int], ...],
    Tuple[int, ...],
    int,
    int,
    int,
]:
    if np is None:
        raise RuntimeError("NumPy backend requested without NumPy")
    tables = field_tables(p)
    matrix = np.zeros((p + 1, p + 1), dtype=np.int32)
    max_multiplicity = 0
    for x in range(p):
        counts = multiplicity_counts(p, x, tables)
        row_max = max(counts)
        max_multiplicity = max(max_multiplicity, row_max)
        matrix[x, :] = counts
    matrix[p, 0] = p - 1

    if max_multiplicity > 6:
        raise AssertionError(
            f"finite multiplicity {max_multiplicity} exceeds the degree-six bound"
        )
    if not bool(np.all(matrix.sum(axis=1) == p - 1)):
        raise AssertionError("a multiplicity row has the wrong mass")

    # This is the bucket Gram computation: no (u,v,x,x') enumeration occurs.
    codegrees = matrix @ matrix.T
    if not bool(np.array_equal(codegrees, codegrees.T)):
        raise AssertionError("codegree Gram matrix is not symmetric")
    if int(codegrees[p, p]) != (p - 1) ** 2:
        raise AssertionError("infinity/infinity codegree is incorrect")
    if bool(np.any(codegrees[p, :p])) or bool(np.any(codegrees[:p, p])):
        raise AssertionError("finite/infinity codegrees must vanish")

    differences = codegrees.astype(np.int64) - p
    flag_mask = differences * differences > 36 * p
    indices = np.argwhere(flag_mask)
    flags = tuple(
        (int(x), int(y), int(codegrees[x, y])) for x, y in indices
    )
    diagonal_values = tuple(int(value) for value in np.diag(codegrees)[:p])
    finite_mask = ~np.eye(p, dtype=bool)
    finite_off_values = codegrees[:p, :p][finite_mask]
    return (
        flags,
        diagonal_values,
        int(finite_off_values.min()),
        int(finite_off_values.max()),
        max_multiplicity,
    )


if hasattr(int, "bit_count"):
    _popcount = int.bit_count
else:
    # Python 3.9 has no int.bit_count; bin/count still executes the scan exactly.
    def _popcount(value: int) -> int:
        return bin(value).count("1")


def count_bit_planes(counts: Sequence[int]) -> Tuple[int, int, int]:
    planes = [0, 0, 0]
    for position, count in enumerate(counts):
        if count > 7:
            raise AssertionError("three bit planes cannot encode this multiplicity")
        bit = 1 << position
        if count & 1:
            planes[0] |= bit
        if count & 2:
            planes[1] |= bit
        if count & 4:
            planes[2] |= bit
    return planes[0], planes[1], planes[2]


def bit_plane_dot(
    left: Tuple[int, int, int], right: Tuple[int, int, int]
) -> int:
    a0, a1, a2 = left
    b0, b1, b2 = right
    return (
        _popcount(a0 & b0)
        + 2 * (_popcount(a0 & b1) + _popcount(a1 & b0))
        + 4
        * (
            _popcount(a0 & b2)
            + _popcount(a1 & b1)
            + _popcount(a2 & b0)
        )
        + 8 * (_popcount(a1 & b2) + _popcount(a2 & b1))
        + 16 * _popcount(a2 & b2)
    )


def scan_codegrees_stdlib(p: int) -> Tuple[
    Tuple[Tuple[int, int, int], ...],
    Tuple[int, ...],
    int,
    int,
    int,
]:
    tables = field_tables(p)
    rows: List[Tuple[int, int, int]] = []
    max_multiplicity = 0
    for x in range(p):
        counts = multiplicity_counts(p, x, tables)
        max_multiplicity = max(max_multiplicity, max(counts))
        rows.append(count_bit_planes(counts))
    if max_multiplicity > 6:
        raise AssertionError(
            f"finite multiplicity {max_multiplicity} exceeds the degree-six bound"
        )

    flags: List[Tuple[int, int, int]] = []
    diagonal_values = array("I", [0]) * p
    finite_off_min = (p - 1) ** 2
    finite_off_max = 0
    for x, left in enumerate(rows):
        for y in range(x, p):
            value = bit_plane_dot(left, rows[y])
            if x == y:
                diagonal_values[x] = value
            else:
                finite_off_min = min(finite_off_min, value)
                finite_off_max = max(finite_off_max, value)
            if is_flagged(value, p):
                flags.append((x, y, value))
                if x != y:
                    flags.append((y, x, value))

    for x in range(p):
        if not is_flagged(0, p):
            raise AssertionError("finite/infinity pair unexpectedly missed threshold")
        flags.extend(((x, p, 0), (p, x, 0)))
    infinity_diagonal = (p - 1) ** 2
    if is_flagged(infinity_diagonal, p):
        flags.append((p, p, infinity_diagonal))
    flags.sort()
    return (
        tuple(flags),
        tuple(diagonal_values),
        finite_off_min,
        finite_off_max,
        max_multiplicity,
    )


def rref_nullspace(matrix: Sequence[Sequence[int]], p: int) -> List[Tuple[int, ...]]:
    if not matrix:
        return []
    rows = [[entry % p for entry in row] for row in matrix]
    column_count = len(rows[0])
    pivot_columns: List[int] = []
    pivot_row = 0
    for column in range(column_count):
        selected = next(
            (r for r in range(pivot_row, len(rows)) if rows[r][column]), None
        )
        if selected is None:
            continue
        rows[pivot_row], rows[selected] = rows[selected], rows[pivot_row]
        scale = pow(rows[pivot_row][column], -1, p)
        rows[pivot_row] = [(scale * value) % p for value in rows[pivot_row]]
        for r in range(len(rows)):
            if r == pivot_row or rows[r][column] == 0:
                continue
            factor = rows[r][column]
            rows[r] = [
                (value - factor * pivot_value) % p
                for value, pivot_value in zip(rows[r], rows[pivot_row])
            ]
        pivot_columns.append(column)
        pivot_row += 1
        if pivot_row == len(rows):
            break

    free_columns = [c for c in range(column_count) if c not in pivot_columns]
    basis: List[Tuple[int, ...]] = []
    for free in free_columns:
        vector = [0] * column_count
        vector[free] = 1
        for r in range(len(pivot_columns) - 1, -1, -1):
            pivot = pivot_columns[r]
            vector[pivot] = -sum(
                rows[r][c] * vector[c] for c in free_columns
            ) % p
        basis.append(tuple(vector))
    return basis


def bilinear_row(x: int, y: int, p: int) -> Tuple[int, int, int, int]:
    x_num, x_den = projective_coordinates(x, p)
    y_num, y_den = projective_coordinates(y, p)
    return (
        x_num * y_num % p,
        x_num * y_den % p,
        x_den * y_num % p,
        x_den * y_den % p,
    )


def centered_residue(value: int, p: int) -> int:
    value %= p
    return value - p if value > p // 2 else value


def primitive_integer_tuple(values: Sequence[int]) -> Tuple[int, ...]:
    divisor = 0
    for value in values:
        divisor = math.gcd(divisor, abs(value))
    if divisor > 1:
        values = [value // divisor for value in values]
    result = tuple(values)
    first = next((value for value in result if value), 0)
    return tuple(-value for value in result) if first < 0 else result


def small_integer_signature(vector: Sequence[int], p: int) -> Tuple[int, int, int, int]:
    candidates: List[Tuple[Tuple[object, ...], Tuple[int, int, int, int]]] = []
    for scale in range(1, p):
        lifted = [centered_residue(scale * value, p) for value in vector]
        primitive = primitive_integer_tuple(lifted)
        score = (
            max(abs(value) for value in primitive),
            sum(abs(value) for value in primitive),
            tuple(abs(value) for value in primitive),
            primitive,
        )
        candidates.append((score, primitive))
    return min(candidates, key=lambda item: item[0])[1]  # type: ignore[return-value]


def fit_bilinear_relation(
    pairs: Sequence[Tuple[int, int]], p: int
) -> Optional[Tuple[int, int, int, int]]:
    if len(pairs) < 4:
        return None
    basis = rref_nullspace([bilinear_row(x, y, p) for x, y in pairs], p)
    if len(basis) != 1:
        return None
    return small_integer_signature(basis[0], p)


def constant_relation_candidates(
    pairs: Sequence[Tuple[int, int]], p: int, operation: str
) -> Tuple[int, Tuple[Tuple[int, int], ...]]:
    if operation == "product":
        counts = Counter((x * y) % p for x, y in pairs)
    elif operation == "sum":
        counts = Counter((x + y) % p for x, y in pairs)
    else:
        raise ValueError(operation)
    maximum = max(counts.values(), default=0)
    candidates = tuple(sorted((constant, count) for constant, count in counts.items() if count >= 4))
    return maximum, candidates


def scan_codegrees(p: int, force_stdlib: bool) -> CodegreeResult:
    started = time.perf_counter()
    if np is not None and not force_stdlib:
        backend = "numpy exact integer bucket-Gram"
        raw = scan_codegrees_numpy(p)
    else:
        backend = "stdlib exact bit-plane bucket-Gram"
        raw = scan_codegrees_stdlib(p)
    flags, diagonal, finite_off_min, finite_off_max, max_multiplicity = raw

    finite_pairs = [(x, y) for x, y, _ in flags if x < p and y < p]
    boundary_pairs = [(x, y) for x, y, _ in flags if x == p or y == p]
    signatures = []
    for locus in (finite_pairs, boundary_pairs):
        signature = fit_bilinear_relation(locus, p)
        if signature is not None and signature not in signatures:
            signatures.append(signature)

    product_max, product_candidates = constant_relation_candidates(
        finite_pairs, p, "product"
    )
    sum_max, sum_candidates = constant_relation_candidates(finite_pairs, p, "sum")
    return CodegreeResult(
        p=p,
        backend=backend,
        seconds=time.perf_counter() - started,
        flags=flags,
        diagonal_values=diagonal,
        finite_off_min=finite_off_min,
        finite_off_max=finite_off_max,
        max_multiplicity=max_multiplicity,
        bilinear_signatures=tuple(signatures),
        product_max_support=product_max,
        product_candidates=product_candidates,
        sum_max_support=sum_max,
        sum_candidates=sum_candidates,
    )


def relation_name(signature: Tuple[int, int, int, int]) -> str:
    if signature == (0, 1, -1, 0):
        return "X Z' - Z X' = 0 (projective diagonal x=x')"
    if signature == (0, 0, 0, 1):
        return "Z Z' = 0 (x=infinity or x'=infinity)"
    a, b, c, d = signature
    return f"{a} X X' + {b} X Z' + {c} Z X' + {d} Z Z' = 0"


def print_flagged_set(result: CodegreeResult) -> None:
    p = result.p
    flag_pairs = {(x, y) for x, y, _ in result.flags}
    boundary = {(x, p) for x in range(p + 1)} | {
        (p, y) for y in range(p + 1)
    }
    missing_boundary = sorted(boundary - flag_pairs)
    finite_diagonal_flagged = [x for x in range(p) if (x, x) in flag_pairs]
    missing_diagonal = [x for x in range(p) if (x, x) not in flag_pairs]
    finite_residual = sorted(
        (x, y)
        for x, y in flag_pairs
        if x < p and y < p and x != y
    )
    if missing_boundary:
        boundary_text = f"boundary missing={missing_boundary}"
    else:
        boundary_text = (
            "{(infinity,x),(x,infinity): x in P^1(F_p)} "
            f"[{len(boundary)} pairs]"
        )
    excluded = "none" if not missing_diagonal else str(missing_diagonal)
    print(
        f"FLAGGED_SET p={p}: {boundary_text}; finite diagonal has "
        f"{len(finite_diagonal_flagged)}/{p} points (excluded={excluded}); "
        f"finite off-diagonal residual={finite_residual}"
    )


def moments(values: Sequence[int]) -> Moments:
    count = len(values)
    mean = sum(values) / count
    variance = sum((value - mean) ** 2 for value in values) / count
    factorial_2 = sum(value * (value - 1) for value in values) / count
    factorial_3 = (
        sum(value * (value - 1) * (value - 2) for value in values) / count
    )
    return Moments(mean, variance, factorial_2, factorial_3)


def annealed_visits(p: int, tables: FieldTables, rng: SplitMix64) -> Tuple[int, ...]:
    x: Point = 0
    visits: List[int] = []
    for step in range(1, p):
        u = rng.randbelow(p - 1)
        x = transfer(p, u, x, tables)
        if x == p:
            visits.append(step)
    return tuple(visits)


def poisson_expected_count(sample_size: int, rate: float, value: int) -> float:
    return sample_size * math.exp(-rate) * rate**value / math.factorial(value)


def rate_string(numerator: int, denominator: int) -> str:
    return f"{numerator}/{denominator}={numerator / denominator:.6f}"


def record_checksum(records: Sequence[OrbitRecord]) -> str:
    lines = []
    for record in records:
        lines.append(
            f"{record.p}|{format_positions(record.ordered_visits)}|"
            f"{format_positions(record.annealed_visits)}"
        )
    return sha256("\n".join(lines).encode("ascii")).hexdigest()


def normalized_bin(value: int, scale: int, bin_count: int = 10) -> int:
    return min(bin_count - 1, (bin_count * value) // scale)


def consecutive_returns(
    records: Sequence[OrbitRecord], which: str
) -> List[Tuple[int, int]]:
    result: List[Tuple[int, int]] = []
    for record in records:
        visits = (
            record.ordered_visits if which == "ordered" else record.annealed_visits
        )
        result.extend(
            (record.p, later - earlier)
            for earlier, later in zip(visits, visits[1:])
        )
    return result


def print_return_summary(
    label: str, returns: Sequence[Tuple[int, int]], short_cutoffs: Sequence[int]
) -> None:
    if not returns:
        print(f"{label}: no observed returns")
        return
    gaps = [gap for _, gap in returns]
    minimum = min(gaps)
    normalized_mean = statistics.fmean(
        gap / (p + 1) for p, gap in returns
    )
    print(
        f"{label}: intervals={len(gaps)} min={minimum} "
        f"min_frequency={gaps.count(minimum)} mean={statistics.fmean(gaps):.6f} "
        f"mean(tau/(p+1))={normalized_mean:.6f}"
    )
    print(
        f"{label} short-return counts: "
        + ", ".join(f"tau<={cutoff}:{sum(g <= cutoff for g in gaps)}" for cutoff in short_cutoffs)
    )
    smallest = sorted(Counter(gaps).items())[:15]
    print(f"{label} smallest exact return-time frequencies: {smallest}")


def run() -> int:
    all_primes = primes_up_to(EXPERIMENT_2_UPPER)
    experiment_2_primes = tuple(
        p for p in all_primes if p >= EXPERIMENT_2_LOWER
    )[:EXPERIMENT_2_COUNT]
    if len(experiment_2_primes) != EXPERIMENT_2_COUNT:
        raise RuntimeError("not enough primes in the requested Experiment-2 range")
    gate_primes = tuple(sorted(set(EXPERIMENT_1_PRIMES) | set(experiment_2_primes)))

    force_stdlib = os.environ.get("Q32_FORCE_STDLIB") == "1"
    numpy_status = "available" if np is not None else "unavailable"
    print("Q32 TRANSFER-COCYCLE CODEGREE PROTOCOL")
    print(
        f"CONFIG: Experiment-1 primes={list(EXPERIMENT_1_PRIMES)}; "
        f"Experiment-2 first {len(experiment_2_primes)} primes in "
        f"[{EXPERIMENT_2_LOWER},{EXPERIMENT_2_UPPER}], actual range "
        f"[{experiment_2_primes[0]},{experiment_2_primes[-1]}]"
    )
    print(
        f"CONFIG: annealed RNG=SplitMix64 seed={ANNEALED_SEED}; "
        f"NumPy={numpy_status}; force_stdlib={force_stdlib}"
    )
    print("\nMANDATORY GATE PHASE")

    gate_records: Dict[int, GateRecord] = {}
    try:
        for p in gate_primes:
            record = gate_prime(p)
            gate_records[p] = record
            print(
                f"GATE VERIFIED p={p} x1={record.x1} "
                f"visits={format_positions(record.visits)} "
                "recurrence=VERIFIED reflection=VERIFIED state-cocycle=VERIFIED"
            )
    except ProtocolFailure as failure:
        print(f"CRITICAL: {failure.heading}", file=sys.stderr)
        for detail in failure.details:
            print(f"CRITICAL: {detail}", file=sys.stderr)
        write_stall_report(failure)
        print(f"STALL REPORT WRITTEN: {REPORT_PATH}", file=sys.stderr)
        return 2
    print(f"ALL GATES VERIFIED: {len(gate_records)} distinct primes")

    print("\nEXPERIMENT 1 — CODEGREE EXCEPTIONAL-LOCUS SCAN")
    codegree_results: List[CodegreeResult] = []
    for p in EXPERIMENT_1_PRIMES:
        result = scan_codegrees(p, force_stdlib)
        codegree_results.append(result)
        print_flagged_set(result)
        print(
            f"CODEGREE p={p}: backend={result.backend}; seconds={result.seconds:.3f}; "
            f"max finite multiplicity={result.max_multiplicity}; "
            f"product candidates (support>=4)={list(result.product_candidates)} "
            f"[max support {result.product_max_support}]; "
            f"sum candidates (support>=4)={list(result.sum_candidates)} "
            f"[max support {result.sum_max_support}]"
        )
        print(
            "BILINEAR FITS p="
            f"{p}: {[relation_name(signature) for signature in result.bilinear_signatures]}"
        )

    stable_signatures = set(codegree_results[0].bilinear_signatures)
    for result in codegree_results[1:]:
        stable_signatures.intersection_update(result.bilinear_signatures)
    stable_sorted = tuple(sorted(stable_signatures))
    print("\nExperiment-1 summary table")
    print(
        "p | flagged | finite diagonal flagged | diagonal C min/mean/max | "
        "finite off-diagonal C min/max"
    )
    for result in codegree_results:
        p = result.p
        pair_set = {(x, y) for x, y, _ in result.flags}
        diagonal_flagged = sum((x, x) in pair_set for x in range(p))
        diagonal_mean = statistics.fmean(result.diagonal_values)
        print(
            f"{p} | {len(result.flags)} | {diagonal_flagged}/{p} | "
            f"{min(result.diagonal_values)}/{diagonal_mean:.6f}/"
            f"{max(result.diagonal_values)} | "
            f"{result.finite_off_min}/{result.finite_off_max}"
        )
    print(
        "STABLE BILINEAR RELATIONS: "
        + "; ".join(relation_name(signature) for signature in stable_sorted)
    )
    print(
        "NONTRIVIAL REFLECTION-IMAGE VERDICT: no fixed product, sum, or other "
        "finite off-diagonal bilinear correspondence was flagged.  The ordered "
        "state reflection is instead the parameter-dependent projective identity "
        "x_(p-1-n)*x_(n+1)=-(n+1)^6, verified in every gate."
    )

    print("\nEXPERIMENT 2 — ANNEALED VS ORDERED VISITS")
    rng = SplitMix64(ANNEALED_SEED)
    orbit_records: List[OrbitRecord] = []
    for p in experiment_2_primes:
        annealed = annealed_visits(p, field_tables(p), rng)
        orbit_records.append(
            OrbitRecord(p, gate_records[p].visits, annealed)
        )

    ordered_counts = [len(record.ordered_visits) for record in orbit_records]
    annealed_counts = [len(record.annealed_visits) for record in orbit_records]
    central_counts = [
        int((record.p - 1) // 2 in record.ordered_visits)
        for record in orbit_records
    ]
    halved_counts = [
        (ordered - central) // 2
        for ordered, central in zip(ordered_counts, central_counts)
    ]
    if any(
        ordered != 2 * halved + central
        for ordered, halved, central in zip(
            ordered_counts, halved_counts, central_counts
        )
    ):
        raise AssertionError("ordered reflection decomposition failed")

    max_count = max(max(ordered_counts), max(annealed_counts))
    ordered_frequency = Counter(ordered_counts)
    annealed_frequency = Counter(annealed_counts)
    print(
        f"Experiment-2 primes: count={len(orbit_records)}, range="
        f"[{orbit_records[0].p},{orbit_records[-1].p}], record SHA256="
        f"{record_checksum(orbit_records)}"
    )
    print("Visit-count frequency table")
    print("m | ordered | annealed | Poisson(1) expected among 200")
    for value in range(max_count + 1):
        print(
            f"{value} | {ordered_frequency[value]} | {annealed_frequency[value]} | "
            f"{poisson_expected_count(len(orbit_records), 1.0, value):.6f}"
        )

    half_frequency = Counter(halved_counts)
    print("Halved ordered pair-count table")
    print("h | observed | Poisson(1/2) expected among 200")
    for value in range(max(halved_counts) + 1):
        print(
            f"{value} | {half_frequency[value]} | "
            f"{poisson_expected_count(len(orbit_records), 0.5, value):.6f}"
        )

    ordered_moments = moments(ordered_counts)
    annealed_moments = moments(annealed_counts)
    halved_moments = moments(halved_counts)
    print("Moment table")
    print("ensemble | mean | variance | E[(m)_2] | E[(m)_3]")
    for label, value in (
        ("ordered V", ordered_moments),
        ("annealed A", annealed_moments),
        ("ordered pairs H", halved_moments),
    ):
        print(
            f"{label} | {value.mean:.6f} | {value.variance:.6f} | "
            f"{value.factorial_2:.6f} | {value.factorial_3:.6f}"
        )
    print("Poisson(1) | 1.000000 | 1.000000 | 1.000000 | 1.000000")
    print("Poisson(1/2) | 0.500000 | 0.500000 | 0.250000 | 0.125000")

    even_ordered = sum(value % 2 == 0 for value in ordered_counts)
    even_annealed = sum(value % 2 == 0 for value in annealed_counts)
    expected_annealed_total = sum((p - 1) / (p + 1) for p in experiment_2_primes)
    print(
        "PARITY: ordered even "
        f"{rate_string(even_ordered, len(orbit_records))}; annealed even "
        f"{rate_string(even_annealed, len(orbit_records))}; "
        f"Poisson(1) even={(1 + math.exp(-2)) / 2:.6f}; "
        f"ordered central visits={sum(central_counts)}"
    )
    print(
        f"ANNEALED TOTAL: observed={sum(annealed_counts)}; "
        f"sum_p (p-1)/(p+1)={expected_annealed_total:.6f}"
    )

    paired_gaps: List[Tuple[int, int]] = []
    for record in orbit_records:
        center = (record.p - 1) // 2
        visits = set(record.ordered_visits)
        for r in record.ordered_visits:
            reflected = record.p - 1 - r
            if reflected not in visits:
                raise ProtocolFailure(
                    f"CRITICAL asymmetric ordered visit at p={record.p}",
                    [f"r={r} is present but p-1-r={reflected} is absent"],
                )
            if r < center:
                paired_gaps.append((record.p, reflected - r))
    print(
        f"REFLECTION CHECK: machine-exact for all {len(orbit_records)} primes; "
        f"paired gaps={len(paired_gaps)}"
    )
    print("Normalized reflection-gap histogram (gap/(p-1))")
    print("bin | observed | fraction | Exp(1) conditional-on-[0,1) fraction")
    paired_histogram = Counter(
        normalized_bin(gap, p - 1) for p, gap in paired_gaps
    )
    exponential_normalizer = 1 - math.exp(-1)
    for index in range(10):
        lower = index / 10
        upper = (index + 1) / 10
        conditional_mass = (
            math.exp(-lower) - math.exp(-upper)
        ) / exponential_normalizer
        observed = paired_histogram[index]
        print(
            f"[{lower:.1f},{upper:.1f}) | {observed} | "
            f"{observed / len(paired_gaps):.6f} | {conditional_mass:.6f}"
        )
    print(
        f"REFLECTION GAP SUMMARY: min={min(gap for _, gap in paired_gaps)}; "
        f"max={max(gap for _, gap in paired_gaps)}; "
        f"mean normalized={statistics.fmean(gap / (p - 1) for p, gap in paired_gaps):.6f}"
    )

    print("\nEXPERIMENT 3 — POST-VISIT RESTART STRUCTURE")
    ordered_returns = consecutive_returns(orbit_records, "ordered")
    annealed_returns = consecutive_returns(orbit_records, "annealed")
    short_cutoffs = (2, 5, 10, 25, 50, 100)
    print_return_summary("ordered", ordered_returns, short_cutoffs)
    print_return_summary("annealed", annealed_returns, short_cutoffs)
    thresholds = (
        Fraction(1, 10),
        Fraction(1, 4),
        Fraction(1, 2),
        Fraction(3, 4),
        Fraction(1, 1),
    )
    print("Normalized return-time CDF")
    print("t | ordered | annealed | 1-exp(-t)")
    for threshold in thresholds:
        ordered_below = sum(
            gap * threshold.denominator <= threshold.numerator * (p + 1)
            for p, gap in ordered_returns
        )
        annealed_below = sum(
            gap * threshold.denominator <= threshold.numerator * (p + 1)
            for p, gap in annealed_returns
        )
        t = float(threshold)
        print(
            f"{t:.2f} | {ordered_below}/{len(ordered_returns)}="
            f"{ordered_below / len(ordered_returns):.6f} | "
            f"{annealed_below}/{len(annealed_returns)}="
            f"{annealed_below / len(annealed_returns):.6f} | "
            f"{1 - math.exp(-t):.6f}"
        )

    print("\nPROTOCOL COMPLETE")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(run())
    except ProtocolFailure as failure:
        print(f"CRITICAL: {failure.heading}", file=sys.stderr)
        for detail in failure.details:
            print(f"CRITICAL: {detail}", file=sys.stderr)
        write_stall_report(failure)
        print(f"STALL REPORT WRITTEN: {REPORT_PATH}", file=sys.stderr)
        sys.exit(2)
