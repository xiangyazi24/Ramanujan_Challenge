#!/usr/bin/env python3
"""Independent exact reference scanner for Q7308.

This deliberately does not use the pair-CRT enumerator or the forthcoming C++
scanner.  It computes exact Apéry zero sets modulo every prime p in (X,2X],
then scatters every event n == r (mod p), 0 <= n < X^2, into a sparse map.
Every map entry of load at least three is certified.

For each retained row n and hit p_i, put

    n = q_i p_i + r_i,   rho_i = min(r_i, p_i - 1 - r_i).

The two cross-zero matrices are

    Q_ij = 1[B_{q_j} == 0 (mod p_i)],
    R_ij = 1[B_{rho_j} == 0 (mod p_i)].

On the deduplicated digit pool S(n) = {q_i,rho_i}, the script computes every
minor

    Delta_i(s,t) = B_s D_t - D_s B_t (mod p_i),  s < t,

where B_0=1,B_1=5 and D_0=0,D_1=1 satisfy the cleared Apéry recurrence.
Only zero minors (projective collisions) are emitted, but every unordered pair
is tested and the exact number of tests is reported.

Classification used throughout:

  D   same canonical zero d=rho_i (the strict R_ii diagonal is a subset);
  F   own reflected zero d=p_i-1-rho_i;
  O   a genuinely off-diagonal cross-zero digit.

For state minors, s=t diagonals are suppressed by construction and counted
separately.  A zero minor with s+t=p_i-1 and both endpoints in the zero fiber
is F (automatic reflection).  Every other collision is off-diagonal, split as
OZ (both endpoints are zeros) or ON (a nonzero projective fiber).

The program uses only the Python standard library.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import DefaultDict, Dict, Iterable, List, Mapping, MutableMapping, Sequence, Tuple


XS = (128, 256, 512)

# Existing repository hm3_exploration.md values.  They are checked only after
# the independent scatter has completed; they are never used to construct rows.
REPO_CHECK = {
    256: {"prime_count": 43, "active_count": 17, "zero_count": 46,
          "canonical_triples": 10, "max_k": 3},
    512: {"prime_count": 75, "active_count": 30, "zero_count": 70,
          "canonical_triples": 25, "max_k": 3},
}


@dataclass(frozen=True)
class PrimeData:
    p: int
    roots: Tuple[int, ...]
    root_set: frozenset[int]
    B: Tuple[int, ...]
    D: Tuple[int, ...]


def poly(n: int) -> int:
    return 34 * n * n * n + 51 * n * n + 27 * n + 5


def primes_upto(n: int) -> List[int]:
    if n < 2:
        return []
    sieve = bytearray(b"\x01") * (n + 1)
    sieve[0:2] = b"\x00\x00"
    for p in range(2, math.isqrt(n) + 1):
        if sieve[p]:
            start = p * p
            sieve[start:n + 1:p] = b"\x00" * (((n - start) // p) + 1)
    return [p for p in range(2, n + 1) if sieve[p]]


def cleared_states(p: int) -> Tuple[Tuple[int, ...], Tuple[int, ...]]:
    """Return B_0..B_{p-1}, D_0..D_{p-1} modulo p.

    Both states obey Y_{m+1}=P(m)Y_m-m^6Y_{m-1}.  The Wronskian identity is
    checked at every index, so an indexing or recurrence error fails closed.
    """

    B = [0] * p
    D = [0] * p
    B[0], D[0] = 1 % p, 0
    if p > 1:
        B[1], D[1] = 5 % p, 1 % p
    for m in range(1, p - 1):
        c = poly(m) % p
        back = pow(m, 6, p)
        B[m + 1] = (c * B[m] - back * B[m - 1]) % p
        D[m + 1] = (c * D[m] - back * D[m - 1]) % p

    fact = 1
    for m in range(0, p - 1):
        lhs = (B[m] * D[m + 1] - B[m + 1] * D[m]) % p
        rhs = pow(fact, 6, p)
        if lhs != rhs:
            raise AssertionError(
                f"Wronskian failure p={p}, m={m}: {lhs} != {rhs}"
            )
        fact = fact * (m + 1) % p
    return tuple(B), tuple(D)


def divided_zero_set(p: int) -> Tuple[int, ...]:
    """Independent zero-set check using the original divided recurrence."""

    values = [0] * p
    values[0] = 1 % p
    if p > 1:
        values[1] = 5 % p
    for m in range(1, p - 1):
        numerator = (poly(m) * values[m] - pow(m, 3, p) * values[m - 1]) % p
        inv = pow(m + 1, -1, p)
        values[m + 1] = numerator * pow(inv, 3, p) % p
    return tuple(i for i, value in enumerate(values) if value == 0)


def build_prime_data(x: int) -> Tuple[List[PrimeData], List[PrimeData]]:
    all_data: List[PrimeData] = []
    for p in (q for q in primes_upto(2 * x) if q > x):
        B, D = cleared_states(p)
        roots = tuple(i for i, value in enumerate(B) if value == 0)
        reference = divided_zero_set(p)
        if roots != reference:
            raise AssertionError(
                f"cleared/divided root mismatch p={p}: {roots} != {reference}"
            )
        root_set = frozenset(roots)
        for r in roots:
            if p - 1 - r not in root_set:
                raise AssertionError(f"reflection failure p={p}, r={r}")
            if D[r] == 0:
                raise AssertionError(f"singular zero state p={p}, r={r}")
        all_data.append(PrimeData(p, roots, root_set, B, D))
    return all_data, [item for item in all_data if item.roots]


def scatter_events(active: Sequence[PrimeData], limit: int) -> Tuple[DefaultDict[int, List[Tuple[int, int]]], int]:
    """Scatter every exact (p,r) hit event directly into its row n."""

    events: DefaultDict[int, List[Tuple[int, int]]] = defaultdict(list)
    event_count = 0
    for item in active:
        for r in item.roots:
            for n in range(r, limit, item.p):
                events[n].append((item.p, r))
                event_count += 1
    for n, row in events.items():
        primes = [p for p, _ in row]
        if len(primes) != len(set(primes)):
            raise AssertionError(f"duplicate characteristic in row n={n}: {row}")
        for p, r in row:
            if n % p != r:
                raise AssertionError(f"bad scattered residue n={n}, p={p}, r={r}")
    return events, event_count


def zero_class(p: int, own_rho: int, digit: int, is_zero: bool) -> str:
    if not is_zero:
        return "."
    if digit == own_rho:
        return "D"
    if digit == p - 1 - own_rho:
        return "F"
    return "O"


def projective_label(b: int, d: int, p: int) -> object:
    """Canonical [B:D] label: B/D if D!=0, otherwise the string 'inf'."""

    if d % p == 0:
        if b % p == 0:
            raise AssertionError(f"zero projective vector modulo {p}")
        return "inf"
    return b * pow(d, -1, p) % p


def edge_types(origins: Mapping[int, frozenset[str]], s: int, t: int) -> Tuple[str, ...]:
    left, right = origins[s], origins[t]
    out: List[str] = []
    if "Q" in left and "Q" in right:
        out.append("QQ")
    if ("Q" in left and "R" in right) or ("R" in left and "Q" in right):
        out.append("QR")
    if "R" in left and "R" in right:
        out.append("RR")
    if not out:
        raise AssertionError(f"untyped digit edge {(s, t)}")
    return tuple(out)


def matrix_and_classes(
    hits: Sequence[Mapping[str, int]],
    prime_map: Mapping[int, PrimeData],
    column_key: str,
) -> Tuple[List[List[int]], List[List[str]], Counter[str], Counter[str]]:
    """Return full matrix plus entry and deduplicated-(p,digit) class counts."""

    columns = [hit[column_key] for hit in hits]
    bits: List[List[int]] = []
    classes: List[List[str]] = []
    entry_counts: Counter[str] = Counter()
    unique_counts: Counter[str] = Counter()
    for hit in hits:
        item = prime_map[hit["p"]]
        row_bits: List[int] = []
        row_classes: List[str] = []
        for digit in columns:
            is_zero = item.B[digit] == 0
            cls = zero_class(item.p, hit["rho"], digit, is_zero)
            row_bits.append(int(is_zero))
            row_classes.append(cls)
            if cls != ".":
                entry_counts[cls] += 1
        for digit in sorted(set(columns)):
            cls = zero_class(item.p, hit["rho"], digit, item.B[digit] == 0)
            if cls != ".":
                unique_counts[cls] += 1
        bits.append(row_bits)
        classes.append(row_classes)
    return bits, classes, entry_counts, unique_counts


def compact_bit_rows(matrix: Sequence[Sequence[int]]) -> List[str]:
    return ["".join(str(value) for value in row) for row in matrix]


def compact_class_rows(matrix: Sequence[Sequence[str]]) -> List[str]:
    return ["".join(row) for row in matrix]


def canonical_digest(value: object) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def certify_row(
    x: int,
    n: int,
    raw_hits: Sequence[Tuple[int, int]],
    prime_map: Mapping[int, PrimeData],
) -> Dict[str, object]:
    hits: List[Dict[str, int]] = []
    for p, scattered_r in sorted(raw_hits):
        r = n % p
        if r != scattered_r or r not in prime_map[p].root_set:
            raise AssertionError(f"row reconstruction failure X={x}, n={n}, p={p}")
        q = n // p
        rho = min(r, p - 1 - r)
        if not (0 <= q < x and 0 <= rho < x):
            raise AssertionError(f"digit bound failure X={x}, n={n}, p={p}")
        if rho not in prime_map[p].root_set:
            raise AssertionError(f"canonical root failure X={x}, n={n}, p={p}")
        hits.append({"p": p, "r": r, "q": q, "rho": rho})

    Q, Qc, q_entries, q_unique = matrix_and_classes(hits, prime_map, "q")
    R, Rc, r_entries, r_unique = matrix_and_classes(hits, prime_map, "rho")

    # The strict R diagonal is the defining hit and must always be D.
    for i in range(len(hits)):
        if R[i][i] != 1 or Rc[i][i] != "D":
            raise AssertionError(f"R diagonal failure X={x}, n={n}, i={i}")

    mutable_origins: MutableMapping[int, set[str]] = defaultdict(set)
    for hit in hits:
        mutable_origins[hit["q"]].add("Q")
        mutable_origins[hit["rho"]].add("R")
    origins: Dict[int, frozenset[str]] = {
        digit: frozenset(tags) for digit, tags in mutable_origins.items()
    }
    digits = sorted(origins)

    collisions: List[Dict[str, object]] = []
    minor_class_counts: Counter[str] = Counter()
    typed_collision_counts: Counter[str] = Counter()
    for hit in hits:
        item = prime_map[hit["p"]]
        for s, t in itertools.combinations(digits, 2):
            delta = (item.B[s] * item.D[t] - item.D[s] * item.B[t]) % item.p
            if delta != 0:
                continue
            zs = item.B[s] == 0
            zt = item.B[t] == 0
            if zs != zt:
                raise AssertionError(
                    f"mixed zero/nonzero projective collision p={item.p}, {(s, t)}"
                )
            if zs and s + t == item.p - 1:
                cls = "F"
            elif zs:
                cls = "OZ"
            else:
                cls = "ON"
            types = edge_types(origins, s, t)
            label_s = projective_label(item.B[s], item.D[s], item.p)
            label_t = projective_label(item.B[t], item.D[t], item.p)
            if label_s != label_t:
                raise AssertionError(
                    f"normalization mismatch p={item.p}, {(s, t)}: {label_s}, {label_t}"
                )
            record = {
                "p": item.p,
                "s": s,
                "t": t,
                "types": list(types),
                "class": cls,
                "label": label_s,
                "state_s": [item.B[s], item.D[s]],
                "state_t": [item.B[t], item.D[t]],
                "minor": 0,
            }
            collisions.append(record)
            minor_class_counts[cls] += 1
            typed_collision_counts.update(types)

    collision_key = lambda c: (c["p"], c["s"], c["t"], c["class"], c["types"])
    collisions.sort(key=collision_key)
    k = len(hits)
    diagonal_suppressed = k * len(digits)
    offdiagonal_tests = k * math.comb(len(digits), 2)

    core: Dict[str, object] = {
        "X": x,
        "n": n,
        "K": k,
        "hits": hits,
        "digits": digits,
        "origins": {str(d): "".join(sorted(origins[d])) for d in digits},
        "Q": compact_bit_rows(Q),
        "Qclass": compact_class_rows(Qc),
        "R": compact_bit_rows(R),
        "Rclass": compact_class_rows(Rc),
        "Q_entry_counts": dict(sorted(q_entries.items())),
        "R_entry_counts": dict(sorted(r_entries.items())),
        "Q_unique_counts": dict(sorted(q_unique.items())),
        "R_unique_counts": dict(sorted(r_unique.items())),
        "minor_diagonals_suppressed": diagonal_suppressed,
        "minor_offdiagonal_tests": offdiagonal_tests,
        "minor_collisions": collisions,
        "minor_class_counts": dict(sorted(minor_class_counts.items())),
        "minor_typed_counts": dict(sorted(typed_collision_counts.items())),
    }
    core["sha256"] = canonical_digest(core)
    return core


def add_counter(target: Counter[str], source: Mapping[str, int]) -> None:
    for key, value in source.items():
        target[key] += value


def analyze_x(x: int) -> Dict[str, object]:
    limit = x * x
    all_data, active = build_prime_data(x)
    prime_map = {item.p: item for item in all_data}
    events, event_count = scatter_events(active, limit)

    load_hist_all = Counter(len(row) for row in events.values())
    load_hist_all[0] = limit - len(events)
    retained_ns = sorted(n for n, row in events.items() if len(row) >= 3)
    rows = [certify_row(x, n, events[n], prime_map) for n in retained_ns]

    q_entry: Counter[str] = Counter()
    r_entry: Counter[str] = Counter()
    q_unique: Counter[str] = Counter()
    r_unique: Counter[str] = Counter()
    minor_class: Counter[str] = Counter()
    minor_typed: Counter[str] = Counter()
    diagonal_suppressed = 0
    offdiagonal_tests = 0
    collision_count = 0
    strict_r_diagonal = 0
    for row in rows:
        add_counter(q_entry, row["Q_entry_counts"])
        add_counter(r_entry, row["R_entry_counts"])
        add_counter(q_unique, row["Q_unique_counts"])
        add_counter(r_unique, row["R_unique_counts"])
        add_counter(minor_class, row["minor_class_counts"])
        add_counter(minor_typed, row["minor_typed_counts"])
        diagonal_suppressed += int(row["minor_diagonals_suppressed"])
        offdiagonal_tests += int(row["minor_offdiagonal_tests"])
        collision_count += len(row["minor_collisions"])
        strict_r_diagonal += int(row["K"])

    canonical_triples = sum(
        math.comb(k, 3) * count for k, count in load_hist_all.items() if k >= 3
    )
    ordered_s3 = 6 * canonical_triples
    direct_s3 = sum(
        k * (k - 1) * (k - 2) * count
        for k, count in load_hist_all.items()
    )
    if ordered_s3 != direct_s3:
        raise AssertionError(f"factorial-moment mismatch X={x}")

    summary: Dict[str, object] = {
        "X": x,
        "limit": limit,
        "prime_count": len(all_data),
        "active_count": len(active),
        "zero_count": sum(len(item.roots) for item in all_data),
        "scattered_events": event_count,
        "load_histogram": {str(k): load_hist_all[k] for k in sorted(load_hist_all)},
        "retained_rows": len(rows),
        "retained_load_histogram": {
            str(k): load_hist_all[k] for k in sorted(load_hist_all) if k >= 3
        },
        "max_k": max(load_hist_all),
        "canonical_triples": canonical_triples,
        "ordered_S3": ordered_s3,
        "strict_R_diagonal": strict_r_diagonal,
        "Q_entry_counts": dict(sorted(q_entry.items())),
        "R_entry_counts": dict(sorted(r_entry.items())),
        "Q_unique_counts": dict(sorted(q_unique.items())),
        "R_unique_counts": dict(sorted(r_unique.items())),
        "minor_diagonals_suppressed": diagonal_suppressed,
        "minor_offdiagonal_tests": offdiagonal_tests,
        "minor_collision_count": collision_count,
        "minor_class_counts": dict(sorted(minor_class.items())),
        "minor_typed_counts": dict(sorted(minor_typed.items())),
        "row_certificate_sha256": canonical_digest(rows),
    }

    expected = REPO_CHECK.get(x)
    if expected is not None:
        for key, value in expected.items():
            if summary[key] != value:
                raise AssertionError(
                    f"repository cross-check failure X={x}, {key}: "
                    f"{summary[key]} != {value}"
                )
        summary["repo_hm3_cross_check"] = "PASS"
    else:
        summary["repo_hm3_cross_check"] = "not previously tabulated"

    return {"summary": summary, "rows": rows}


def count_triplet(mapping: Mapping[str, int], labels: Sequence[str]) -> str:
    return "/".join(str(mapping.get(label, 0)) for label in labels)


def markdown_report(document: Mapping[str, object]) -> str:
    results = document["results"]
    lines: List[str] = []
    lines.append("# Q7308 independent scatter verification")
    lines.append("")
    lines.append("The run is a pure Python, standard-library reference computation.  It does not use pair CRT or any C++ output: every exact residue event is scattered directly into `n < X^2`, and only then are rows with `K_X(n) >= 3` retained.")
    lines.append("")
    lines.append("## Classification")
    lines.append("")
    lines.append("For cross-zero entries, `D` means the tested digit equals the row characteristic's canonical zero `rho_i`; `F` means its reflected mate `p_i-1-rho_i`; and `O` is genuinely off-diagonal.  Thus the strict `R_ii` diagonal is automatic, equal-digit duplicate columns remain `D`, and only `O` contributes new cross-zero information.")
    lines.append("")
    lines.append("For projective minors, all `s=t` diagonals are omitted and counted under `minor_diagonals_suppressed`.  Among tested `s<t`, `F` is an automatic zero/reflection edge (`s+t=p-1`), `OZ` is a nonreflection collision in the zero fiber, and `ON` is a collision in a nonzero projective fiber.  `OZ` and `ON` are the genuinely off-diagonal records.  Edge types are all applicable labels among `QQ`, `QR`, and `RR`; a digit belonging to both pools can give more than one type.")
    lines.append("")
    lines.append("## Exact per-X summary")
    lines.append("")
    headers = [
        "X", "#p", "active", "sum |Zp|", "events", "K>=3 rows", "K>=3 hist",
        "max K", "Q D/F/O", "R D/F/O", "minor tests", "diag suppressed",
        "minor F/OZ/ON", "row SHA-256"
    ]
    lines.append("| " + " | ".join(headers) + " |")
    lines.append("|" + "|".join("---" for _ in headers) + "|")
    for result in results:
        s = result["summary"]
        row = [
            s["X"], s["prime_count"], s["active_count"], s["zero_count"],
            s["scattered_events"], s["retained_rows"],
            json.dumps(s["retained_load_histogram"], separators=(",", ":")),
            s["max_k"],
            count_triplet(s["Q_entry_counts"], ("D", "F", "O")),
            count_triplet(s["R_entry_counts"], ("D", "F", "O")),
            s["minor_offdiagonal_tests"], s["minor_diagonals_suppressed"],
            count_triplet(s["minor_class_counts"], ("F", "OZ", "ON")),
            "`" + s["row_certificate_sha256"] + "`",
        ]
        lines.append("| " + " | ".join(str(v) for v in row) + " |")
    lines.append("")
    lines.append("Entry counts above retain all matrix columns.  The following deduplicated counts treat each `(characteristic,digit)` once, as required when equal q or rho digits occur.")
    lines.append("")
    lines.append("| X | Q unique D/F/O | R unique D/F/O | strict R diagonal | typed minor incidences QQ/QR/RR | canonical triples | ordered S3 | prior hm3 check |")
    lines.append("|---|---|---|---|---|---|---|---|")
    for result in results:
        s = result["summary"]
        row = [
            s["X"],
            count_triplet(s["Q_unique_counts"], ("D", "F", "O")),
            count_triplet(s["R_unique_counts"], ("D", "F", "O")),
            s["strict_R_diagonal"],
            count_triplet(s["minor_typed_counts"], ("QQ", "QR", "RR")),
            s["canonical_triples"], s["ordered_S3"], s["repo_hm3_cross_check"],
        ]
        lines.append("| " + " | ".join(str(v) for v in row) + " |")
    lines.append("")
    lines.append("## Row certificates")
    lines.append("")
    lines.append("Each compact JSON object is canonical and complete for the requested checks. `hits` entries are `[p,r,q,rho]`; `Q`, `R` and their class matrices are row strings; each minor collision is `[p,s,t,types,class,label,state_s,state_t]`.  The per-row SHA-256 is computed before adding the digest field.")
    lines.append("")
    for result in results:
        s = result["summary"]
        lines.append(f"### X = {s['X']}")
        lines.append("")
        if not result["rows"]:
            lines.append("No row has `K_X(n) >= 3`.")
            lines.append("")
            continue
        for cert in result["rows"]:
            compact = {
                "n": cert["n"],
                "K": cert["K"],
                "hits": [[h["p"], h["r"], h["q"], h["rho"]] for h in cert["hits"]],
                "S": cert["digits"],
                "orig": cert["origins"],
                "Q": cert["Q"],
                "Qc": cert["Qclass"],
                "R": cert["R"],
                "Rc": cert["Rclass"],
                "diag": cert["minor_diagonals_suppressed"],
                "tests": cert["minor_offdiagonal_tests"],
                "M": [
                    [m["p"], m["s"], m["t"], "+".join(m["types"]), m["class"],
                     m["label"], m["state_s"], m["state_t"]]
                    for m in cert["minor_collisions"]
                ],
                "sha256": cert["sha256"],
            }
            lines.append("`" + json.dumps(compact, sort_keys=True, separators=(",", ":")) + "`")
            lines.append("")
    lines.append("## Machine-check digest")
    lines.append("")
    lines.append("The JSON document digest (with `document_sha256` omitted while hashing) is `" + document["document_sha256"] + "`.")
    lines.append("")
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--markdown", default="q7308-output.md")
    parser.add_argument("--json", default="q7308-output.json")
    args = parser.parse_args()

    results = [analyze_x(x) for x in XS]
    document: Dict[str, object] = {
        "schema": "Q7308-independent-scatter-v1",
        "X_values": list(XS),
        "results": results,
    }
    document["document_sha256"] = canonical_digest(document)

    Path(args.json).write_text(
        json.dumps(document, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    Path(args.markdown).write_text(markdown_report(document), encoding="utf-8")

    for result in results:
        s = result["summary"]
        print(
            f"X={s['X']} rows={s['retained_rows']} hist={s['retained_load_histogram']} "
            f"minor(F/OZ/ON)={count_triplet(s['minor_class_counts'], ('F','OZ','ON'))} "
            f"sha256={s['row_certificate_sha256']}"
        )
    print(f"document_sha256={document['document_sha256']}")


if __name__ == "__main__":
    main()
