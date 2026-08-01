#!/usr/bin/env -S sage -python
"""Exact mirror-quotient Stepanov syzygy test from Q6523 section 2.

Run with Sage's Python::

    sage -python CRON_lowrank_alg.py

The main map is

    F_p[u,y]_{deg_u<=A, deg_y<=B_eff}
      -> direct_sum_{1<=h<=H} F_p[u]/(G_h(u)^M),

where ``B_eff=min(B,H-1)`` and ``G_h`` is the squarefree mirror quotient of
the Apéry gap continuant.  All ranks are exact finite-field ranks.  The
script runs both the Q6523 A-window and the largest pressure window with
Q<T, together with the three controls requested there.
"""

from __future__ import annotations

import argparse
import json
import random
import sys
import time
from dataclasses import dataclass
from typing import Iterable

try:
    from sage.all import GF, PolynomialRing, binomial, matrix, set_random_seed
except ImportError as exc:  # pragma: no cover - exercised only outside Sage
    raise SystemExit(
        "CRON_lowrank_alg.py requires SageMath; run "
        "`sage -python CRON_lowrank_alg.py`."
    ) from exc


# Q6523 section 2.7.  A_pressure is determined below by its defining rule:
# it is the largest A for which Q=(A+1)(B_eff+1)<T.
CASES = (
    (1009, 12, 1, 21, 3),
    (2003, 20, 1, 41, 5),
    (5003, 16, 3, 94, 4),
    (9001, 32, 2, 146, 8),
)

CONTROL_KINDS = ("degree_random", "mirror_random", "clock_shuffle")


@dataclass(frozen=True)
class Case:
    p: int
    H: int
    M: int
    A: int
    B: int


def expected_delta(h: int) -> int:
    """The theoretical degree of the mirror quotient G_h."""
    return (3 * (h - 1) - (1 if h % 2 == 0 else 0)) // 2


def apery_P(t):
    return 34 * t**3 + 51 * t**2 + 27 * t + 5


def gap_polynomials(F, H: int):
    """Return N[0],...,N[H] using the verified h-direction recurrence."""
    Rx = PolynomialRing(F, "x")
    x = Rx.gen()
    N = [Rx.zero(), Rx.one()]
    for h in range(1, H):
        N.append(apery_P(x + h) * N[h] - (x + h) ** 6 * N[h - 1])
    return Rx, N


def squarefree_monic(f):
    if f.is_zero():
        raise ValueError("zero gap polynomial")
    g = f.gcd(f.derivative())
    sf = f // g
    return sf.monic()


def mirror_quotient(Nh, h: int, F, Rz, Ru):
    """Remove the forced even-h mirror root and descend z^2 to u."""
    z = Rz.gen()
    u = Ru.gen()
    inv2 = F(2) ** (-1)
    fz = sum(
        (F(Nh[i]) * ((z - F(h + 1)) * inv2) ** i for i in range(Nh.degree() + 1)),
        Rz.zero(),
    )
    if h % 2 == 0:
        if fz[0] != 0:
            raise AssertionError(("missing forced mirror factor", h))
        fz = fz // z
    for e in range(1, fz.degree() + 1, 2):
        if fz[e] != 0:
            raise AssertionError(("mirror parity failure", h, e))
    gu = sum(
        (Ru(fz[2 * j]) * u**j for j in range(fz.degree() // 2 + 1)),
        Ru.zero(),
    )
    return squarefree_monic(gu)


def real_mirror_quotients(p: int, H: int):
    F = GF(p)
    _, N = gap_polynomials(F, H)
    Rz = PolynomialRing(F, "z")
    Ru = PolynomialRing(F, "u")
    G = {}
    anomalies = []
    center_roots = []
    for h in range(1, H + 1):
        gh = mirror_quotient(N[h], h, F, Rz, Ru)
        expected = expected_delta(h)
        if gh.degree() != expected:
            anomalies.append(
                {"h": h, "expected_degree": expected, "actual_degree": int(gh.degree())}
            )
        if h % 2 == 1 and gh.degree() > 0 and gh[0] == 0:
            center_roots.append(h)
        G[h] = gh
    return F, Ru, G, anomalies, center_roots


def random_squarefree_monic(Ru, degree: int, rng: random.Random):
    if degree == 0:
        return Ru.one()
    u = Ru.gen()
    p = Ru.base_ring().characteristic()
    while True:
        f = u**degree + sum(
            (Ru(rng.randrange(p)) * u**i for i in range(degree)), Ru.zero()
        )
        if f.gcd(f.derivative()).degree() == 0:
            return f


def controlled_family(kind: str, real_G: dict, Ru, H: int, seed: int):
    """Construct one of Q6523's three controls and its clock labels."""
    rng = random.Random(seed)
    labels = {h: h for h in range(1, H + 1)}
    if kind == "clock_shuffle":
        shuffled = list(range(1, H + 1))
        rng.shuffle(shuffled)
        labels = {h: shuffled[h - 1] for h in range(1, H + 1)}
        return dict(real_G), labels

    G = {
        h: random_squarefree_monic(Ru, real_G[h].degree(), rng)
        for h in range(1, H + 1)
    }
    if kind == "mirror_random":
        # Q6523 asks that this control be lifted back to z with the full
        # even/odd mirror parity.  The main test subsequently quotients by
        # that parity, so the round-trip must recover exactly the sampled G.
        F = Ru.base_ring()
        Rz = PolynomialRing(F, "z_control")
        z = Rz.gen()
        u = Ru.gen()
        for h, gh in G.items():
            lifted = sum(
                (Rz(gh[j]) * z ** (2 * j) for j in range(gh.degree() + 1)),
                Rz.zero(),
            )
            if h % 2 == 0:
                lifted *= z
                lifted //= z
            roundtrip = sum(
                (Ru(lifted[2 * j]) * u**j for j in range(lifted.degree() // 2 + 1)),
                Ru.zero(),
            )
            if roundtrip != gh:
                raise AssertionError(("mirror-control round-trip failure", h))
    elif kind != "degree_random":
        raise ValueError(f"unknown control kind: {kind}")
    return G, labels


def build_matrix(F, Ru, G: dict, labels: dict, H: int, M: int, A: int, B: int):
    """Build the exact direct-mod-G_h^M matrix of Q6523 section 2.3."""
    u = Ru.gen()
    Beff = min(B, H - 1)
    cols = [(i, j) for i in range(A + 1) for j in range(Beff + 1)]
    rows = []
    for h in range(1, H + 1):
        modulus = G[h] ** M
        d = modulus.degree()
        if d == 0:
            continue
        rem = [(u**i) % modulus for i in range(A + 1)]
        hpows = [F(labels[h]) ** j for j in range(Beff + 1)]
        for k in range(d):
            rows.append([hpows[j] * rem[i][k] for i, j in cols])
    return matrix(F, rows), cols


def build_hasse_matrix(F, Ru, G: dict, labels: dict, H: int, M: int, A: int, B: int):
    """Equivalent squarefree Hasse matrix, used only by the self-test."""
    u = Ru.gen()
    Beff = min(B, H - 1)
    cols = [(i, j) for i in range(A + 1) for j in range(Beff + 1)]
    rows = []
    for h in range(1, H + 1):
        d = G[h].degree()
        if d == 0:
            continue
        hpows = [F(labels[h]) ** j for j in range(Beff + 1)]
        for q in range(M):
            for k in range(d):
                row = []
                for i, j in cols:
                    if i < q:
                        row.append(F.zero())
                    else:
                        remainder = (u ** (i - q)) % G[h]
                        row.append(hpows[j] * F(binomial(i, q)) * remainder[k])
                rows.append(row)
    return matrix(F, rows)


def rank_record(Mat, cols, case: Case, A: int, variant: str, kind: str, seed, anomalies, center_roots):
    started = time.time()
    rank = int(Mat.rank())
    elapsed = time.time() - started
    Q = int(Mat.ncols())
    T = int(Mat.nrows())
    nullity = Q - rank
    nullity0 = max(0, Q - T)
    surplus = nullity - nullity0
    record = {
        "record": "algebra_rank",
        "p": case.p,
        "H": case.H,
        "M": case.M,
        "A": A,
        "B": case.B,
        "Beff": min(case.B, case.H - 1),
        "variant": variant,
        "control": kind,
        "seed": seed,
        "rows": T,
        "cols": Q,
        "rank": rank,
        "nullity": nullity,
        "nullity0": nullity0,
        "surplus": surplus,
        "candidate_anomaly": surplus >= 1,
        "degree_anomalies": anomalies,
        "odd_center_roots": center_roots,
        "rank_seconds": round(elapsed, 6),
    }

    if surplus > 0 and kind == "real":
        basis = Mat.right_kernel_matrix(basis="echelon")
        sparse_basis = []
        for vector in basis.rows():
            support = []
            for index, coefficient in enumerate(vector):
                if coefficient:
                    i, j = cols[index]
                    support.append([i, j, int(coefficient)])
            sparse_basis.append(support)
        record["kernel_rref_sparse"] = sparse_basis
        record["lowest_support"] = min(sparse_basis, key=len) if sparse_basis else []
    return record


def pressure_A(G: dict, H: int, M: int, B: int) -> int:
    """Largest A with (A+1)(B_eff+1)<M*sum_h deg(G_h)."""
    T = int(M * sum(G[h].degree() for h in range(1, H + 1)))
    width_y = min(B, H - 1) + 1
    value = int((T - 1) // width_y - 1)
    if value < 0:
        raise ValueError("no positive pressure window with Q<T")
    assert (value + 1) * width_y < T
    assert (value + 2) * width_y >= T
    return value


def orbit_points(p: int):
    """Small-instance projective orbit for the continuant cross-check."""
    def key(x, y):
        if x % p:
            return (1, y * pow(x, -1, p) % p)
        return (0, 1)

    b = [0] * (p - 1)
    c = [0] * (p - 1)
    b[0], b[1] = 1, 5 % p
    c[0], c[1] = 0, 1
    for n in range(2, p - 1):
        A = (34 * n**3 - 51 * n**2 + 27 * n - 5) % p
        B = (n - 1) ** 3 % p
        inv = pow(n**3 % p, -1, p)
        b[n] = (A * b[n - 1] - B * b[n - 2]) * inv % p
        c[n] = (A * c[n - 1] - B * c[n - 2]) * inv % p
    return [key(b[n], c[n]) for n in range(p - 1)]


def self_test():
    """Machine-check the recurrence, collision roots, parity, and jet matrix."""
    p, H = 101, 8
    F, Ru, G, anomalies, _ = real_mirror_quotients(p, H)
    assert not anomalies

    # Coefficient-by-coefficient agreement with CRON_avggcd.py's verified
    # low-to-high polynomial implementation.
    from CRON_avggcd import gap_polys as list_gap_polys

    _, sage_N = gap_polynomials(F, H)
    list_N = list_gap_polys(p, H)
    for h in range(1, H + 1):
        coeffs = [int(sage_N[h][i]) for i in range(sage_N[h].degree() + 1)]
        assert coeffs == list_N[h], (h, coeffs, list_N[h])
        assert G[h].degree() == expected_delta(h)

    # First-principles root/collision equivalence on the whole valid interval.
    pts = orbit_points(p)
    Rx = sage_N[1].parent()
    for h in range(1, H + 1):
        for r in range(0, p - 1 - h):
            assert (sage_N[h](F(r)) == 0) == (pts[r] == pts[r + h]), (h, r)

    labels = {h: h for h in range(1, H + 1)}
    direct, _ = build_matrix(F, Ru, G, labels, H=5, M=2, A=9, B=2)
    hasse = build_hasse_matrix(F, Ru, G, labels, H=5, M=2, A=9, B=2)
    assert direct.rank() == hasse.rank()

    # A hand-checkable one-row case: H=2 has only the linear G_2.
    one_row, _ = build_matrix(F, Ru, G, labels, H=2, M=1, A=0, B=0)
    assert (one_row.nrows(), one_row.ncols(), one_row.rank()) == (1, 1, 1)
    print(json.dumps({"record": "self_test", "test": "algebra", "status": "OK"}))


def selected_cases(indices: Iterable[int]) -> list[Case]:
    all_cases = [Case(*values) for values in CASES]
    chosen = []
    for index in indices:
        if not 1 <= index <= len(all_cases):
            raise SystemExit(f"case index must lie in 1..{len(all_cases)}")
        chosen.append(all_cases[index - 1])
    return chosen


def run_case(case: Case, controls: int):
    F, Ru, real_G, anomalies, center_roots = real_mirror_quotients(case.p, case.H)
    A_pressure = pressure_A(real_G, case.H, case.M, case.B)
    if A_pressure < case.A:
        raise AssertionError(("pressure A is below Q6523 A", case, A_pressure))
    print(
        json.dumps(
            {
                "record": "case_start",
                "p": case.p,
                "H": case.H,
                "M": case.M,
                "A": case.A,
                "A_pressure": A_pressure,
                "B": case.B,
                "controls_per_kind": controls,
            },
            sort_keys=True,
        ),
        flush=True,
    )

    datasets = [("real", None, real_G, {h: h for h in range(1, case.H + 1)})]
    for kind_index, kind in enumerate(CONTROL_KINDS):
        for trial in range(controls):
            seed = 1000 + 100 * kind_index + trial
            G, labels = controlled_family(kind, real_G, Ru, case.H, seed)
            datasets.append((kind, seed, G, labels))

    records = []
    for kind, seed, G, labels in datasets:
        # Build the pressure matrix once.  The A-window is its initial column
        # block because columns are ordered first by u-degree, then y-degree.
        pressure_matrix, pressure_cols = build_matrix(
            F, Ru, G, labels, case.H, case.M, A_pressure, case.B
        )
        q_main = (case.A + 1) * (min(case.B, case.H - 1) + 1)
        main_matrix = pressure_matrix.matrix_from_columns(range(q_main))
        main_cols = pressure_cols[:q_main]
        for Mat, cols, A, variant in (
            (main_matrix, main_cols, case.A, "A"),
            (pressure_matrix, pressure_cols, A_pressure, "A_pressure"),
        ):
            record = rank_record(
                Mat,
                cols,
                case,
                A,
                variant,
                kind,
                seed,
                anomalies,
                center_roots,
            )
            records.append(record)
            print(json.dumps(record, sort_keys=True), flush=True)

    for variant in ("A", "A_pressure"):
        real = next(r for r in records if r["control"] == "real" and r["variant"] == variant)
        controls_here = [r for r in records if r["control"] != "real" and r["variant"] == variant]
        same_or_larger = sum(r["surplus"] >= real["surplus"] for r in controls_here)
        summary = {
            "record": "case_summary",
            "p": case.p,
            "H": case.H,
            "M": case.M,
            "variant": variant,
            "real_surplus": real["surplus"],
            "control_trials": len(controls_here),
            "controls_with_surplus_at_least_real": same_or_larger,
            "empirical_control_frequency": (
                same_or_larger / len(controls_here) if controls_here else None
            ),
            "verdict": "candidate_anomaly" if real["surplus"] >= 1 else "no_surplus",
        }
        print(json.dumps(summary, sort_keys=True), flush=True)
    return records


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--case",
        type=int,
        action="append",
        dest="cases",
        help="run only this 1-based Q6523 case (repeatable)",
    )
    parser.add_argument(
        "--controls",
        type=int,
        default=5,
        help="trials for each of the three controls (default: 5, as in Q6523 skeleton)",
    )
    parser.add_argument("--self-test-only", action="store_true")
    return parser.parse_args()


def main():
    args = parse_args()
    if args.controls < 1:
        raise SystemExit("--controls must be positive")
    started = time.time()
    self_test()
    if args.self_test_only:
        return
    cases = selected_cases(args.cases or range(1, len(CASES) + 1))
    all_records = []
    for case in cases:
        all_records.extend(run_case(case, args.controls))
    real_records = [r for r in all_records if r["control"] == "real"]
    positive = [r for r in real_records if r["surplus"] > 0]
    final = {
        "record": "algebra_final",
        "cases": len(cases),
        "real_windows": len(real_records),
        "positive_windows": len(positive),
        "all_degree_checks_ok": all(not r["degree_anomalies"] for r in real_records),
        "verdict": "candidate_anomaly" if positive else "no_syzygy_surplus",
        "elapsed_seconds": round(time.time() - started, 3),
    }
    print(json.dumps(final, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
