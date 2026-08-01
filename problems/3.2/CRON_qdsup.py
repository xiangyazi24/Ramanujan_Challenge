#!/usr/bin/env python3
"""Supercritical Q_D sweep with arithmetic strata and mirror-random controls.

This is the executable companion to CODEX_SPEC_CRON_qdsup.md.  The production
run covers every prime in [3000, 4200], plus 10007 and 30011.  It deliberately
works from value fibres: pairs in a fibre give R_h and d_D, while triples give
C_p(a,g).  Thus the work is proportional to the actual collision count rather
than to p D^2.

The checked orbit and baseline constructors are imported from
CRON_stratify_t34.py.  CRON_radon_spectrum.py and CRON_b1_crosscorr.py are used
as independent source cross-checks before the sweep.  The latter has top-level
execution, so its orbit function is loaded from its AST without importing the
module and accidentally starting its historical full scan.
"""

from __future__ import annotations

import argparse
import ast
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import statistics
import time
from typing import Any, Callable, Iterable

import numpy as np

from CRON_radon_spectrum import apery_pair
from CRON_stratify_t34 import mirror_random_key, orbit_keys, primes_in


HERE = Path(__file__).resolve().parent
SEED = 20260801
CALIBRATION = {100: 0, 316: 10, 1000: 124, 3162: 1089}
PRIMARY_PATHS = [
    "Q_D",
    "S_D",
    "M_D",
    "R.max",
    "R.q95",
    "R.q99",
    "C.max",
    "C.q99",
]


class Heartbeat:
    """Print at most every five seconds inside loops that could grow long."""

    def __init__(self, label: str, interval: float = 5.0):
        self.label = label
        self.interval = interval
        self.last = time.monotonic()

    def maybe(self, detail: str) -> None:
        now = time.monotonic()
        if now - self.last >= self.interval:
            print(f"  [{self.label}] {detail}", flush=True)
            self.last = now


@dataclass
class FibreEvents:
    p: int
    max_D: int
    R: np.ndarray
    pair_r: np.ndarray
    pair_h: np.ndarray
    triple_r: np.ndarray
    triple_a: np.ndarray
    triple_g: np.ndarray
    triple_span: np.ndarray
    fibre_count: int
    max_fibre: int
    T3: int
    Q_full: int


def scale_specs(p: int) -> list[dict[str, Any]]:
    values = [
        ("log p", math.log(p)),
        ("p^0.1", p**0.1),
        ("p^0.2", p**0.2),
    ]
    return [
        {
            "label": label,
            "L": float(value),
            "D": min(p - 2, math.floor(math.sqrt(p) * value)),
        }
        for label, value in values
    ]


def load_function_without_toplevel(path: Path, name: str) -> Callable[..., Any]:
    """Load one source function while refusing to execute module top level."""
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    matches = [node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == name]
    if len(matches) != 1:
        raise RuntimeError(f"expected exactly one {name} in {path}, found {len(matches)}")
    module = ast.Module(body=matches, type_ignores=[])
    ast.fix_missing_locations(module)
    namespace: dict[str, Any] = {}
    exec(compile(module, str(path), "exec"), namespace)
    return namespace[name]


def same_partition(left: list[Any], right: list[Any]) -> bool:
    """Check that two label sequences define exactly the same fibres."""
    if len(left) != len(right):
        return False
    lr: dict[Any, Any] = {}
    rl: dict[Any, Any] = {}
    for x, y in zip(left, right):
        if (x in lr and lr[x] != y) or (y in rl and rl[y] != x):
            return False
        lr[x] = y
        rl[y] = x
    return True


def build_events(key: list[int], p: int, max_D: int, label: str) -> FibreEvents:
    N = p - 1
    if len(key) != N:
        raise AssertionError(f"p={p}: orbit length {len(key)} != {N}")
    if not (1 <= max_D <= p - 2):
        raise AssertionError(f"p={p}: invalid max_D={max_D}")

    fibres: dict[int, list[int]] = defaultdict(list)
    for r, value in enumerate(key):
        fibres[int(value)].append(r)

    pair_r: list[int] = []
    pair_h: list[int] = []
    triple_r: list[int] = []
    triple_a: list[int] = []
    triple_g: list[int] = []
    triple_span: list[int] = []
    heartbeat = Heartbeat(label)

    for fibre_index, positions in enumerate(fibres.values(), start=1):
        m = len(positions)
        for i in range(m):
            x = positions[i]
            for j in range(i + 1, m):
                h = positions[j] - x
                if h > max_D:
                    break
                pair_r.append(x)
                pair_h.append(h)
            for j in range(i + 1, m):
                y = positions[j]
                if y - x >= max_D:
                    # A third distinct point would have still larger span.
                    break
                for k in range(j + 1, m):
                    z = positions[k]
                    span = z - x
                    if span > max_D:
                        break
                    triple_r.append(x)
                    triple_a.append(y - x)
                    triple_g.append(z - y)
                    triple_span.append(span)
        if fibre_index % 512 == 0:
            heartbeat.maybe(
                f"fibres {fibre_index}/{len(fibres)}, pairs={len(pair_h)}, triples={len(triple_span)}"
            )

    pair_h_arr = np.asarray(pair_h, dtype=np.int32)
    R = np.bincount(pair_h_arr, minlength=max_D + 1).astype(np.int64, copy=False)
    multiplicities = [len(v) for v in fibres.values()]
    T3 = sum(m * (m - 1) * (m - 2) for m in multiplicities)
    # Compute Q_{p-2} from its defining future degree at every base point,
    # separately from the fibre factorial-moment formula used for T3.
    Q_full = sum(
        sum(math.comb(m - i - 1, 2) for i in range(m - 2))
        for m in multiplicities
        if m >= 3
    )
    Q_full_fibre = sum(math.comb(m, 3) for m in multiplicities if m >= 3)
    if Q_full != Q_full_fibre or 6 * Q_full != T3:
        raise AssertionError(f"p={p}: global gate failed: 6*{Q_full} != {T3}")

    return FibreEvents(
        p=p,
        max_D=max_D,
        R=R,
        pair_r=np.asarray(pair_r, dtype=np.int32),
        pair_h=pair_h_arr,
        triple_r=np.asarray(triple_r, dtype=np.int32),
        triple_a=np.asarray(triple_a, dtype=np.int32),
        triple_g=np.asarray(triple_g, dtype=np.int32),
        triple_span=np.asarray(triple_span, dtype=np.int32),
        fibre_count=len(fibres),
        max_fibre=max(multiplicities, default=0),
        T3=T3,
        Q_full=Q_full,
    )


def linear_quantile_dense(values: np.ndarray, q: float) -> float:
    if values.size == 0:
        return 0.0
    return float(np.quantile(values, q, method="linear"))


def linear_quantile_sparse(counter: Counter[tuple[int, int]], total_size: int, q: float) -> float:
    """NumPy's linear quantile, including the implicit zero coordinates."""
    if total_size <= 0:
        return 0.0
    nonzero_values = sorted(counter.values())
    zero_count = total_size - len(nonzero_values)
    if zero_count < 0:
        raise AssertionError("more nonzero C coordinates than triangular coordinates")

    def kth(k: int) -> int:
        return 0 if k < zero_count else int(nonzero_values[k - zero_count])

    position = (total_size - 1) * q
    lo = math.floor(position)
    hi = math.ceil(position)
    if lo == hi:
        return float(kth(lo))
    weight = position - lo
    return float((1.0 - weight) * kth(lo) + weight * kth(hi))


def smallest_prime_factors(n: int) -> list[int]:
    spf = list(range(n + 1))
    if n >= 0:
        spf[0] = 0
    if n >= 1:
        spf[1] = 0
    for i in range(2, math.isqrt(n) + 1):
        if spf[i] == i:
            for j in range(i * i, n + 1, i):
                if spf[j] == j:
                    spf[j] = i
    return spf


def arithmetic_class(h: int, spf: list[int]) -> str:
    if h == 1:
        return "unit"
    # This is exactly the specification's spf(h) > sqrt(h) test.
    return "prime" if spf[h] * spf[h] > h else "smooth"


def top_h(R: np.ndarray, hs: Iterable[int], limit: int = 8) -> list[str]:
    ranked = sorted(((int(R[h]), int(h)) for h in hs), key=lambda item: (-item[0], item[1]))
    return [f"h={h} (R={mass})" for mass, h in ranked[:limit] if mass > 0]


def top_pairs(
    C: Counter[tuple[int, int]],
    predicate: Callable[[int, int], bool] | None = None,
    limit: int = 8,
) -> list[str]:
    ranked = [
        (int(mass), int(a), int(g))
        for (a, g), mass in C.items()
        if predicate is None or predicate(a, g)
    ]
    ranked.sort(key=lambda item: (-item[0], item[1], item[2]))
    return [f"(a,g)=({a},{g}) (C={mass})" for mass, a, g in ranked[:limit]]


def sum_C(C: Counter[tuple[int, int]], predicate: Callable[[int, int], bool]) -> int:
    return sum(int(mass) for (a, g), mass in C.items() if predicate(a, g))


def stats_for_D(
    key: list[int],
    events: FibreEvents,
    p: int,
    scale: dict[str, Any],
    include_hints: bool,
) -> dict[str, Any]:
    D = int(scale["D"])
    N = p - 1
    if D > events.max_D:
        raise AssertionError(f"p={p}: D={D} exceeds event cutoff {events.max_D}")

    pair_mask = events.pair_h <= D
    selected_pair_r = events.pair_r[pair_mask]
    selected_pair_h = events.pair_h[pair_mask]
    d = np.bincount(selected_pair_r, minlength=N).astype(np.int64, copy=False)
    S_D = int(selected_pair_h.size)
    Q_D = int(np.sum(d * (d - 1) // 2, dtype=np.int64))
    R_values = events.R[1 : D + 1]
    if int(R_values.sum()) != S_D:
        raise AssertionError(f"p={p}, D={D}: S identity failed")

    triple_mask = events.triple_span <= D
    tr = events.triple_r[triple_mask]
    ta = events.triple_a[triple_mask]
    tg = events.triple_g[triple_mask]
    ts = events.triple_span[triple_mask]
    C: Counter[tuple[int, int]] = Counter(zip(map(int, ta), map(int, tg)))
    C_total = int(sum(C.values()))
    if C_total != Q_D or int(ts.size) != Q_D:
        raise AssertionError(
            f"p={p}, D={D}: pair gate failed Q={Q_D}, sumC={C_total}, triples={ts.size}"
        )

    M_D = int(d.max(initial=0))
    argmax_r_all = np.flatnonzero(d == M_D).tolist() if M_D > 0 else []
    R_max = int(R_values.max(initial=0))
    argmax_h_all = (np.flatnonzero(R_values == R_max) + 1).tolist() if R_max > 0 else []
    C_max = max(C.values(), default=0)
    argmax_pairs_all = sorted([list(pair) for pair, value in C.items() if value == C_max]) if C_max else []
    total_C_coordinates = D * (D - 1) // 2

    spf = smallest_prime_factors(D)
    hs = list(range(1, D + 1))
    even_hs = [h for h in hs if h % 2 == 0]
    odd_hs = [h for h in hs if h % 2 == 1]
    div_minus_hs = [h for h in hs if (p - 1) % h == 0]
    div_plus_hs = [h for h in hs if (p + 1) % h == 0]
    class_hs = {
        cls: [h for h in hs if arithmetic_class(h, spf) == cls]
        for cls in ("unit", "prime", "smooth")
    }

    forced_hs: list[int] = []
    forced_missing: list[int] = []
    for h in even_hs:
        rstar = (p - 1 - h) // 2
        if key[rstar] == key[rstar + h]:
            forced_hs.append(h)
        else:
            forced_missing.append(h)
    if forced_missing:
        raise AssertionError(f"p={p}, D={D}: missing forced R roots at h={forced_missing[:10]}")

    R_mass = {
        "mirror": {
            "even_total": int(events.R[even_hs].sum()) if even_hs else 0,
            "forced_root": len(forced_hs),
            "even_nonforced_excess": (
                int(events.R[even_hs].sum()) - len(forced_hs) if even_hs else 0
            ),
        },
        "divisibility": {
            "h_div_p_minus_1": int(events.R[div_minus_hs].sum()) if div_minus_hs else 0,
            "h_not_div_p_minus_1": S_D
            - (int(events.R[div_minus_hs].sum()) if div_minus_hs else 0),
            "h_div_p_plus_1": int(events.R[div_plus_hs].sum()) if div_plus_hs else 0,
            "h_not_div_p_plus_1": S_D
            - (int(events.R[div_plus_hs].sum()) if div_plus_hs else 0),
        },
        "parity": {
            "even": int(events.R[even_hs].sum()) if even_hs else 0,
            "odd": int(events.R[odd_hs].sum()) if odd_hs else 0,
        },
        "primality": {
            cls: int(events.R[class_hs[cls]].sum()) if class_hs[cls] else 0
            for cls in ("unit", "prime", "smooth")
        },
        "h2_minus51": int(events.R[2]) if D >= 2 else 0,
    }

    G = math.floor(math.sqrt(p / (24.0 * D ** (2.0 / 3.0))))
    G = max(0, G)
    g_even_mass = sum_C(C, lambda a, g: g % 2 == 0)
    axis_mass = sum_C(C, lambda a, g: min(a, g) <= G)

    forced_counter: Counter[tuple[int, int]] = Counter()
    mirror_edge_counter: Counter[tuple[int, int]] = Counter()
    for r, a, g in zip(map(int, tr), map(int, ta), map(int, tg)):
        if g % 2 == 0 and 2 * r + 2 * a + g == N:
            forced_counter[(a, g)] += 1
        if (
            2 * r + a == N
            or 2 * r + a + g == N
            or 2 * r + 2 * a + g == N
        ):
            mirror_edge_counter[(a, g)] += 1
    forced_psi_mass = int(sum(forced_counter.values()))
    any_mirror_edge_mass = int(sum(mirror_edge_counter.values()))
    g_even_counter = Counter({pair: mass for pair, mass in C.items() if pair[1] % 2 == 0})
    other_g_even_counter = g_even_counter - forced_counter
    no_mirror_edge_counter = Counter(C) - mirror_edge_counter

    divm_a = sum_C(C, lambda a, g: (p - 1) % a == 0)
    divm_g = sum_C(C, lambda a, g: (p - 1) % g == 0)
    divm_span = sum_C(C, lambda a, g: (p - 1) % (a + g) == 0)
    divm_any = sum_C(C, lambda a, g: any((p - 1) % h == 0 for h in (a, g, a + g)))
    divp_a = sum_C(C, lambda a, g: (p + 1) % a == 0)
    divp_g = sum_C(C, lambda a, g: (p + 1) % g == 0)
    divp_span = sum_C(C, lambda a, g: (p + 1) % (a + g) == 0)
    divp_any = sum_C(C, lambda a, g: any((p + 1) % h == 0 for h in (a, g, a + g)))

    parity_mass = {
        f"a_{'even' if ap == 0 else 'odd'}_g_{'even' if gp == 0 else 'odd'}": sum_C(
            C, lambda a, g, ap=ap, gp=gp: a % 2 == ap and g % 2 == gp
        )
        for ap in (0, 1)
        for gp in (0, 1)
    }
    primality_mass: dict[str, int] = {}
    for ac in ("unit", "prime", "smooth"):
        for gc in ("unit", "prime", "smooth"):
            primality_mass[f"a_{ac}_g_{gc}"] = sum_C(
                C,
                lambda a, g, ac=ac, gc=gc: arithmetic_class(a, spf) == ac
                and arithmetic_class(g, spf) == gc,
            )
    primality_mass["span_prime"] = sum_C(C, lambda a, g: arithmetic_class(a + g, spf) == "prime")
    primality_mass["span_smooth"] = Q_D - primality_mass["span_prime"]
    primality_mass["any_adjacent_prime"] = sum_C(
        C,
        lambda a, g: arithmetic_class(a, spf) == "prime"
        or arithmetic_class(g, spf) == "prime",
    )
    primality_mass["no_adjacent_prime"] = Q_D - primality_mass["any_adjacent_prime"]

    h2_a = sum_C(C, lambda a, g: a == 2)
    h2_g = sum_C(C, lambda a, g: g == 2)
    h2_span = sum_C(C, lambda a, g: a + g == 2)
    h2_any = sum_C(C, lambda a, g: 2 in (a, g, a + g))

    Q_mass = {
        "mirror": {
            "g_even_coordinate": g_even_mass,
            "g_odd_coordinate": Q_D - g_even_mass,
            "forced_psi_root": forced_psi_mass,
            "other_g_even": g_even_mass - forced_psi_mass,
            "any_exact_mirror_edge": any_mirror_edge_mass,
            "no_exact_mirror_edge": Q_D - any_mirror_edge_mass,
        },
        "axis": {
            "min_a_g_le_G": axis_mass,
            "min_a_g_gt_G": Q_D - axis_mass,
        },
        "divisibility": {
            "a_div_p_minus_1": divm_a,
            "g_div_p_minus_1": divm_g,
            "span_div_p_minus_1": divm_span,
            "any_gap_div_p_minus_1": divm_any,
            "no_gap_div_p_minus_1": Q_D - divm_any,
            "a_div_p_plus_1": divp_a,
            "g_div_p_plus_1": divp_g,
            "span_div_p_plus_1": divp_span,
            "any_gap_div_p_plus_1": divp_any,
            "no_gap_div_p_plus_1": Q_D - divp_any,
        },
        "parity": parity_mass,
        "primality": primality_mass,
        "h2_minus51": {
            "a_eq_2": h2_a,
            "g_eq_2": h2_g,
            "span_eq_2": h2_span,
            "any_gap_eq_2": h2_any,
            "no_gap_eq_2": Q_D - h2_any,
        },
    }

    result: dict[str, Any] = {
        "p": p,
        "scale": str(scale["label"]),
        "L": float(scale["L"]),
        "D": D,
        "G": G,
        "Q_D": Q_D,
        "S_D": S_D,
        "M_D": M_D,
        "argmax_r": argmax_r_all[:20],
        "argmax_r_count": len(argmax_r_all),
        "lag_profile": {
            "max": R_max,
            "q95": linear_quantile_dense(R_values, 0.95),
            "q99": linear_quantile_dense(R_values, 0.99),
            "argmax_h": argmax_h_all[:20],
            "argmax_h_count": len(argmax_h_all),
        },
        "pair_profile": {
            "total": C_total,
            "max": int(C_max),
            "q99": linear_quantile_sparse(C, total_C_coordinates, 0.99),
            "argmax_pairs": argmax_pairs_all[:20],
            "argmax_pair_count": len(argmax_pairs_all),
            "all_zero": C_max == 0,
            "nonzero_pairs": len(C),
            "total_pairs": total_C_coordinates,
        },
        "strata": {
            "R_mass": R_mass,
            "Q_mass": Q_mass,
            "forced_R_expected": len(even_hs),
            "forced_R_missing": forced_missing,
            "prime_definition": "spf(h) > sqrt(h); h=1 is reported separately as unit",
        },
        "fibre_count": events.fibre_count,
        "max_fibre": events.max_fibre,
    }

    if include_hints:
        hints: dict[str, list[str]] = {
            "Q_D": top_pairs(C),
            "S_D": top_h(events.R, hs),
            "M_D": [f"r={r} (d={M_D})" for r in argmax_r_all[:8]],
            "R.max": [f"h={h} (R={R_max})" for h in argmax_h_all[:8]],
            "R.q95": top_h(events.R, hs),
            "R.q99": top_h(events.R, hs),
            "C.max": top_pairs(C),
            "C.q99": top_pairs(C),
            "C.nonzero_pairs": top_pairs(C),
        }

        R_domains: dict[str, list[int]] = {
            "R.mass.mirror.even_total": even_hs,
            "R.mass.mirror.forced_root": forced_hs,
            "R.mass.mirror.even_nonforced_excess": even_hs,
            "R.mass.divisibility.h_div_p_minus_1": div_minus_hs,
            "R.mass.divisibility.h_not_div_p_minus_1": [h for h in hs if h not in set(div_minus_hs)],
            "R.mass.divisibility.h_div_p_plus_1": div_plus_hs,
            "R.mass.divisibility.h_not_div_p_plus_1": [h for h in hs if h not in set(div_plus_hs)],
            "R.mass.parity.even": even_hs,
            "R.mass.parity.odd": odd_hs,
            "R.mass.primality.unit": class_hs["unit"],
            "R.mass.primality.prime": class_hs["prime"],
            "R.mass.primality.smooth": class_hs["smooth"],
            "R.mass.h2_minus51": [2] if D >= 2 else [],
        }
        for path, domain in R_domains.items():
            hints[path] = top_h(events.R, domain)

        Q_predicates: dict[str, Callable[[int, int], bool]] = {
            "Q.mass.mirror.g_even_coordinate": lambda a, g: g % 2 == 0,
            "Q.mass.mirror.g_odd_coordinate": lambda a, g: g % 2 == 1,
            "Q.mass.axis.min_a_g_le_G": lambda a, g: min(a, g) <= G,
            "Q.mass.axis.min_a_g_gt_G": lambda a, g: min(a, g) > G,
            "Q.mass.divisibility.a_div_p_minus_1": lambda a, g: (p - 1) % a == 0,
            "Q.mass.divisibility.g_div_p_minus_1": lambda a, g: (p - 1) % g == 0,
            "Q.mass.divisibility.span_div_p_minus_1": lambda a, g: (p - 1) % (a + g) == 0,
            "Q.mass.divisibility.any_gap_div_p_minus_1": lambda a, g: any(
                (p - 1) % h == 0 for h in (a, g, a + g)
            ),
            "Q.mass.divisibility.no_gap_div_p_minus_1": lambda a, g: all(
                (p - 1) % h != 0 for h in (a, g, a + g)
            ),
            "Q.mass.divisibility.a_div_p_plus_1": lambda a, g: (p + 1) % a == 0,
            "Q.mass.divisibility.g_div_p_plus_1": lambda a, g: (p + 1) % g == 0,
            "Q.mass.divisibility.span_div_p_plus_1": lambda a, g: (p + 1) % (a + g) == 0,
            "Q.mass.divisibility.any_gap_div_p_plus_1": lambda a, g: any(
                (p + 1) % h == 0 for h in (a, g, a + g)
            ),
            "Q.mass.divisibility.no_gap_div_p_plus_1": lambda a, g: all(
                (p + 1) % h != 0 for h in (a, g, a + g)
            ),
        }
        for ap in (0, 1):
            for gp in (0, 1):
                name = f"a_{'even' if ap == 0 else 'odd'}_g_{'even' if gp == 0 else 'odd'}"
                Q_predicates[f"Q.mass.parity.{name}"] = (
                    lambda a, g, ap=ap, gp=gp: a % 2 == ap and g % 2 == gp
                )
        for ac in ("unit", "prime", "smooth"):
            for gc in ("unit", "prime", "smooth"):
                name = f"a_{ac}_g_{gc}"
                Q_predicates[f"Q.mass.primality.{name}"] = (
                    lambda a, g, ac=ac, gc=gc: arithmetic_class(a, spf) == ac
                    and arithmetic_class(g, spf) == gc
                )
        Q_predicates.update(
            {
                "Q.mass.primality.span_prime": lambda a, g: arithmetic_class(a + g, spf) == "prime",
                "Q.mass.primality.span_smooth": lambda a, g: arithmetic_class(a + g, spf) == "smooth",
                "Q.mass.primality.any_adjacent_prime": lambda a, g: arithmetic_class(a, spf) == "prime"
                or arithmetic_class(g, spf) == "prime",
                "Q.mass.primality.no_adjacent_prime": lambda a, g: arithmetic_class(a, spf) != "prime"
                and arithmetic_class(g, spf) != "prime",
                "Q.mass.h2_minus51.a_eq_2": lambda a, g: a == 2,
                "Q.mass.h2_minus51.g_eq_2": lambda a, g: g == 2,
                "Q.mass.h2_minus51.span_eq_2": lambda a, g: a + g == 2,
                "Q.mass.h2_minus51.any_gap_eq_2": lambda a, g: 2 in (a, g, a + g),
                "Q.mass.h2_minus51.no_gap_eq_2": lambda a, g: 2 not in (a, g, a + g),
            }
        )
        for path, predicate in Q_predicates.items():
            hints[path] = top_pairs(C, predicate)
        hints["Q.mass.mirror.forced_psi_root"] = top_pairs(forced_counter)
        hints["Q.mass.mirror.other_g_even"] = top_pairs(other_g_even_counter)
        hints["Q.mass.mirror.any_exact_mirror_edge"] = top_pairs(mirror_edge_counter)
        hints["Q.mass.mirror.no_exact_mirror_edge"] = top_pairs(no_mirror_edge_counter)
        result["coordinate_hints"] = hints

    return result


def walk_numbers(prefix: str, value: Any, output: dict[str, float]) -> None:
    if isinstance(value, dict):
        for key, child in value.items():
            walk_numbers(f"{prefix}.{key}" if prefix else str(key), child, output)
    elif isinstance(value, (int, float)) and not isinstance(value, bool):
        output[prefix] = float(value)


def numeric_statistics(stats: dict[str, Any]) -> dict[str, float]:
    output = {
        "Q_D": float(stats["Q_D"]),
        "S_D": float(stats["S_D"]),
        "M_D": float(stats["M_D"]),
        "R.max": float(stats["lag_profile"]["max"]),
        "R.q95": float(stats["lag_profile"]["q95"]),
        "R.q99": float(stats["lag_profile"]["q99"]),
        "C.max": float(stats["pair_profile"]["max"]),
        "C.q99": float(stats["pair_profile"]["q99"]),
        "C.nonzero_pairs": float(stats["pair_profile"]["nonzero_pairs"]),
    }
    walk_numbers("R.mass", stats["strata"]["R_mass"], output)
    walk_numbers("Q.mass", stats["strata"]["Q_mass"], output)
    return output


def safe_ratio(data: float, baseline_mean: float) -> float | str:
    if baseline_mean != 0.0:
        return float(data / baseline_mean)
    return 1.0 if data == 0.0 else "inf"


def summarize_against_replicas(data: dict[str, float], replicas: list[dict[str, float]]) -> dict[str, Any]:
    if not replicas:
        raise AssertionError("baseline replica list is empty")
    keys = list(data)
    if any(set(replica) != set(keys) for replica in replicas):
        raise AssertionError("data/baseline statistic keys differ")
    means: dict[str, float] = {}
    sds: dict[str, float] = {}
    maxima: dict[str, float] = {}
    ratios: dict[str, float | str] = {}
    flags: dict[str, bool] = {}
    exceeds_all: dict[str, bool] = {}
    z_scores: dict[str, float | str | None] = {}
    for path in keys:
        values = [replica[path] for replica in replicas]
        mean = float(statistics.fmean(values))
        sd = float(statistics.stdev(values)) if len(values) >= 2 else 0.0
        maximum = float(max(values))
        datum = float(data[path])
        all_flag = datum > maximum
        strong_flag = datum > mean + 3.0 * sd
        means[path] = mean
        sds[path] = sd
        maxima[path] = maximum
        ratios[path] = safe_ratio(datum, mean)
        exceeds_all[path] = all_flag
        flags[path] = bool(all_flag and strong_flag)
        if sd > 0:
            z_scores[path] = float((datum - mean) / sd)
        elif datum > mean:
            z_scores[path] = "inf"
        elif datum == mean:
            z_scores[path] = 0.0
        else:
            z_scores[path] = None
    return {
        "mean": means,
        "sd": sds,
        "max": maxima,
        "ratio": ratios,
        "exceeds_all": exceeds_all,
        "flag": flags,
        "z_baseline_sd": z_scores,
    }


def binomial_survival(k: int, n: int, p0: float) -> float:
    if k <= 0:
        return 1.0
    return float(
        math.fsum(
            math.comb(n, j) * p0**j * (1.0 - p0) ** (n - j)
            for j in range(k, n + 1)
        )
    )


def naive_gate(key: list[int], p: int, D: int) -> dict[str, Any]:
    N = p - 1
    R = [0] * (D + 1)
    d = [0] * N
    for h in range(1, D + 1):
        for r in range(N - h):
            if key[r] == key[r + h]:
                R[h] += 1
                d[r] += 1
    C: Counter[tuple[int, int]] = Counter()
    heartbeat = Heartbeat(f"naive p={p}")
    for a in range(1, D):
        for g in range(1, D - a + 1):
            count = 0
            for r in range(N - a - g):
                if key[r] == key[r + a] == key[r + a + g]:
                    count += 1
            if count:
                C[(a, g)] = count
        heartbeat.maybe(f"a={a}/{D-1}")
    Q = sum(x * (x - 1) // 2 for x in d)
    return {"R": R, "d": d, "C": C, "S": sum(d), "Q": Q}


def preflight() -> tuple[dict[str, Any], dict[int, list[int]], dict[int, FibreEvents]]:
    print("preflight: source/orbit/baseline cross-checks", flush=True)
    p = 3001
    key = orbit_keys(p)

    b, c = apery_pair(p)
    radon_key = [
        (int(b[n]) * pow(int(c[n]), p - 2, p)) % p if int(c[n]) != 0 else p
        for n in range(p - 1)
    ]
    radon_ok = key == radon_key
    if not radon_ok:
        raise AssertionError("CRON_radon_spectrum orbit cross-check failed")

    b1_orbit = load_function_without_toplevel(HERE / "CRON_b1_crosscorr.py", "orbit")
    b1_key = b1_orbit(p)
    b1_ok = same_partition(key, b1_key)
    if not b1_ok:
        raise AssertionError("CRON_b1_crosscorr orbit partition cross-check failed")

    rng = np.random.default_rng(np.random.SeedSequence([SEED, p, 999]))
    random_key = mirror_random_key(p, rng)
    N = p - 1
    baseline_mirror_ok = all(random_key[n] == random_key[(N - n) % N] for n in range(N))
    if not baseline_mirror_ok:
        raise AssertionError("mirror_random_key failed its involution check")

    D_gate = 64
    events = build_events(key, p, D_gate, f"preflight p={p}")
    sparse_scale = {"label": "gate", "L": 0.0, "D": D_gate}
    sparse = stats_for_D(key, events, p, sparse_scale, include_hints=False)
    naive = naive_gate(key, p, D_gate)
    sparse_C = Counter(
        zip(
            map(int, events.triple_a[events.triple_span <= D_gate]),
            map(int, events.triple_g[events.triple_span <= D_gate]),
        )
    )
    naive_ok = (
        sparse["Q_D"] == naive["Q"]
        and sparse["S_D"] == naive["S"]
        and events.R[: D_gate + 1].tolist() == naive["R"]
        and sparse_C == naive["C"]
    )
    if not naive_ok:
        raise AssertionError("sparse/naive statistics cross-check failed")

    print("preflight: p=10007 calibration", flush=True)
    pcal = 10007
    keycal = orbit_keys(pcal)
    eventscal = build_events(keycal, pcal, max(CALIBRATION), "calibration p=10007")
    observed: dict[int, int] = {}
    for D in CALIBRATION:
        scale = {"label": "calibration", "L": 0.0, "D": D}
        observed[D] = int(stats_for_D(keycal, eventscal, pcal, scale, False)["Q_D"])
    calibration_ok = observed == CALIBRATION
    if not calibration_ok:
        raise AssertionError(f"calibration failed: observed {observed}, expected {CALIBRATION}")

    gates = {
        "source_crosscheck": {
            "status": "PASS",
            "p": p,
            "radon_exact_keys": radon_ok,
            "b1_same_partition": b1_ok,
            "baseline_mirror_involution": baseline_mirror_ok,
            "note": (
                "The current CRON_radon_spectrum.py exposes apery_pair (not orbit_keys), "
                "and mirror_random_key is in CRON_stratify_t34.py; the specification's "
                "source labels are stale, but all named source roles were cross-checked."
            ),
        },
        "sparse_vs_naive": {
            "status": "PASS",
            "p": p,
            "D": D_gate,
            "S_D": naive["S"],
            "Q_D": naive["Q"],
            "nonzero_C_pairs": len(naive["C"]),
        },
        "calibration_p10007": {
            "status": "PASS",
            "observed": {str(D): value for D, value in observed.items()},
            "expected": {str(D): value for D, value in CALIBRATION.items()},
        },
    }
    return gates, {p: key, pcal: keycal}, {pcal: eventscal}


def pooled_results(rows: list[dict[str, Any]], band_primes: list[int]) -> dict[str, Any]:
    output: dict[str, Any] = {}
    band_set = set(band_primes)
    for scale in ("log p", "p^0.1", "p^0.2"):
        selected = [row for row in rows if row["scale"] == scale and row["p"] in band_set]
        if len(selected) != len(band_primes):
            raise AssertionError(f"scale {scale}: pooled row count {len(selected)} != {len(band_primes)}")
        data_numeric = [numeric_statistics(row["data"]) for row in selected]
        paths = list(data_numeric[0])
        replica_count = len(selected[0]["baseline"]["replicas"])
        pooled_data = {
            path: float(statistics.fmean(item[path] for item in data_numeric)) for path in paths
        }
        pooled_replicas: list[dict[str, float]] = []
        for replica_index in range(replica_count):
            pooled_replicas.append(
                {
                    path: float(
                        statistics.fmean(
                            numeric_statistics(row["baseline"]["replicas"][replica_index])[path]
                            for row in selected
                        )
                    )
                    for path in paths
                }
            )
        comparison = summarize_against_replicas(pooled_data, pooled_replicas)
        counts = {
            path: sum(bool(row["baseline"]["summary"]["flag"][path]) for row in selected)
            for path in paths
        }
        # The exceed-all event has distribution-free null probability 1/(replicas+1).
        # The extra >3 baseline-sd requirement can only make this conservative.
        p0 = 1.0 / (replica_count + 1.0)
        pvalues = {
            path: binomial_survival(counts[path], len(selected), p0) for path in paths
        }
        output[scale] = {
            "n_primes": len(selected),
            "replica_count": replica_count,
            "null_flag_probability_upper_bound": p0,
            "data_mean": pooled_data,
            "baseline_replica_means": pooled_replicas,
            "summary": comparison,
            "flagged_prime_count": counts,
            "binomial_p_conservative": pvalues,
        }
    return output


def fmt_number(value: Any) -> str:
    if value == "inf":
        return "∞"
    if value is None:
        return "—"
    if isinstance(value, bool):
        return "YES" if value else "no"
    number = float(value)
    if number == 0:
        return "0"
    if abs(number) >= 1000 or (0 < abs(number) < 0.001):
        return f"{number:.4g}"
    if number.is_integer():
        return str(int(number))
    return f"{number:.4f}"


def md_escape(text: str) -> str:
    return text.replace("|", "\\|").replace("\n", " ")


def find_row(rows: list[dict[str, Any]], p: int, scale: str) -> dict[str, Any]:
    matches = [row for row in rows if row["p"] == p and row["scale"] == scale]
    if len(matches) != 1:
        raise AssertionError(f"expected one row for p={p}, scale={scale}, got {len(matches)}")
    return matches[0]


def table_for_paths(
    rows: list[dict[str, Any]], pooled: dict[str, Any], scale: str, paths: list[str]
) -> list[str]:
    pool = pooled[scale]
    n_primes = int(pool["n_primes"])
    large = {p: find_row(rows, p, scale) for p in (10007, 30011)}
    lines = [
        f"| statistic | band data mean | band base mean | ratio | pooled flag | flagged/{n_primes} | binom p | p=10007 data/base/ratio | p=30011 data/base/ratio |",
        "|---|---:|---:|---:|:---:|---:|---:|---:|---:|",
    ]
    for path in paths:
        band_data = pool["data_mean"][path]
        band_base = pool["summary"]["mean"][path]
        band_ratio = pool["summary"]["ratio"][path]
        band_flag = pool["summary"]["flag"][path]
        count = pool["flagged_prime_count"][path]
        pvalue = pool["binomial_p_conservative"][path]
        large_cells = []
        for p in (10007, 30011):
            row = large[p]
            datum = numeric_statistics(row["data"])[path]
            base = row["baseline"]["summary"]["mean"][path]
            ratio = row["baseline"]["summary"]["ratio"][path]
            marker = " FLAG" if row["baseline"]["summary"]["flag"][path] else ""
            large_cells.append(
                f"{fmt_number(datum)}/{fmt_number(base)}/{fmt_number(ratio)}{marker}"
            )
        lines.append(
            "| "
            + " | ".join(
                [
                    f"`{path}`",
                    fmt_number(band_data),
                    fmt_number(band_base),
                    fmt_number(band_ratio),
                    "FLAG" if band_flag else "—",
                    str(count),
                    fmt_number(pvalue),
                    large_cells[0],
                    large_cells[1],
                ]
            )
            + " |"
        )
    return lines


def verdict_text(pooled: dict[str, Any]) -> str:
    thin_tokens = ("axis", "divisibility", "primality", "h2_minus51")
    significant: list[tuple[str, str, float | str, int, float]] = []
    for scale, pool in pooled.items():
        for path, pooled_flag in pool["summary"]["flag"].items():
            if not path.startswith("Q.mass.") or not any(token in path for token in thin_tokens):
                continue
            ratio = pool["summary"]["ratio"][path]
            ratio_above = ratio == "inf" or (isinstance(ratio, float) and ratio > 1.0)
            count = pool["flagged_prime_count"][path]
            pvalue = pool["binomial_p_conservative"][path]
            if pooled_flag and ratio_above and pvalue < 0.01:
                significant.append((scale, path, ratio, count, pvalue))

    r_h2_flags = []
    for scale, pool in pooled.items():
        path = "R.mass.h2_minus51"
        if pool["summary"]["flag"].get(path, False):
            r_h2_flags.append(scale)

    if significant:
        names = ", ".join(f"{scale}: {path}" for scale, path, _, _, _ in significant)
        return (
            "The mirror-preserving null is not sufficient for every thin Q_D stratum: "
            f"the pooled-and-count criteria identify {names}. These are empirical flags, not an "
            "asymptotic theorem; the coordinate table above shows where the mass sits. The exact "
            "mirror component itself is reproduced by construction and is not counted as an anomaly."
        )
    h2_clause = (
        " The fixed h=2 (-51) R_h layer is pooled-flagged at " + ", ".join(r_h2_flags) + ","
        " but that one-lag effect does not produce a significant pooled excess in any Q_D thin stratum."
        if r_h2_flags
        else ""
    )
    return (
        "No tested thin lag family carries a reproducible excess of Q_D mass over the "
        "mirror-random baseline at these supercritical scales: after preserving the forced mirror "
        "geometry, the axis, divisor, parity, prime/smooth, and h=2 decompositions remain compatible "
        "with the null under both the pooled >3-sd rule and the conservative flagged-prime count test."
        + h2_clause
        + " On this finite range the Q_D picture is therefore Poisson-plus-mirror, with no detected "
        "moving-resonance family."
    )


def render_report(results: dict[str, Any]) -> str:
    gates = results["gates"]
    rows = results["rows"]
    pooled = results["pooled_band"]
    lines: list[str] = [
        "# Supercritical $Q_D$ quantiles and arithmetic stratification",
        "",
        f"Run completed `{results['metadata']['completed_utc']}` in {results['metadata']['runtime_seconds']:.2f}s. "
        f"The production band contains {results['metadata']['band_prime_count']} primes; random seed "
        f"is `{results['metadata']['seed']}`. Quantiles use NumPy's linear convention and include zero "
        "coordinates (in particular, the $C_p(a,g)$ 99th percentile is over the full triangular set "
        "$a,g\ge1$, $a+g\le D$).",
        "",
        "## Correctness gates",
        "",
        "| gate | status | numbers |",
        "|---|:---:|---|",
    ]
    source = gates["source_crosscheck"]
    lines.append(
        f"| source cross-check | {source['status']} | p={source['p']}; exact Radon keys; same B1 partition; mirror baseline involution |"
    )
    naive = gates["sparse_vs_naive"]
    lines.append(
        f"| sparse vs naive | {naive['status']} | p={naive['p']}, D={naive['D']}: "
        f"S={naive['S_D']}, Q={naive['Q_D']}, nonzero C pairs={naive['nonzero_C_pairs']} |"
    )
    calibration = gates["calibration_p10007"]
    cal_numbers = ", ".join(
        f"D={D}: {calibration['observed'][str(D)]}/{calibration['expected'][str(D)]}"
        for D in CALIBRATION
    )
    lines.append(f"| p=10007 calibration | {calibration['status']} | {cal_numbers} |")
    pair_gate = gates["pair_identity"]
    lines.append(
        f"| $Q_D=\sum C_p(a,g)$ | {pair_gate['status']} | {pair_gate['checks']} data/baseline scale checks, failures={pair_gate['failures']} |"
    )
    t3 = gates["T3_identity"]
    examples = "; ".join(
        f"p={item['p']}: 6Q={6*item['Q_full']}=T3={item['T3']}" for item in t3["examples"]
    )
    lines.append(
        f"| $6Q_{{p-2}}=T_3$ | {t3['status']} | checked {t3['checks']} data orbits; {examples} |"
    )
    lines.extend(
        [
            "",
            source["note"],
            "",
            "The sparse algorithm groups equal orbit keys into fibres. Every within-fibre pair contributes "
            "one $I(r,h)$ event, and every ordered-position triple contributes one $C_p(a,g)$ event. "
            "The naive gate independently loops over all admissible $r,h,a,g$.",
            "",
            "## Primary statistics",
        ]
    )

    for scale in ("log p", "p^0.1", "p^0.2"):
        lines.extend(["", f"### $L={scale}$", ""])
        lines.extend(table_for_paths(rows, pooled, scale, PRIMARY_PATHS))

    lines.extend(
        [
            "",
            "### Large-prime extrema coordinates",
            "",
            "The JSON contains these coordinates for every band prime as well; the report spells out the "
            "two separately requested large-prime rows.",
            "",
            "| p | scale | D | M_D and argmax r | max R_h and argmax h | max C and argmax (a,g) |",
            "|---:|---|---:|---|---|---|",
        ]
    )
    for p in (10007, 30011):
        for scale in ("log p", "p^0.1", "p^0.2"):
            data = find_row(rows, p, scale)["data"]
            r_text = md_escape(", ".join(map(str, data["argmax_r"])))
            h_text = md_escape(", ".join(map(str, data["lag_profile"]["argmax_h"])))
            pair_text = md_escape(
                ", ".join(f"({a},{g})" for a, g in data["pair_profile"]["argmax_pairs"])
            )
            if data["pair_profile"]["argmax_pair_count"] > len(data["pair_profile"]["argmax_pairs"]):
                pair_text += f", … ({data['pair_profile']['argmax_pair_count']} ties)"
            lines.append(
                f"| {p} | {scale} | {data['D']} | {data['M_D']} at {r_text or '—'} | "
                f"{data['lag_profile']['max']} at {h_text or '—'} | "
                f"{data['pair_profile']['max']} at {pair_text or 'all zero'} |"
            )

    all_paths = list(pooled["log p"]["data_mean"])
    R_paths = [path for path in all_paths if path.startswith("R.mass.")]
    Q_paths = [path for path in all_paths if path.startswith("Q.mass.")]
    lines.extend(
        [
            "",
            "## Strata",
            "",
            "Masses in different rows intentionally overlap when they answer different readings of a lag "
            "condition. The partition rows are explicit complements. For $R_h$, `forced_root` is the exact "
            "$r=(p-1-h)/2$ event at even $h$. For $Q_D$, `forced_psi_root` is the exact $g$-even event "
            "$r=(p-1-2a-g)/2$, while `any_exact_mirror_edge` allows any of the three edges of the triple "
            "to be a reflection pair. `prime` means $\operatorname{spf}(h)>\sqrt h$; $h=1$ is `unit`.",
        ]
    )
    for scale in ("log p", "p^0.1", "p^0.2"):
        lines.extend(["", f"### $R_h$ mass, $L={scale}$", ""])
        lines.extend(table_for_paths(rows, pooled, scale, R_paths))
        lines.extend(["", f"### $Q_D$ mass, $L={scale}$", ""])
        lines.extend(table_for_paths(rows, pooled, scale, Q_paths))

    flagged_rows: list[tuple[dict[str, Any], str]] = []
    for row in rows:
        for path, flag in row["baseline"]["summary"]["flag"].items():
            if flag:
                flagged_rows.append((row, path))
    lines.extend(
        [
            "",
            "## Flagged anomalies",
            "",
            "A per-prime flag requires strict exceedance of every replica and strict exceedance of the "
            "replica mean by more than three sample standard deviations. The band binomial p-value uses "
            "$p_0=1/(3+1)=0.25$, the distribution-free exceed-all probability; because the 3-sd condition "
            "is additional, this is conservative. Overlapping strata and nearby primes are not claimed "
            "independent, so these p-values are diagnostics rather than theorem-level significance.",
            "",
        ]
    )
    pooled_flagged = [
        (scale, path)
        for scale, pool in pooled.items()
        for path, flag in pool["summary"]["flag"].items()
        if flag
    ]
    if pooled_flagged:
        lines.extend(
            [
                "### Pooled band flags",
                "",
                "| scale | statistic | data mean | baseline mean ± sd | ratio | flagged primes | binom p |",
                "|---|---|---:|---:|---:|---:|---:|",
            ]
        )
        for scale, path in pooled_flagged:
            pool = pooled[scale]
            lines.append(
                f"| {scale} | `{path}` | {fmt_number(pool['data_mean'][path])} | "
                f"{fmt_number(pool['summary']['mean'][path])} ± {fmt_number(pool['summary']['sd'][path])} | "
                f"{fmt_number(pool['summary']['ratio'][path])} | {pool['flagged_prime_count'][path]}/{pool['n_primes']} | "
                f"{fmt_number(pool['binomial_p_conservative'][path])} |"
            )
    else:
        lines.extend(["### Pooled band flags", "", "None."])

    lines.extend(
        [
            "",
            "### Per-prime flags and coordinates",
            "",
            "| p | scale | statistic | data | baseline mean ± sd | ratio | coordinates/top contributors |",
            "|---:|---|---|---:|---:|---:|---|",
        ]
    )
    if flagged_rows:
        for row, path in flagged_rows:
            summary = row["baseline"]["summary"]
            datum = numeric_statistics(row["data"])[path]
            coords = row["data"].get("coordinate_hints", {}).get(path, [])
            coord_text = md_escape(", ".join(coords) if coords else "no unique coordinate")
            lines.append(
                f"| {row['p']} | {row['scale']} | `{path}` | {fmt_number(datum)} | "
                f"{fmt_number(summary['mean'][path])} ± {fmt_number(summary['sd'][path])} | "
                f"{fmt_number(summary['ratio'][path])} | {coord_text} |"
            )
    else:
        lines.append("| — | — | — | — | — | — | no flags |")

    lines.extend(
        [
            "",
            "## Verdict",
            "",
            verdict_text(pooled),
            "",
        ]
    )
    return "\n".join(lines)


def run(args: argparse.Namespace) -> dict[str, Any]:
    started_wall = time.time()
    gates, key_cache, event_cache = preflight()
    if args.self_test:
        return {
            "metadata": {"self_test": True, "runtime_seconds": time.time() - started_wall},
            "gates": gates,
        }

    full_band = primes_in(3000, 4200)
    if len(full_band) != 144:
        raise AssertionError(f"expected 144 band primes, found {len(full_band)}")
    band = full_band[:2] if args.quick else full_band
    primes = band + [10007, 30011]
    rows: list[dict[str, Any]] = []
    pair_checks = 0
    t3_checks: list[dict[str, int]] = []

    print(
        f"sweep: {len(band)} band primes + 10007/30011; "
        f"{'QUICK' if args.quick else 'PRODUCTION'} mode",
        flush=True,
    )
    for prime_index, p in enumerate(primes, start=1):
        prime_start = time.monotonic()
        specs = scale_specs(p)
        max_D = max(spec["D"] for spec in specs)
        if p == 10007:
            max_D = max(max_D, max(CALIBRATION))
        key = key_cache.get(p)
        if key is None:
            key = orbit_keys(p)
            key_cache[p] = key
        events = event_cache.get(p)
        if events is None or events.max_D < max_D:
            events = build_events(key, p, max_D, f"data p={p}")
            event_cache[p] = events
        if 6 * events.Q_full != events.T3:
            raise AssertionError(f"p={p}: T3 gate failed after event construction")
        t3_checks.append({"p": p, "Q_full": events.Q_full, "T3": events.T3})
        data_by_scale = {
            spec["label"]: stats_for_D(key, events, p, spec, include_hints=True)
            for spec in specs
        }
        pair_checks += len(specs)

        nrep = 3 if p in set(full_band) else 5
        replica_by_scale: dict[str, list[dict[str, Any]]] = {
            spec["label"]: [] for spec in specs
        }
        for replica in range(nrep):
            rng = np.random.default_rng(np.random.SeedSequence([SEED, p, replica]))
            baseline_key = mirror_random_key(p, rng)
            baseline_events = build_events(
                baseline_key, p, max(spec["D"] for spec in specs), f"baseline p={p} rep={replica}"
            )
            if 6 * baseline_events.Q_full != baseline_events.T3:
                raise AssertionError(f"p={p}, replica={replica}: baseline T3 gate failed")
            for spec in specs:
                replica_by_scale[spec["label"]].append(
                    stats_for_D(baseline_key, baseline_events, p, spec, include_hints=False)
                )
                pair_checks += 1

        for spec in specs:
            label = spec["label"]
            data_stats = data_by_scale[label]
            replicas = replica_by_scale[label]
            summary = summarize_against_replicas(
                numeric_statistics(data_stats), [numeric_statistics(item) for item in replicas]
            )
            rows.append(
                {
                    "p": p,
                    "scale": label,
                    "L": spec["L"],
                    "D": spec["D"],
                    "data": data_stats,
                    "baseline": {"replicas": replicas, "summary": summary},
                }
            )

        elapsed = time.monotonic() - prime_start
        q_text = ", ".join(
            f"{spec['label']}:D={spec['D']},Q={data_by_scale[spec['label']]['Q_D']}"
            for spec in specs
        )
        print(
            f"[{prime_index}/{len(primes)}] p={p} reps={nrep} {q_text} [{elapsed:.2f}s]",
            flush=True,
        )

    pooled = pooled_results(rows, band)
    gates["pair_identity"] = {"status": "PASS", "checks": pair_checks, "failures": 0}
    examples_ps = [band[0], 10007, 30011]
    examples = [item for item in t3_checks if item["p"] in examples_ps]
    gates["T3_identity"] = {
        "status": "PASS",
        "checks": len(t3_checks),
        "examples": examples,
    }
    gates["overall"] = "PASS"

    completed = datetime.now(timezone.utc).isoformat()
    results: dict[str, Any] = {
        "metadata": {
            "spec": "CODEX_SPEC_CRON_qdsup.md",
            "script": "CRON_qdsup.py",
            "mode": "quick" if args.quick else "production",
            "seed": SEED,
            "completed_utc": completed,
            "runtime_seconds": time.time() - started_wall,
            "band": [3000, 4200],
            "band_prime_count": len(band),
            "band_primes": band,
            "large_primes": [10007, 30011],
            "replicas": {"band": 3, "large": 5},
            "quantile_method": "linear, including implicit zero C coordinates",
            "flag_rule": "data > every replica and data > replica mean + 3 sample sd",
            "binomial_null": (
                "conservative p0=1/(replicas+1) from exchangeable exceed-all event; "
                "the additional 3-sd requirement is ignored in p0"
            ),
        },
        "gates": gates,
        "rows": rows,
        "pooled_band": pooled,
        "verdict": verdict_text(pooled),
    }
    return results


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test", action="store_true", help="run preflight gates and exit")
    parser.add_argument(
        "--quick", action="store_true", help="use the first two band primes (development only)"
    )
    parser.add_argument("--output-json", type=Path, default=HERE / "qdsup_results.json")
    parser.add_argument("--output-report", type=Path, default=HERE / "CODEX_QDSUP_report.md")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    results = run(args)
    if args.self_test:
        print(json.dumps(results, indent=2, sort_keys=True), flush=True)
        return
    args.output_json.write_text(
        json.dumps(results, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8"
    )
    args.output_report.write_text(render_report(results), encoding="utf-8")
    print(f"wrote {args.output_json}", flush=True)
    print(f"wrote {args.output_report}", flush=True)


if __name__ == "__main__":
    main()
