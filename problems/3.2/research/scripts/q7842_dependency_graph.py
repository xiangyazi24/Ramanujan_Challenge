#!/usr/bin/env python3
"""Exact X=512 overlap/dependence graph audit for Q7842.

This imports the repository's independently cross-checked Apéry zero-set
implementation, then computes pair intersections by CRT on several intervals.
Only Python's standard library is used.
"""

from __future__ import annotations

import importlib.util
import itertools
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
HM3_PATH = ROOT / "problems" / "3.2" / "hm3_explore.py"
spec = importlib.util.spec_from_file_location("hm3_explore", HM3_PATH)
assert spec is not None and spec.loader is not None
hm3 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(hm3)

X = 512
M = X * X


def ceil_div(a: int, b: int) -> int:
    return -((-a) // b)


def pair_count(left, right, lo: int, hi: int) -> int:
    """Count m in [lo,hi) hitting both prime zero sets exactly."""
    p, q = left.p, right.p
    inv_p_mod_q = pow(p, -1, q)
    pq = p * q
    count = 0
    for rp in left.zeros:
        for rq in right.zeros:
            a0 = rp + p * (((rq - rp) * inv_p_mod_q) % q)
            k = ceil_div(lo - a0, pq)
            a = a0 + k * pq
            if a < hi:
                count += 1
    return count


def marginal_count(item, lo: int, hi: int) -> int:
    p = item.p
    total = 0
    for r in item.zeros:
        k0 = ceil_div(lo - r, p)
        a = r + k0 * p
        if a < hi:
            total += (hi - 1 - a) // p + 1
    return total


def audit_interval(all_data, active, lo: int, hi: int, label: str) -> None:
    adjacency = {item.p: set() for item in all_data}
    exact_dependence = {item.p: set() for item in all_data}
    pair_hist = Counter()
    active_marginals = {item.p: marginal_count(item, lo, hi) for item in active}
    length = hi - lo

    for left, right in itertools.combinations(active, 2):
        j = pair_count(left, right, lo, hi)
        pair_hist[j] += 1
        if j > 0:
            adjacency[left.p].add(right.p)
            adjacency[right.p].add(left.p)
        # For uniform m in the finite interval, exact independence is
        # J/length = A_p A_q/length^2.
        if j * length != active_marginals[left.p] * active_marginals[right.p]:
            exact_dependence[left.p].add(right.p)
            exact_dependence[right.p].add(left.p)

    overlap_degrees = {p: len(qs) for p, qs in adjacency.items()}
    dep_degrees = {p: len(qs) for p, qs in exact_dependence.items()}
    max_overlap = max(overlap_degrees.values(), default=0)
    max_dep = max(dep_degrees.values(), default=0)
    overlap_max_vertices = sorted(p for p, d in overlap_degrees.items() if d == max_overlap)
    dep_max_vertices = sorted(p for p, d in dep_degrees.items() if d == max_dep)
    active_overlap_degrees = [overlap_degrees[item.p] for item in active]
    active_dep_degrees = [dep_degrees[item.p] for item in active]
    edges = sum(overlap_degrees.values()) // 2
    dep_edges = sum(dep_degrees.values()) // 2
    possible_active_edges = len(active) * (len(active) - 1) // 2

    print(f"INTERVAL {label} lo={lo} hi={hi} length={length}")
    print(f"  overlap_edges={edges} possible_active_edges={possible_active_edges}")
    print(f"  overlap_density_active={edges/possible_active_edges if possible_active_edges else 0:.12f}")
    print(f"  max_overlap_degree_all={max_overlap}")
    print(f"  max_overlap_vertices={overlap_max_vertices}")
    print(f"  active_overlap_degree_hist={sorted(Counter(active_overlap_degrees).items())}")
    print(f"  pair_intersection_hist={sorted(pair_hist.items())}")
    print(f"  exact_dependence_edges={dep_edges}")
    print(f"  exact_dependence_density_active={dep_edges/possible_active_edges if possible_active_edges else 0:.12f}")
    print(f"  max_exact_dependence_degree_all={max_dep}")
    print(f"  max_exact_dependence_vertices={dep_max_vertices}")
    print(f"  active_exact_dependence_degree_hist={sorted(Counter(active_dep_degrees).items())}")


all_data, checked = hm3.build_zero_sets(X, deep_check=True)
active = [item for item in all_data if item.zeros]
print("Q7842_DEPENDENCY_GRAPH_AUDIT")
print(f"X={X} M={M}")
print(f"prime_count={len(all_data)}")
print(f"active_prime_count={len(active)}")
print(f"inactive_prime_count={len(all_data)-len(active)}")
print(f"sum_Z={sum(len(item.zeros) for item in all_data)}")
print(f"deep_recurrence_checks={checked}")
print("active_zero_sizes=" + repr([(item.p, len(item.zeros)) for item in active]))

audit_interval(all_data, active, 0, M, "[0,X^2)")
audit_interval(all_data, active, M, 2*M, "[X^2,2X^2)")

# Independent m-scatter check for the standard HM3 interval.
k_values = hm3.scatter_k(active, M)
print(f"max_K_standard={max(k_values)}")
print(f"K_hist_standard={sorted(Counter(k_values).items())}")
print("Q7842_AUDIT_SUCCESS")
