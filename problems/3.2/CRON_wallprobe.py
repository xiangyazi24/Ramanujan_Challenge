#!/usr/bin/env python3
"""Translated-box and gauge-scan falsification probe for Q6573 sections 17--18, 21.

The authoritative domain is the non-wrapping triangle

    Omega_p = {(r,h): 0 <= r < r+h <= p-2}.

For every requested translated I x J box this program computes the zero-fiber
row/column variances, the all-level coincidence defect in three gauges, and
(when the effective matrix is at most 2000 by 2000) a largest-singular-value
comparison with three Bernoulli replicas.  It writes both machine-readable JSON
and a complete Markdown report.

Default run:

    python3 CRON_wallprobe.py

Optional fourth prime, resuming an already completed default run:

    python3 CRON_wallprobe.py --include-optional --resume
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import platform
import threading
import time
import warnings
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import scipy
from scipy.linalg import svdvals
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import ArpackNoConvergence, LinearOperator, svds


REQUIRED_PRIMES = (1009, 3001, 10007)
OPTIONAL_PRIME = 30011
BASE_RANDOM_SEED = 65731821
SCALE_SPECS = (
    ("p^(1/3)", 1.0 / 3.0, 16),
    ("p^(1/2)", 1.0 / 2.0, 12),
    ("p^(2/3)", 2.0 / 3.0, 8),
)
SVD_DIMENSION_LIMIT = 2000
RANDOM_REPLICAS = 3
DEFAULT_MAX_CHUNK_CELLS = 2_000_000


class Heartbeat:
    """Print a progress line at least every eight seconds while work is active."""

    def __init__(self, interval: float = 8.0) -> None:
        self.interval = interval
        self._status = "initializing"
        self._stop = threading.Event()
        self._lock = threading.Lock()
        self._thread = threading.Thread(target=self._run, daemon=True)

    def __enter__(self) -> "Heartbeat":
        self._thread.start()
        return self

    def __exit__(self, exc_type: Any, exc: Any, tb: Any) -> None:
        self._stop.set()
        self._thread.join(timeout=1.0)

    def set(self, status: str, *, emit: bool = False) -> None:
        with self._lock:
            self._status = status
        if emit:
            print(f"[wallprobe] {status}", flush=True)

    def _run(self) -> None:
        while not self._stop.wait(self.interval):
            with self._lock:
                status = self._status
            print(f"[wallprobe heartbeat] {status}", flush=True)


def round_half_up(x: float) -> int:
    return int(math.floor(x + 0.5))


def evenly_spaced_starts(max_start: int, count: int) -> list[int]:
    """Integer scan locations spanning both endpoints without duplicates."""
    if max_start < 0:
        raise ValueError(f"negative max_start={max_start}")
    if max_start == 0 or count <= 1:
        return [0]
    count = min(count, max_start + 1)
    raw = [round_half_up(i * max_start / (count - 1)) for i in range(count)]
    starts = sorted(set(raw))
    if len(starts) != count:
        # This is only relevant for tiny custom test primes.
        starts = [round_half_up(x) for x in np.linspace(0, max_start, count)]
        starts = sorted(set(starts))
    return starts


def modular_inverse_table(p: int) -> np.ndarray:
    inv = np.zeros(p, dtype=np.int64)
    inv[1] = 1
    for x in range(2, p):
        inv[x] = (p - (p // x) * int(inv[p % x])) % p
    return inv


def apery_pair(p: int) -> tuple[np.ndarray, np.ndarray]:
    """Verified Apéry solution and companion on indices 0,...,p-2.

    This is the recurrence used in CRON_radon_spectrum.py and
    CRON_b1_crosscorr.py:

      (n+1)^3 y_(n+1) = P(n)y_n - n^3 y_(n-1),
      P(n)=34n^3+51n^2+27n+5.
    """
    npoints = p - 1
    b = np.zeros(npoints, dtype=np.int64)
    c = np.zeros(npoints, dtype=np.int64)
    b[0], b[1] = 1, 5 % p
    c[0], c[1] = 0, 1
    inv = modular_inverse_table(p)
    for n in range(1, p - 2):
        pn = (34 * n**3 + 51 * n * n + 27 * n + 5) % p
        inv_cube = pow(int(inv[n + 1]), 3, p)
        b[n + 1] = ((pn * int(b[n]) - pow(n, 3, p) * int(b[n - 1])) * inv_cube) % p
        c[n + 1] = ((pn * int(c[n]) - pow(n, 3, p) * int(c[n - 1])) * inv_cube) % p
    return b, c


def projective_energy(b: np.ndarray, c: np.ndarray, p: int) -> dict[str, Any]:
    inv = modular_inverse_table(p)
    colors = np.empty(p - 1, dtype=np.int64)
    nonzero = b != 0
    colors[~nonzero] = p  # point at infinity in the b != 0 chart
    colors[nonzero] = (c[nonzero] * inv[b[nonzero]]) % p
    multiplicities = np.bincount(colors, minlength=p + 1)
    energy = int(np.dot(multiplicities, multiplicities))
    return {
        "energy": energy,
        "energy_over_p": energy / p,
        "max_multiplicity": int(multiplicities.max()),
        "zero_b_count": int((~nonzero).sum()),
    }


def lc_sequence_mod(p: int, max_h: int) -> np.ndarray:
    """lc(N_h) modulo p: ell_1=1, ell_2=34, ell_(h+1)=34ell_h-ell_(h-1)."""
    lc = np.zeros(max_h + 1, dtype=np.int64)
    if max_h >= 1:
        lc[1] = 1
    if max_h >= 2:
        lc[2] = 34 % p
    for h in range(2, max_h):
        lc[h + 1] = (34 * int(lc[h]) - int(lc[h - 1])) % p
    return lc


def _poly_mul_mod(a: list[int], b: list[int], modulus: int) -> list[int]:
    out = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            out[i + j] = (out[i + j] + ai * bj) % modulus
    return out


def _poly_sub_mod(a: list[int], b: list[int], modulus: int) -> list[int]:
    out = [0] * max(len(a), len(b))
    for i in range(len(a)):
        out[i] = a[i]
    for i in range(len(b)):
        out[i] = (out[i] - b[i]) % modulus
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def _shifted_p_coeffs(d: int, modulus: int) -> list[int]:
    # P(X+d), low coefficient first.
    return [
        (34 * d**3 + 51 * d * d + 27 * d + 5) % modulus,
        (102 * d * d + 102 * d + 27) % modulus,
        (102 * d + 51) % modulus,
        34 % modulus,
    ]


def _shifted_sixth_coeffs(d: int, modulus: int) -> list[int]:
    return [(math.comb(6, k) * d ** (6 - k)) % modulus for k in range(7)]


def verify_lc_recurrence(max_h: int = 12) -> dict[str, Any]:
    """Direct polynomial recurrence check, independent of the scalar lc code."""
    modulus = 1_000_003
    polys: dict[int, list[int]] = {1: [1], 2: _shifted_p_coeffs(1, modulus)}
    scalar = lc_sequence_mod(modulus, max_h)
    checked: list[dict[str, int]] = []
    for h in range(1, max_h + 1):
        if h >= 3 and h not in polys:
            d = h - 1
            polys[h] = _poly_sub_mod(
                _poly_mul_mod(_shifted_p_coeffs(d, modulus), polys[h - 1], modulus),
                _poly_mul_mod(_shifted_sixth_coeffs(d, modulus), polys[h - 2], modulus),
                modulus,
            )
        degree = len(polys[h]) - 1
        leading = polys[h][-1]
        expected_degree = 3 * (h - 1)
        expected_lc = int(scalar[h])
        if degree != expected_degree or leading != expected_lc:
            raise AssertionError(
                f"lc recurrence check failed at h={h}: degree={degree}, lc={leading}, "
                f"expected degree={expected_degree}, lc={expected_lc}"
            )
        checked.append({"h": h, "degree": degree, "lc_mod_1000003": leading})
    return {"modulus": modulus, "max_h": max_h, "checks": checked, "passed": True}


def verify_collision_identity(p: int = 181) -> dict[str, Any]:
    """Exhaustively compare determinant zeros to the N_h continuant at one prime."""
    b, c = apery_pair(p)
    tested = 0
    collisions = 0
    for r in range(p - 2):
        w_prev, w = 0, 1  # N_1(r)=1
        for h in range(1, p - 1 - r):
            if h >= 2:
                d = h - 1
                t = r + d
                pt = (34 * t**3 + 51 * t * t + 27 * t + 5) % p
                w_prev, w = w, (pt * w - pow(t, 6, p) * w_prev) % p
            det = (int(b[r]) * int(c[r + h]) - int(c[r]) * int(b[r + h])) % p
            if (det == 0) != (w == 0):
                raise AssertionError(f"collision cross-check failed at p={p}, r={r}, h={h}")
            tested += 1
            collisions += int(det == 0)
    return {"prime": p, "pairs_tested": tested, "collisions": collisions, "passed": True}


def run_validation() -> dict[str, Any]:
    identity = verify_collision_identity(181)
    lc_check = verify_lc_recurrence(12)
    b, c = apery_pair(1009)
    energy = projective_energy(b, c, 1009)
    if energy["energy"] != 3030:
        raise AssertionError(
            f"banked orbit-energy check changed: p=1009, E={energy['energy']} (expected 3030)"
        )
    return {
        "collision_identity": identity,
        "leading_coefficient_recurrence": lc_check,
        "banked_orbit_energy": {
            "prime": 1009,
            **energy,
            "expected_exact_energy": 3030,
            "passed": True,
        },
    }


@dataclass(frozen=True)
class JSpec:
    label: str
    kind: str
    start: int
    stop: int

    @property
    def size(self) -> int:
        return self.stop - self.start


@dataclass(frozen=True)
class HWindow:
    scale_label: str
    scale_exponent: float
    H: int
    translate_index: int
    translate_count: int
    H0: int

    @property
    def h_lo(self) -> int:
        return self.H0 + 1

    @property
    def h_hi(self) -> int:
        return self.H0 + self.H

    @property
    def h_values(self) -> np.ndarray:
        return np.arange(self.h_lo, self.h_hi + 1, dtype=np.int64)


def make_h_windows(p: int, translate_counts: tuple[int, int, int]) -> tuple[int, list[HWindow]]:
    scan_limit = min(p - 2, round_half_up(4.0 * p ** (2.0 / 3.0)))
    windows: list[HWindow] = []
    for (label, exponent, default_count), requested_count in zip(SCALE_SPECS, translate_counts):
        H = max(1, round_half_up(p**exponent))
        count = requested_count if requested_count > 0 else default_count
        starts = evenly_spaced_starts(scan_limit - H, count)
        for i, h0 in enumerate(starts):
            windows.append(HWindow(label, exponent, H, i, len(starts), h0))
    return scan_limit, windows


def make_j_specs(p: int, translated_count: int) -> list[JSpec]:
    x_size = p - 1
    local_size = min(x_size, round_half_up(p / 4.0))
    starts = evenly_spaced_starts(x_size - local_size, translated_count)
    specs = [JSpec("full", "full", 0, x_size)]
    specs.extend(
        JSpec(f"t{i:02d}", "translated", start, start + local_size)
        for i, start in enumerate(starts)
    )
    return specs


def deterministic_seed(box_id: str, replica: int) -> int:
    payload = f"{BASE_RANDOM_SEED}:{box_id}:{replica}".encode("ascii")
    return int.from_bytes(hashlib.sha256(payload).digest()[:4], "little")


class CenteredFerrersOperator(LinearOperator):
    """Sparse zero incidences minus (1/p) times a triangular/Ferrers mask.

    The boxes are intersections of intervals with r+h<=p-2.  After inactive
    outer rows and columns are removed, every row of the validity mask is a
    prefix.  Prefix/suffix sums therefore apply the dense centering term in
    linear time, while the 1[Delta=0] term remains a very sparse CSR matrix.
    """

    def __init__(self, zero_matrix: csr_matrix, row_cutoffs: np.ndarray, p: int) -> None:
        self.zero_matrix = zero_matrix
        self.row_cutoffs = np.asarray(row_cutoffs, dtype=np.int64)
        self.p = p
        super().__init__(dtype=np.dtype(np.float64), shape=zero_matrix.shape)

    def _matvec(self, vector: np.ndarray) -> np.ndarray:
        vector = np.asarray(vector).reshape(-1)
        prefix = np.cumsum(vector)
        mask_product = np.zeros(self.shape[0], dtype=np.float64)
        nonempty = self.row_cutoffs > 0
        mask_product[nonempty] = prefix[self.row_cutoffs[nonempty] - 1]
        sparse_product = np.asarray(self.zero_matrix @ vector).reshape(-1)
        return sparse_product - mask_product / self.p

    def _rmatvec(self, vector: np.ndarray) -> np.ndarray:
        vector = np.asarray(vector).reshape(-1)
        buckets = np.bincount(
            self.row_cutoffs, weights=vector, minlength=self.shape[1] + 1
        )
        # Column j is valid in precisely the rows whose cutoff is >= j+1.
        mask_product = np.cumsum(buckets[::-1])[::-1][1 : self.shape[1] + 1]
        sparse_product = np.asarray(self.zero_matrix.T @ vector).reshape(-1)
        return sparse_product - mask_product / self.p


def largest_singular_value(operator: CenteredFerrersOperator, seed: int) -> tuple[float, str]:
    if min(operator.shape) == 1:
        if operator.shape[1] == 1:
            value = np.linalg.norm(operator.matvec(np.ones(1)))
        else:
            value = np.linalg.norm(operator.rmatvec(np.ones(1)))
        return float(value), "norm"
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            value = svds(
                operator,
                k=1,
                which="LM",
                return_singular_vectors=False,
                solver="arpack",
                tol=1e-7,
                maxiter=max(1000, 20 * min(operator.shape)),
                ncv=min(10, min(operator.shape) - 1),
                random_state=seed,
            )
        return float(value[-1]), "scipy.sparse.linalg.svds(LinearOperator)"
    except (ArpackNoConvergence, ValueError, np.linalg.LinAlgError):
        identity = np.eye(operator.shape[1], dtype=np.float64)
        dense = operator.matmat(identity)
        return float(svdvals(dense, check_finite=False)[0]), "scipy.linalg.svdvals-fallback"


class BoxAccumulator:
    def __init__(
        self,
        p: int,
        window: HWindow,
        jspec: JSpec,
        h_values: np.ndarray,
        lc_nonzero: np.ndarray,
    ) -> None:
        self.p = p
        self.window = window
        self.jspec = jspec
        self.h_values = h_values
        self.lc_nonzero = lc_nonzero
        self.box_id = (
            f"p{p}-{window.scale_label.replace('^', '').replace('/', '_').replace('(', '').replace(')', '')}"
            f"-H{window.H}-i{window.translate_index:02d}-J{jspec.label}"
        )
        m = len(h_values)
        self.v0 = 0.0
        self.col_zero = np.zeros(m, dtype=np.int64)
        self.col_q = np.zeros(m, dtype=np.int64)
        self.hist = {
            "g1_raw": np.zeros(p, dtype=np.int64),
            "g2_lc": np.zeros(p, dtype=np.int64),
            "g3_chart": np.zeros(p, dtype=np.int64),
        }
        self.omega = {"g1_raw": 0, "g2_lc": 0, "g3_chart": 0}
        self.zero_count = 0
        self.max_abs_row_residual = -1.0
        self.max_row_r: int | None = None
        self.max_row_zero_count = 0
        self.max_row_q = 0

        # Rows/columns outside the triangular support are exactly zero and can be
        # removed without changing singular values.
        active_r_stop = min(jspec.stop, p - 1 - int(h_values.min()))
        self.svd_r_start = jspec.start
        self.svd_r_stop = max(jspec.start, active_r_stop)
        self.svd_h_mask = h_values <= (p - 2 - jspec.start)
        active_shape = (self.svd_r_stop - self.svd_r_start, int(self.svd_h_mask.sum()))
        self.svd_active_shape = active_shape
        self.svd_eligible = (
            active_shape[0] > 0
            and active_shape[1] > 0
            and active_shape[0] <= SVD_DIMENSION_LIMIT
            and active_shape[1] <= SVD_DIMENSION_LIMIT
        )
        self.svd_zero = np.zeros(active_shape, dtype=bool) if self.svd_eligible else None
        self.svd_valid = np.zeros(active_shape, dtype=bool) if self.svd_eligible else None

    def update(
        self,
        chunk_start: int,
        chunk_stop: int,
        delta: np.ndarray,
        valid: np.ndarray,
        delta_g2: np.ndarray,
        valid_g2: np.ndarray,
        delta_g3: np.ndarray,
        valid_g3: np.ndarray,
    ) -> None:
        lo = max(chunk_start, self.jspec.start)
        hi = min(chunk_stop, self.jspec.stop)
        if lo >= hi:
            return
        left, right = lo - chunk_start, hi - chunk_start
        d = delta[:, left:right]
        v = valid[:, left:right]
        zeros = (d == 0) & v

        row_zero = zeros.sum(axis=0, dtype=np.int64)
        row_q = v.sum(axis=0, dtype=np.int64)
        residual = row_zero.astype(np.float64) - row_q.astype(np.float64) / self.p
        self.v0 += float(np.dot(residual, residual))
        if residual.size:
            idx = int(np.argmax(np.abs(residual)))
            abs_resid = abs(float(residual[idx]))
            if abs_resid > self.max_abs_row_residual:
                self.max_abs_row_residual = abs_resid
                self.max_row_r = lo + idx
                self.max_row_zero_count = int(row_zero[idx])
                self.max_row_q = int(row_q[idx])

        self.col_zero += zeros.sum(axis=1, dtype=np.int64)
        self.col_q += v.sum(axis=1, dtype=np.int64)
        nvalid = int(v.sum())
        self.omega["g1_raw"] += nvalid
        self.zero_count += int(zeros.sum())
        self.hist["g1_raw"] += np.bincount(d[v], minlength=self.p)

        vg2 = valid_g2[:, left:right]
        dg2 = delta_g2[:, left:right]
        self.omega["g2_lc"] += int(vg2.sum())
        self.hist["g2_lc"] += np.bincount(dg2[vg2], minlength=self.p)

        vg3 = valid_g3[:, left:right]
        dg3 = delta_g3[:, left:right]
        self.omega["g3_chart"] += int(vg3.sum())
        self.hist["g3_chart"] += np.bincount(dg3[vg3], minlength=self.p)

        if self.svd_eligible and self.svd_zero is not None and self.svd_valid is not None:
            slo = max(lo, self.svd_r_start)
            shi = min(hi, self.svd_r_stop)
            if slo < shi:
                source_left = slo - chunk_start
                source_right = shi - chunk_start
                target_left = slo - self.svd_r_start
                target_right = shi - self.svd_r_start
                hmask = self.svd_h_mask
                sv = valid[hmask, source_left:source_right]
                sz = (delta[hmask, source_left:source_right] == 0) & sv
                self.svd_zero[target_left:target_right, :] = sz.T
                self.svd_valid[target_left:target_right, :] = sv.T

    def _gauge_result(self, name: str) -> dict[str, Any]:
        omega = int(self.omega[name])
        if omega == 0:
            return {
                "omega": 0,
                "coverage_of_raw": 0.0,
                "n_coinc": 0,
                "p_times_defect": 0,
                "defect": 0.0,
                "defect_over_omega": None,
            }
        histogram = self.hist[name]
        n_coinc = int(np.dot(histogram, histogram))
        p_times_defect = self.p * n_coinc - omega * omega
        # Cauchy guarantees nonnegative; fail loudly if arithmetic/overflow changed.
        if p_times_defect < 0:
            raise AssertionError(f"negative all-level defect in {self.box_id}/{name}")
        return {
            "omega": omega,
            "coverage_of_raw": omega / max(1, self.omega["g1_raw"]),
            "n_coinc": n_coinc,
            "p_times_defect": int(p_times_defect),
            "defect": p_times_defect / self.p,
            "defect_over_omega": p_times_defect / (self.p * omega),
        }

    def _svd_result(self) -> dict[str, Any]:
        if not self.svd_eligible or self.svd_zero is None or self.svd_valid is None:
            return {
                "computed": False,
                "dimension_limit": SVD_DIMENSION_LIMIT,
                "nominal_shape": [self.jspec.size, len(self.h_values)],
                "active_shape": list(self.svd_active_shape),
                "reason": "effective support exceeds the 2000-by-2000 protocol limit",
            }
        observed_density = self.zero_count / max(1, self.omega["g1_raw"])
        # B is centered at the wall's reference density 1/p, so the matching
        # iid benchmark is Bernoulli(1/p).  Using the empirical density makes
        # three-replica comparisons unstable in one-hit boxes: all replicas can
        # be empty and create an artificial near-zero denominator.
        benchmark_density = 1.0 / self.p
        row_cutoffs = self.svd_valid.sum(axis=1, dtype=np.int64)
        expected_mask = np.arange(self.svd_valid.shape[1])[None, :] < row_cutoffs[:, None]
        if not np.array_equal(self.svd_valid, expected_mask):
            raise AssertionError(f"non-Ferrers SVD mask in {self.box_id}")
        data_operator = CenteredFerrersOperator(
            csr_matrix(self.svd_zero.astype(np.float64)), row_cutoffs, self.p
        )
        data_seed = deterministic_seed(self.box_id, -1)
        data_value, data_method = largest_singular_value(data_operator, data_seed)
        random_values: list[float] = []
        random_methods: list[str] = []
        for replica in range(RANDOM_REPLICAS):
            seed = deterministic_seed(self.box_id, replica)
            rng = np.random.default_rng(seed)
            draw = (rng.random(self.svd_zero.shape) < benchmark_density) & self.svd_valid
            # The LinearOperator applies exactly the same 1/p centering as in
            # the data matrix.  Thus even a zero-density box has a fair baseline.
            benchmark_operator = CenteredFerrersOperator(
                csr_matrix(draw.astype(np.float64)), row_cutoffs, self.p
            )
            value, method = largest_singular_value(benchmark_operator, seed)
            random_values.append(value)
            random_methods.append(method)
        benchmark_mean = float(np.mean(random_values))
        ratio = data_value / benchmark_mean if benchmark_mean > 0 else None
        return {
            "computed": True,
            "dimension_limit": SVD_DIMENSION_LIMIT,
            "nominal_shape": [self.jspec.size, len(self.h_values)],
            "active_shape": list(self.svd_active_shape),
            "observed_zero_density": observed_density,
            "benchmark_bernoulli_probability": benchmark_density,
            "centering": "1/p for data and random replicas",
            "data_singular_value": data_value,
            "data_method": data_method,
            "random_singular_values": random_values,
            "random_methods": random_methods,
            "random_mean": benchmark_mean,
            "ratio_data_over_random_mean": ratio,
        }

    def finalize(self) -> dict[str, Any]:
        h_size = len(self.h_values)
        j_size = self.jspec.size
        col_residual = self.col_zero.astype(np.float64) - self.col_q.astype(np.float64) / self.p
        vcol = float(np.dot(col_residual, col_residual))
        raw_omega = self.omega["g1_raw"]
        return {
            "box_id": self.box_id,
            "p": self.p,
            "scale": self.window.scale_label,
            "scale_exponent": self.window.scale_exponent,
            "H": self.window.H,
            "translate_index": self.window.translate_index,
            "translate_count": self.window.translate_count,
            "H0": self.window.H0,
            "I": [self.window.h_lo, self.window.h_hi],
            "J_label": self.jspec.label,
            "J_kind": self.jspec.kind,
            "J": [self.jspec.start, self.jspec.stop - 1],
            "J_size": j_size,
            "omega": int(raw_omega),
            "zero_count": int(self.zero_count),
            "zero_density": self.zero_count / max(1, raw_omega),
            "V0": self.v0,
            "V0_over_I": self.v0 / h_size,
            "V_col": vcol,
            "V_col_over_J": vcol / j_size,
            "row_extreme": {
                "r": self.max_row_r,
                "zero_count": self.max_row_zero_count,
                "q": self.max_row_q,
                "absolute_residual": self.max_abs_row_residual,
            },
            "lc_zero_h_count": int((~self.lc_nonzero).sum()),
            "gauges": {
                "g1_raw": self._gauge_result("g1_raw"),
                "g2_lc": self._gauge_result("g2_lc"),
                "g3_chart": self._gauge_result("g3_chart"),
            },
            "svd": self._svd_result(),
        }


def compute_window_boxes(
    p: int,
    b: np.ndarray,
    c: np.ndarray,
    inv_residue: np.ndarray,
    inv_b: np.ndarray,
    lc: np.ndarray,
    window: HWindow,
    jspecs: list[JSpec],
    max_chunk_cells: int,
    heartbeat: Heartbeat,
) -> list[dict[str, Any]]:
    h_values = window.h_values
    lc_values = lc[h_values]
    lc_nonzero = lc_values != 0
    inv_lc = np.zeros_like(lc_values)
    inv_lc[lc_nonzero] = inv_residue[lc_values[lc_nonzero]]
    accumulators = [BoxAccumulator(p, window, jspec, h_values, lc_nonzero) for jspec in jspecs]

    chunk_width = max(1, max_chunk_cells // len(h_values))
    x_stop = p - 1
    for chunk_start in range(0, x_stop, chunk_width):
        chunk_stop = min(x_stop, chunk_start + chunk_width)
        heartbeat.set(
            f"p={p}, I=[{window.h_lo},{window.h_hi}], r-chunk=[{chunk_start},{chunk_stop})"
        )
        r = np.arange(chunk_start, chunk_stop, dtype=np.int64)
        s = h_values[:, None] + r[None, :]
        valid = s <= p - 2
        safe_s = np.minimum(s, p - 2)

        br = b[r][None, :]
        cr = c[r][None, :]
        bs = b[safe_s]
        cs = c[safe_s]
        delta = (br * cs - cr * bs) % p

        delta_g2 = (delta * inv_lc[:, None]) % p
        valid_g2 = valid & lc_nonzero[:, None]

        chart = valid & (br != 0) & (bs != 0)
        delta_g3 = (delta * inv_b[r][None, :]) % p
        delta_g3 = (delta_g3 * inv_b[safe_s]) % p

        for accumulator in accumulators:
            accumulator.update(
                chunk_start,
                chunk_stop,
                delta,
                valid,
                delta_g2,
                valid_g2,
                delta_g3,
                chart,
            )

    results: list[dict[str, Any]] = []
    for i, accumulator in enumerate(accumulators, start=1):
        heartbeat.set(
            f"p={p}, I=[{window.h_lo},{window.h_hi}], finalizing box {i}/{len(accumulators)}",
            emit=False,
        )
        results.append(accumulator.finalize())
    return results


def metric_summary(boxes: list[dict[str, Any]], getter: Any) -> dict[str, Any] | None:
    pairs = [(float(value), box["box_id"]) for box in boxes if (value := getter(box)) is not None]
    if not pairs:
        return None
    values = np.array([value for value, _ in pairs], dtype=np.float64)
    max_value, max_box = max(pairs, key=lambda pair: pair[0])
    return {
        "min": float(values.min()),
        "median": float(np.median(values)),
        "max": max_value,
        "argmax_box": max_box,
    }


def summarize_prime(p_result: dict[str, Any]) -> list[dict[str, Any]]:
    summaries: list[dict[str, Any]] = []
    boxes = p_result["boxes"]
    for scale in [spec[0] for spec in SCALE_SPECS]:
        selected = [box for box in boxes if box["scale"] == scale]
        if not selected:
            continue
        summaries.append(
            {
                "scale": scale,
                "H": selected[0]["H"],
                "h_translates": selected[0]["translate_count"],
                "box_count": len(selected),
                "V0_over_I": metric_summary(selected, lambda b: b["V0_over_I"]),
                "V_col_over_J": metric_summary(selected, lambda b: b["V_col_over_J"]),
                "g1_raw": metric_summary(
                    selected, lambda b: b["gauges"]["g1_raw"]["defect_over_omega"]
                ),
                "g2_lc": metric_summary(
                    selected, lambda b: b["gauges"]["g2_lc"]["defect_over_omega"]
                ),
                "g3_chart": metric_summary(
                    selected, lambda b: b["gauges"]["g3_chart"]["defect_over_omega"]
                ),
                "svd_ratio": metric_summary(
                    selected,
                    lambda b: b["svd"].get("ratio_data_over_random_mean")
                    if b["svd"]["computed"]
                    else None,
                ),
                "svd_boxes_computed": sum(box["svd"]["computed"] for box in selected),
            }
        )
    return summaries


def compute_prime(
    p: int,
    translate_counts: tuple[int, int, int],
    j_translates: int,
    max_chunk_cells: int,
    heartbeat: Heartbeat,
) -> dict[str, Any]:
    started = time.perf_counter()
    heartbeat.set(f"p={p}: generating Apéry orbit", emit=True)
    b, c = apery_pair(p)
    orbit = projective_energy(b, c, p)
    inv_residue = modular_inverse_table(p)
    inv_b = np.zeros_like(b)
    nonzero_b = b != 0
    inv_b[nonzero_b] = inv_residue[b[nonzero_b]]

    scan_limit, windows = make_h_windows(p, translate_counts)
    jspecs = make_j_specs(p, j_translates)
    lc = lc_sequence_mod(p, scan_limit)
    boxes: list[dict[str, Any]] = []
    for index, window in enumerate(windows, start=1):
        heartbeat.set(
            f"p={p}: h-window {index}/{len(windows)} I=[{window.h_lo},{window.h_hi}]",
            emit=True,
        )
        boxes.extend(
            compute_window_boxes(
                p,
                b,
                c,
                inv_residue,
                inv_b,
                lc,
                window,
                jspecs,
                max_chunk_cells,
                heartbeat,
            )
        )
    elapsed = time.perf_counter() - started
    result = {
        "p": p,
        "scan_limit": scan_limit,
        "orbit": orbit,
        "lc_zero_h": [h for h in range(1, scan_limit + 1) if int(lc[h]) == 0],
        "h_windows": [
            {
                "scale": window.scale_label,
                "H": window.H,
                "translate_index": window.translate_index,
                "translate_count": window.translate_count,
                "H0": window.H0,
                "I": [window.h_lo, window.h_hi],
            }
            for window in windows
        ],
        "j_windows": [
            {
                "label": jspec.label,
                "kind": jspec.kind,
                "J": [jspec.start, jspec.stop - 1],
                "size": jspec.size,
            }
            for jspec in jspecs
        ],
        "boxes": boxes,
        "elapsed_seconds": elapsed,
    }
    result["scale_summaries"] = summarize_prime(result)
    heartbeat.set(f"p={p}: complete in {elapsed:.1f}s ({len(boxes)} boxes)", emit=True)
    return result


def fmt(x: Any, digits: int = 4) -> str:
    if x is None:
        return "--"
    if isinstance(x, int):
        return str(x)
    return f"{float(x):.{digits}g}"


def short_box(box_id: str | None) -> str:
    return "--" if box_id is None else f"`{box_id}`"


def all_boxes(results: dict[str, Any]) -> list[dict[str, Any]]:
    return [box for p_result in results["primes"] for box in p_result["boxes"]]


def compute_verdicts(results: dict[str, Any]) -> list[dict[str, str]]:
    boxes = all_boxes(results)
    max_v0 = max(boxes, key=lambda b: b["V0_over_I"])
    svd_boxes = [
        box
        for box in boxes
        if box["svd"]["computed"]
        and box["svd"].get("ratio_data_over_random_mean") is not None
    ]
    max_svd = max(
        svd_boxes,
        key=lambda b: b["svd"]["ratio_data_over_random_mean"],
        default=None,
    )
    gauge_max: dict[str, tuple[float, dict[str, Any]]] = {}
    for gauge in ("g1_raw", "g2_lc", "g3_chart"):
        candidates = [
            (box["gauges"][gauge]["defect_over_omega"], box)
            for box in boxes
            if box["gauges"][gauge]["defect_over_omega"] is not None
        ]
        gauge_max[gauge] = max(candidates, key=lambda item: item[0])

    weighted_detected = (
        max_svd is not None and max_svd["svd"]["ratio_data_over_random_mean"] > 3.0
    )
    all_level_detected = any(value > 3.0 for value, _ in gauge_max.values())
    localized_detected = max_v0["V0_over_I"] > 3.0
    return [
        {
            "mode": "localized heavy row",
            "verdict": "detected at a scanned box" if localized_detected else "not detected at these scales",
            "evidence": (
                f"max V0/|I|={max_v0['V0_over_I']:.4g} in {max_v0['box_id']}; "
                f"extreme row r={max_v0['row_extreme']['r']} has "
                f"d={max_v0['row_extreme']['zero_count']} of q={max_v0['row_extreme']['q']}"
            ),
        },
        {
            "mode": "thin resonance class",
            "verdict": "detected spectrally" if weighted_detected else "not detected at these scales",
            "evidence": (
                "no eligible SVD boxes"
                if max_svd is None
                else f"largest data/random singular ratio="
                f"{max_svd['svd']['ratio_data_over_random_mean']:.4g} in {max_svd['box_id']}"
            ),
        },
        {
            "mode": "nonzero-level concentration / gauge dependence",
            "verdict": "detected at a scanned box" if all_level_detected else "not detected at these scales",
            "evidence": "; ".join(
                f"{gauge} max={value:.4g} ({box['box_id']})"
                for gauge, (value, box) in gauge_max.items()
            ),
        },
        {
            "mode": "weighted singular vector",
            "verdict": "detected (ratio > 3)" if weighted_detected else "not detected (no ratio > 3)",
            "evidence": (
                "no eligible SVD boxes"
                if max_svd is None
                else f"max ratio={max_svd['svd']['ratio_data_over_random_mean']:.4g}; "
                f"data={max_svd['svd']['data_singular_value']:.4g}, "
                f"benchmark={max_svd['svd']['random_mean']:.4g}"
            ),
        },
    ]


def render_report(results: dict[str, Any]) -> str:
    lines: list[str] = []
    lines.append("# Translated-box wall falsification probe")
    lines.append("")
    lines.append(
        "This report executes the translated-window, all-level gauge, and operator-norm "
        "protocol of Q6573 sections 17--18 and 21 on the exact non-wrapping triangle "
        "`0 <= r < r+h <= p-2`."
    )
    lines.append("")
    lines.append("## Conventions and validation")
    lines.append("")
    lines.append(
        "For `I=(H0,H0+H]`, `d_I(r)` counts zero determinants among admissible h and "
        "`q_I(r)` counts the corresponding admissible pairs. For a column h, `J_h` is "
        "the subset of J remaining in the triangle. The 8 translated J-windows have "
        "nearest-integer length p/4 and evenly spaced starts; h-window starts span the "
        "full interval from 1 through the rounded `4 p^(2/3)` endpoint."
    )
    lines.append("")
    validation = results["validation"]
    energy = validation["banked_orbit_energy"]
    collision = validation["collision_identity"]
    lines.append(
        f"- Orbit code cross-check: all {collision['pairs_tested']} pairs at p={collision['prime']} "
        f"satisfy `Delta=0 iff N_h(r)=0` ({collision['collisions']} collisions)."
    )
    lines.append(
        f"- Banked statistic: at p={energy['prime']}, `E^pi={energy['energy']}` and "
        f"`E^pi/p={energy['energy_over_p']:.6f}` (exact expected energy {energy['expected_exact_energy']})."
    )
    lines.append(
        "- The direct polynomial recurrence through h=12 agrees with "
        "`lc(N_1)=1`, `lc(N_2)=34`, and `lc(N_(h+1))=34 lc(N_h)-lc(N_(h-1))`."
    )
    lines.append(
        "- Gauge g2 omits an h-column when its leading coefficient is 0 modulo p; gauge g3 "
        "uses only pairs where both b-coordinates are nonzero. Each gauge is normalized by "
        "its own domain size, and the retained coverage is reported."
    )
    lines.append(
        "- SVD matrices set entries outside the triangular domain to zero and discard all-zero "
        "outer rows/columns. Eligible effective shapes are at most 2000 by 2000. Each random "
        "replica is iid Bernoulli with the wall reference mean `1/p` on the same mask and "
        "uses the same `1/p` centering."
    )
    lines.append("")

    lines.append("## Maxima by prime and scale")
    lines.append("")
    lines.append(
        "| p | scale | H | h scans | boxes | max V0/H | max Vcol/|J| | "
        "max g1 | max g2 | max g3 | SVD boxes | max SVD ratio |"
    )
    lines.append("|---:|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for p_result in results["primes"]:
        for summary in p_result["scale_summaries"]:
            lines.append(
                f"| {p_result['p']} | {summary['scale']} | {summary['H']} | "
                f"{summary['h_translates']} | {summary['box_count']} | "
                f"{fmt(summary['V0_over_I']['max'])} | {fmt(summary['V_col_over_J']['max'])} | "
                f"{fmt(summary['g1_raw']['max'])} | {fmt(summary['g2_lc']['max'])} | "
                f"{fmt(summary['g3_chart']['max'])} | {summary['svd_boxes_computed']} | "
                f"{fmt(summary['svd_ratio']['max'] if summary['svd_ratio'] else None)} |"
            )
    lines.append("")

    lines.append("## Gauge flatness across all translations")
    lines.append("")
    lines.append(
        "A random uniform level map has expected normalized defect close to 1. "
        "The table excludes empty triangular intersections; minimum coverage is the "
        "smallest retained gauge domain relative to the raw domain."
    )
    lines.append("")
    lines.append("| p | gauge | boxes | min | median | max | minimum coverage |")
    lines.append("|---:|:---|---:|---:|---:|---:|---:|")
    for p_result in results["primes"]:
        for gauge in ("g1_raw", "g2_lc", "g3_chart"):
            entries = [
                box["gauges"][gauge]
                for box in p_result["boxes"]
                if box["gauges"][gauge]["defect_over_omega"] is not None
            ]
            values = np.array([entry["defect_over_omega"] for entry in entries])
            minimum_coverage = min(entry["coverage_of_raw"] for entry in entries)
            lines.append(
                f"| {p_result['p']} | {gauge} | {len(entries)} | {values.min():.4g} | "
                f"{np.median(values):.4g} | {values.max():.4g} | {minimum_coverage:.4f} |"
            )
    lines.append("")

    empty_count = sum(box["omega"] == 0 for box in all_boxes(results))
    if empty_count:
        lines.append(
            f"The fixed translation grid contains {empty_count} empty boxes: for these, the "
            "rightmost J-window and a large translated gap-window do not intersect the "
            "non-wrapping triangle. They are retained in the complete table with undefined "
            "all-level and SVD entries and are excluded from extrema."
        )
        lines.append("")

    lines.append("## Argmax boxes")
    lines.append("")
    for p_result in results["primes"]:
        lines.append(f"### p={p_result['p']}")
        lines.append("")
        lines.append(
            f"Orbit `E^pi/p={p_result['orbit']['energy_over_p']:.6f}`; "
            f"scan endpoint h={p_result['scan_limit']}; leading-coefficient-zero gaps: "
            f"{p_result['lc_zero_h'] if p_result['lc_zero_h'] else 'none'}."
        )
        lines.append("")
        lines.append("| scale | metric | maximum | argmax box |")
        lines.append("|:---|:---|---:|:---|")
        for summary in p_result["scale_summaries"]:
            for label, key in (
                ("V0/H", "V0_over_I"),
                ("Vcol/|J|", "V_col_over_J"),
                ("all-level g1", "g1_raw"),
                ("all-level g2", "g2_lc"),
                ("all-level g3", "g3_chart"),
                ("SVD ratio", "svd_ratio"),
            ):
                item = summary[key]
                lines.append(
                    f"| {summary['scale']} | {label} | "
                    f"{fmt(item['max'] if item else None)} | "
                    f"{short_box(item['argmax_box'] if item else None)} |"
                )
        lines.append("")

    lines.append("## Failure-mode verdicts")
    lines.append("")
    lines.append(
        "The finite-scan detection thresholds used here are `V0/|I| > 3` for a localized "
        "row anomaly, all-level defect/|Omega| > 3 for level concentration, and the "
        "protocol's explicit data/random SVD ratio > 3 for a weighted mode. Values below "
        "these thresholds are evidence only at the scanned primes and translations."
    )
    lines.append("")
    lines.append("| failure mode | verdict | numerical evidence |")
    lines.append("|:---|:---|:---|")
    for verdict in results["verdicts"]:
        lines.append(
            f"| {verdict['mode']} | {verdict['verdict']} | {verdict['evidence']} |"
        )
    lines.append("")

    lines.append("## Complete per-box table")
    lines.append("")
    lines.append(
        "`cov2` and `cov3` are the g2/g3 domain fractions relative to raw Omega. "
        "The SVD cell is `data/random-mean (ratio)`; `--` means the effective support "
        "exceeded the 2000-by-2000 limit."
    )
    lines.append("")
    for p_result in results["primes"]:
        lines.append(f"### p={p_result['p']}: every scanned box")
        lines.append("")
        lines.append(
            "| scale | I | J | |Omega| | V0/H | Vcol/|J| | g1 | g2 | cov2 | g3 | cov3 | SVD |"
        )
        lines.append("|:---|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|:---|")
        for box in p_result["boxes"]:
            gauges = box["gauges"]
            svd = box["svd"]
            if svd["computed"]:
                svd_cell = (
                    f"{svd['data_singular_value']:.4g}/{svd['random_mean']:.4g} "
                    f"({svd['ratio_data_over_random_mean']:.3g})"
                )
            else:
                svd_cell = "--"
            lines.append(
                f"| {box['scale']} | {box['I'][0]}--{box['I'][1]} | "
                f"{box['J_label']}:{box['J'][0]}--{box['J'][1]} | {box['omega']} | "
                f"{box['V0_over_I']:.5g} | {box['V_col_over_J']:.5g} | "
                f"{fmt(gauges['g1_raw']['defect_over_omega'], 5)} | "
                f"{fmt(gauges['g2_lc']['defect_over_omega'], 5)} | "
                f"{gauges['g2_lc']['coverage_of_raw']:.3f} | "
                f"{fmt(gauges['g3_chart']['defect_over_omega'], 5)} | "
                f"{gauges['g3_chart']['coverage_of_raw']:.3f} | {svd_cell} |"
            )
        lines.append("")

    lines.append("## Reproduction")
    lines.append("")
    primes = ",".join(str(p_result["p"]) for p_result in results["primes"])
    lines.append("```bash")
    lines.append(f"python3 CRON_wallprobe.py --primes {primes}")
    lines.append("```")
    lines.append("")
    lines.append(
        f"Runtime recorded in JSON: "
        f"{sum(p_result['elapsed_seconds'] for p_result in results['primes']):.1f} seconds "
        "of per-prime wall time on this machine."
    )
    lines.append("")
    return "\n".join(lines)


def write_json(path: Path, payload: dict[str, Any]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=False) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def write_report(path: Path, report: str) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(report, encoding="utf-8")
    os.replace(temporary, path)


def parse_primes(text: str | None, include_optional: bool) -> tuple[int, ...]:
    if text:
        primes = tuple(int(token.strip()) for token in text.split(",") if token.strip())
    else:
        primes = REQUIRED_PRIMES + ((OPTIONAL_PRIME,) if include_optional else ())
    if not primes:
        raise ValueError("empty prime list")
    return primes


def parse_translate_counts(text: str) -> tuple[int, int, int]:
    values = tuple(int(token.strip()) for token in text.split(","))
    if len(values) != 3 or any(value < 1 for value in values):
        raise ValueError("--h-translate-counts requires three positive integers")
    return values  # type: ignore[return-value]


def base_payload(args: argparse.Namespace, validation: dict[str, Any]) -> dict[str, Any]:
    return {
        "schema_version": 2,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "protocol": "Q6573 sections 17--18 and 21 / CODEX_SPEC_CRON_wallprobe.md",
        "configuration": {
            "required_primes": list(REQUIRED_PRIMES),
            "optional_prime": OPTIONAL_PRIME,
            "random_seed": BASE_RANDOM_SEED,
            "random_replicas": RANDOM_REPLICAS,
            "svd_dimension_limit": SVD_DIMENSION_LIMIT,
            "svd_bernoulli_probability": "1/p",
            "h_translate_counts": list(args.h_translate_counts),
            "j_translates": args.j_translates,
            "max_chunk_cells": args.max_chunk_cells,
            "rounding": "nearest integer, half up",
            "numpy_version": np.__version__,
            "scipy_version": scipy.__version__,
            "python_version": platform.python_version(),
            "platform": platform.platform(),
        },
        "validation": validation,
        "primes": [],
        "verdicts": [],
    }


def compatible_resume(existing: dict[str, Any], args: argparse.Namespace) -> bool:
    config = existing.get("configuration", {})
    return (
        config.get("random_seed") == BASE_RANDOM_SEED
        and config.get("random_replicas") == RANDOM_REPLICAS
        and config.get("svd_dimension_limit") == SVD_DIMENSION_LIMIT
        and config.get("svd_bernoulli_probability") == "1/p"
        and config.get("h_translate_counts") == list(args.h_translate_counts)
        and config.get("j_translates") == args.j_translates
        and config.get("max_chunk_cells") == args.max_chunk_cells
    )


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primes", help="comma-separated prime list (default: required primes)")
    parser.add_argument("--include-optional", action="store_true", help="add p=30011")
    parser.add_argument(
        "--h-translate-counts",
        default="16,12,8",
        help="numbers of h-window scan locations at the three scales",
    )
    parser.add_argument("--j-translates", type=int, default=8)
    parser.add_argument("--max-chunk-cells", type=int, default=DEFAULT_MAX_CHUNK_CELLS)
    parser.add_argument("--results", type=Path, default=Path("wallprobe_results.json"))
    parser.add_argument("--report", type=Path, default=Path("CODEX_WALLPROBE_report.md"))
    parser.add_argument("--resume", action="store_true", help="reuse matching primes in existing JSON")
    args = parser.parse_args(list(argv) if argv is not None else None)
    args.h_translate_counts = parse_translate_counts(args.h_translate_counts)
    if args.j_translates < 1:
        parser.error("--j-translates must be positive")
    if args.max_chunk_cells < 1:
        parser.error("--max-chunk-cells must be positive")
    primes = parse_primes(args.primes, args.include_optional)

    for path in (args.results, args.report):
        if path.parent != Path(".") and not path.parent.exists():
            parser.error(f"output directory does not exist: {path.parent}")

    print("[wallprobe] validating recurrence conventions and banked statistics", flush=True)
    validation = run_validation()
    if args.resume and args.results.exists():
        existing = json.loads(args.results.read_text(encoding="utf-8"))
        if not compatible_resume(existing, args):
            parser.error("existing result JSON is incompatible with the requested scan configuration")
        payload = existing
        payload["validation"] = validation
    else:
        payload = base_payload(args, validation)

    completed = {item["p"] for item in payload["primes"]}
    with Heartbeat() as heartbeat:
        for p in primes:
            if p in completed:
                heartbeat.set(f"p={p}: retained from compatible resume file", emit=True)
                continue
            result = compute_prime(
                p,
                args.h_translate_counts,
                args.j_translates,
                args.max_chunk_cells,
                heartbeat,
            )
            payload["primes"].append(result)
            payload["primes"].sort(key=lambda item: item["p"])
            payload["verdicts"] = compute_verdicts(payload)
            payload["generated_at_utc"] = datetime.now(timezone.utc).isoformat()
            write_json(args.results, payload)
            write_report(args.report, render_report(payload))
            heartbeat.set(f"p={p}: checkpoint written", emit=True)

    payload["primes"] = [item for item in payload["primes"] if item["p"] in primes]
    payload["primes"].sort(key=lambda item: item["p"])
    payload["verdicts"] = compute_verdicts(payload)
    payload["generated_at_utc"] = datetime.now(timezone.utc).isoformat()
    write_json(args.results, payload)
    write_report(args.report, render_report(payload))
    print(
        f"[wallprobe] complete: {args.results} and {args.report} "
        f"({sum(len(item['boxes']) for item in payload['primes'])} boxes)",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
