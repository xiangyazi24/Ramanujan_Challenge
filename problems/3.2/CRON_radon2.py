#!/usr/bin/env python3
"""Q6556 Sections 7--8: refined cyclic Radon-spectrum experiment.

The script implements four tables:

* the Apéry Casoratian table on the full cyclic (r,h)-plane;
* the central-binomial negative control from Q6556 (7.2)--(7.4);
* a reflected random projective table with the same diagonal/mirror zeros;
* that same random table on the nonwrapping triangular mask.

For the Apéry and random tables, ``circ`` removes the union of the diagonal,
the forced reflection line s=-1-r, and the restart row/column r=-1 or s=-1.
Q6556 did not give a numerical formula for any additional constant/main-term
subtraction.  The literal primary zero-push therefore uses exactly that cell
removal and no extra center.  An ``active_cells(h)/p`` model-centered value is
reported separately as a sensitivity check.  The binomial control uses the
uncentered collision count exactly as Section 7 specifies.

Q6556's diagnostics are

    F(t,xi) = sum_{r,h} e_p(t D(r,r+h) + xi h),
    B_p     = p^-4 sum_xi |sum_t F^circ(t,xi)|^2,
    B_p_max = p^-3/2 max_xi |sum_t F^circ(t,xi)|.

With C^circ(h) the centered zero count, orthogonality and Parseval give the
computationally exact and more stable formulas

    B_p     = (1/p) sum_h |C^circ(h)|^2,
    B_p_max = max_xi |FFT(C^circ)(xi)| / sqrt(p).

NumPy's forward FFT has the negative sign.  Reported Apéry/random coordinates
are mapped back to Q6556's positive-sign (t,xi) convention.
"""

from __future__ import annotations

import argparse
import gc
import json
import math
import time
from pathlib import Path
from typing import Callable, Iterable

import numpy as np


DEFAULT_PRIMES = (101, 211, 401, 809, 1601, 3001)
DEFAULT_THRESHOLDS = (4.0, 5.0, 6.0, 8.0, 10.0)
DEFAULT_TOP_K = 12
DEFAULT_SEED = 20260801

# Q6556 asks for a fixed p-independent list xi=a*t+b but does not prescribe
# the coefficients.  Lock a modest list before looking at the output.
FIXED_LINE_SLOPES = tuple(range(-3, 4))
FIXED_LINE_INTERCEPTS = tuple(range(-2, 3))


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2
    return True


def parse_int_list(value: str) -> tuple[int, ...]:
    out = tuple(int(x.strip()) for x in value.split(",") if x.strip())
    if not out:
        raise argparse.ArgumentTypeError("expected a nonempty comma-separated list")
    if any(not is_prime(p) for p in out):
        bad = [p for p in out if not is_prime(p)]
        raise argparse.ArgumentTypeError(f"all entries must be prime; bad entries: {bad}")
    return out


def parse_float_list(value: str) -> tuple[float, ...]:
    out = tuple(float(x.strip()) for x in value.split(",") if x.strip())
    if not out:
        raise argparse.ArgumentTypeError("expected a nonempty comma-separated list")
    return out


def apery_pair(p: int) -> tuple[np.ndarray, np.ndarray]:
    """Return the two Apéry recurrence solutions on indices 0,...,p-1."""
    b = np.zeros(p, dtype=np.int64)
    c = np.zeros(p, dtype=np.int64)
    b[0], b[1] = 1, 5 % p
    c[0], c[1] = 0, 1
    for n in range(1, p - 1):
        pn = (34 * n**3 + 51 * n * n + 27 * n + 5) % p
        n3 = n**3 % p
        inv = pow((n + 1) ** 3 % p, p - 2, p)
        b[n + 1] = ((pn * int(b[n]) - n3 * int(b[n - 1])) % p) * inv % p
        c[n + 1] = ((pn * int(c[n]) - n3 * int(c[n - 1])) % p) * inv % p
    return b, c


def central_binom_table(p: int) -> np.ndarray:
    """a(r)=binom(2r,r) mod p, computed by Q6556's first-order recurrence."""
    a = np.zeros(p, dtype=np.int64)
    a[0] = 1
    for r in range(p - 1):
        a[r + 1] = (
            int(a[r]) * 2 * (2 * r + 1) * pow(r + 1, -1, p)
        ) % p
    return a


def reflected_random_projective_pair(
    p: int, seed: int
) -> tuple[np.ndarray, np.ndarray]:
    """Random P^1 table with pi(p-1-r)=pi(r), in canonical scalar lifts."""
    rng = np.random.default_rng(seed + 1_000_003 * p)
    values = np.full(p, -1, dtype=np.int64)
    for r in range(p):
        if values[r] >= 0:
            continue
        s = p - 1 - r
        value = int(rng.integers(0, p + 1))  # p denotes infinity
        values[r] = value
        values[s] = value
    infinity = values == p
    b = np.where(infinity, 1, values).astype(np.int64)
    c = np.where(infinity, 0, 1).astype(np.int64)
    return b, c


def full_domain(r: np.ndarray, s: np.ndarray, p: int) -> np.ndarray:
    del s, p
    return np.ones(r.shape, dtype=bool)


def nonwrapping_triangle(r: np.ndarray, s: np.ndarray, p: int) -> np.ndarray:
    """The ordered mask 0 <= r < s <= p-2 from Q6556 Section 1.1."""
    return (r < s) & (s <= p - 2)


def known_mode_parts(
    r: np.ndarray, s: np.ndarray, p: int, domain: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Disjoint diagonal, reflection, restart pieces and their union."""
    diagonal = domain & (s == r)
    reflection = domain & ~diagonal & (s == (p - 1 - r))
    restart = domain & ~(diagonal | reflection) & ((r == p - 1) | (s == p - 1))
    known = diagonal | reflection | restart
    return diagonal, reflection, restart, known


def build_histogram(
    p: int,
    b: np.ndarray,
    c: np.ndarray,
    domain_fn: Callable[[np.ndarray, np.ndarray, int], np.ndarray],
    remove_known: bool = True,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, int]]:
    """Build hist[h,d]=#{r active: D(r,r+h)=d}, with explicit mode removal."""
    r = np.arange(p, dtype=np.int64)
    hist = np.zeros((p, p), dtype=np.float64)
    zero_counts = np.zeros(p, dtype=np.int64)
    active_counts = np.zeros(p, dtype=np.int64)
    mode_counts = {"diagonal": 0, "reflection": 0, "restart": 0}

    for h in range(p):
        s = (r + h) % p
        domain = domain_fn(r, s, p)
        diagonal, reflection, restart, known = known_mode_parts(r, s, p, domain)
        mode_counts["diagonal"] += int(diagonal.sum())
        mode_counts["reflection"] += int(reflection.sum())
        mode_counts["restart"] += int(restart.sum())
        active = domain & ~known if remove_known else domain
        dvals = (b * c[s] - b[s] * c) % p
        row = np.bincount(dvals[active], minlength=p).astype(np.float64)
        hist[h, :] = row
        zero_counts[h] = int(row[0])
        active_counts[h] = int(active.sum())

    mode_counts["total"] = sum(mode_counts.values())
    mode_counts["domain"] = int(active_counts.sum()) + mode_counts["total"]
    mode_counts["active"] = int(active_counts.sum())
    return hist, zero_counts, active_counts, mode_counts


def transform_histogram(hist: np.ndarray) -> np.ndarray:
    """Return F[xi,t] using NumPy's negative-sign FFT convention."""
    a = np.fft.fft(hist, axis=1)
    f = np.fft.fft(a, axis=0)
    del a
    return f


def q_coordinate(index: int, p: int, sign: int) -> int:
    return (sign * index) % p


def line_labels(t: int, xi: int, p: int) -> list[str]:
    labels: list[str] = []
    if t == 0:
        labels.append("t=0")
    if xi == 0:
        labels.append("xi=0")
    for a in FIXED_LINE_SLOPES:
        for b in FIXED_LINE_INTERCEPTS:
            if xi == (a * t + b) % p:
                labels.append(f"xi={a}*t{b:+d}")
    return labels


def distribution_stats(values: np.ndarray, thresholds: Iterable[float]) -> dict:
    finite = values.ravel()
    n = int(finite.size)
    exceed = {}
    for threshold in thresholds:
        count = int(np.count_nonzero(finite > threshold))
        exceed[f"{threshold:g}"] = {"count": count, "fraction": count / n}
    return {
        "points": n,
        "max": float(finite.max()),
        "median": float(np.median(finite)),
        "q99": float(np.quantile(finite, 0.99)),
        "mean": float(finite.mean()),
        "thresholds": exceed,
    }


def exceptional_locus_stats(
    values: np.ndarray,
    p: int,
    t_sign: int,
    xi_sign: int,
    top_k: int,
) -> dict:
    """Analyze values[xi_index,t_index-1] against a fixed affine line list."""
    if values.shape != (p, p - 1):
        raise ValueError(f"unexpected spectrum shape {values.shape}")

    flat = values.ravel()
    k = min(top_k, flat.size)
    chosen = np.argpartition(flat, flat.size - k)[-k:]
    chosen = chosen[np.argsort(flat[chosen])[::-1]]
    top = []
    top_hits = 0
    for flat_index in chosen:
        xi_index, t_offset = divmod(int(flat_index), p - 1)
        t = q_coordinate(t_offset + 1, p, t_sign)
        xi = q_coordinate(xi_index, p, xi_sign)
        labels = line_labels(t, xi, p)
        if labels:
            top_hits += 1
        top.append(
            {
                "value": float(flat[flat_index]),
                "t": t,
                "xi": xi,
                "lines": labels,
            }
        )

    t_indices = np.arange(1, p, dtype=np.int64)
    t_coords = (t_sign * t_indices) % p
    candidate = np.zeros(values.shape, dtype=bool)
    for xi_index in range(p):
        xi = q_coordinate(xi_index, p, xi_sign)
        row = xi == 0
        if row:
            candidate[xi_index, :] = True
        for a in FIXED_LINE_SLOPES:
            for b in FIXED_LINE_INTERCEPTS:
                candidate[xi_index, :] |= xi == ((a * t_coords + b) % p)

    candidate_max = float(values[candidate].max())
    off_max = 0.0
    for xi_index in range(p):
        keep = ~candidate[xi_index]
        if np.any(keep):
            off_max = max(off_max, float(values[xi_index, keep].max()))

    return {
        "fixed_family": {
            "axes": ["t=0", "xi=0"],
            "slopes": list(FIXED_LINE_SLOPES),
            "intercepts": list(FIXED_LINE_INTERCEPTS),
            "formula": "xi = a*t + b (mod p)",
        },
        "candidate_points": int(candidate.sum()),
        "candidate_max": candidate_max,
        "off_candidate_max": off_max,
        "top_k": top,
        "top_k_on_fixed_family": top_hits,
    }


def spectrum_stats(
    values: np.ndarray,
    p: int,
    thresholds: Iterable[float],
    top_k: int,
    t_sign: int,
    xi_sign: int,
) -> dict:
    """Statistics for values=|F^circ|/p, shaped [xi,t!=0]."""
    result = {
        "t_nonzero_all_xi": distribution_stats(values, thresholds),
        "t_nonzero_xi_nonzero": distribution_stats(values[1:, :], thresholds),
        "xi_zero_max": float(values[0, :].max()),
        "exceptional_locus": exceptional_locus_stats(
            values, p, t_sign=t_sign, xi_sign=xi_sign, top_k=top_k
        ),
    }
    return result


def zero_push_stats(
    zero_counts: np.ndarray,
    active_counts: np.ndarray,
    p: int,
    center_denominator: int,
    primary_center: str = "none",
) -> dict:
    """Compute Q6556 B_p and B_p^max from the exact zero-count transform."""
    model_main = active_counts.astype(np.float64) / center_denominator
    uncentered = zero_counts.astype(np.float64)
    model_centered = uncentered - model_main
    empirical = uncentered - uncentered.mean()

    if primary_center == "none":
        primary = uncentered
        primary_description = "none beyond the explicit geometric cell-mode removal"
    elif primary_center == "model":
        primary = model_centered
        primary_description = f"active_cells(h)/{center_denominator}"
    else:
        raise ValueError(f"unknown primary center {primary_center!r}")

    primary_hat = np.fft.fft(primary)

    b_value = float(np.mean(primary * primary))
    b_max_values = np.abs(primary_hat) / math.sqrt(p)
    b_max_index = int(np.argmax(b_max_values))
    b_max = float(b_max_values[b_max_index])

    # Literal Q6556 formula: T=p*FFT(C^circ), then p^-4 sum |T|^2.
    t_push = p * primary_hat
    b_from_t_push = float(np.sum(np.abs(t_push) ** 2) / p**4)

    uncentered_hat = np.fft.fft(uncentered)
    model_centered_hat = np.fft.fft(model_centered)
    empirical_hat = np.fft.fft(empirical)

    top_indices = np.argsort(b_max_values)[-8:][::-1]
    top_frequencies = [
        {
            "xi": int((-int(index)) % p),
            "normalized_abs": float(b_max_values[index]),
        }
        for index in top_indices
    ]

    return {
        "primary_center": primary_description,
        "B_p": b_value,
        "B_p_from_literal_t_push": b_from_t_push,
        "parseval_abs_error": abs(b_value - b_from_t_push),
        "B_p_max": b_max,
        "B_p_max_xi": int((-b_max_index) % p),
        "top_frequencies": top_frequencies,
        "mean_zero_count": float(uncentered.mean()),
        "mean_model_main": float(model_main.mean()),
        "max_zero_count": int(zero_counts.max()),
        "uncentered_B": float(np.mean(uncentered * uncentered)),
        "uncentered_B_max": float(np.abs(uncentered_hat).max() / math.sqrt(p)),
        "uncentered_B_max_off_xi_zero": float(
            np.abs(uncentered_hat[1:]).max() / math.sqrt(p)
        ),
        "model_center": f"active_cells(h)/{center_denominator}",
        "model_centered_B": float(np.mean(model_centered * model_centered)),
        "model_centered_B_max": float(
            np.abs(model_centered_hat).max() / math.sqrt(p)
        ),
        "model_centered_B_max_off_xi_zero": float(
            np.abs(model_centered_hat[1:]).max() / math.sqrt(p)
        ),
        "empirical_mean_centered_B": float(np.mean(empirical * empirical)),
        "empirical_mean_centered_B_max": float(
            np.abs(empirical_hat).max() / math.sqrt(p)
        ),
    }


def analyze_histogram_model(
    name: str,
    p: int,
    b: np.ndarray,
    c: np.ndarray,
    domain_fn: Callable[[np.ndarray, np.ndarray, int], np.ndarray],
    thresholds: Iterable[float],
    top_k: int,
    center_denominator: int,
) -> dict:
    start = time.time()
    hist, zero_counts, active_counts, mode_counts = build_histogram(
        p, b, c, domain_fn=domain_fn, remove_known=True
    )
    f_grid = transform_histogram(hist)
    del hist
    normalized = np.abs(f_grid[:, 1:]) / p
    del f_grid
    raw = spectrum_stats(
        normalized,
        p,
        thresholds=thresholds,
        top_k=top_k,
        t_sign=-1,
        xi_sign=-1,
    )
    del normalized
    zero_push = zero_push_stats(
        zero_counts,
        active_counts,
        p,
        center_denominator=center_denominator,
        primary_center="none",
    )
    result = {
        "model": name,
        "p": p,
        "mode_removal": mode_counts,
        "raw_transform": raw,
        "zero_push": zero_push,
        "elapsed_seconds": time.time() - start,
    }
    gc.collect()
    return result


def analyze_apery(
    p: int, thresholds: Iterable[float], top_k: int
) -> dict:
    b, c = apery_pair(p)
    result = analyze_histogram_model(
        "apery_cyclic",
        p,
        b,
        c,
        domain_fn=full_domain,
        thresholds=thresholds,
        top_k=top_k,
        center_denominator=p,
    )
    result["table_checks"] = {
        "b_endpoint": int(b[-1]),
        "c_endpoint": int(c[-1]),
        "b_reflection": bool(np.array_equal(b, b[::-1])),
        "c_reflection": bool(np.array_equal(c, c[::-1])),
    }
    return result


def binomial_zero_counts(a: np.ndarray, p: int) -> tuple[np.ndarray, np.ndarray]:
    r = np.arange(p)
    counts = np.zeros(p, dtype=np.int64)
    active = np.zeros(p, dtype=np.int64)
    for h in range(p):
        s = (r + h) % p
        equal = a == a[s]
        # Q6556 Section 7 uses the raw, uncentered collision function, including
        # h=0.  Do not import the Apéry ``circ`` subtraction into this control.
        counts[h] = int(equal.sum())
        active[h] = p
    return counts, active


def analyze_binomial(
    p: int, thresholds: Iterable[float], top_k: int
) -> dict:
    start = time.time()
    a = central_binom_table(p)
    # values[xi,t-1] = F_bin(t,xi)/p.  This is the raw factorized spectrum
    # requested in Q6556 (7.3), so the diagonal is not subtracted here.
    values = np.empty((p, p - 1), dtype=np.float64)
    phase_scale = 2j * np.pi / p
    for t in range(1, p):
        v = np.exp(phase_scale * t * a)
        s_hat = np.fft.fft(v)
        values[:, t - 1] = np.abs(s_hat) ** 2 / p

    raw = spectrum_stats(
        values,
        p,
        thresholds=thresholds,
        top_k=top_k,
        t_sign=1,
        xi_sign=1,
    )
    xi_one_max = float(values[1, :].max())
    xi_one_t = int(np.argmax(values[1, :])) + 1

    signed_low_frequencies = list(range(0, 9)) + list(range(-8, 0))
    ridge = {
        str(xi): float(values[xi % p, :].max()) for xi in signed_low_frequencies
    }

    zero_counts, active_counts = binomial_zero_counts(a, p)
    zero_push = zero_push_stats(
        zero_counts,
        active_counts,
        p,
        center_denominator=p,
        primary_center="none",
    )
    zero_plateau = np.flatnonzero(a == 0)
    result = {
        "model": "central_binomial",
        "p": p,
        "zero_plateau_start": int(zero_plateau[0]),
        "zero_plateau_length": int(zero_plateau.size),
        "raw_transform": raw,
        "xi_one_control": {"max_F_over_p": xi_one_max, "t": xi_one_t},
        "low_frequency_ridge": ridge,
        "zero_push": zero_push,
        "elapsed_seconds": time.time() - start,
    }
    del values
    gc.collect()
    return result


def analyze_random_controls(
    p: int, thresholds: Iterable[float], top_k: int, seed: int
) -> tuple[dict, dict]:
    b, c = reflected_random_projective_pair(p, seed)
    cyclic = analyze_histogram_model(
        "random_reflected_cyclic",
        p,
        b,
        c,
        domain_fn=full_domain,
        thresholds=thresholds,
        top_k=top_k,
        center_denominator=p + 1,
    )
    cyclic["seed"] = seed
    masked = analyze_histogram_model(
        "random_reflected_nonwrapping_mask",
        p,
        b,
        c,
        domain_fn=nonwrapping_triangle,
        thresholds=thresholds,
        top_k=top_k,
        center_denominator=p + 1,
    )
    masked["seed"] = seed
    return cyclic, masked


def direct_transform(hist: np.ndarray, p: int) -> np.ndarray:
    """Slow negative-sign transform used only by the small-prime validator."""
    out = np.zeros((p, p), dtype=np.complex128)
    zeta = np.exp(-2j * np.pi / p)
    for xi in range(p):
        for t in range(p):
            total = 0j
            for h in range(p):
                for d in range(p):
                    total += hist[h, d] * zeta ** ((t * d + xi * h) % p)
            out[xi, t] = total
    return out


def validate_small_primes(verbose: bool = True) -> dict:
    """Exhaustive correctness gates before the large sweep."""
    checks = 0
    max_fft_error = 0.0
    max_bin_factor_error = 0.0
    max_push_error = 0.0

    for p in (5, 7, 11):
        a = central_binom_table(p)
        expected = np.array([math.comb(2 * r, r) % p for r in range(p)])
        assert np.array_equal(a, expected)
        checks += 1
        expected_zero = np.arange((p + 1) // 2, p)
        assert np.array_equal(np.flatnonzero(a == 0), expected_zero)
        checks += 1

        # Factorization (7.4) versus an independent double sum.
        zeta = np.exp(2j * np.pi / p)
        for t in range(p):
            for xi in range(p):
                factored = abs(
                    sum(zeta ** ((t * int(a[r]) - xi * r) % p) for r in range(p))
                ) ** 2
                doubled = sum(
                    zeta
                    ** (
                        t * (int(a[r]) - int(a[s])) + xi * ((s - r) % p)
                    )
                    for r in range(p)
                    for s in range(p)
                )
                error = abs(factored - doubled)
                max_bin_factor_error = max(max_bin_factor_error, error)
                assert error < 1e-8 * p * p
        checks += 1

        b, c = apery_pair(p)
        assert np.array_equal(b, b[::-1])
        assert np.array_equal(c, c[::-1])
        assert (int(b[-1]), int(c[-1])) == (1, 0)
        checks += 1

        hist, zero_counts, active_counts, modes = build_histogram(
            p, b, c, full_domain, remove_known=True
        )
        assert np.array_equal(hist.sum(axis=1).astype(np.int64), active_counts)
        assert modes["active"] + modes["total"] == p * p
        checks += 1
        fast = transform_histogram(hist)
        slow = direct_transform(hist, p)
        error = float(np.max(np.abs(fast - slow)))
        max_fft_error = max(max_fft_error, error)
        assert error < 1e-8 * p * p
        checks += 1

        # Exact orthogonality: sum_t F(xi,t)=p*FFT_h(C)(xi).
        pushed = fast.sum(axis=1)
        expected_push = p * np.fft.fft(zero_counts.astype(np.float64))
        error = float(np.max(np.abs(pushed - expected_push)))
        max_push_error = max(max_push_error, error)
        assert error < 1e-8 * p * p
        checks += 1

        zp = zero_push_stats(zero_counts, active_counts, p, p)
        assert zp["parseval_abs_error"] < 1e-10
        checks += 1

    result = {
        "status": "PASS",
        "primes": [5, 7, 11],
        "checks": checks,
        "max_fft_error": max_fft_error,
        "max_binomial_factorization_error": max_bin_factor_error,
        "max_zero_push_identity_error": max_push_error,
    }
    if verbose:
        print(
            "VALIDATION PASS "
            f"checks={checks} max_fft_err={max_fft_error:.3e} "
            f"max_bin_factor_err={max_bin_factor_error:.3e} "
            f"max_push_err={max_push_error:.3e}",
            flush=True,
        )
    return result


def loglog_slope(points: list[tuple[int, float]]) -> dict:
    usable = [(p, y) for p, y in points if y > 0]
    if len(usable) < 2:
        return {"slope": None, "intercept": None, "r_squared": None}
    x = np.log([p for p, _ in usable])
    y = np.log([value for _, value in usable])
    slope, intercept = np.polyfit(x, y, 1)
    fitted = intercept + slope * x
    ss_res = float(np.sum((y - fitted) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot else 1.0
    return {"slope": float(slope), "intercept": float(intercept), "r_squared": r2}


def summarize_scaling(rows: list[dict]) -> dict:
    apery_a = []
    apery_a_off_lines = []
    apery_b = []
    bin_xi1 = []
    bin_b = []
    for row in rows:
        p = row["p"]
        apery = row["apery"]
        binomial = row["binomial"]
        apery_a.append(
            (p, apery["raw_transform"]["t_nonzero_xi_nonzero"]["max"])
        )
        apery_a_off_lines.append(
            (p, apery["raw_transform"]["exceptional_locus"]["off_candidate_max"])
        )
        apery_b.append((p, apery["zero_push"]["B_p"]))
        bin_xi1.append((p, binomial["xi_one_control"]["max_F_over_p"]))
        bin_b.append((p, binomial["zero_push"]["B_p"]))
    return {
        "apery_raw_off_axis": loglog_slope(apery_a),
        "apery_raw_off_fixed_lines": loglog_slope(apery_a_off_lines),
        "apery_zero_push_B": loglog_slope(apery_b),
        "binomial_xi_one": loglog_slope(bin_xi1),
        "binomial_zero_push_B": loglog_slope(bin_b),
    }


def print_model_summary(model: dict) -> None:
    raw = model["raw_transform"]
    all_stats = raw["t_nonzero_all_xi"]
    off_stats = raw["t_nonzero_xi_nonzero"]
    exc = raw["exceptional_locus"]
    zp = model["zero_push"]
    maximum = exc["top_k"][0]
    label_text = ",".join(maximum["lines"]) if maximum["lines"] else "off-fixed"
    print(
        f"  {model['model']}: max={all_stats['max']:.6g} "
        f"offaxis={off_stats['max']:.6g} median={all_stats['median']:.6g} "
        f"q99={all_stats['q99']:.6g} B={zp['B_p']:.6g} "
        f"Bmax={zp['B_p_max']:.6g} "
        f"argmax=({maximum['t']},{maximum['xi']})[{label_text}] "
        f"offlines={exc['off_candidate_max']:.6g} "
        f"elapsed={model['elapsed_seconds']:.2f}s",
        flush=True,
    )


def make_plot(results: dict, output: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rows = results["rows"]
    p = np.array([row["p"] for row in rows])
    apery_raw = np.array(
        [
            row["apery"]["raw_transform"]["t_nonzero_xi_nonzero"]["max"]
            for row in rows
        ]
    )
    apery_off_lines = np.array(
        [
            row["apery"]["raw_transform"]["exceptional_locus"][
                "off_candidate_max"
            ]
            for row in rows
        ]
    )
    bin_xi1 = np.array(
        [row["binomial"]["xi_one_control"]["max_F_over_p"] for row in rows]
    )
    apery_b = np.array([row["apery"]["zero_push"]["B_p"] for row in rows])
    bin_b = np.array([row["binomial"]["zero_push"]["B_p"] for row in rows])
    apery_bmax = np.array(
        [row["apery"]["zero_push"]["B_p_max"] for row in rows]
    )
    bin_bmax = np.array(
        [row["binomial"]["zero_push"]["B_p_max"] for row in rows]
    )

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    ax = axes[0]
    ax.loglog(p, apery_raw, "o-", label="Apéry max, t≠0, ξ≠0")
    ax.loglog(p, apery_off_lines, "s--", label="Apéry max off fixed lines")
    ax.loglog(p, bin_xi1, "^-", label="binomial max at ξ=1")
    ax.loglog(p, p / p[0] * bin_xi1[0], ":", color="0.4", label="slope 1")
    ax.set_xlabel("prime p")
    ax.set_ylabel("normalized raw spectrum |F|/p")
    ax.set_title("Raw Radon transform")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=8)

    ax = axes[1]
    ax.loglog(p, apery_b, "o-", label="Apéry B_p")
    ax.loglog(p, bin_b, "^-", label="binomial B_p")
    ax.loglog(p, apery_bmax, "s--", label="Apéry B_p^max")
    ax.loglog(p, bin_bmax, "v--", label="binomial B_p^max")
    ax.axhline(1.0, color="0.5", linestyle=":", linewidth=1)
    ax.set_xlabel("prime p")
    ax.set_ylabel("zero-push diagnostic")
    ax.set_title("GPRV zero-push (Q6556 (0.1))")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=8)

    fig.suptitle("Apéry versus central-binomial negative control")
    fig.tight_layout()
    fig.savefig(output, dpi=180)
    plt.close(fig)


def run(args: argparse.Namespace) -> dict:
    validation = validate_small_primes(verbose=True)
    if args.validate_only:
        return {"validation": validation}

    thresholds = tuple(args.thresholds)
    rows = []
    for p in args.primes:
        print(f"p={p}", flush=True)
        apery = analyze_apery(p, thresholds=thresholds, top_k=args.top_k)
        print_model_summary(apery)
        binomial = analyze_binomial(p, thresholds=thresholds, top_k=args.top_k)
        print_model_summary(binomial)
        print(
            f"    binomial xi=1 max/p={binomial['xi_one_control']['max_F_over_p']:.6g} "
            f"at t={binomial['xi_one_control']['t']}",
            flush=True,
        )
        row = {"p": p, "apery": apery, "binomial": binomial}
        if not args.skip_random_controls:
            random_cyclic, random_masked = analyze_random_controls(
                p,
                thresholds=thresholds,
                top_k=args.top_k,
                seed=args.seed,
            )
            print_model_summary(random_cyclic)
            print_model_summary(random_masked)
            row["random_cyclic"] = random_cyclic
            row["random_masked"] = random_masked
        rows.append(row)

    results = {
        "metadata": {
            "protocol": "Q6556 Sections 7-8",
            "primes": list(args.primes),
            "thresholds": list(thresholds),
            "threshold_note": "Q6556 did not prescribe values; fixed before sweep",
            "top_k": args.top_k,
            "seed": args.seed,
            "random_controls": not args.skip_random_controls,
            "fixed_line_slopes": list(FIXED_LINE_SLOPES),
            "fixed_line_intercepts": list(FIXED_LINE_INTERCEPTS),
            "mode_removal": (
                "diagonal union reflection s=-1-r union restart r=-1 or s=-1; "
                "primary zero-push has no further center; active/p (random P1: /p+1) "
                "is reported only as a sensitivity statistic"
            ),
            "coordinate_convention": "reported coordinates use positive e_p(tD+xi*h)",
        },
        "validation": validation,
        "rows": rows,
    }
    results["scaling_fits"] = summarize_scaling(rows)
    return results


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--primes",
        type=parse_int_list,
        default=DEFAULT_PRIMES,
        help="comma-separated prime list (default: 101,211,401,809,1601,3001)",
    )
    parser.add_argument(
        "--thresholds",
        type=parse_float_list,
        default=DEFAULT_THRESHOLDS,
        help="fixed |F|/p thresholds (default: 4,5,6,8,10)",
    )
    parser.add_argument("--top-k", type=int, default=DEFAULT_TOP_K)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument(
        "--skip-random-controls",
        action="store_true",
        help="skip Q6556 Section 8.4 random cyclic/masked controls",
    )
    parser.add_argument("--validate-only", action="store_true")
    parser.add_argument(
        "--output-json", type=Path, default=Path("CRON_RADON2_RESULTS.json")
    )
    parser.add_argument(
        "--output-plot", type=Path, default=Path("CRON_RADON2_COMPARISON.png")
    )
    parser.add_argument("--no-plot", action="store_true")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    results = run(args)
    if args.validate_only:
        return
    args.output_json.write_text(
        json.dumps(results, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"RESULT_JSON {args.output_json}", flush=True)
    if not args.no_plot:
        make_plot(results, args.output_plot)
        print(f"PLOT {args.output_plot}", flush=True)
    print("SCALING_FITS " + json.dumps(results["scaling_fits"], sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
