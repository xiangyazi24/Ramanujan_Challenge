#!/usr/bin/env python3
"""Q7807: exact multiplicative and E2-weighted ratio energy for Apéry zeros.

This script is standalone (Python 3 + NumPy).  It computes:

1. The intrinsic multiplicative energy of Z_p for all primes p <= 1000.
2. The Q7798 Section 2 E2-weighted Gram matrix W_p for every p <= 500,
   placing p in its unique power-of-two dyadic block X < p <= 2X and using
   L=X^2 exactly as in Q7798.
3. The centered trivial-character off-diagonal contraction and the exact
   diagonal scale from (ARD_E2).

Outputs:
  q7807_energy_p1000.csv
  q7807_e2_p500.csv
  q7807_summary.md
"""

from __future__ import annotations

import csv
import math
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np


OUT = Path("q7807_results")
OUT.mkdir(exist_ok=True)


# ---------------------------------------------------------------------------
# Prime generation and Apéry zeros
# ---------------------------------------------------------------------------


def primes_upto(n: int) -> List[int]:
    sieve = bytearray(b"\x01") * (n + 1)
    if n >= 0:
        sieve[0] = 0
    if n >= 1:
        sieve[1] = 0
    for p in range(2, math.isqrt(n) + 1):
        if sieve[p]:
            sieve[p * p : n + 1 : p] = b"\x00" * (((n - p * p) // p) + 1)
    return [i for i in range(2, n + 1) if sieve[i]]


def apery_values_mod_p(p: int) -> Tuple[int, ...]:
    """Return (b_0,...,b_{p-1}) modulo p in O(p) operations."""
    if p < 3:
        raise ValueError("odd prime expected")
    b = [0] * p
    b[0] = 1
    if p > 1:
        b[1] = 5 % p

    inv = [0] * p
    inv[1] = 1
    for a in range(2, p):
        inv[a] = (-(p // a) * inv[p % a]) % p

    for n in range(1, p - 1):
        middle = ((2 * n + 1) * (17 * n * n + 17 * n + 5)) % p
        rhs = (middle * b[n] - (n * n * n % p) * b[n - 1]) % p
        u = inv[n + 1]
        b[n + 1] = rhs * (u * u % p) * u % p

    return tuple(b)


def apery_zero_set(p: int) -> Tuple[int, ...]:
    b = apery_values_mod_p(p)
    z = tuple(i for i, x in enumerate(b) if x == 0)

    # Banked integrity checks; failures indicate a normalization/code error.
    if p >= 5:
        if b[p - 1] != 1:
            raise AssertionError(("endpoint", p, b[p - 1]))
        if any(b[r] != b[p - 1 - r] for r in range(p)):
            raise AssertionError(("reflection", p))
        if any(b[r] == 0 and b[r + 1] == 0 for r in range(p - 1)):
            raise AssertionError(("consecutive zeros", p))
    if 0 in z:
        raise AssertionError(("b_0 should be a unit", p))
    return z


# ---------------------------------------------------------------------------
# Intrinsic multiplicative energy
# ---------------------------------------------------------------------------


def multiplicative_ratio_counts(p: int, zeros: Sequence[int]) -> np.ndarray:
    """M[h] = #{(x,y) in Z^2: x = h y mod p}, h in F_p."""
    counts = np.zeros(p, dtype=np.int64)
    for y0 in zeros:
        y = int(y0)
        yi = pow(y, -1, p)
        for x0 in zeros:
            counts[(int(x0) * yi) % p] += 1
    if counts[0] != 0:
        raise AssertionError(("zero ratio from nonzero zero-set", p))
    if int(np.sum(counts)) != len(zeros) ** 2:
        raise AssertionError(("ratio mass", p))
    return counts


def intrinsic_energy_row(p: int, zeros: Sequence[int]) -> dict:
    z = len(zeros)
    M = multiplicative_ratio_counts(p, zeros)
    energy = int(np.dot(M, M))
    sidon_floor = 2 * z * z - z
    off = M.copy()
    off[1] = 0
    off_ratio_size = int(np.count_nonzero(off))
    max_off_mult = int(off.max(initial=0))
    collision_excess = energy - sidon_floor
    if collision_excess < 0:
        raise AssertionError(("below Sidon floor", p, energy, sidon_floor))
    return {
        "p": p,
        "z": z,
        "zeros": " ".join(str(x) for x in zeros),
        "energy": energy,
        "sidon_floor_2z2_minus_z": sidon_floor,
        "energy_minus_floor": collision_excess,
        "sidon": int(energy == sidon_floor),
        "z_cubed": z ** 3,
        "energy_over_z_cubed": energy / (z ** 3),
        "offdiag_ordered_pairs": z * (z - 1),
        "offdiag_ratio_set_size": off_ratio_size,
        "max_offdiag_ratio_multiplicity": max_off_mult,
    }


# ---------------------------------------------------------------------------
# Q7798 E2 completion
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class BlockData:
    X: int
    L: int
    primes: Tuple[int, ...]
    zeros: Dict[int, Tuple[int, ...]]
    masks: Tuple[np.ndarray, ...]
    pair_hits: Tuple[Tuple[int, int, int, int, np.ndarray], ...]


def unique_power_two_block(p: int) -> int:
    """Unique power of two X with X < p <= 2X (p odd)."""
    return 1 << (p.bit_length() - 1)


def build_block(X: int, zero_cache: Dict[int, Tuple[int, ...]]) -> BlockData:
    L = X * X
    primes = tuple(p for p in primes_upto(2 * X) if X < p <= 2 * X)
    masks: List[np.ndarray] = []
    for p in primes:
        mask = np.zeros(L, dtype=np.bool_)
        for r in zero_cache[p]:
            mask[r:L:p] = True
        masks.append(mask)

    pair_hits: List[Tuple[int, int, int, int, np.ndarray]] = []
    for i in range(len(primes)):
        for j in range(i + 1, len(primes)):
            pos = np.flatnonzero(masks[i] & masks[j]).astype(np.int32)
            if pos.size:
                pair_hits.append((i, j, primes[i], primes[j], pos))

    return BlockData(
        X=X,
        L=L,
        primes=primes,
        zeros={p: zero_cache[p] for p in primes},
        masks=tuple(masks),
        pair_hits=tuple(pair_hits),
    )


def beta_matrix_physical(block: BlockData, p: int) -> Tuple[np.ndarray, np.ndarray]:
    """Compute beta_{p,a}(u) from the exact physical-space identity.

    Q7798 (1.3), after CRT Fourier inversion in q and ell, becomes

      beta_{p,a}(u) = (1/p) sum_{0<=m<L} H_{p,u}(m) e_p(-a u m),

    where H_{p,u}(m) counts ordered q != ell, both distinct from p,
    hitting m and satisfying (q ell)^(-1) = u mod p.

    Returns:
      B[u-1,a-1] = beta_{p,a}(u), 1<=u,a<p,
      C[u,r]      = sum_{m == r mod p} H_{p,u}(m).
    """
    if p not in block.primes:
        raise ValueError((p, block.X))

    p_index = block.primes.index(p)
    C = np.zeros((p, p), dtype=np.float64)

    # pair_hits stores unordered q<ell. Q7798 uses ordered q,ell, so each
    # physical pair contributes twice. The semiprime residue q*ell is symmetric.
    for i, j, q, ell, pos in block.pair_hits:
        if i == p_index or j == p_index:
            continue
        u = pow((q * ell) % p, -1, p)
        residues = pos % p
        np.add.at(C[u], residues, 2.0)

    # NumPy FFT uses exp(-2*pi*i*k*r/p), matching e_p(-a*u*r).
    F = np.fft.fft(C, axis=1) / p
    uvals = np.arange(1, p, dtype=np.int64)
    avals = np.arange(1, p, dtype=np.int64)
    freq = (uvals[:, None] * avals[None, :]) % p
    B = F[uvals[:, None], freq]
    return B, C


def e2_weighted_row(block: BlockData, p: int) -> dict:
    zeros = block.zeros[p]
    z = len(zeros)
    M = multiplicative_ratio_counts(p, zeros)
    T = p * M - z * z
    T1 = int(T[1])
    if T1 != p * z - z * z:
        raise AssertionError(("T1", p))

    B, C = beta_matrix_physical(block, p)
    ordered_e2_incidence_mass = int(round(float(C.sum())))

    if not np.any(B):
        return {
            "X": block.X,
            "L": block.L,
            "p": p,
            "z": z,
            "ordered_E2_pair_incidence_mass": ordered_e2_incidence_mass,
            "trace_W": 0.0,
            "frobenius_W_sq": 0.0,
            "max_abs_offdiag_W": 0.0,
            "T1": T1,
            "tau_offdiag": p * z * (z - 1) / (p - 2) - z * z,
            "diagonal_scale": 0.0,
            "offdiag_raw": 0.0,
            "offdiag_centered": 0.0,
            "abs_centered_over_diagonal": float("nan"),
            "signed_centered_over_diagonal": float("nan"),
            "raw_offdiag_over_diagonal": float("nan"),
            "hermitian_error": 0.0,
            "imaginary_contraction_error": 0.0,
            "status": "zero E2 completion",
        }

    # Q7798 (2.5): W(a,a')=(p-1)^(-1) sum_u beta_a(u) conj(beta_a'(u)).
    W = (B.T @ B.conj()) / (p - 1)
    herm_err = float(np.max(np.abs(W - W.conj().T)))

    avals = np.arange(1, p, dtype=np.int64)
    inva = np.asarray([pow(int(a), -1, p) for a in avals], dtype=np.int64)
    hmat = (inva[:, None] * avals[None, :]) % p
    Tmat = T[hmat].astype(np.float64)
    off_mask = ~np.eye(p - 1, dtype=np.bool_)

    tau = p * z * (z - 1) / (p - 2) - z * z
    traceW = float(np.trace(W).real)
    diag_scale = float(T1 * traceW)
    off_raw_c = np.sum(W[off_mask] * Tmat[off_mask])
    off_ctr_c = np.sum(W[off_mask] * (Tmat[off_mask] - tau))
    imag_err = max(abs(float(off_raw_c.imag)), abs(float(off_ctr_c.imag)))
    off_raw = float(off_raw_c.real)
    off_ctr = float(off_ctr_c.real)

    offW = W.copy()
    np.fill_diagonal(offW, 0.0)

    return {
        "X": block.X,
        "L": block.L,
        "p": p,
        "z": z,
        "ordered_E2_pair_incidence_mass": ordered_e2_incidence_mass,
        "trace_W": traceW,
        "frobenius_W_sq": float(np.sum(np.abs(W) ** 2)),
        "max_abs_offdiag_W": float(np.max(np.abs(offW), initial=0.0)),
        "T1": T1,
        "tau_offdiag": tau,
        "diagonal_scale": diag_scale,
        "offdiag_raw": off_raw,
        "offdiag_centered": off_ctr,
        "abs_centered_over_diagonal": abs(off_ctr) / diag_scale if diag_scale else float("nan"),
        "signed_centered_over_diagonal": off_ctr / diag_scale if diag_scale else float("nan"),
        "raw_offdiag_over_diagonal": off_raw / diag_scale if diag_scale else float("nan"),
        "hermitian_error": herm_err,
        "imaginary_contraction_error": imag_err,
        "status": "ok" if diag_scale > 0 else "zero diagonal",
    }


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------


def write_csv(path: Path, rows: Sequence[dict]) -> None:
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def md_float(x: float, digits: int = 6) -> str:
    if not math.isfinite(x):
        return "NA"
    return f"{x:.{digits}g}"


def main() -> None:
    t0 = time.perf_counter()
    primes1000 = [p for p in primes_upto(1000) if p >= 7]
    zero_cache: Dict[int, Tuple[int, ...]] = {
        p: apery_zero_set(p) for p in primes_upto(1024) if p >= 5
    }

    energy_rows = [
        intrinsic_energy_row(p, zero_cache[p])
        for p in primes1000
        if len(zero_cache[p]) >= 2
    ]
    write_csv(OUT / "q7807_energy_p1000.csv", energy_rows)

    # Compute Q7798's W_p in each prime's unique power-of-two block.
    e2_rows: List[dict] = []
    block_summaries: List[dict] = []
    for X in (4, 8, 16, 32, 64, 128, 256):
        block = build_block(X, zero_cache)
        rows_this_block: List[dict] = []
        for p in block.primes:
            if p > 500 or p < 7 or len(zero_cache[p]) < 2:
                continue
            row = e2_weighted_row(block, p)
            rows_this_block.append(row)
            e2_rows.append(row)

        diag = sum(float(r["diagonal_scale"]) for r in rows_this_block)
        abs_ctr = sum(abs(float(r["offdiag_centered"])) for r in rows_this_block)
        signed_ctr = sum(float(r["offdiag_centered"]) for r in rows_this_block)
        block_summaries.append(
            {
                "X": X,
                "number_of_window_primes": len(block.primes),
                "number_of_output_primes_z_ge_2_p_le_500": len(rows_this_block),
                "sum_diagonal_scale": diag,
                "sum_abs_centered_offdiag": abs_ctr,
                "sum_signed_centered_offdiag": signed_ctr,
                "ARD_abs_ratio": abs_ctr / diag if diag else float("nan"),
                "ARD_signed_ratio": signed_ctr / diag if diag else float("nan"),
            }
        )

    write_csv(OUT / "q7807_e2_p500.csv", e2_rows)
    write_csv(OUT / "q7807_e2_block_summary.csv", block_summaries)

    sidon_rows = [r for r in energy_rows if r["sidon"] == 1]
    nonsidon_rows = [r for r in energy_rows if r["sidon"] == 0]

    # Aggregate weighted ratio exactly in the form of Q7798 (ARD_E2), blockwise.
    valid_e2 = [r for r in e2_rows if float(r["diagonal_scale"]) > 0]
    total_diag = sum(float(r["diagonal_scale"]) for r in valid_e2)
    total_abs_ctr = sum(abs(float(r["offdiag_centered"])) for r in valid_e2)
    total_signed_ctr = sum(float(r["offdiag_centered"]) for r in valid_e2)

    lines: List[str] = []
    lines.append("# Q7807 exact ratio-energy computation")
    lines.append("")
    lines.append(f"Runtime: {time.perf_counter() - t0:.3f} seconds (before Markdown write).")
    lines.append("")
    lines.append("## Intrinsic multiplicative energy, p <= 1000")
    lines.append("")
    lines.append(f"Primes with z_p >= 2: **{len(energy_rows)}**.")
    lines.append(f"Sidon-floor primes: **{len(sidon_rows)}**.")
    lines.append(f"Strictly above the Sidon floor: **{len(nonsidon_rows)}**.")
    lines.append("")
    lines.append("The exact floor used is `2 z^2 - z`, not `z^3`.")
    lines.append("")
    lines.append("| p | z | E_x(Z_p) | 2z^2-z | excess | |R_off| | max off mult | zeros |")
    lines.append("|---:|---:|---:|---:|---:|---:|---:|:---|")
    for r in energy_rows:
        lines.append(
            f"| {r['p']} | {r['z']} | {r['energy']} | "
            f"{r['sidon_floor_2z2_minus_z']} | {r['energy_minus_floor']} | "
            f"{r['offdiag_ratio_set_size']} | {r['max_offdiag_ratio_multiplicity']} | "
            f"{r['zeros']} |"
        )

    lines.append("")
    lines.append("## Q7798 E2-weighted trivial-character ratio contraction, p <= 500")
    lines.append("")
    lines.append(
        "Each p is computed in its unique power-of-two block `X < p <= 2X`, "
        "with `L=X^2`."
    )
    lines.append("")
    lines.append(
        "| X | p | z | E2 incidence mass | diagonal | centered offdiag | "
        "signed/diag | abs/diag | status |"
    )
    lines.append("|---:|---:|---:|---:|---:|---:|---:|---:|:---|")
    for r in e2_rows:
        lines.append(
            f"| {r['X']} | {r['p']} | {r['z']} | "
            f"{r['ordered_E2_pair_incidence_mass']} | "
            f"{md_float(float(r['diagonal_scale']))} | "
            f"{md_float(float(r['offdiag_centered']))} | "
            f"{md_float(float(r['signed_centered_over_diagonal']))} | "
            f"{md_float(float(r['abs_centered_over_diagonal']))} | "
            f"{r['status']} |"
        )

    lines.append("")
    lines.append("### Blockwise ARD_E2 ledger")
    lines.append("")
    lines.append("| X | outputs | sum diagonal | sum |centered| | abs ratio | signed ratio |")
    lines.append("|---:|---:|---:|---:|---:|---:|")
    for r in block_summaries:
        lines.append(
            f"| {r['X']} | {r['number_of_output_primes_z_ge_2_p_le_500']} | "
            f"{md_float(float(r['sum_diagonal_scale']))} | "
            f"{md_float(float(r['sum_abs_centered_offdiag']))} | "
            f"{md_float(float(r['ARD_abs_ratio']))} | "
            f"{md_float(float(r['ARD_signed_ratio']))} |"
        )

    lines.append("")
    lines.append(
        f"Across all valid blocks: sum diagonal = {md_float(total_diag, 10)}, "
        f"sum absolute centered = {md_float(total_abs_ctr, 10)}, "
        f"ratio = {md_float(total_abs_ctr / total_diag if total_diag else float('nan'), 10)}, "
        f"signed ratio = {md_float(total_signed_ctr / total_diag if total_diag else float('nan'), 10)}."
    )
    lines.append("")
    lines.append("## Numerical integrity checks")
    lines.append("")
    max_herm = max((float(r["hermitian_error"]) for r in e2_rows), default=0.0)
    max_imag = max((float(r["imaginary_contraction_error"]) for r in e2_rows), default=0.0)
    lines.append(f"Maximum Hermitian error in W_p: `{max_herm:.3e}`.")
    lines.append(f"Maximum imaginary contraction residue: `{max_imag:.3e}`.")
    lines.append("")
    lines.append("CSV files contain all floating-point Gram diagnostics.")

    summary = "\n".join(lines) + "\n"
    (OUT / "q7807_summary.md").write_text(summary, encoding="utf-8")
    print(summary)
    print(f"FINAL_RUNTIME_SECONDS={time.perf_counter() - t0:.6f}")


if __name__ == "__main__":
    main()
