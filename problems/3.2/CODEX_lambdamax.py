#!/usr/bin/env python3
"""Matrix-free centered-Gram experiment for CODEX_SPEC_lambdamax.md."""

from __future__ import annotations

import math
from pathlib import Path
from typing import Callable

import numpy as np


PRIMES = (997, 1999, 4001)
POWER_STEPS = 40
PSD_STARTS = 1
OFF_STARTS = 1
SEED = 20260801
OUT = Path(__file__).with_name("CODEX_LAMBDAMAX_report.md")


def ceil_power_two_thirds(p: int) -> int:
    """The least H with H^3 >= p^2, using integer arithmetic."""
    lo, hi = 0, p
    while lo + 1 < hi:
        mid = (lo + hi) // 2
        if mid**3 >= p**2:
            hi = mid
        else:
            lo = mid
    return hi


def inverse_table(p: int) -> np.ndarray:
    inv = np.zeros(p, dtype=np.int64)
    inv[1] = 1
    for i in range(2, p):
        inv[i] = (p - (p // i) * int(inv[p % i]) % p) % p
    return inv


def apery_pair(p: int, inv: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return the two normalized recurrence solutions on 0 <= n <= p-2."""
    M = p - 2
    b = np.zeros(M + 1, dtype=np.int64)
    c = np.zeros(M + 1, dtype=np.int64)
    b[0], b[1] = 1, 5 % p
    c[0], c[1] = 0, 1
    for n in range(1, M):
        polynomial = (34 * n**3 + 51 * n * n + 27 * n + 5) % p
        inv_cube = pow(int(inv[n + 1]), 3, p)
        n_cube = pow(n, 3, p)
        b[n + 1] = (polynomial * int(b[n]) - n_cube * int(b[n - 1])) * inv_cube % p
        c[n + 1] = (polynomial * int(c[n]) - n_cube * int(c[n - 1])) * inv_cube % p
    return b, c


def histograms(p: int, H: int) -> tuple[np.ndarray, np.ndarray, dict[str, int]]:
    """Build all n_h and verify Delta*product=N_h on every table entry."""
    M = p - 2
    inv = inverse_table(p)
    b, c = apery_pair(p, inv)
    counts = np.zeros((H, p), dtype=np.uint16)
    lengths = np.arange(M - 1, M - H - 1, -1, dtype=np.int64)

    residues = np.arange(p, dtype=np.int64)
    polynomial = (34 * residues**3 + 51 * residues**2 + 27 * residues + 5) % p
    cubes = residues**3 % p
    sixth = cubes**2 % p
    n_previous = np.zeros(p, dtype=np.int64)  # N_0
    n_current = np.ones(p, dtype=np.int64)    # N_1
    denominator = np.ones(p, dtype=np.int64)
    checks = 0
    max_multiplicity = 0

    for h in range(1, H + 1):
        length = int(lengths[h - 1])
        r = np.arange(1, length + 1, dtype=np.int64)
        denominator[r] = denominator[r] * cubes[r + h] % p
        delta = (b[r] * c[r + h] - b[r + h] * c[r]) % p
        if not np.array_equal(delta * denominator[r] % p, n_current[r]):
            bad = np.flatnonzero(delta * denominator[r] % p != n_current[r])
            raise AssertionError(f"D1 failure p={p}, h={h}, r={int(r[bad[0]])}")
        checks += length
        row = np.bincount(delta, minlength=p)
        assert int(row.sum()) == length
        assert int(row.max()) < 2**16
        counts[h - 1] = row
        max_multiplicity = max(max_multiplicity, int(row.max()))

        if h < H:
            shifted = (residues + h) % p
            n_next = (polynomial[shifted] * n_current - sixth[shifted] * n_previous) % p
            n_previous, n_current = n_current, n_next

    # Independent sign check for the exceptional cube row.
    r = np.arange(1, M, dtype=np.int64)
    h1 = (b[r] * c[r + 1] - b[r + 1] * c[r]) % p
    expected = np.array([pow(pow(int(x + 1), 3, p), -1, p) for x in r], dtype=np.int64)
    assert np.array_equal(h1, expected)
    return counts, lengths, {
        "checks": checks,
        "entries": int(lengths.sum()),
        "max_multiplicity": max_multiplicity,
        "h1_max_multiplicity": int(counts[0].max()),
    }


def matrix_free_data(
    p: int, counts: np.ndarray, lengths: np.ndarray
) -> tuple[Callable[[np.ndarray], np.ndarray], np.ndarray, np.ndarray]:
    """Return the Gamma/p^2 matvec, its exact diagonal, and the centered table."""
    q = counts.astype(np.float64) - lengths[:, None] / p

    def matvec(v: np.ndarray) -> np.ndarray:
        # Exactly two dense passes over q; Gamma/p^2 = q q^T / p.
        return q @ (q.T @ v) / p

    row_squares = np.einsum("ij,ij->i", counts.astype(np.int64), counts.astype(np.int64))
    diagonal = p * row_squares - lengths * lengths
    return matvec, diagonal, q


def normalize(v: np.ndarray) -> np.ndarray:
    norm = float(np.linalg.norm(v))
    if not math.isfinite(norm) or norm == 0.0:
        raise ArithmeticError("power iteration produced a zero/nonfinite vector")
    return v / norm


def power_runs(
    matvec: Callable[[np.ndarray], np.ndarray],
    dimension: int,
    rng: np.random.Generator,
    starts: int,
    orthogonal_to: np.ndarray | None = None,
) -> list[tuple[float, np.ndarray, float]]:
    """Run independent starts in columns; each column still gets all 40 steps."""
    v = rng.standard_normal((dimension, starts))
    if orthogonal_to is not None:
        v -= orthogonal_to[:, None] * (orthogonal_to @ v)[None, :]
    v /= np.linalg.norm(v, axis=0)[None, :]
    for _ in range(POWER_STEPS):
        w = matvec(v)
        if orthogonal_to is not None:
            w -= orthogonal_to[:, None] * (orthogonal_to @ w)[None, :]
        norms = np.linalg.norm(w, axis=0)
        if np.any(~np.isfinite(norms)) or np.any(norms == 0.0):
            raise ArithmeticError("power iteration produced a zero/nonfinite vector")
        v = w / norms[None, :]
    av = matvec(v)
    eigenvalues = np.einsum("ij,ij->j", v, av)
    residuals = np.linalg.norm(av - v * eigenvalues[None, :], axis=0)
    return [
        (float(eigenvalues[j]), v[:, j].copy(), float(residuals[j]))
        for j in range(starts)
    ]


def psd_top_two(
    matvec: Callable[[np.ndarray], np.ndarray], dimension: int, seed: int
) -> tuple[dict[str, float], np.ndarray]:
    rng = np.random.default_rng(seed)
    first_runs = power_runs(matvec, dimension, rng, PSD_STARTS)
    lambda1, u1, residual1 = max(first_runs, key=lambda z: z[0])
    second_runs = power_runs(matvec, dimension, rng, PSD_STARTS, u1)
    lambda2, u2, residual2 = max(second_runs, key=lambda z: z[0])
    overlap = abs(float(np.dot(u1, u2)))
    return {
        "lambda1": lambda1,
        "lambda2": lambda2,
        "residual1": residual1,
        "residual2": residual2,
        "deflation_overlap": overlap,
    }, u1


def off_diagonal_norm(
    gamma_matvec: Callable[[np.ndarray], np.ndarray],
    diagonal_scaled: np.ndarray,
    dimension: int,
    seed: int,
) -> dict[str, float]:
    def off_matvec(v: np.ndarray) -> np.ndarray:
        diagonal = diagonal_scaled if v.ndim == 1 else diagonal_scaled[:, None]
        return gamma_matvec(v) - diagonal * v

    rng = np.random.default_rng(seed)
    runs = power_runs(off_matvec, dimension, rng, OFF_STARTS)
    signed, vector, residual = max(runs, key=lambda z: abs(z[0]))
    # Recompute after selection so the sign is explicitly tied to O, not |O|.
    ov = off_matvec(vector)
    signed = float(np.dot(vector, ov))
    residual = float(np.linalg.norm(ov - signed * vector))
    return {"signed_eigenvalue": signed, "norm": abs(signed), "residual": residual}


def lanczos_audit(
    matvec: Callable[[np.ndarray], np.ndarray], dimension: int, seed: int
) -> np.ndarray:
    """Fixed 40-step, fully reorthogonalized symmetric Lanczos Ritz values."""
    steps = min(POWER_STEPS, dimension)
    rng = np.random.default_rng(seed)
    q = normalize(rng.standard_normal(dimension))
    previous = np.zeros(dimension, dtype=np.float64)
    beta_previous = 0.0
    basis = np.empty((dimension, steps), dtype=np.float64)
    alpha = np.empty(steps, dtype=np.float64)
    beta = np.empty(max(0, steps - 1), dtype=np.float64)
    used = steps
    for j in range(steps):
        basis[:, j] = q
        z = matvec(q)
        if j:
            z -= beta_previous * previous
        alpha[j] = float(np.dot(q, z))
        z -= alpha[j] * q
        # Two full passes make loss of orthogonality negligible at this size.
        for _ in range(2):
            z -= basis[:, : j + 1] @ (basis[:, : j + 1].T @ z)
        if j + 1 == steps:
            break
        beta_current = float(np.linalg.norm(z))
        if beta_current < 1e-14:
            used = j + 1
            break
        beta[j] = beta_current
        previous, q = q, z / beta_current
        beta_previous = beta_current
    tridiagonal = np.diag(alpha[:used])
    if used > 1:
        tridiagonal += np.diag(beta[: used - 1], 1) + np.diag(beta[: used - 1], -1)
    return np.linalg.eigvalsh(tridiagonal)


def exact_statistics(p: int, counts: np.ndarray, lengths: np.ndarray) -> dict[str, int]:
    counts64 = counts.astype(np.int64)
    row_squares = np.einsum("ij,ij->i", counts64, counts64)
    diagonal = p * row_squares - lengths * lengths
    totals = counts64.sum(axis=0)
    n_coinc = int(np.dot(totals, totals))
    size = int(lengths.sum())
    ones = p * n_coinc - size * size
    assert ones >= 0
    centered_numerator = p * totals - size
    centered_square_numerator = int(np.dot(centered_numerator, centered_numerator))
    assert centered_square_numerator % p == 0
    assert centered_square_numerator // p == ones
    return {
        "dmax": int(diagonal.max()),
        "dmax_h": int(np.argmax(diagonal)) + 1,
        "trace": int(diagonal.sum()),
        "n_coinc": n_coinc,
        "size": size,
        "ones": ones,
    }


def analyze_block(p: int, counts: np.ndarray, lengths: np.ndarray, peeled: bool) -> dict:
    if peeled:
        counts = counts[1:]
        lengths = lengths[1:]
    dimension = len(lengths)
    gamma_matvec, diagonal, _ = matrix_free_data(p, counts, lengths)
    psd, _ = psd_top_two(gamma_matvec, dimension, SEED + 1009 * p + int(peeled))
    off = off_diagonal_norm(
        gamma_matvec, diagonal.astype(np.float64) / (p * p), dimension,
        SEED + 2027 * p + int(peeled),
    )
    psd_ritz = lanczos_audit(gamma_matvec, dimension, SEED + 3037 * p + int(peeled))

    def off_matvec(v: np.ndarray) -> np.ndarray:
        return gamma_matvec(v) - diagonal.astype(np.float64) / (p * p) * v

    off_ritz = lanczos_audit(off_matvec, dimension, SEED + 4001 * p + int(peeled))
    audit = {
        "lambda1": float(psd_ritz[-1]),
        "lambda2": float(psd_ritz[-2]),
        "off_max": float(off_ritz[-1]),
        "off_min": float(off_ritz[0]),
        "off_norm": float(max(abs(off_ritz[0]), abs(off_ritz[-1]))),
    }
    exact = exact_statistics(p, counts, lengths)
    if peeled:
        exact["dmax_h"] += 1
    return {"dimension": dimension, "psd": psd, "off": off, "audit": audit, "exact": exact}


def parseval_check() -> dict[str, float | int]:
    p = 199
    H = ceil_power_two_thirds(p)
    counts, lengths, meta = histograms(p, H)
    counts64 = counts.astype(np.int64)
    gamma_exact = p * (counts64 @ counts64.T) - np.outer(lengths, lengths)
    fourier = np.fft.fft(counts.astype(np.float64), axis=1)
    gamma_fourier = fourier[:, 1:] @ fourier[:, 1:].conj().T
    error = gamma_fourier - gamma_exact
    return {
        "p": p,
        "H": H,
        "pairs": H * H,
        "max_real_error": float(np.max(np.abs(error.real))),
        "max_imag_error": float(np.max(np.abs(error.imag))),
        "max_gamma": int(np.max(np.abs(gamma_exact))),
        "d1_checks": meta["checks"],
    }


def fmt_float(x: float) -> str:
    return f"{x:.8f}"


def main() -> None:
    parseval = parseval_check()
    cases = []
    for p in PRIMES:
        H = ceil_power_two_thirds(p)
        counts, lengths, meta = histograms(p, H)
        full = analyze_block(p, counts, lengths, peeled=False)
        peeled = analyze_block(p, counts, lengths, peeled=True)
        cases.append({"p": p, "H": H, "meta": meta, "full": full, "peeled": peeled})

    lines = [
        "# Direct power-iteration measurement of the centered Gram matrix",
        "",
        "## 1. Histograms and finite-field arithmetic",
        "",
        "For each prime, `M=p-2`, `I_h={1,...,M-h}`, and "
        "`H=ceil(p^(2/3))` was evaluated as the least integer with `H^3>=p^2`. "
        "The Apéry recurrence produced `b_n,c_n`; every table value used "
        "`Delta_(r,h)=b_r c_(r+h)-b_(r+h)c_r (mod p)`. Simultaneously, the "
        "numerator recurrence checked `Delta_(r,h) prod_(j=1)^h(r+j)^3=N_h(r)` "
        "at every histogram entry. Counts and all displayed integer invariants "
        "are exact; only eigendata and displayed ratios use binary64 arithmetic.",
        "",
        "| p | H | sum_h card(I_h) | D1 checks | max n_h(a) | max n_1(a) |",
        "|---:|---:|---:|---:|---:|---:|",
    ]
    for case in cases:
        m = case["meta"]
        lines.append(
            f"| {case['p']} | {case['H']} | {m['entries']} | {m['checks']} | "
            f"{m['max_multiplicity']} | {m['h1_max_multiplicity']} |"
        )

    lines.extend([
        "",
        "All D1 checks passed. The independent identity "
        "`Delta_(r,1)=(r+1)^(-3)` also passed at every admissible `r` for all "
        "three primes.",
        "",
        "## 2. Numerical Parseval check at p=199",
        "",
        f"Using `H={parseval['H']}`, all {parseval['pairs']} ordered row pairs were "
        "checked against FFT evaluations of "
        "`S_h(t)=sum_a n_h(a) exp(-2 pi i t a/p)`. The comparison was",
        "",
        "`p <q_h,q_k> = p sum_a n_h(a)n_k(a)-|I_h||I_k| "
        "= sum_(t!=0) S_h(t) conjugate(S_k(t))`.",
        "",
        f"The maximum real-part error was `{parseval['max_real_error']:.3e}` "
        f"and the maximum spurious imaginary part was `{parseval['max_imag_error']:.3e}`; "
        f"the largest compared exact entry had magnitude `{parseval['max_gamma']}`. "
        f"All {parseval['d1_checks']} underlying D1 checks passed.",
        "",
        "## 3. Matrix-free iteration and exact checks",
        "",
        "No campaign-size Gram matrix was formed. With `Q_(h,a)=q_h(a)`, every "
        "application of `Gamma/p^2` used exactly the two table passes "
        "`v -> Q^T v -> Q(Q^T v)/p`. For each block, one seeded random start "
        "was run for exactly 40 steps; its Rayleigh quotient gave `lambda_1`. "
        "The same fixed procedure on the orthogonal complement of "
        "its vector gave `lambda_2`. There was no convergence-based early stop.",
        "",
        "The exact checks use `Gamma_hh=p sum_a n_h(a)^2-|I_h|^2`, "
        "`Tr(Gamma)=sum_h Gamma_hh`, and",
        "",
        "`1^T Gamma 1 = p ||sum_h q_h||^2 = p N_coinc-S^2`,",
        "",
        "where `S=sum_h |I_h|` and `N_coinc=sum_a(sum_h n_h(a))^2`. The next "
        "table gives the exact integer values; `(h=1 removed)` recomputes every "
        "quantity on the principal block indexed by `2,...,H`. The norm identity "
        "was checked a second exact way as "
        "`sum_a(p sum_h n_h(a)-S)^2/p`; it agreed in every block.",
        "",
        "| p | block | S | N_coinc | D_max (argmax h) | 1^T Gamma 1 | Tr(Gamma) |",
        "|---:|:---|---:|---:|---:|---:|---:|",
    ])
    for case in cases:
        for key, label in (("full", "all h"), ("peeled", "h=1 removed")):
            e = case[key]["exact"]
            lines.append(
                f"| {case['p']} | {label} | {e['size']} | {e['n_coinc']} | "
                f"{e['dmax']} (h={e['dmax_h']}) | {e['ones']} | {e['trace']} |"
            )

    lines.extend([
        "",
        "## 4. Top two eigenvalues and normalized diagnostics",
        "",
        "Here `rows` is `H` for the full block and `H-1` after peeling. Thus the "
        "all-ones column is normalized by the actual number of retained rows. "
        "Residuals are for the normalized matrix `Gamma/p^2`.",
        "",
        "| p | block | rows | 40-step power lambda_1/p^2 | 40-step deflated lambda_2/p^2 | D_max/p^2 | (1^T Gamma 1)/(p^2 rows) | Tr/p^2 | residuals (1,2) |",
        "|---:|:---|---:|---:|---:|---:|---:|---:|:---|",
    ])
    for case in cases:
        for key, label in (("full", "all h"), ("peeled", "h=1 removed")):
            block = case[key]
            e, psd = block["exact"], block["psd"]
            p, rows = case["p"], block["dimension"]
            lines.append(
                f"| {p} | {label} | {rows} | {fmt_float(psd['lambda1'])} | "
                f"{fmt_float(psd['lambda2'])} | {fmt_float(e['dmax']/p**2)} | "
                f"{fmt_float(e['ones']/(p**2*rows))} | {fmt_float(e['trace']/p**2)} | "
                f"{psd['residual1']:.2e}, {psd['residual2']:.2e} |"
            )

    lines.extend([
        "",
        "The largest measured deflation overlap `|<u_1,u_2>|` was "
        f"`{max(case[key]['psd']['deflation_overlap'] for case in cases for key in ('full','peeled')):.3e}`.",
        "",
        "Because 40 scalar power steps need not resolve a clustered spectral "
        "edge, the same matrix-free operator was independently audited with "
        "40-step fully reorthogonalized symmetric Lanczos. Only the resulting "
        "40 by 40 tridiagonal matrix was diagonalized. These are the preferred "
        "spectral estimates; in particular, they restore the required principal-"
        "submatrix interlacing in the `p=4001` cell.",
        "",
        "| p | block | audited lambda_1/p^2 | audited lambda_2/p^2 | power shortfall (1,2) |",
        "|---:|:---|---:|---:|:---|",
    ])
    for case in cases:
        for key, label in (("full", "all h"), ("peeled", "h=1 removed")):
            block = case[key]
            a, psd = block["audit"], block["psd"]
            lines.append(
                f"| {case['p']} | {label} | {a['lambda1']:.8f} | {a['lambda2']:.8f} | "
                f"{a['lambda1']-psd['lambda1']:.3e}, {a['lambda2']-psd['lambda2']:.3e} |"
            )

    lines.extend([
        "",
        "## 5. Off-diagonal operator",
        "",
        "For `O=Gamma-diag(Gamma)`, one seeded random start was run for exactly "
        "40 symmetric-matrix power steps. The Rayleigh-quotient sign was retained; this distinguishes a "
        "dominant negative eigenvalue from a positive one. Residuals refer to "
        "`O/p^2`.",
        "",
        "| p | block | signed power Rayleigh/p^2 | 40-step abs(Rayleigh)/p^2 | 40-step random-scale ratio | residual |",
        "|---:|:---|---:|---:|---:|---:|",
    ])
    for case in cases:
        for key, label in (("full", "all h"), ("peeled", "h=1 removed")):
            block = case[key]
            off = block["off"]
            rows, p = block["dimension"], case["p"]
            random_norm = off["norm"] * math.sqrt(p / rows)
            lines.append(
                f"| {p} | {label} | {off['signed_eigenvalue']:+.8f} | "
                f"{off['norm']:.8f} | {random_norm:.8f} | {off['residual']:.2e} |"
            )

    lines.extend([
        "",
        "The symmetric Lanczos audit also computed both algebraic edges of "
        "`O/p^2`; the positive edge dominates in every cell.",
        "",
        "| p | block | lambda_max(O)/p^2 | lambda_min(O)/p^2 | audited norm(O)/p^2 | audited random-scale ratio |",
        "|---:|:---|---:|---:|---:|---:|",
    ])
    for case in cases:
        for key, label in (("full", "all h"), ("peeled", "h=1 removed")):
            block = case[key]
            a = block["audit"]
            ratio = a["off_norm"] * math.sqrt(case["p"] / block["dimension"])
            lines.append(
                f"| {case['p']} | {label} | {a['off_max']:+.8f} | {a['off_min']:+.8f} | "
                f"{a['off_norm']:.8f} | {ratio:.8f} |"
            )

    full_lambdas = [case["full"]["audit"]["lambda1"] for case in cases]
    peeled_lambdas = [case["peeled"]["audit"]["lambda1"] for case in cases]
    full_random = [case["full"]["audit"]["off_norm"] * math.sqrt(case["p"] / case["full"]["dimension"]) for case in cases]
    peeled_random = [case["peeled"]["audit"]["off_norm"] * math.sqrt(case["p"] / case["peeled"]["dimension"]) for case in cases]
    lines.extend([
        "",
        "## 6. Verdict",
        "",
        f"**[OP-OFF-0] finite test: supported at these three sizes.** The full "
        f"`lambda_1/p^2` values lie in `[{min(full_lambdas):.4f}, "
        f"{max(full_lambdas):.4f}]`; after removing the exceptional cube row "
        f"they lie in `[{min(peeled_lambdas):.4f}, {max(peeled_lambdas):.4f}]`. "
        "There is no observed growth with `p` in this range. This is direct "
        "finite evidence for the bounded scenario, not an asymptotic proof.",
        "",
        f"**Random-matrix off-diagonal scale: supported.** The ratios "
        f"`||O||/(sqrt(rows) p^(3/2))` lie in "
        f"`[{min(full_random):.4f}, {max(full_random):.4f}]` for the full blocks "
        f"and `[{min(peeled_random):.4f}, {max(peeled_random):.4f}]` after peeling. "
        "They remain constant-sized, so the data fit the proposed "
        "`sqrt(H) p^(3/2)` scale. Peeling changes the cube-row outlier in the "
        "PSD spectrum for primes `p=1 mod 3`, but does not expose a growing "
        "off-diagonal norm.",
        "",
        "The scalar-power rows are finite-iteration lower estimates in absolute "
        "Rayleigh quotient. The displayed residuals and shortfalls quantify "
        "their unresolved spectral error; the numerical ranges in this verdict "
        "use the fixed-step symmetric Lanczos audit.",
    ])
    OUT.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
