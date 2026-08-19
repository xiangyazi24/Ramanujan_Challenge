#!/usr/bin/env python3
"""Q547: exact degrees and finite-field root counts for Apéry transfer T_12.

This is dependency-free and uses the same continuant normalization as
problems/3.2/fiber_verify.py.

For P(n)=34n^3+51n^2+27n+5 and
  M(n) = [[P(n)/(n+1)^3, -n^3/(n+1)^3], [1, 0]],
write
  T_h(r) = M(r+h-1)...M(r).
Then
  T_h(r)[0,1] = -r^3 N_h(r) / prod_{j=1}^h (r+j)^3,
where
  N_0=0, N_1=1,
  N_{h+1}(r)=P(r+h)N_h(r)-(r+h)^6N_{h-1}(r).

The script counts roots only where the rational transfer matrix is defined:
r=0,...,p-h-1.  It also reports the project-relevant nontrivial count on
r=1,...,p-h-1, after removing the universal r^3 factor.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from statistics import mean, median
from math import isqrt

Poly = tuple[int, ...]  # coefficients in increasing powers
ZERO: Poly = (0,)
ONE: Poly = (1,)
P_BASE: Poly = (5, 27, 51, 34)
MAX_H = 20
LOW_P = 100
HIGH_P = 500


def trim(values: list[int] | tuple[int, ...]) -> Poly:
    out = list(values)
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return tuple(out)


def poly_add(f: Poly, g: Poly) -> Poly:
    out = [0] * max(len(f), len(g))
    for i, a in enumerate(f):
        out[i] += a
    for i, b in enumerate(g):
        out[i] += b
    return trim(out)


def poly_neg(f: Poly) -> Poly:
    return trim([-a for a in f])


def poly_sub(f: Poly, g: Poly) -> Poly:
    return poly_add(f, poly_neg(g))


def poly_mul(f: Poly, g: Poly) -> Poly:
    if f == ZERO or g == ZERO:
        return ZERO
    out = [0] * (len(f) + len(g) - 1)
    for i, a in enumerate(f):
        for j, b in enumerate(g):
            out[i + j] += a * b
    return trim(out)


def poly_pow(f: Poly, exponent: int) -> Poly:
    result = ONE
    base = f
    n = exponent
    while n:
        if n & 1:
            result = poly_mul(result, base)
        base = poly_mul(base, base)
        n >>= 1
    return result


def linear(shift: int) -> Poly:
    return (shift, 1)


def poly_shift(f: Poly, shift: int) -> Poly:
    """Return f(x+shift) over Z."""
    result = ZERO
    for coefficient in reversed(f):
        result = poly_add(poly_mul(result, linear(shift)), (coefficient,))
    return result


def poly_eval_mod(f: Poly, value: int, modulus: int) -> int:
    result = 0
    for coefficient in reversed(f):
        result = (result * value + coefficient) % modulus
    return result


def degree_mod(f: Poly, modulus: int) -> int:
    for i in range(len(f) - 1, -1, -1):
        if f[i] % modulus:
            return i
    return -1


def build_n_polynomials(max_h: int) -> list[Poly]:
    n_polys = [ZERO for _ in range(max_h + 1)]
    n_polys[0] = ZERO
    n_polys[1] = ONE
    for h in range(1, max_h):
        n_polys[h + 1] = poly_sub(
            poly_mul(poly_shift(P_BASE, h), n_polys[h]),
            poly_mul(poly_pow(linear(h), 6), n_polys[h - 1]),
        )
    return n_polys


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    d = 3
    while d <= isqrt(n):
        if n % d == 0:
            return False
        d += 2
    return True


def primes_between(low: int, high: int) -> list[int]:
    return [n for n in range(low, high + 1) if is_prime(n)]


def p_value(n: int, modulus: int) -> int:
    return (34 * n**3 + 51 * n**2 + 27 * n + 5) % modulus


def matmul2(left: tuple[tuple[int, int], tuple[int, int]],
            right: tuple[tuple[int, int], tuple[int, int]],
            p: int) -> tuple[tuple[int, int], tuple[int, int]]:
    return (
        (
            (left[0][0] * right[0][0] + left[0][1] * right[1][0]) % p,
            (left[0][0] * right[0][1] + left[0][1] * right[1][1]) % p,
        ),
        (
            (left[1][0] * right[0][0] + left[1][1] * right[1][0]) % p,
            (left[1][0] * right[0][1] + left[1][1] * right[1][1]) % p,
        ),
    )


def verify_transfer_identity(n_polys: list[Poly], primes: list[int]) -> int:
    """Directly verify D_h*T_12=-r^3*N_h for every tested p,r,h."""
    checks = 0
    identity = ((1, 0), (0, 1))
    for p in primes:
        for r in range(p):
            transfer = identity
            denominator_product = 1
            max_h_here = min(MAX_H, p - 1 - r)
            for h in range(1, max_h_here + 1):
                n = r + h - 1
                denominator = pow(n + 1, 3, p)
                inv_denominator = pow(denominator, -1, p)
                matrix = (
                    (
                        p_value(n, p) * inv_denominator % p,
                        (-pow(n, 3, p) * inv_denominator) % p,
                    ),
                    (1, 0),
                )
                transfer = matmul2(matrix, transfer, p)
                denominator_product = denominator_product * denominator % p
                lhs = transfer[0][1] * denominator_product % p
                rhs = (-pow(r, 3, p) * poly_eval_mod(n_polys[h], r, p)) % p
                assert lhs == rhs, (p, r, h, lhs, rhs)
                checks += 1
    return checks


def format_float(value: float) -> str:
    return f"{value:.3f}"


def markdown_csv_matrix(primes: list[int], by_prime: dict[int, list[int]]) -> str:
    lines = ["p," + ",".join(f"h{h}" for h in range(1, MAX_H + 1))]
    for p in primes:
        lines.append(str(p) + "," + ",".join(str(v) for v in by_prime[p]))
    return "\n".join(lines)


def run(output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    primes = primes_between(LOW_P, HIGH_P)
    assert len(primes) == 70, len(primes)

    n_polys = build_n_polynomials(MAX_H)
    leading = [0] * (MAX_H + 1)
    leading[0] = 0
    leading[1] = 1
    for h in range(1, MAX_H + 1):
        assert len(n_polys[h]) - 1 == 3 * (h - 1)
        leading[h] = n_polys[h][-1]
        if h >= 2:
            assert leading[h] == 34 * leading[h - 1] - leading[h - 2]
        assert leading[h] > 0

    transfer_checks = verify_transfer_identity(n_polys, primes)

    rows: list[dict[str, int | float | bool]] = []
    gap_matrix: dict[int, list[int]] = {p: [] for p in primes}
    t12_matrix: dict[int, list[int]] = {p: [] for p in primes}
    exceptions: list[dict[str, int]] = []

    for h in range(1, MAX_H + 1):
        n_poly = n_polys[h]
        generic_deg_n = 3 * (h - 1)
        generic_deg_q = 3 * h
        for p in primes:
            deg_n_p = degree_mod(n_poly, p)
            assert deg_n_p >= 0
            deg_q_p = deg_n_p + 3
            if deg_n_p != generic_deg_n:
                exceptions.append(
                    {
                        "h": h,
                        "p": p,
                        "generic_deg_N": generic_deg_n,
                        "deg_N_mod_p": deg_n_p,
                        "generic_deg_Q": generic_deg_q,
                        "deg_Q_mod_p": deg_q_p,
                    }
                )

            # Admissible values are precisely r=0,...,p-h-1; the h residues
            # -1,...,-h are poles of the rational transfer product.
            gap_roots = sum(
                poly_eval_mod(n_poly, r, p) == 0 for r in range(1, p - h)
            )
            t12_roots = 1 + gap_roots  # universal r=0 root from r^3
            direct_t12_roots = sum(
                (r == 0 or poly_eval_mod(n_poly, r, p) == 0)
                for r in range(0, p - h)
            )
            assert direct_t12_roots == t12_roots

            n_roots_including_zero = sum(
                poly_eval_mod(n_poly, r, p) == 0 for r in range(0, p - h)
            )
            pole_roots = sum(
                poly_eval_mod(n_poly, (-j) % p, p) == 0
                for j in range(1, h + 1)
            )

            gap_matrix[p].append(gap_roots)
            t12_matrix[p].append(t12_roots)
            rows.append(
                {
                    "h": h,
                    "p": p,
                    "generic_deg_N": generic_deg_n,
                    "deg_N_mod_p": deg_n_p,
                    "generic_deg_Q": generic_deg_q,
                    "deg_Q_mod_p": deg_q_p,
                    "gap_roots_r_1_to_p_h_minus_1": gap_roots,
                    "T12_roots_admissible_including_r0": t12_roots,
                    "N_roots_admissible_including_r0": n_roots_including_zero,
                    "N_roots_at_excluded_poles": pole_roots,
                    "factor3_vs_N_degree": (h == 1) or (3 * gap_roots <= deg_n_p),
                    "factor3_vs_Q_degree": 3 * t12_roots <= deg_q_p,
                    "gap_root_fraction_of_N_degree": (
                        0.0 if deg_n_p == 0 else gap_roots / deg_n_p
                    ),
                    "T12_root_fraction_of_Q_degree": t12_roots / deg_q_p,
                }
            )

    # Long-form CSV.
    csv_path = output_dir / "q547-root-counts-long.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    # Prime-by-h matrices.
    gap_path = output_dir / "q547-gap-root-matrix.csv"
    t12_path = output_dir / "q547-t12-root-matrix.csv"
    for path, matrix in ((gap_path, gap_matrix), (t12_path, t12_matrix)):
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(["p"] + [f"h{h}" for h in range(1, MAX_H + 1)])
            for p in primes:
                writer.writerow([p] + matrix[p])

    summaries: list[dict[str, object]] = []
    for h in range(1, MAX_H + 1):
        h_rows = [row for row in rows if row["h"] == h]
        gaps = [int(row["gap_roots_r_1_to_p_h_minus_1"]) for row in h_rows]
        t_roots = [int(row["T12_roots_admissible_including_r0"]) for row in h_rows]
        max_gap = max(gaps)
        max_t = max(t_roots)
        summaries.append(
            {
                "h": h,
                "deg_N": 3 * (h - 1),
                "deg_Q": 3 * h,
                "leading_coefficient_N": leading[h],
                "mean_gap_roots": mean(gaps),
                "median_gap_roots": median(gaps),
                "min_gap_roots": min(gaps),
                "max_gap_roots": max_gap,
                "primes_at_max_gap": [
                    int(row["p"])
                    for row in h_rows
                    if int(row["gap_roots_r_1_to_p_h_minus_1"]) == max_gap
                ],
                "mean_T12_roots": mean(t_roots),
                "max_T12_roots": max_t,
                "factor3_N_count": sum(bool(row["factor3_vs_N_degree"]) for row in h_rows),
                "factor3_Q_count": sum(bool(row["factor3_vs_Q_degree"]) for row in h_rows),
                "mean_gap_fraction_of_N_degree": mean(
                    float(row["gap_root_fraction_of_N_degree"]) for row in h_rows
                ),
                "mean_T12_fraction_of_Q_degree": mean(
                    float(row["T12_root_fraction_of_Q_degree"]) for row in h_rows
                ),
            }
        )

    machine = {
        "parameters": {
            "h_min": 1,
            "h_max": MAX_H,
            "prime_min": LOW_P,
            "prime_max": HIGH_P,
            "prime_count": len(primes),
            "primes": primes,
            "admissible_r": "0 <= r <= p-h-1",
            "gap_r": "1 <= r <= p-h-1",
        },
        "identity": "T12=-r^3*N_h/prod_{j=1}^h(r+j)^3",
        "transfer_identity_checks": transfer_checks,
        "degree_exceptions": exceptions,
        "summaries": summaries,
    }
    (output_dir / "q547-summary.json").write_text(
        json.dumps(machine, indent=2) + "\n", encoding="utf-8"
    )

    degree_lines = []
    for h in range(1, MAX_H + 1):
        exc = [e for e in exceptions if e["h"] == h]
        exc_text = "none" if not exc else "; ".join(
            f"p={e['p']}: {e['generic_deg_N']}->{e['deg_N_mod_p']}"
            for e in exc
        )
        degree_lines.append(
            f"| {h} | {3*(h-1)} | {3*h} | {leading[h]} | {exc_text} |"
        )

    summary_lines = []
    for item in summaries:
        max_primes = ",".join(str(p) for p in item["primes_at_max_gap"])
        summary_lines.append(
            "| {h} | {deg_N} | {mean_gap} | {median_gap} | {max_gap} | {max_primes} | "
            "{factor_n}/70 | {mean_t} | {factor_q}/70 |".format(
                h=item["h"],
                deg_N=item["deg_N"],
                mean_gap=format_float(float(item["mean_gap_roots"])),
                median_gap=format_float(float(item["median_gap_roots"])),
                max_gap=item["max_gap_roots"],
                max_primes=max_primes,
                factor_n=item["factor3_N_count"],
                mean_t=format_float(float(item["mean_T12_roots"])),
                factor_q=item["factor3_Q_count"],
            )
        )

    nontrivial_rows = [row for row in rows if int(row["h"]) >= 2]
    global_factor3_n = sum(bool(row["factor3_vs_N_degree"]) for row in nontrivial_rows)
    global_pairs_n = len(nontrivial_rows)
    global_factor3_q = sum(bool(row["factor3_vs_Q_degree"]) for row in rows)
    global_pairs_q = len(rows)
    all_gap_values = [int(row["gap_roots_r_1_to_p_h_minus_1"]) for row in nontrivial_rows]

    source = Path(__file__).read_text(encoding="utf-8")
    report = f"""# Q547 Apéry transfer-root computation

## Correction and normalization

The normalized matrix entry is a rational function, not literally a polynomial:

`T12(r,h) = -r^3 N_h(r) / product_(j=1)^h (r+j)^3`.

Thus there are two exact polynomial degrees:

- full cleared numerator `Q_h(r)=-r^3 N_h(r)`: `deg Q_h=3h`;
- nontrivial gap/continuant polynomial after removing the universal `r^3` factor: `deg N_h=3(h-1)`.

The recurrence is

`N_0=0`, `N_1=1`, and
`N_(h+1)(r)=P(r+h)N_h(r)-(r+h)^6N_(h-1)(r)`,
where `P(x)=34x^3+51x^2+27x+5`.

The leading coefficient `ell_h` obeys `ell_0=0`, `ell_1=1`,
`ell_(h+1)=34 ell_h-ell_(h-1)`, so it is positive and nonzero over Z.
The direct normalized-matrix identity was checked at every admissible `(p,r,h)`:
**{transfer_checks:,} exact modular checks**.

## Counting convention

For a prime `p` and `h<p`, the rational transfer product is defined at exactly
`r=0,...,p-h-1`; the omitted residues `-1,...,-h` are poles.  The table reports:

- `gap roots`: roots of `N_h` on `r=1,...,p-h-1` (the project-relevant, nontrivial count);
- `T12 roots`: roots of the rational `T12` on `r=0,...,p-h-1`.

Because `r=0` is always a root of the universal `r^3` factor,
`T12 roots = 1 + gap roots` under this convention.

## Exact degree table

The last column records every degree drop of `N_h mod p` for primes in `[100,500]`.
The corresponding full-numerator degree is always three larger.

| h | deg N_h over Z | deg Q_h over Z | leading coefficient ell_h | mod-p degree drops |
|---:|---:|---:|---:|:---|
{chr(10).join(degree_lines)}

## Root-count summary over all 70 primes in [100,500]

`factor >=3` means the number of distinct roots is at most one third of the
actual polynomial degree over that field.

| h | deg N_h | mean gap roots | median | max | primes attaining max | factor >=3 vs N_h | mean T12 roots | factor >=3 vs Q_h |
|---:|---:|---:|---:|---:|:---|---:|---:|---:|
{chr(10).join(summary_lines)}

Across the {global_pairs_n} nontrivial `(h,p)` pairs (`h=2,...,20`),
{global_factor3_n} satisfy `#gap roots <= deg_Fp(N_h)/3`.
Across all {global_pairs_q} pairs, {global_factor3_q} satisfy
`#T12 roots <= deg_Fp(Q_h)/3`.
The mean nontrivial root count over all `h=2,...,20` and all tested primes is
{mean(all_gap_values):.3f}; the maximum is {max(all_gap_values)}.

## Full per-prime gap-root matrix

Each entry is `#{{1 <= r <= p-h-1 : N_h(r)=0 mod p}}`.

```csv
{markdown_csv_matrix(primes, gap_matrix)}
```

## Full per-prime T12-root matrix

Each entry is `#{{0 <= r <= p-h-1 : T12(r,h)=0 mod p}}`.

```csv
{markdown_csv_matrix(primes, t12_matrix)}
```

## Reproducible computation

```python
{source}
```
"""
    (output_dir / "q547-report.md").write_text(report, encoding="utf-8")

    print(f"primes={len(primes)} rows={len(rows)} transfer_checks={transfer_checks}")
    print(f"degree_exceptions={len(exceptions)}")
    print(f"global_factor3_N={global_factor3_n}/{global_pairs_n}")
    print(f"global_factor3_Q={global_factor3_q}/{global_pairs_q}")
    print(f"mean_gap={mean(all_gap_values):.6f} max_gap={max(all_gap_values)}")
    print(f"report={output_dir / 'q547-report.md'}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("problems/3.2/research/computations/q547"),
    )
    args = parser.parse_args()
    run(args.output_dir)


if __name__ == "__main__":
    main()
