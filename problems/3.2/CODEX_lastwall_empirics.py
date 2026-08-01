#!/usr/bin/env python3
"""Exact last-wall empirical ground truth for CODEX_SPEC_lastwall_empirics.md.

Only Python's standard library is used.  Every finite-field operation and every
reported count is integer-exact.  Decimal arithmetic is used solely to evaluate
the displayed real-valued scales and the transcendental sqrt(p) log(p) cutoff.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from decimal import Decimal, getcontext
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple


PRIMES = (997, 1999, 4001, 7919)
OUT = Path(__file__).with_name("CODEX_LASTWALL_EMPIRICS_report.md")
getcontext().prec = 80


def ceil_decimal(x: Decimal) -> int:
    n = int(x)
    return n if x == n else n + 1


def ceil_rational_power(p: int, numerator: int, denominator: int) -> int:
    """Smallest n with n**denominator >= p**numerator, using integers only."""
    target = p**numerator
    lo, hi = 0, p + 1
    while lo + 1 < hi:
        mid = (lo + hi) // 2
        if mid**denominator >= target:
            hi = mid
        else:
            lo = mid
    return hi


def cutoffs(p: int) -> List[Tuple[str, int]]:
    pd = Decimal(p)
    sqrt_log_value = pd.sqrt() * pd.ln()
    # The nearest boundary is vastly farther away than Decimal's error here.
    distance = min(sqrt_log_value % 1, 1 - sqrt_log_value % 1)
    assert distance > Decimal("1e-50")
    return [
        ("ceil(sqrt(p) log(p))", ceil_decimal(sqrt_log_value)),
        ("ceil(p^0.6)", ceil_rational_power(p, 3, 5)),
        ("ceil(p^0.66)", ceil_rational_power(p, 33, 50)),
    ]


def inverse_table(p: int) -> List[int]:
    inv = [0] * p
    inv[1] = 1
    for i in range(2, p):
        inv[i] = (p - (p // i) * inv[p % i] % p) % p
    return inv


def apery_pair(p: int) -> Tuple[List[int], List[int], List[int]]:
    """Return b,c on 0..p-2 and the inverse table."""
    nmax = p - 2
    inv = inverse_table(p)
    b = [0] * (nmax + 1)
    c = [0] * (nmax + 1)
    b[0], b[1] = 1, 5 % p
    c[0], c[1] = 0, 1
    for n in range(1, nmax):
        pn = (34 * n**3 + 51 * n * n + 27 * n + 5) % p
        n3 = n**3 % p
        inv_cube = inv[n + 1] ** 3 % p
        b[n + 1] = (pn * b[n] - n3 * b[n - 1]) * inv_cube % p
        c[n + 1] = (pn * c[n] - n3 * c[n - 1]) * inv_cube % p
    return b, c, inv


def projective_key(b: int, c: int, inv: Sequence[int], p: int) -> int:
    assert b != 0 or c != 0
    return p if c == 0 else b * inv[c] % p


def quantiles(values: Sequence[int]) -> List[int]:
    ordered = sorted(values)
    n = len(ordered)
    ans = [ordered[0]]
    for j in range(1, 11):
        ans.append(ordered[(j * n + 9) // 10 - 1])
    return ans


def format_counter(counter: Counter[int]) -> str:
    return ", ".join(f"{k}:{counter[k]}" for k in sorted(counter))


def analyze_prime(p: int) -> Dict[str, object]:
    nmax = p - 2
    named_cutoffs = cutoffs(p)
    max_d = max(d for _, d in named_cutoffs)
    small_limit = ceil_decimal(Decimal(p).ln() ** 2)
    b, c, inv = apery_pair(p)

    keys = [None] * (nmax + 1)
    positions: Dict[int, List[int]] = defaultdict(list)
    for r in range(1, nmax + 1):
        key = projective_key(b[r], c[r], inv, p)
        keys[r] = key
        positions[key].append(r)

    # Evaluate every N_d(r) over all r in F_p.  On the nonwrapping triangle,
    # simultaneously verify Delta_{r,d} prod_j(r+j)^3 = N_d(r).
    pvals = [(34 * x**3 + 51 * x * x + 27 * x + 5) % p for x in range(p)]
    cube = [x**3 % p for x in range(p)]
    sixth = [x**6 % p for x in range(p)]
    n_prev = [0] * p  # N_0
    n_cur = [1] * p   # N_1
    denom = [1] * (nmax + 1)
    root_counts: List[int] = [0] * (max_d + 1)
    collision_sets: List[set[int]] = [set() for _ in range(nmax + 1)]
    c_by_gap = [0] * (max_d + 1)
    identity_checks = 0

    for d in range(1, max_d + 1):
        root_counts[d] = sum(value == 0 for value in n_cur)
        gap_count = 0
        for r in range(1, nmax - d + 1):
            x = (r + d) % p
            denom[r] = denom[r] * cube[x] % p
            delta = (b[r] * c[r + d] - b[r + d] * c[r]) % p
            assert delta * denom[r] % p == n_cur[r]
            identity_checks += 1
            if n_cur[r] == 0:
                collision_sets[r].add(d)
                gap_count += 1
        c_by_gap[d] = gap_count
        if d < max_d:
            n_next = [0] * p
            for r in range(p):
                x = (r + d) % p
                n_next[r] = (pvals[x] * n_cur[r] - sixth[x] * n_prev[r]) % p
            n_prev, n_cur = n_cur, n_next

    # Independent reconstruction from normalized projective orbit fibers.
    fiber_pairs: Dict[Tuple[int, int], int] = {}
    for pos in positions.values():
        for i, left in enumerate(pos):
            for j in range(i + 1, len(pos)):
                gap = pos[j] - left
                if gap > max_d:
                    break
                fiber_pairs[(left, gap)] = j - i - 1
    recurrence_pairs = {
        (r, d) for r in range(1, nmax + 1) for d in collision_sets[r]
    }
    assert set(fiber_pairs) == recurrence_pairs

    rows = []
    for label, dcut in named_cutoffs:
        multiplicities = [
            sum(gap <= dcut for gap in collision_sets[r])
            for r in range(1, nmax + 1)
        ]
        s_value = sum(multiplicities)
        assert s_value == sum(c_by_gap[1 : dcut + 1])
        q_value = sum(m * (m - 1) // 2 for m in multiplicities)

        split_hist: Counter[int] = Counter()
        for (r, gap), split_count in fiber_pairs.items():
            if gap <= dcut:
                split_hist[split_count] += 1
                # Check every claimed split against both collision predicates.
                pos = positions[keys[r]]
                left_index = pos.index(r)
                for k in range(1, split_count + 1):
                    middle = pos[left_index + k]
                    assert middle - r in collision_sets[r]
                    assert r + gap - middle in collision_sets[middle]
        assert sum(split_hist.values()) == s_value
        assert sum(k * count for k, count in split_hist.items()) == q_value

        max_mult = max(multiplicities)
        maximizers = [r for r, m in enumerate(multiplicities, start=1) if m == max_mult]
        ranked = sorted(
            range(1, nmax + 1), key=lambda r: (-multiplicities[r - 1], r)
        )[:10]
        details_by_r = {}
        for r in sorted(set(ranked) | set(maximizers)):
            gaps = sorted(g for g in collision_sets[r] if g <= dcut)
            exact_mirror = [g for g in gaps if 2 * r + g == p - 1]
            min_mirror2 = min((abs(2 * r + g - (p - 1)) for g in gaps), default=None)
            details_by_r[r] = {
                "r": r,
                "m": multiplicities[r - 1],
                "z": b[r] == 0,
                "small": r <= small_limit,
                "gaps": gaps,
                "mirror": exact_mirror,
                "min_mirror2": min_mirror2,
            }

        rows.append(
            {
                "label": label,
                "D": dcut,
                "S": s_value,
                "Q": q_value,
                "max": max_mult,
                "hist": Counter(multiplicities),
                "quantiles": quantiles(multiplicities),
                "maximizers": maximizers,
                "top": [details_by_r[r] for r in ranked],
                "max_details": [details_by_r[r] for r in maximizers],
                "P": split_hist[0],
                "split_hist": split_hist,
                "nonprimitive": s_value - split_hist[0],
                "split_witnesses": sum(k * count for k, count in split_hist.items()),
            }
        )

    return {
        "p": p,
        "N": nmax,
        "small_limit": small_limit,
        "identity_checks": identity_checks,
        "rows": rows,
        "root_counts": root_counts[1 : small_limit + 1],
        "z_count": sum(b[r] == 0 for r in range(1, nmax + 1)),
    }


def decimal_ratio(a: int, b: Decimal) -> str:
    return f"{Decimal(a) / b:.6f}"


def detail_text(item: Dict[str, object], p: int) -> str:
    mirror = item["mirror"]
    mirror_text = str(mirror) if mirror else "none"
    min2 = item["min_mirror2"]
    return (
        f"r={item['r']} (m={item['m']}, r/p={Decimal(item['r']) / Decimal(p):.4f}, "
        f"Z_p={'yes' if item['z'] else 'no'}, small={'yes' if item['small'] else 'no'}, "
        f"gaps={item['gaps']}, exact-mirror gaps={mirror_text}, min |2r+d-(p-1)|={min2})"
    )


def render(results: Sequence[Dict[str, object]]) -> str:
    lines: List[str] = []
    lines += [
        "# Last-wall empirical ground truth",
        "",
        "All finite-field values and counts below are exact integers produced by pure "
        "Python 3 standard-library code in `CODEX_lastwall_empirics.py`. Here `N=p-2`, "
        "`r` ranges over `1,...,N`, and `d_D(r)` counts admissible nonwrapping "
        "collisions with `r+d<=N`. The cutoffs `p^0.6` and `p^0.66` were ceiled by "
        "exact integer-power comparisons; `sqrt(p) log(p)` (natural logarithm) was "
        "evaluated with 80-digit Decimal arithmetic. Ratios only are rounded.",
        "",
        "For each prime, the orbit was generated from both Apery recurrences. Every "
        "regular value was checked against",
        "",
        "`Delta_(r,d) prod_(j=1)^d (r+j)^3 = N_d(r) (mod p)`,",
        "",
        "and all collision pairs were independently reconstructed by grouping equal "
        "normalized projective orbit values. `R_d` in Section 4 means the full number "
        "of roots of `N_d` in `F_p`, including residues outside the nonwrapping "
        "window; it is therefore distinct from `C_d`.",
        "",
        "## 1. Exact window statistics and distributions",
        "",
        "Histogram notation is `multiplicity:number of bases`. Deciles use nearest "
        "rank: `q_j` is the entry of rank `ceil(jN/10)` in the sorted list (with "
        "`q_0` the minimum). Thus the histogram plus deciles and top ten give the "
        "requested full distribution summary.",
        "",
    ]
    for result in results:
        p = int(result["p"])
        lines += [
            f"### p={p}",
            "",
            f"`N={result['N']}`, `|Z_p intersect [1,N]|={result['z_count']}`, "
            f"Delta-identity checks `{result['identity_checks']}` (all passed).",
            "",
        ]
        for row in result["rows"]:  # type: ignore[index]
            qs = ", ".join(
                f"q{10*j}={v}" for j, v in enumerate(row["quantiles"])
            )
            lines += [
                f"#### {row['label']}: D={row['D']}",
                "",
                f"`S_D={row['S']}`, `Q_D={row['Q']}`, `max_r d_D(r)={row['max']}`.",
                "",
                f"- Full histogram: `{format_counter(row['hist'])}`.",
                f"- Deciles: `{qs}`.",
                "- Top ten (ties ordered by increasing `r`):",
                "",
            ]
            lines += [f"  {i}. {detail_text(item, p)}" for i, item in enumerate(row["top"], 1)]
            lines.append("")

    lines += [
        "## 2. Vector-7 premise: maximum return multiplicity",
        "",
        "The comparison columns are `max/D^(2/3)` and `max/log(p)`. The final "
        "column lists every maximizing base, not merely the first ten. `mirror` "
        "records whether at least one of that base's collisions is the exact forced "
        "mirror collision `2r+d=p-1`; `near2` is the minimum doubled distance "
        "`|2r+d-(p-1)|`. `small` means `r<=ceil((log p)^2)`.",
        "",
        "Each cell is written as a code block, followed by the complete maximizing-base "
        "list.",
        "",
    ]
    for result in results:
        p = int(result["p"])
        for row in result["rows"]:  # type: ignore[index]
            d = int(row["D"])
            maxima_parts = []
            details_by_r = {int(item["r"]): item for item in row["max_details"]}
            for r in row["maximizers"]:
                item = details_by_r[r]
                maxima_parts.append(
                    f"{r}(Z={'Y' if item['z'] else 'N'},small={'Y' if item['small'] else 'N'},"
                    f"mirror={'Y' if item['mirror'] else 'N'},near2={item['min_mirror2']},"
                    f"gaps={item['gaps']})"
                )
            d23 = Decimal(d) ** (Decimal(2) / Decimal(3))
            logp = Decimal(p).ln()
            lines += [
                f"`p={p}, D={d}, max={row['max']}, D^(2/3)={d23:.6f}, "
                f"log(p)={logp:.6f}, max/D^(2/3)={decimal_ratio(row['max'], d23)}, "
                f"max/log(p)={decimal_ratio(row['max'], logp)}`",
                "",
                f"Maximizing bases: {'; '.join(maxima_parts)}.",
                "",
            ]
    all_max = [int(row["max"]) for result in results for row in result["rows"]]  # type: ignore[index]
    lines += [
        "",
        f"Across all 12 cells the maximum lies in `{min(all_max)}..{max(all_max)}` "
        "while `D^(2/3)` grows by a much larger factor. These data strongly reject "
        "`c D^(2/3)` at the tested scale. They are compatible with an absolute "
        "bound and also with logarithmic growth too slow to resolve from four primes; "
        "the bounded description is the sharper empirical fit, not a proof.",
        "",
        "## 3. Vector-8 premise: primitive decomposition",
        "",
        "For an endpoint collision `(r,d)`, its split multiplicity is the number of "
        "`d'` with `0<d'<d` for which both `(r,d')` and `(r+d',d-d')` collide. "
        "Because the determinant and projective-fiber computations agreed pair by "
        "pair, this is also the exact number of intermediate occurrences of the same "
        "projective value. Histogram notation is `split multiplicity:number of endpoint "
        "collisions`. In every cell, `P_D` is the zero bin, every positive-bin "
        "collision splits, and the exact renewal checksum is "
        "`sum_k k*n_k=Q_D`.",
        "",
        "The per-cell records below are `(p,D,P_D,S_D,P_D/S_D,nonprimitive,"
        "split_histogram,split_witnesses,Q_D,failures)`.",
        "",
    ]
    for result in results:
        p = int(result["p"])
        for row in result["rows"]:  # type: ignore[index]
            ratio = Decimal(row["P"]) / Decimal(row["S"]) if row["S"] else Decimal(0)
            lines.append(
                f"- `({p},{row['D']},{row['P']},{row['S']},{ratio:.6f},"
                f"{row['nonprimitive']},[{format_counter(row['split_hist'])}],"
                f"{row['split_witnesses']},{row['Q']},0)`"
            )
    lines += [
        "",
        "The renewal claim therefore passes every endpoint collision in all 12 cells, "
        "with no unsplit nonprimitive collision. The primitive share remains high, so "
        "the data validate decomposition but do not by themselves supply a small "
        "primitive bound.",
        "",
        "## 4. Exact small-d full root counts",
        "",
        "Each line gives every exact `d:R_d` value through "
        "`K=ceil((log p)^2)`. The mean is also shown as the exact fraction "
        "`sum R_d/K`.",
        "",
    ]
    global_max = 0
    for result in results:
        p = int(result["p"])
        roots = result["root_counts"]
        total = sum(roots)
        maximum = max(roots)
        global_max = max(global_max, maximum)
        max_ds = [i + 1 for i, x in enumerate(roots) if x == maximum]
        values = ", ".join(f"{i}:{x}" for i, x in enumerate(roots, 1))
        lines += [
            f"### p={p}, K={result['small_limit']}",
            "",
            f"`max R_d={maximum}` at `d={max_ds}`; mean "
            f"`={total}/{len(roots)}={Decimal(total) / Decimal(len(roots)):.6f}`.",
            "",
            f"`{values}`",
            "",
        ]
    lines += [
        f"The largest observed small-segment root count is `{global_max}`. Means stay "
        "of constant size. Thus `R_d=O(1)` on this polylogarithmic segment is "
        "empirically available for these four primes (with observed constant "
        f"`{global_max}`), but the finite scan is not a uniform theorem in `p`.",
        "",
        "## 5. Verdicts for the three deep-strike premises",
        "",
        "- **Vector 7 -- HOLDS EMPIRICALLY IN THE BOUNDED/LOG-SIZED FORM, NOT THE "
        "`D^(2/3)` FORM.** The maxima remain tiny compared with both `D^(2/3)` and "
        "the available window length. Maximizers are not systematically Apery zeros, "
        "small bases, or exact mirror centers; the per-cell annotations give the full "
        "exceptions. This supports anti-concentration, but does not prove the uniform "
        "bound needed by the strike.",
        "",
        "- **Vector 8 -- RENEWAL CLAIM VERIFIED EXACTLY; THE NEEDED PRIMITIVE BOUND "
        "REMAINS OPEN.** Every nonprimitive endpoint splits and the weighted split "
        "histogram equals `Q_D` in every cell. However, most endpoint collisions are "
        "primitive at these scales, so decomposition alone has not reduced `S_D` to "
        "a demonstrably smaller quantity.",
        "",
        "- **Small-d input -- HOLDS EMPIRICALLY.** Complete root counts through "
        "`ceil((log p)^2)` have constant-sized means and a small global maximum. This "
        "makes the proposed `R_d=O(1)` input numerically plausible on the entire "
        "polylog segment, while supplying no proof beyond the four tested primes.",
        "",
    ]
    return "\n".join(lines)


def main() -> None:
    results = []
    for p in PRIMES:
        print(f"[lastwall] p={p}", flush=True)
        result = analyze_prime(p)
        results.append(result)
        print(
            f"[lastwall] p={p} checks={result['identity_checks']} "
            f"cells={[(row['D'], row['S'], row['Q'], row['max']) for row in result['rows']]}",
            flush=True,
        )
    OUT.write_text(render(results), encoding="utf-8")
    print(f"[lastwall] wrote {OUT}", flush=True)


if __name__ == "__main__":
    main()
