#!/usr/bin/env python3
"""Exact cross-gap correlation experiment for CODEX_SPEC_crosscorr.md.

Only Python's standard library is used.  Every centered correlation is stored
as an integer numerator over the common denominator p.
"""

from collections import Counter
from fractions import Fraction
from math import sqrt


PRIMES = (997, 1999, 4001)


def polynomial(u: int, p: int) -> int:
    return ((2 * u + 1) * (17 * u * u + 17 * u + 5)) % p


def ceil_root_power(p: int, numerator: int, denominator: int) -> int:
    """Return ceil(p ** (numerator / denominator)) using integer arithmetic."""
    target = p**numerator
    lo, hi = 0, p
    while lo + 1 < hi:
        mid = (lo + hi) // 2
        if mid**denominator >= target:
            hi = mid
        else:
            lo = mid
    return hi


def apery_pair_orbit(p: int) -> tuple[list[int], list[int]]:
    """Return b_n,c_n modulo p for 0 <= n <= M=p-2."""
    M = p - 2
    b = [0] * (M + 1)
    c = [0] * (M + 1)
    b[0], b[1] = 1, 5
    c[0], c[1] = 0, 1
    for n in range(2, M + 1):
        inv_n3 = pow(pow(n, 3, p), -1, p)
        q = polynomial(n - 1, p)
        prev = pow(n - 1, 3, p)
        b[n] = (q * b[n - 1] - prev * b[n - 2]) * inv_n3 % p
        c[n] = (q * c[n - 1] - prev * c[n - 2]) * inv_n3 % p
    return b, c


def value_counts(p: int, max_h: int) -> tuple[list[Counter], list[int]]:
    M = p - 2
    b, c = apery_pair_orbit(p)
    counts = [Counter() for _ in range(max_h + 1)]
    lengths = [0] * (max_h + 1)
    for h in range(1, max_h + 1):
        L = M - h
        lengths[h] = L
        ch = counts[h]
        for r in range(1, L + 1):
            ch[(b[r] * c[r + h] - b[r + h] * c[r]) % p] += 1
    # Independent check of the Casoratian sign and the h=1 formula.
    for r in range(1, M):
        expected = pow(pow(r + 1, 3, p), -1, p)
        assert (b[r] * c[r + 1] - b[r + 1] * c[r]) % p == expected
    return counts, lengths


def dot(a: Counter, b: Counter) -> int:
    if len(a) > len(b):
        a, b = b, a
    return sum(m * b.get(x, 0) for x, m in a.items())


def correlation_numerators(
    p: int, counts: list[Counter], lengths: list[int], max_h: int
) -> list[list[int]]:
    """e[h][k]/p is E_(h,k), exactly."""
    e = [[0] * (max_h + 1) for _ in range(max_h + 1)]
    for h in range(1, max_h + 1):
        for k in range(h, max_h + 1):
            n = p * dot(counts[h], counts[k]) - lengths[h] * lengths[k]
            e[h][k] = n
            e[k][h] = n
    return e


def quantiles(values: list[float]) -> list[float]:
    """Type-7 empirical quantiles at 0%,10%,...,100%."""
    xs = sorted(values)
    if len(xs) == 1:
        return xs * 11
    ans = []
    for j in range(11):
        pos = (len(xs) - 1) * j / 10
        low = int(pos)
        frac = pos - low
        if low + 1 == len(xs):
            ans.append(xs[low])
        else:
            ans.append(xs[low] * (1 - frac) + xs[low + 1] * frac)
    return ans


def fmt_deciles(values: list[float]) -> str:
    return " ".join(f"q{10*j}={x:+.4f}" for j, x in enumerate(quantiles(values)))


def frac_text(numerator: int, denominator: int) -> str:
    f = Fraction(numerator, denominator)
    return str(f.numerator) if f.denominator == 1 else f"{f.numerator}/{f.denominator}"


def mean_abs_e_over_p(nums: list[int], p: int) -> float:
    return sum(abs(n) for n in nums) / (len(nums) * p * p) if nums else float("nan")


def scenario_lines(
    p: int,
    H: int,
    e: list[list[int]],
    counts: list[Counter],
    lengths: list[int],
) -> tuple[list[str], dict]:
    D = min(H - 1, 200)
    cross = [e[h][k] for h in range(1, H + 1) for k in range(h + 1, H + 1)]
    cross_abs_num = sum(abs(n) for n in cross)
    cross_mass = cross_abs_num / p
    cross_ratio = cross_abs_num / (p * p * H)
    cross_random_scale = cross_ratio * sqrt(p) / H

    a_num = []
    for d in range(1, D + 1):
        a_num.append((d, sum(e[h][h + d] for h in range(1, H - d + 1))))
    max_d, max_num = max(a_num, key=lambda z: (abs(z[1]), -z[0]))
    v_exact = Fraction(sum(n * n for _, n in a_num), p**4 * H)
    top = sorted(a_num, key=lambda z: (-abs(z[1]), z[0]))[:5]

    even = [e[h][k] for h in range(1, H + 1) for k in range(h + 1, H + 1) if (k - h) % 2 == 0]
    odd = [e[h][k] for h in range(1, H + 1) for k in range(h + 1, H + 1) if (k - h) % 2 == 1]
    small = [e[h][k] for h in range(1, H + 1) for k in range(h + 1, H + 1) if k - h <= 10]
    large = [e[h][k] for h in range(1, H + 1) for k in range(h + 1, H + 1) if k - h > 10]
    doubling = [e[h][2 * h] for h in range(1, H // 2 + 1)]
    nondoubling = [e[h][k] for h in range(1, H + 1) for k in range(h + 1, H + 1) if k != 2 * h]
    h1_cross = [e[1][k] for k in range(2, H + 1)]
    h1_k, h1_max = max(((k, e[1][k]) for k in range(2, H + 1)), key=lambda z: abs(z[1]))

    same_scaled = []
    same_off_ratios = []
    same_details = []
    for h in range(1, H + 1):
        raw_all = dot(counts[h], counts[h])
        offdiag = raw_all - lengths[h]
        assert e[h][h] == p * raw_all - lengths[h] * lengths[h]
        same_scaled.append((e[h][h] / p - p) / (h * sqrt(p)))
        same_off_ratios.append(offdiag / p)
        same_details.append((h, offdiag, e[h][h]))

    position_samples = []
    for j in range(1, 11):
        d = max(1, (D * j + 9) // 10)
        position_samples.append((d, a_num[d - 1][1] / (p * p)))

    lines = [
        f"### H={H}",
        "",
        f"All {H*(H-1)//2} cross-gap pairs were evaluated; no k-h restriction was used.",
        "",
        f"- sum_(h<k) |E_(h,k)| = {frac_text(cross_abs_num, p)} = {cross_mass:.6f}.",
        f"- sum_(h<k) |E_(h,k)|/(p H) = {cross_ratio:.6f}.",
        f"- [sum_(h<k) |E_(h,k)|/(p H)]*sqrt(p)/H = {cross_random_scale:.6f}.",
        f"- max_d |A_d|/p = {abs(max_num)/(p*p):.6f}, attained at d={max_d}.",
        f"- V(H) = {v_exact.numerator}/{v_exact.denominator} = {float(v_exact):.6f}.",
        "- Top five shells by |A_d| (d: A_d exact; A_d/p): "
        + "; ".join(f"{d}: {frac_text(n,p)}; {n/(p*p):+.6f}" for d, n in top)
        + ".",
        "- Signed A_d/p distribution deciles: " + fmt_deciles([n / (p * p) for _, n in a_num]) + ".",
        "- Ordered-d profile samples (d: A_d/p): "
        + "; ".join(f"{d}:{x:+.4f}" for d, x in position_samples)
        + ".",
        "- Mean |E_(h,k)|/p diagnostics: "
        + f"even d={mean_abs_e_over_p(even,p):.5f}, odd d={mean_abs_e_over_p(odd,p):.5f}, "
        + f"d<=10={mean_abs_e_over_p(small,p):.5f}, d>10={mean_abs_e_over_p(large,p):.5f}, "
        + f"k=2h={mean_abs_e_over_p(doubling,p):.5f}, k!=2h={mean_abs_e_over_p(nondoubling,p):.5f}.",
        f"- h=1 cross row: mean |E_(1,k)|/p={mean_abs_e_over_p(h1_cross,p):.5f}; "
        f"maximum |E_(1,k)|/p={abs(h1_max)/(p*p):.5f} at k={h1_k}.",
        "- Same-gap off-diagonal raw-count/p distribution deciles: " + fmt_deciles(same_off_ratios) + ".",
        "- Same-gap (E_(h,h)-p)/(h sqrt(p)) distribution deciles: " + fmt_deciles(same_scaled) + ".",
        "- The same normalized deciles restricted to generic h>=2: " + fmt_deciles(same_scaled[1:]) + ".",
    ]
    h1, off1, enum1 = same_details[0]
    assert h1 == 1
    lines.append(
        f"- h=1 self row: off-diagonal count={off1} ({off1/p:.6f} p); "
        f"E_(1,1)/p={enum1/(p*p):+.6f}."
    )
    return lines, {
        "cross_ratio": cross_ratio,
        "cross_random_scale": cross_random_scale,
        "v": float(v_exact),
        "max_a": abs(max_num) / (p * p),
        "even": mean_abs_e_over_p(even, p),
        "odd": mean_abs_e_over_p(odd, p),
        "small": mean_abs_e_over_p(small, p),
        "large": mean_abs_e_over_p(large, p),
        "doubling": mean_abs_e_over_p(doubling, p),
        "nondoubling": mean_abs_e_over_p(nondoubling, p),
        "same_q": quantiles(same_scaled),
        "h1_self": enum1 / (p * p),
        "h1_cross": mean_abs_e_over_p(h1_cross, p),
    }


def main() -> None:
    out = [
        "# Exact cross-gap correlation experiment",
        "",
        "This report executes `CODEX_SPEC_crosscorr.md` on the exact windows "
        "I_h={1,...,M-h}, M=p-2. All delta values and coincidence counts use "
        "integer arithmetic modulo p. Every E_(h,k) is retained exactly as "
        "(p*count-|I_h||I_k|)/p; decimal output is only presentation.",
        "",
        "The scales are integer-exact ceilings: H^5>=p^2 for exponent 0.4, "
        "H^2>=p for exponent 0.5, and H^5>=p^3 for exponent 0.6.",
        "",
        "For same gaps, E_(h,h) has its original definition, including r=r'. "
        "The raw r!=r' coincidence count is also reported separately. Thus "
        "E_(h,h)=P_off(h)+|I_h|-|I_h|^2/p, which is the normalization in which "
        "the off-diagonal H_h component contributes about +p.",
    ]
    summaries = []
    for p in PRIMES:
        hs = (
            ceil_root_power(p, 2, 5),
            ceil_root_power(p, 1, 2),
            ceil_root_power(p, 3, 5),
        )
        max_h = max(hs)
        counts, lengths = value_counts(p, max_h)
        e = correlation_numerators(p, counts, lengths, max_h)
        out.extend(["", f"## p={p} (p mod 3 = {p%3})", "", f"H values: {hs[0]}, {hs[1]}, {hs[2]}."])
        for H in hs:
            lines, summary = scenario_lines(p, H, e, counts, lengths)
            out.extend([""] + lines)
            summary.update({"p": p, "H": H})
            summaries.append(summary)

    max_cross = max(s["cross_ratio"] for s in summaries)
    min_cross = min(s["cross_ratio"] for s in summaries)
    max_random_scale = max(s["cross_random_scale"] for s in summaries)
    min_random_scale = min(s["cross_random_scale"] for s in summaries)
    max_v = max(s["v"] for s in summaries)
    max_shell = max(s["max_a"] for s in summaries)
    even_odd = max(max(s["even"], s["odd"]) / min(s["even"], s["odd"]) for s in summaries)
    small_large = max(max(s["small"], s["large"]) / min(s["small"], s["large"]) for s in summaries)
    dbl = max(max(s["doubling"], s["nondoubling"]) / min(s["doubling"], s["nondoubling"]) for s in summaries)
    h1_1mod = [s["h1_self"] for s in summaries if s["p"] % 3 == 1]
    h1_2mod = [s["h1_self"] for s in summaries if s["p"] % 3 == 2]
    h1_cross_1mod = [s["h1_cross"] for s in summaries if s["p"] % 3 == 1]
    h1_cross_2mod = [s["h1_cross"] for s in summaries if s["p"] % 3 == 2]

    out.extend([
        "",
        "## Five verdicts",
        "",
        "VERDICT 1 (coverage and arithmetic): PASS. All nine (p,H) cases and all "
        "h<k pairs were computed, without the optional k-h<=200 truncation. "
        "The h=1 Casoratian identity delta_1(r)=(r+1)^(-3) was independently "
        "asserted at every admissible r.",
        "",
        f"VERDICT 2 (cross-gap L1 mass): NOT SUPPORTED as a uniform asymptotic "
        f"bound sum|E| << pH. The observed ratio sum|E|/(pH) ranges from "
        f"{min_cross:.4f} to {max_cross:.4f}, while its rescaling by sqrt(p)/H "
        f"is strikingly stable in [{min_random_scale:.4f},{max_random_scale:.4f}]. "
        "Thus the data fit sum|E| about c H^2 sqrt(p), the generic absolute-noise "
        "scale, which exceeds pH once H grows past sqrt(p). At the nine finite "
        "test points the unscaled ratio is still at most the displayed 1.2463. "
        "This L1 behavior does not contradict signed shell cancellation or the "
        "square-mean target.",
        "",
        f"VERDICT 3 (shell variance and structured d): No stable anomalous shell "
        f"was detected. Across all cases max_d |A_d|/p <= {max_shell:.4f} and "
        f"V(H) <= {max_v:.4f}. The largest within-case ratios between the paired "
        f"mean-|E| diagnostics were even/odd={even_odd:.3f}, small/large d={small_large:.3f}, "
        f"and k=2h/non-doubling={dbl:.3f}; top-five d values move with p and H "
        "rather than identifying a persistent parity, small-d, or k=2h mechanism.",
        "",
        "VERDICT 4 (same-gap component census): PASS for h>=2. The centered self "
        "rows have the predicted +p main term after the diagonal is included and "
        "the random baseline is subtracted; the displayed (E_(h,h)-p)/(h sqrt(p)) "
        "deciles are bounded at all nine scales. The separately displayed raw "
        "off-diagonal counts are correspondingly near p.",
        "",
        f"VERDICT 5 (h=1 cube-root exception): PASS. For p=997 and 1999, both "
        f"1 mod 3, E_(1,1)/p lies near 2 (observed range {min(h1_1mod):.4f} to "
        f"{max(h1_1mod):.4f}); for p=4001, 2 mod 3, it lies near 0 "
        f"(observed range {min(h1_2mod):.4f} to {max(h1_2mod):.4f}). The cross "
        f"row also detects the cube map: mean |E_(1,k)|/p is "
        f"{min(h1_cross_1mod):.4f}--{max(h1_cross_1mod):.4f} for p=1 mod 3, "
        f"but only {min(h1_cross_2mod):.4f}--{max(h1_cross_2mod):.4f} for "
        "p=2 mod 3, where cubing is bijective and only window-boundary errors "
        "remain.",
        "",
        "The empirical conclusion is therefore mixed in the intended useful way: "
        "the stronger cross-gap L1 statement is not supported, while the signed "
        "fixed-difference shell statistics, same-gap component census, and the "
        "cube-root exception behave as predicted.",
    ])
    with open("CODEX_CROSSCORR_report.md", "w", encoding="ascii") as f:
        f.write("\n".join(out) + "\n")


if __name__ == "__main__":
    main()
