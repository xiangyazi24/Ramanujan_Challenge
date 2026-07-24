#!/usr/bin/env python3
"""Round 2, V1/V2: audit the proposed trace identity and reversal.

This script deliberately uses only the definitions recorded in
problems/3.2/proof.tex:

    b_n = sum_k binom(n,k)^2 binom(n+k,k)^2,
    H_p(t) = sum_{j=0}^{p-1} b_j t^j  in F_p[t].

For each requested prime it computes every coefficient twice (by the Apéry
recurrence and by the defining binomial sum), evaluates H_p at the two
conifold parameters when they lie in F_p, and checks every reversal ratio.

The default report is /tmp/round2_v1_v2.txt.
"""

from __future__ import annotations

import argparse
from math import comb, isqrt
from pathlib import Path
from typing import Iterable, List, Optional, Sequence, Tuple


PRIMES: Tuple[int, ...] = (
    13,
    17,
    19,
    23,
    29,
    31,
    37,
    41,
    43,
    47,
    53,
    59,
    61,
    67,
    71,
    73,
    79,
    89,
    97,
    101,
)


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    for q in range(2, isqrt(n) + 1):
        if n % q == 0:
            return False
    return True


def legendre(a: int, p: int) -> int:
    """Return the Legendre symbol (a/p), for odd prime p."""
    value = pow(a % p, (p - 1) // 2, p)
    if value == p - 1:
        return -1
    return value


def apery_binomial_mod_p(n: int, p: int) -> int:
    return sum(comb(n, k) ** 2 * comb(n + k, k) ** 2 for k in range(n + 1)) % p


def apery_recurrence_mod_p(p: int) -> List[int]:
    """Return b_0,...,b_{p-1} modulo p using the Apéry recurrence."""
    b = [0] * p
    b[0] = 1
    b[1] = 5 % p
    for n in range(1, p - 1):
        n2 = n * n
        n3 = n2 * n
        middle = (34 * n3 + 51 * n2 + 27 * n + 5) % p
        denominator = pow(n + 1, 3, p)
        numerator = (middle * b[n] - (n3 % p) * b[n - 1]) % p
        b[n + 1] = numerator * pow(denominator, -1, p) % p
    return b


def evaluate_polynomial(coefficients: Sequence[int], t: int, p: int) -> int:
    value = 0
    for coefficient in reversed(coefficients):
        value = (value * t + coefficient) % p
    return value


def least_square_root(a: int, p: int) -> Optional[int]:
    roots = [x for x in range(p) if x * x % p == a % p]
    return min(roots) if roots else None


def conifold_parameters(p: int) -> Optional[Tuple[int, int, int]]:
    """Return (chosen sqrt(2), t1, t2), or None if 2 is nonsquare.

    There is no canonical square root in F_p.  We choose the least integer
    representative r and label t1=17-12r, t2=17+12r.  Replacing r by -r
    merely swaps t1 and t2, leaving every tested sum unchanged.
    """
    root = least_square_root(2, p)
    if root is None:
        return None
    t1 = (17 - 12 * root) % p
    t2 = (17 + 12 * root) % p
    assert (t1 * t2) % p == 1
    assert (t1 * t1 - 34 * t1 + 1) % p == 0
    assert (t2 * t2 - 34 * t2 + 1) % p == 0
    return root, t1, t2


def chunks(items: Sequence[str], size: int = 12) -> Iterable[str]:
    for start in range(0, len(items), size):
        yield " ".join(items[start : start + size])


def format_vector(name: str, values: Sequence[int], width: int = 12) -> List[str]:
    entries = [f"{j}:{value}" for j, value in enumerate(values)]
    result = [f"{name} (j:value):"]
    result.extend(f"  {line}" for line in chunks(entries, width))
    return result


def analyze_prime(p: int) -> dict:
    assert is_prime(p)
    b = apery_recurrence_mod_p(p)
    b_direct = [apery_binomial_mod_p(j, p) for j in range(p)]
    assert b == b_direct, f"recurrence/binomial disagreement at p={p}"

    palindromic_failures = [j for j in range(p) if b[p - 1 - j] != b[j]]
    ratios = [b[p - 1 - j] * pow(b[j], -1, p) % p for j in range(p) if b[j]]
    ratio_values = sorted(set(ratios))

    parameters = conifold_parameters(p)
    correction = [0] * p
    c1: Optional[int] = None
    c2: Optional[int] = None
    sqrt2: Optional[int] = None
    t1: Optional[int] = None
    t2: Optional[int] = None
    if parameters is not None:
        sqrt2, t1, t2 = parameters
        c1 = evaluate_polynomial(b, t1, p)
        c2 = evaluate_polynomial(b, t2, p)
        # Reduction of omega(t)^(-j) modulo p is t^(-j).
        correction = [
            (c1 * pow(t1, -j, p) + c2 * pow(t2, -j, p)) % p
            for j in range(p)
        ]

    residual = [(b[j] + correction[j]) % p for j in range(p)]
    # The displayed star identity has b_j = -(alpha(j) + correction(j)),
    # hence alpha(j) = -b_j-correction(j) = -residual(j).
    inferred_alpha = [(-value) % p for value in residual]
    zero_indices = [j for j, value in enumerate(b) if value == 0]
    residual_at_b_zeros = [(j, residual[j]) for j in zero_indices]
    zero_implication = all(value == 0 for _, value in residual_at_b_zeros)

    defect = "Q=t^2-34t+1" if p % 24 in {13, 17, 19, 23} else "1"
    if parameters is not None and defect.startswith("Q"):
        assert c1 == 0 and c2 == 0
    if parameters is not None:
        # Palindromy and t2=t1^-1 force equal evaluations.
        assert c1 == c2

    return {
        "p": p,
        "b": b,
        "binomial_verified": b == b_direct,
        "legendre2": legendre(2, p),
        "sqrt2": sqrt2,
        "t1": t1,
        "t2": t2,
        "c1": c1,
        "c2": c2,
        "defect": defect,
        "correction": correction,
        "residual": residual,
        "inferred_alpha": inferred_alpha,
        "zero_indices": zero_indices,
        "residual_at_b_zeros": residual_at_b_zeros,
        "zero_implication": zero_implication,
        "palindromic_failures": palindromic_failures,
        "ratio_values": ratio_values,
        "nonzero_ratio_count": len(ratios),
    }


def render_report(rows: Sequence[dict]) -> str:
    split_rows = [row for row in rows if row["legendre2"] == 1]
    nonsplit_rows = [row for row in rows if row["legendre2"] == -1]
    nonzero_boundary_rows = [row for row in split_rows if row["c1"] or row["c2"]]
    zero_rows = [row for row in rows if row["zero_indices"]]
    zero_test_failures = [
        (row["p"], row["residual_at_b_zeros"])
        for row in rows
        if not row["zero_implication"]
    ]
    reversal_failures = [
        (row["p"], row["palindromic_failures"], row["ratio_values"])
        for row in rows
        if row["palindromic_failures"] or row["ratio_values"] != [1]
    ]

    lines = [
        "Round 2: V1 trace-identity audit and V2 reversal test",
        "=" * 72,
        "",
        "Verified definitions (problems/3.2/proof.tex)",
        "------------------------------------------------",
        "b_n = sum_{k=0}^n binom(n,k)^2 binom(n+k,k)^2.",
        "H_p(t) = sum_{j=0}^{p-1} b_j t^j in F_p[t].",
        "Every coefficient below was computed by the recurrence and independently",
        "checked against the defining binomial sum.",
        "",
        "V1 conventions and logical audit",
        "--------------------------------",
        "For (2/p)=1, choose the least integer r with r^2=2 mod p and put",
        "t1=17-12r, t2=17+12r.  The other choice swaps t1,t2.",
        "c1=H_p(t1), c2=H_p(t2), exactly as stated in V1.",
        "Modulo p, omega(t)^(-j) reduces to t^(-j).  Thus correction(j) is",
        "c1*t1^(-j)+c2*t2^(-j), and is zero when (2/p)=-1 (epsilon_p=0).",
        "residual(j)=b_j+correction(j).",
        "The finite Mellin inversion b_j=-sum_t H_p(t)t^(-j) is valid as",
        "stated only for 1<=j<=p-2.  At j=0 or p-1, the k=0 and k=p-1",
        "coefficients have the same character.  Endpoint table entries are kept",
        "because V1 requested all b_j and V2 uses them, but they are not treated",
        "as an interior Mellin decomposition test.",
        "The displayed star identity implies alpha(j)=-residual(j), not",
        "alpha(j)=residual(j).  Accordingly inferred_alpha below is -residual.",
        "",
        "Crucial identifiability issue: no independent definition or computation of",
        "the pure object M_j or its unit root alpha(j) is supplied.  Therefore star",
        "can always be made true by defining alpha(j):=-residual(j).  The phrase",
        "'depends smoothly on j' has no stated finite-field criterion.  V1 is thus",
        "not verifiable as a trace decomposition from these data alone.  Only the",
        "separately requested implication b_j=0 => residual(j)=0 is falsifiable here.",
        "There is also a terminology obstruction: a literal Frobenius unit root is",
        "a p-adic unit and therefore has nonzero reduction.  Yet inferred_alpha=0",
        "at zero positions with no boundary term (for example p=19, j=8,10).",
        "Thus alpha cannot simultaneously mean a genuine unit root at every j and",
        "extend by zero on the nonordinary locus.  On that locus it is either",
        "undefined or is a Hasse-type scalar rather than a unit root.",
        "",
        "V1 summary",
        "----------",
        f"Split primes (2/p)=1 ({len(split_rows)}): "
        + ", ".join(str(row["p"]) for row in split_rows),
        f"Nonsplit primes (2/p)=-1 ({len(nonsplit_rows)}): "
        + ", ".join(str(row["p"]) for row in nonsplit_rows),
        "Split primes with nonzero c1=c2: "
        + (", ".join(str(row["p"]) for row in nonzero_boundary_rows) or "none"),
        "Primes having at least one zero b_j: "
        + (", ".join(str(row["p"]) for row in zero_rows) or "none"),
        "Failures of b_j=0 => residual(j)=0: "
        + (str(zero_test_failures) if zero_test_failures else "none"),
        "Hence the requested zero-position residual test: "
        + ("FAIL" if zero_test_failures else "PASS (including vacuous cases)"),
        "",
        "Structural simplification: t1*t2=1.  Palindromy of H_p gives",
        "H_p(t1)=H_p(t2), hence c1=c2 for every split prime.  In every split",
        "case with defect S_p=t^2-34t+1, both c-values vanish automatically.",
        "",
        "V2 summary",
        "----------",
        "For every j with b_j nonzero, ratio(j)=b_{p-1-j}/b_j was computed.",
        "Failures: " + (str(reversal_failures) if reversal_failures else "none"),
        "All 20 primes have ratio set {1}; the sign/unit is +1 exactly.",
        "Thus t^(p-1) H_p(1/t)=H_p(t), not merely an unspecified +/- form,",
        "for every tested prime.",
        "In fact the computation reflects an exact binomial proof.  With",
        "n=p-1-j, for 0<=k<=n one has modulo p",
        "  binom(n,k)=(-1)^k binom(j+k,k),",
        "  binom(n+k,k)=(-1)^k binom(j,k).",
        "The second factor is zero when k>j, so summing the squared products",
        "gives b_(p-1-j)=b_j for every odd prime, not just the sample.",
        "",
        "Compact per-prime audit",
        "-----------------------",
        "Columns: p (2/p) defect sqrt2 t1 t2 c1 c2 zeros residual_at_zeros",
        "         zero_implication reversal_ratio_set nonzero_ratios",
    ]

    for row in rows:
        split_data = (
            f"sqrt2={row['sqrt2']} t1={row['t1']} t2={row['t2']} "
            f"c1={row['c1']} c2={row['c2']}"
            if row["sqrt2"] is not None
            else "sqrt2=- t1=- t2=- c1=- c2=-"
        )
        lines.append(
            f"p={row['p']:3d} (2/p)={row['legendre2']:+d} defect={row['defect']:<16} "
            f"{split_data} zeros={row['zero_indices']} "
            f"res@zeros={row['residual_at_b_zeros']} "
            f"zero_test={'PASS' if row['zero_implication'] else 'FAIL'} "
            f"ratios={row['ratio_values']} nratios={row['nonzero_ratio_count']}"
        )

    lines.extend(["", "Complete per-prime tables", "-------------------------"])
    for row in rows:
        lines.extend(
            [
                "",
                f"p={row['p']}",
                f"(2/p)={row['legendre2']:+d}; defect={row['defect']}; "
                f"sqrt2={row['sqrt2']}; t1={row['t1']}; t2={row['t2']}; "
                f"c1={row['c1']}; c2={row['c2']}",
                f"binomial cross-check={row['binomial_verified']}; "
                f"palindromic failures={row['palindromic_failures']}; "
                f"reversal ratio set={row['ratio_values']}",
                f"b-zero indices={row['zero_indices']}; "
                f"residual at b-zeros={row['residual_at_b_zeros']}; "
                f"zero implication={row['zero_implication']}",
                *format_vector("b", row["b"]),
                *format_vector("correction", row["correction"]),
                *format_vector("residual", row["residual"]),
                *format_vector("inferred_alpha=-residual", row["inferred_alpha"]),
            ]
        )

    lines.extend(
        [
            "",
            "Final verdict",
            "-------------",
            "V2 is confirmed exactly on all requested data: ratio=1 at every",
            "defined index.  The independently stated V1 zero-residual test is",
            "summarized above.  The star identity itself remains tautological until",
            "alpha(j), or the pure part M_j whose unit root it is claimed to be, is",
            "constructed independently with a normalization fixing the sign.",
        ]
    )
    return "\n".join(lines) + "\n"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("/tmp/round2_v1_v2.txt"),
        help="report path (default: /tmp/round2_v1_v2.txt)",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    rows = [analyze_prime(p) for p in PRIMES]
    report = render_report(rows)
    args.output.write_text(report, encoding="utf-8")
    print(report, end="")


if __name__ == "__main__":
    main()
