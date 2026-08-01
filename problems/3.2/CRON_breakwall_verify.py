#!/usr/bin/env python3
"""Exact checks for CODEX_BREAKWALL_report.md.

The symbolic gates certify the new characteristic-zero identities used in
[NO-RIGHT-2-3] and the sharp degree-four comparison cover.  The finite-field
gate extracts (without importing the module-level experiment) the exact
``orbit(p)`` function from CRON_b1_crosscorr.py and checks [ZERO-TAIL-2] at
the three primes required by CODEX_SPEC_CRON_breakwall.md.
"""

from __future__ import annotations

import ast
from collections import Counter
from fractions import Fraction
from math import ceil, isqrt, log
from pathlib import Path

from sympy import Poly, QQ, cancel, diff, expand, fraction, gcd, prod, symbols


HERE = Path(__file__).resolve().parent
X, Z = symbols("X Z")


def ceil_sqrt(n: int) -> int:
    """Return the least integer t with t^2 >= n."""
    t = isqrt(n)
    return t if t * t == n else t + 1


def symbolic_gates() -> None:
    """Verify every polynomial identity used by the report."""
    P = lambda t: 34 * t**3 + 51 * t**2 + 27 * t + 5
    N0 = Poly(0, X, domain=QQ)
    N1 = Poly(1, X, domain=QQ)
    N2 = Poly(P(X + 1), X, domain=QQ)
    N3 = Poly(P(X + 2), X, domain=QQ) * N2 - Poly(
        (X + 2) ** 6, X, domain=QQ
    ) * N1
    expected_N3 = Poly(
        1155 * X**6
        + 13860 * X**5
        + 68535 * X**4
        + 178680 * X**3
        + 259059 * X**2
        + 198156 * X
        + 62531,
        X,
        domain=QQ,
    )
    assert N0.is_zero and N3 == expected_N3

    q3 = Poly(prod(X + j for j in range(1, 4)), X, domain=QQ)
    q3_sq = q3 * q3
    # If N3=A*q3^2+B*q3+C, its X^6 coefficient forces A=1155;
    # then its X^4 coefficient would be 1155*58, which is false.
    assert q3_sq.nth(4) == 58
    assert N3.nth(6) == 1155
    assert N3.nth(4) == 68535
    assert N3.nth(4) != N3.nth(6) * q3_sq.nth(4)

    # The degree-four comparison cover.  Its two pole blocks form a
    # degree-two Prouhet partition of {1,...,8}.
    block_1 = (1, 4, 6, 7)
    block_2 = (2, 3, 5, 8)
    P1 = expand(prod(X + j for j in block_1))
    P2 = expand(prod(X + j for j in block_2))
    A = expand(P1 - P2)
    assert A == -16 * X - 72
    assert sum(block_1) == sum(block_2) == 18
    assert sum(j * j for j in block_1) == sum(j * j for j in block_2) == 102

    u = cancel(A / P1)
    phi = cancel(-Z / (1 - Z))
    iota = -9 - X
    assert cancel(P1.subs(X, iota) - P2) == 0
    assert cancel(P2.subs(X, iota) - P1) == 0
    assert cancel(u.subs(X, iota) - phi.subs(Z, u)) == 0

    g = cancel(
        Z
        * (Z - 4)
        * (Z - 3)
        * (Z - 2)
        * (2 * Z - 3)
        * (3 * Z - 4)
        / (16 * (Z - 1) ** 3)
    )
    assert cancel(g.subs(Z, phi) + g) == 0
    comparison = cancel(g.subs(Z, u))
    numerator, denominator = map(expand, fraction(comparison))
    q8 = expand(prod(X + j for j in range(1, 9)))
    assert Poly(denominator, X, domain=QQ).monic() == Poly(
        q8**3, X, domain=QQ
    ).monic()
    assert Poly(numerator, X, domain=QQ).degree() == 21
    assert Poly(denominator, X, domain=QQ).degree() == 24
    assert gcd(Poly(numerator, X), Poly(denominator, X)).degree() == 0
    assert gcd(Poly(numerator, X), Poly(diff(numerator, X), X)).degree() == 0
    assert cancel(comparison.subs(X, iota) + comparison) == 0

    # Moment calculation in the h == 0 (mod 3) case of [NO-RIGHT-3].
    h = symbols("h")
    sigma = 3 * (h + 1) / 2
    tau = (h + 1) * (2 * h + 1) / 2
    endpoint_discriminant = expand(2 * (tau - 1) - (sigma - 1) ** 2)
    assert cancel(endpoint_discriminant + (h - 1) * (h - 5) / 4) == 0

    print("SYMBOLIC_GATES PASS")
    print(f"N3_X4_MISMATCH {N3.nth(4)} != {N3.nth(6)}*{q3_sq.nth(4)}")
    print(f"DEGREE4_PROUHET P1-P2={A}")
    print(
        "DEGREE4_COMPARISON "
        f"deg_num={Poly(numerator, X).degree()} "
        f"deg_den={Poly(denominator, X).degree()} "
        "squarefree_num=yes reflection_odd=yes"
    )


def load_orbit_function():
    """Compile only orbit(p) from CRON_b1_crosscorr.py.

    That source has a campaign experiment at module scope, so importing it
    would run unrelated work.  AST extraction ensures this verifier uses the
    stipulated orbit code verbatim while avoiding those side effects.
    """
    source_path = HERE / "CRON_b1_crosscorr.py"
    source = source_path.read_text(encoding="utf-8")
    tree = ast.parse(source, filename=str(source_path))
    functions = [
        node
        for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name == "orbit"
    ]
    assert len(functions) == 1
    module = ast.Module(body=functions, type_ignores=[])
    ast.fix_missing_locations(module)
    namespace: dict[str, object] = {}
    exec(compile(module, str(source_path), "exec"), namespace)
    return namespace["orbit"]


def root_counts(points: list[tuple[int, int]], D: int) -> list[int]:
    N = len(points)
    return [
        sum(points[r] == points[r + h] for r in range(N - h))
        for h in range(1, D + 1)
    ]


def zero_tail_data(counts: list[int], N: int) -> dict[str, object]:
    """Return and exactly check the layer-cake [ZERO-TAIL-2] inequality."""
    D = len(counts)
    T0 = 1
    while T0 * T0 * D < N:
        T0 += 1
    maximum = max(counts, default=0)
    tails = {t: sum(value >= t for value in counts) for t in range(1, maximum + 1)}
    ZD = max((t * t * tails[t] for t in range(T0 + 1, maximum + 1)), default=0)
    S = sum(counts)
    layer_cake = sum(tails.values())
    exact_bound = D * T0 + ZD * sum(
        (Fraction(1, t * t) for t in range(T0 + 1, maximum + 1)),
        Fraction(0),
    )
    coarse_bound = Fraction(D * T0, 1) + Fraction(ZD, T0)
    assert S == layer_cake
    assert Fraction(S, 1) <= exact_bound <= coarse_bound
    return {
        "T0": T0,
        "S": S,
        "max_R": maximum,
        "ZD": ZD,
        "exact_bound": exact_bound,
        "coarse_bound": coarse_bound,
        "histogram": dict(sorted(Counter(counts).items())),
    }


def marginal_extremizer(N: int, D: int) -> dict[str, object]:
    """Build the sharp abstract per-row L2 countermodel from the report."""
    counts = []
    collision_sums = []
    for h in range(1, D + 1):
        mass = N - h
        zero_mass = min(3 * (h - 1), isqrt(mass))
        # Histogram: zero_mass at 0, and mass-zero_mass singleton cells.
        # There are p-1=N nonzero cells, so this is always realizable.
        collision_sum = zero_mass * zero_mass + mass - zero_mass
        assert 0 <= zero_mass <= 3 * (h - 1)
        assert mass - zero_mass <= N
        assert collision_sum <= 2 * mass
        counts.append(zero_mass)
        collision_sums.append(collision_sum)
    tail = zero_tail_data(counts, N)
    return {
        "S": sum(counts),
        "max_R": max(counts, default=0),
        "max_collision_ratio": max(
            (Fraction(collision_sums[h - 1], N - h) for h in range(1, D + 1)),
            default=Fraction(0),
        ),
        "ZD": tail["ZD"],
    }


def finite_field_gates() -> None:
    orbit = load_orbit_function()
    print("FINITE_FIELD_GATES")
    print(
        "p N L D T0 S_D maxR Z_D Z_D/N coarse_bound "
        "model_S model_maxR model_Z_D model_maxCollision/mass"
    )
    for p in (1009, 3001, 10007):
        N = p - 1
        L = ceil(log(log(p)))
        D = ceil_sqrt(N * L * L)
        assert D < N
        points = orbit(p)
        assert len(points) == N
        counts = root_counts(points, D)
        assert all(value <= 3 * h - 3 for h, value in enumerate(counts, start=1))
        tail = zero_tail_data(counts, N)
        model = marginal_extremizer(N, D)
        print(
            f"{p} {N} {L} {D} {tail['T0']} {tail['S']} {tail['max_R']} "
            f"{tail['ZD']} {float(Fraction(tail['ZD'], N)):.9f} "
            f"{float(tail['coarse_bound']):.6f} "
            f"{model['S']} {model['max_R']} {model['ZD']} "
            f"{float(model['max_collision_ratio']):.9f}"
        )
        print(f"  R_HIST {tail['histogram']}")
    print("FINITE_FIELD_GATES PASS")


def main() -> None:
    symbolic_gates()
    finite_field_gates()
    print("ALL_BREAKWALL_GATES_PASS")


if __name__ == "__main__":
    main()
