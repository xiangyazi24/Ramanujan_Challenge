#!/usr/bin/env python3
"""Exact verification gates for ``CODEX_BREAKWALL2_report.md``.

The symbolic part independently rebuilds the Apéry gap continuants, checks
the renewal identity, the reflection factor, and small-pair characteristic-
zero coprimality.  The finite-field part extracts ``orbit(p)`` verbatim from
``CRON_b1_crosscorr.py`` and evaluates two residuals strictly below the old
truncated weak-L2 wall at the three primes required by the specification.

For

    W32 = max_{t >= 1} t^(3/2) * #{h : R_h >= t},

all comparisons are kept integral by storing

    W32_SQ = max_{t >= 1} t^3 * #{h : R_h >= t}^2.

Thus ``W32 <= N^(3/2)/sqrt(D)`` is exactly ``D*W32_SQ <= N^3``.

For the declared integrated-tail residual, put

    E = sum_h (R_h - T0)_+ = sum_{t > T0} #{h : R_h >= t}.

The concrete face ``[ZERO-EXCESS-1/2]`` is checked without floating point as
``E^2*T0 <= N^2``.
"""

from __future__ import annotations

import ast
from itertools import product
from math import ceil, isqrt, log
from pathlib import Path

from sympy import Poly, ZZ, gcd, symbols


HERE = Path(__file__).resolve().parent
X = symbols("X")


def apery_polynomial(z):
    """The cubic coefficient in the Apéry recurrence."""
    return 34 * z**3 + 51 * z**2 + 27 * z + 5


def ceil_sqrt(value: int) -> int:
    """Return the least integer whose square is at least ``value``."""
    root = isqrt(value)
    return root if root * root == value else root + 1


def build_continuants(limit: int) -> list[Poly]:
    """Build K_0,...,K_limit over Z[X] in the gap variable."""
    values = [
        Poly(1, X, domain=ZZ),
        Poly(apery_polynomial(X), X, domain=ZZ),
    ]
    for m in range(1, limit):
        values.append(
            Poly(apery_polynomial(X + m), X, domain=ZZ) * values[m]
            - Poly((X + m) ** 6, X, domain=ZZ) * values[m - 1]
        )
    return values


def build_gap_numerators(limit: int) -> list[Poly | None]:
    """Build N_h, where N_1=1 and deg N_h=3h-3."""
    values: list[Poly | None] = [None, Poly(1, X, domain=ZZ)]
    if limit >= 2:
        values.append(Poly(apery_polynomial(X + 1), X, domain=ZZ))
    for h in range(2, limit):
        assert values[h] is not None and values[h - 1] is not None
        values.append(
            Poly(apery_polynomial(X + h), X, domain=ZZ) * values[h]
            - Poly((X + h) ** 6, X, domain=ZZ) * values[h - 1]
        )
    return values


def symbolic_gates() -> None:
    """Check the exact algebra used in the route audit."""
    continuants = build_continuants(11)
    renewal_checks = 0
    for m in range(0, 5):
        for g in range(1, 5):
            lhs = continuants[m + g + 1]
            rhs = Poly(
                continuants[m + 1].as_expr()
                * continuants[g].as_expr().subs(X, X + m + 1)
                - (X + m + 1) ** 6
                * continuants[m].as_expr()
                * continuants[g - 1].as_expr().subs(X, X + m + 2),
                X,
                domain=ZZ,
            )
            assert lhs == rhs
            renewal_checks += 1

    numerators = build_gap_numerators(10)
    reflection_checks = 0
    for h in range(2, 11):
        numerator = numerators[h]
        assert numerator is not None
        reflected = Poly(
            numerator.as_expr().subs(X, -h - 1 - X), X, domain=ZZ
        )
        assert reflected == ((-1) ** (h - 1)) * numerator
        if h % 2 == 0:
            central = Poly(2 * X + h + 1, X, domain=ZZ)
            assert (numerator % central).is_zero
        reflection_checks += 1

    gcd_checks = 0
    resultant_checks = 0
    for left in range(2, 10):
        for right in range(left + 1, 11):
            assert numerators[left] is not None and numerators[right] is not None
            assert gcd(numerators[left], numerators[right]).degree() == 0
            gcd_checks += 1
            if right <= 7:
                assert numerators[left].resultant(numerators[right]) != 0
                resultant_checks += 1

    print(
        "SYMBOLIC_GATES PASS "
        f"renewal={renewal_checks} reflection={reflection_checks} "
        f"cross_gcd={gcd_checks} cross_resultant={resultant_checks}"
    )


def load_orbit_function():
    """Compile only ``orbit(p)`` from the campaign experiment."""
    source_path = HERE / "CRON_b1_crosscorr.py"
    tree = ast.parse(source_path.read_text(encoding="utf-8"), str(source_path))
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


def root_counts(points: list[tuple[int, int]], diameter: int) -> list[int]:
    """Return the exact nonwrapping lag profile R_1,...,R_D."""
    length = len(points)
    return [
        sum(points[r] == points[r + h] for r in range(length - h))
        for h in range(1, diameter + 1)
    ]


def tail_count(counts: list[int], threshold: int) -> int:
    """Return #{h : R_h >= threshold}."""
    return sum(value >= threshold for value in counts)


def weak_two(counts: list[int]) -> int:
    """Return max_t t^2 A(t)."""
    maximum = max(counts, default=0)
    return max(
        (t * t * tail_count(counts, t) for t in range(1, maximum + 1)),
        default=0,
    )


def lorentz_three_halves_sq(counts: list[int]) -> int:
    """Return the exact square of max_t t^(3/2) A(t)."""
    maximum = max(counts, default=0)
    return max(
        (
            t**3 * tail_count(counts, t) ** 2
            for t in range(1, maximum + 1)
        ),
        default=0,
    )


def zero_tail_two(counts: list[int], length: int) -> tuple[int, int]:
    """Return T0 and the old truncated weak-L2 quantity."""
    diameter = len(counts)
    threshold = 1
    while threshold * threshold * diameter < length:
        threshold += 1
    maximum = max(counts, default=0)
    value = max(
        (
            t * t * tail_count(counts, t)
            for t in range(threshold + 1, maximum + 1)
        ),
        default=0,
    )
    return threshold, value


def zero_excess(counts: list[int], threshold: int) -> int:
    """Return sum_h (R_h-T0)_+, the integrated tail above T0."""
    return sum(max(value - threshold, 0) for value in counts)


def mirror_skeleton(diameter: int) -> list[int]:
    """The one forced central collision at every even lag."""
    return [1 if h % 2 == 0 else 0 for h in range(1, diameter + 1)]


def spike_profile(diameter: int) -> list[int]:
    """A mirror-compatible profile separating the new and old residuals."""
    counts = mirror_skeleton(diameter)
    heavy_lag = diameter if diameter % 2 == 0 else diameter - 1
    assert heavy_lag >= 2
    counts[heavy_lag - 1] = heavy_lag + 1
    for h, value in enumerate(counts, start=1):
        assert value <= 3 * (h - 1)
        assert value % 2 == (1 if h % 2 == 0 else 0)
    return counts


def w1_not_lorentz_profile(length: int, diameter: int) -> list[int]:
    """Build an O(N)-mass profile which fails the Lorentz residual.

    Heavy entries are placed only in the upper half of the lag interval, so
    the degree cap and the mirror parity can both be respected.
    """
    counts = mirror_skeleton(diameter)
    remaining = length - sum(counts)
    base = max(2, diameter // 2)
    for h in range((diameter + 1) // 2, diameter + 1):
        required_parity = 1 if h % 2 == 0 else 0
        heavy = base if base % 2 == required_parity else base + 1
        delta = heavy - counts[h - 1]
        if delta <= remaining and heavy <= 3 * (h - 1):
            counts[h - 1] = heavy
            remaining -= delta
    assert sum(counts) <= length
    return counts


def exhaustive_lorentz_gate() -> None:
    """Exhaustively test S^3 <= 27 D W32^2 on small profiles."""
    checked = 0
    for diameter in range(1, 7):
        for counts_tuple in product(range(5), repeat=diameter):
            counts = list(counts_tuple)
            total = sum(counts)
            square = lorentz_three_halves_sq(counts)
            assert total**3 <= 27 * diameter * square
            checked += 1
    print(f"LORENTZ_FINITE_GATE PASS profiles={checked}")


def finite_field_gates() -> None:
    """Evaluate all exact quantities at the three mandated primes."""
    orbit = load_orbit_function()
    print(
        "FINITE_FIELD_GATES\n"
        "p N L D T0 S maxR W2 Z2 E E^2*T0/N^2 W32sq "
        "D*W32sq/N^3 maxR*D/N eps spike_old/excess/lorentz "
        "W1model_S/fails_excess/fails_lorentz"
    )
    for prime in (1009, 3001, 10007):
        length = prime - 1
        slow = ceil(log(log(prime)))
        diameter = ceil_sqrt(length * slow * slow)
        points = orbit(prime)
        assert len(points) == length
        counts = root_counts(points, diameter)
        assert all(
            value <= 3 * (h - 1)
            for h, value in enumerate(counts, start=1)
        )

        # On the truncated nonwrapping orbit, reflection pairs all roots
        # except r=0.  Hence R_h=kappa_h+2q_h+epsilon_h exactly, where
        # kappa_h=1_{2|h} and epsilon_h=1_{pi(h)=pi(0)}.
        endpoint_count = 0
        for h, value in enumerate(counts, start=1):
            central = 1 if h % 2 == 0 else 0
            endpoint = 1 if points[0] == points[h] else 0
            endpoint_count += endpoint
            assert value >= central + endpoint
            assert (value - central - endpoint) % 2 == 0
            if central:
                base = (prime - 1 - h) // 2
                assert 0 <= base < base + h < length
                assert points[base] == points[base + h]

        threshold, old_tail = zero_tail_two(counts, length)
        weak2_value = weak_two(counts)
        lorentz_sq = lorentz_three_halves_sq(counts)
        excess = zero_excess(counts, threshold)
        total = sum(counts)

        # Exact forms of both new residuals and the Lorentz embedding.
        assert excess * excess * threshold <= length**2
        assert total**3 <= 27 * diameter * lorentz_sq
        assert diameter * lorentz_sq <= length**3

        spike = spike_profile(diameter)
        spike_threshold, spike_old = zero_tail_two(spike, length)
        assert spike_threshold == threshold
        spike_excess = zero_excess(spike, threshold)
        spike_lorentz_sq = lorentz_three_halves_sq(spike)
        assert spike_old > length
        assert spike_excess * spike_excess * threshold <= length**2
        assert diameter * spike_lorentz_sq <= length**3
        assert sum(value * value for value in spike) > length

        boundary = w1_not_lorentz_profile(length, diameter)
        boundary_excess = zero_excess(boundary, threshold)
        boundary_lorentz_sq = lorentz_three_halves_sq(boundary)
        assert sum(boundary) <= length
        assert boundary_excess * boundary_excess * threshold > length**2
        assert diameter * boundary_lorentz_sq > length**3

        print(
            f"{prime} {length} {slow} {diameter} {threshold} {total} "
            f"{max(counts)} {weak2_value} {old_tail} {excess} "
            f"{excess * excess * threshold}/{length**2} {lorentz_sq} "
            f"{diameter * lorentz_sq}/{length**3} "
            f"{max(counts) * diameter}/{length} "
            f"{endpoint_count} "
            f"{spike_old}/{spike_excess * spike_excess * threshold <= length**2}"
            f"/{diameter * spike_lorentz_sq <= length**3} "
            f"{sum(boundary)}/{boundary_excess * boundary_excess * threshold > length**2}"
            f"/{diameter * boundary_lorentz_sq > length**3}"
        )
    print("FINITE_FIELD_GATES PASS")


def main() -> None:
    symbolic_gates()
    exhaustive_lorentz_gate()
    finite_field_gates()
    print("ALL_BREAKWALL2_GATES_PASS")


if __name__ == "__main__":
    main()
