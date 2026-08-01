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
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import product
from math import ceil, factorial, isqrt, log
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


def apery_numbers(limit: int) -> list[int]:
    """Return b_0,...,b_limit from the integral Apéry recurrence."""
    values = [1, 5]
    for n in range(2, limit + 1):
        numerator = (
            apery_polynomial(n - 1) * values[-1]
            - (n - 1) ** 3 * values[-2]
        )
        denominator = n**3
        assert numerator % denominator == 0
        values.append(numerator // denominator)
    return values[: limit + 1]


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

    numerators = build_gap_numerators(11)
    zero = Poly(0, X, domain=ZZ)

    def numerator(index: int) -> Poly:
        if index == 0:
            return zero
        value = numerators[index]
        assert value is not None
        return value

    apery = apery_numbers(6)
    assert apery[:6] == [1, 5, 73, 1445, 33001, 819005]
    reflection_checks = 0
    for h in range(2, 11):
        numerator_h = numerator(h)
        reflected = Poly(
            numerator_h.as_expr().subs(X, -h - 1 - X), X, domain=ZZ
        )
        assert reflected == ((-1) ** (h - 1)) * numerator_h
        if h % 2 == 0:
            central = Poly(2 * X + h + 1, X, domain=ZZ)
            assert (numerator_h % central).is_zero
        else:
            j = (h + 1) // 2
            pole_value = abs(int(numerator_h.eval(-j)))
            assert pole_value == factorial(j - 1) ** 6 * apery[j - 1] ** 2
        reflection_checks += 1

    renewal_n_checks = 0
    shifted_gcd_checks = 0
    for a in range(1, 7):
        for d in range(1, 7):
            if a + d > 10:
                continue
            shifted_d = Poly(
                numerator(d).as_expr().subs(X, X + a), X, domain=ZZ
            )
            shifted_previous = Poly(
                numerator(d - 1).as_expr().subs(X, X + a + 1),
                X,
                domain=ZZ,
            )
            rhs = (
                numerator(a + 1) * shifted_d
                - Poly((X + a + 1) ** 6, X, domain=ZZ)
                * numerator(a)
                * shifted_previous
            )
            assert numerator(a + d) == rhs
            assert gcd(numerator(a), shifted_d).degree() == 0
            renewal_n_checks += 1
            shifted_gcd_checks += 1

    backward_checks = 0
    for h in range(2, 9):
        for g in range(1, h):
            carrier = Poly(1, X, domain=ZZ)
            for j in range(h - g + 1, h):
                carrier *= Poly((X + j) ** 6, X, domain=ZZ)
            shifted_g = Poly(
                numerator(g).as_expr().subs(X, X + h - g), X, domain=ZZ
            )
            shifted_previous = Poly(
                numerator(g - 1).as_expr().subs(X, X + h - g),
                X,
                domain=ZZ,
            )
            assert carrier * numerator(h - g) == (
                numerator(h - 1) * shifted_g
                - numerator(h) * shifted_previous
            )
            backward_checks += 1

    gcd_checks = 0
    resultant_checks = 0
    for left in range(2, 10):
        for right in range(left + 1, 11):
            assert gcd(numerator(left), numerator(right)).degree() == 0
            gcd_checks += 1
            if right <= 7:
                assert numerator(left).resultant(numerator(right)) != 0
                resultant_checks += 1

    print(
        "SYMBOLIC_GATES PASS "
        f"renewal={renewal_checks} reflection={reflection_checks} "
        f"gap_renewal={renewal_n_checks} backward={backward_checks} "
        f"shifted_gcd={shifted_gcd_checks} cross_gcd={gcd_checks} "
        f"cross_resultant={resultant_checks}"
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


def mirror_components(
    points: list[tuple[int, int]], counts: list[int]
) -> tuple[list[int], list[int], list[int], int]:
    """Return kappa, endpoint epsilon, paired multiplicity mu, and checks.

    The exact decomposition is R_h=kappa_h+2*mu_h+epsilon_h.  Reflection
    pairs the interior bases by r -> N-h-r; r=0 is the sole truncated
    endpoint not preserved by this involution.
    """
    length = len(points)
    kappas: list[int] = []
    endpoints: list[int] = []
    pairs: list[int] = []
    partner_checks = 0
    for h, value in enumerate(counts, start=1):
        roots = {
            r
            for r in range(length - h)
            if points[r] == points[r + h]
        }
        assert len(roots) == value
        central = 1 if h % 2 == 0 else 0
        endpoint = 1 if 0 in roots else 0
        interior = roots - {0}
        for r in interior:
            partner = length - h - r
            assert 1 <= partner <= length - h - 1
            assert partner in interior
            partner_checks += 1
        if central:
            base = (length - h) // 2
            assert base in interior
            assert points[base] == points[base + h]
        assert len(interior) >= central
        assert (len(interior) - central) % 2 == 0
        pair_count = (len(interior) - central) // 2
        assert value == central + 2 * pair_count + endpoint
        kappas.append(central)
        endpoints.append(endpoint)
        pairs.append(pair_count)
    return kappas, endpoints, pairs, partner_checks


def mirror_weak_one(
    pairs: list[int], threshold: int
) -> tuple[int, Fraction]:
    """Return W_mu and the fixed harmonic factor in [MIRROR-WL1]."""
    diameter = len(pairs)
    upper = max(0, 3 * (diameter - 1) // 2)
    weak_one = max(
        (
            t * tail_count(pairs, t)
            for t in range(threshold + 1, upper + 1)
        ),
        default=0,
    )
    harmonic = sum(
        (Fraction(1, t) for t in range(threshold + 1, upper + 1)),
        start=Fraction(0, 1),
    )
    return weak_one, harmonic


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


class DisjointSet:
    """Tiny union-find used to materialize the reflection word model."""

    def __init__(self, size: int) -> None:
        self.parent = list(range(size))

    def find(self, value: int) -> int:
        while self.parent[value] != value:
            self.parent[value] = self.parent[self.parent[value]]
            value = self.parent[value]
        return value

    def union(self, left: int, right: int) -> None:
        left_root = self.find(left)
        right_root = self.find(right)
        if left_root != right_root:
            self.parent[right_root] = left_root


def mirror_word_countermodel_gate() -> None:
    """Verify the strict ZERO-TAIL-2/MIRROR-WL1 separation word."""
    length = 200_000
    diameter = 2_000
    heavy_lag = diameter - 1  # odd
    merged_orbits = 3 * (heavy_lag - 1) // 2

    sources = list(range(1, heavy_lag + 1))
    sources += list(
        range(
            2 * heavy_lag + 1,
            2 * heavy_lag + 1 + merged_orbits - heavy_lag,
        )
    )
    assert len(sources) == merged_orbits
    used_left = sources + [base + heavy_lag for base in sources]
    assert len(set(used_left)) == 2 * merged_orbits
    assert max(used_left) < 4 * heavy_lag
    assert length - 2 * max(used_left) > diameter

    partition = DisjointSet(length + 1)
    for index in range(length + 1):
        partition.union(index, length - index)
    for base in sources:
        partition.union(base, base + heavy_lag)

    colors: dict[int, list[int]] = defaultdict(list)
    for index in range(length):  # endpoint length is outside the orbit
        colors[partition.find(index)].append(index)

    counts_by_lag: Counter[int] = Counter()
    forward_degree: Counter[int] = Counter()
    for color_class in colors.values():
        assert len(color_class) <= 4
        for left_index, left in enumerate(color_class):
            for right in color_class[left_index + 1 :]:
                lag = right - left
                if lag <= diameter:
                    counts_by_lag[lag] += 1
                    forward_degree[left] += 1

    counts = [counts_by_lag[h] for h in range(1, diameter + 1)]
    kappas = mirror_skeleton(diameter)
    endpoints = [
        int(partition.find(0) == partition.find(h))
        for h in range(1, diameter + 1)
    ]
    assert sum(endpoints) == 0
    pairs = []
    for h, (value, central) in enumerate(zip(counts, kappas), start=1):
        assert value <= 3 * (h - 1)
        assert value >= central and (value - central) % 2 == 0
        pairs.append((value - central) // 2)
    assert counts[heavy_lag - 1] == 3 * (heavy_lag - 1)
    assert max(forward_degree.values(), default=0) <= 1
    pair_energy = sum(
        degree * (degree - 1) // 2 for degree in forward_degree.values()
    )
    assert pair_energy == 0
    color_energy = sum(len(color_class) ** 2 for color_class in colors.values())
    assert color_energy < 2 * length + 16 * diameter

    threshold, old_tail = zero_tail_two(counts, length)
    weak_one, harmonic = mirror_weak_one(pairs, threshold)
    stripped_square = sum(
        (value - central) ** 2
        for value, central in zip(counts, kappas)
    )
    assert old_tail > length
    assert weak_one * harmonic < length
    assert stripped_square > length
    assert stripped_square * stripped_square > diameter * length
    assert sum(pairs) <= diameter * threshold + weak_one * harmonic

    # A profile at the W1 boundary which fails [MIRROR-WL1], showing that
    # the new residual is not merely a relabeling of W1.
    root = 400
    boundary_length = root * root
    boundary_diameter = 5 * root
    boundary_pairs = [0] * boundary_diameter
    active = 0
    for h in range(boundary_diameter // 2 + 1, boundary_diameter + 1):
        if active == root:
            break
        boundary_pairs[h - 1] = root // 2
        active += 1
    boundary_counts = [
        central + 2 * pair_count
        for central, pair_count in zip(
            mirror_skeleton(boundary_diameter), boundary_pairs
        )
    ]
    for h, value in enumerate(boundary_counts, start=1):
        assert value <= 3 * (h - 1)
    boundary_threshold, _ = zero_tail_two(
        boundary_counts, boundary_length
    )
    boundary_weak_one, boundary_harmonic = mirror_weak_one(
        boundary_pairs, boundary_threshold
    )
    assert sum(boundary_counts) <= 2 * boundary_length
    assert boundary_weak_one * boundary_harmonic > boundary_length

    print(
        "MIRROR_WORD_GATE PASS "
        f"N={length} D={diameter} d={heavy_lag} "
        f"R_d={counts[heavy_lag - 1]} Q_D={pair_energy} "
        f"energy={color_energy} stripped2={stripped_square} "
        f"Z2/N={old_tail}/{length} "
        f"Wmu*H/N={float(weak_one * harmonic / length):.6f} "
        f"W1model_S={sum(boundary_counts)} "
        f"W1model_Wmu*H/N="
        f"{float(boundary_weak_one * boundary_harmonic / boundary_length):.6f}"
    )


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
        "D*W32sq/N^3 B sumMu2 Wmu WmuH/N maxR*D/N "
        "spike_old/excess/lorentz "
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

        kappas, endpoints, pairs, partner_checks = mirror_components(
            points, counts
        )
        endpoint_count = sum(endpoints)
        assert endpoint_count**3 <= 64 * (diameter + 1) ** 2
        pair_square = sum(pair_count**2 for pair_count in pairs)
        stripped_square = sum(
            (value - central) ** 2
            for value, central in zip(counts, kappas)
        )
        assert 4 * pair_square <= stripped_square
        assert stripped_square <= 8 * pair_square + 2 * endpoint_count
        assert sum(value**2 for value in counts) <= (
            3 * sum(kappas) + 12 * pair_square + 3 * endpoint_count
        )

        threshold, old_tail = zero_tail_two(counts, length)
        weak2_value = weak_two(counts)
        lorentz_sq = lorentz_three_halves_sq(counts)
        excess = zero_excess(counts, threshold)
        total = sum(counts)
        mirror_weak, mirror_harmonic = mirror_weak_one(pairs, threshold)

        # Exact forms of all new residuals and their elementary bridges.
        assert excess * excess * threshold <= length**2
        assert total**3 <= 27 * diameter * lorentz_sq
        assert diameter * lorentz_sq <= length**3
        assert sum(pairs) <= (
            diameter * threshold + mirror_weak * mirror_harmonic
        )
        assert mirror_weak * mirror_harmonic <= length
        assert 4 * (threshold + 1) * mirror_weak <= old_tail

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
            f"{endpoint_count} {pair_square} {mirror_weak} "
            f"{float(mirror_weak * mirror_harmonic / length):.9f} "
            f"{max(counts) * diameter}/{length} "
            f"{spike_old}/{spike_excess * spike_excess * threshold <= length**2}"
            f"/{diameter * spike_lorentz_sq <= length**3} "
            f"{sum(boundary)}/{boundary_excess * boundary_excess * threshold > length**2}"
            f"/{diameter * boundary_lorentz_sq > length**3} "
            f"partners={partner_checks}"
        )
    print("FINITE_FIELD_GATES PASS")


def main() -> None:
    symbolic_gates()
    exhaustive_lorentz_gate()
    mirror_word_countermodel_gate()
    finite_field_gates()
    print("ALL_BREAKWALL2_GATES_PASS")


if __name__ == "__main__":
    main()
