#!/usr/bin/env python3
"""Exact gates for the last-wall restart, primitive, and small-gap audits.

The Apéry convention is the one in ``CODEX_SPEC_lastwall.md``:

    N = p - 2,
    1 <= r <= N - d,
    C_d = #{r : pi(r) = pi(r+d)}.

All finite-field and counting checks use integer arithmetic.  Every cutoff
through floor(p^(3/5)) is checked.  Fractions are kept exact; the only
floating-point values printed are human-readable normalizations of quantities
already computed exactly.

The second half checks the explicit abstract countermodel requested in the
last-wall audit.  Its ``a_n`` are positions occupied by one repeated color;
every other position of the ambient word has a private color.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from fractions import Fraction
import hashlib
import json
import math


APERY_PRIMES = (199, 499, 997)
COUNTERMODEL_PRIMES = (3, 5, 7, 11, 13, 17, 19)


def apery_P(n: int) -> int:
    return 34 * n**3 + 51 * n**2 + 27 * n + 5


def legendre_symbol(a: int, p: int) -> int:
    value = pow(a % p, (p - 1) // 2, p)
    if value == 0:
        return 0
    return 1 if value == 1 else -1


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    for divisor in range(3, math.isqrt(n) + 1, 2):
        if n % divisor == 0:
            return False
    return True


def floor_three_fifths(p: int) -> int:
    """Return floor(p^(3/5)) using integer comparisons only."""

    lower, upper = 0, p
    while lower + 1 < upper:
        middle = (lower + upper) // 2
        if middle**5 <= p**3:
            lower = middle
        else:
            upper = middle
    assert lower**5 <= p**3 < (lower + 1) ** 5
    return lower


def restart_exact_bound(span: int) -> tuple[int, int]:
    """Best short/long adjacent-gap bound and an optimizing cutoff K."""

    if span == 0:
        return (0, 1)
    return min(
        (
            span // (cutoff + 1) + 3 * cutoff * (cutoff - 1) // 2,
            cutoff,
        )
        for cutoff in range(1, span + 1)
    )


def projective_key(x: int, y: int, p: int) -> tuple[int, int]:
    x %= p
    y %= p
    assert x != 0 or y != 0
    if x == 0:
        return (0, 1)
    return (1, y * pow(x, -1, p) % p)


def apery_orbit(p: int) -> tuple[list[int], list[int], list[tuple[int, int]]]:
    """Return b, c, and pi on indices 0,...,p-2."""

    N = p - 2
    b = [0] * (N + 1)
    c = [0] * (N + 1)
    b[0], b[1] = 1, 5 % p
    c[0], c[1] = 0, 1
    for n in range(1, N):
        inverse = pow((n + 1) ** 3 % p, -1, p)
        b[n + 1] = (
            (apery_P(n) * b[n] - n**3 * b[n - 1]) * inverse
        ) % p
        c[n + 1] = (
            (apery_P(n) * c[n] - n**3 * c[n - 1]) * inverse
        ) % p
    keys = [projective_key(bn, cn, p) for bn, cn in zip(b, c)]
    return b, c, keys


def trim(poly: list[int], p: int) -> list[int]:
    result = [coefficient % p for coefficient in poly]
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return result


def poly_sub(left: list[int], right: list[int], p: int) -> list[int]:
    size = max(len(left), len(right))
    result = [0] * size
    for i in range(size):
        result[i] = (
            (left[i] if i < len(left) else 0)
            - (right[i] if i < len(right) else 0)
        ) % p
    return trim(result, p)


def poly_mul(left: list[int], right: list[int], p: int) -> list[int]:
    result = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            result[i + j] = (result[i + j] + a * b) % p
    return trim(result, p)


def linear_power(shift: int, exponent: int, p: int) -> list[int]:
    """Coefficients of (X+shift)^exponent, low degree first."""

    return [
        math.comb(exponent, degree)
        * pow(shift, exponent - degree, p)
        % p
        for degree in range(exponent + 1)
    ]


def shifted_P(shift: int, p: int) -> list[int]:
    result = [0] * 4
    for coefficient, exponent in ((34, 3), (51, 2), (27, 1), (5, 0)):
        power = linear_power(shift, exponent, p)
        for degree, value in enumerate(power):
            result[degree] = (result[degree] + coefficient * value) % p
    return trim(result, p)


def poly_eval(poly: list[int], x: int, p: int) -> int:
    value = 0
    for coefficient in reversed(poly):
        value = (value * x + coefficient) % p
    return value


def gap_polynomials(p: int, D: int) -> list[list[int] | None]:
    """Build N_1,...,N_D directly in F_p[X]."""

    polynomials: list[list[int] | None] = [None] * (D + 1)
    polynomials[1] = [1]
    if D >= 2:
        polynomials[2] = shifted_P(1, p)
    for h in range(2, D):
        assert polynomials[h] is not None and polynomials[h - 1] is not None
        first = poly_mul(shifted_P(h, p), polynomials[h], p)
        second = poly_mul(linear_power(h, 6, p), polynomials[h - 1], p)
        polynomials[h + 1] = poly_sub(first, second, p)
    return polynomials


def lc_apparition_data(p: int, D: int) -> dict[str, object]:
    """Check the Lucas apparition law for lc(N_h) modulo p."""

    discriminant = 34**2 - 4
    ambient = p - legendre_symbol(discriminant, p)
    ell = [0] * (ambient + 1)
    ell[0] = 0
    ell[1] = 1
    for h in range(1, ambient):
        ell[h + 1] = (34 * ell[h] - ell[h - 1]) % p
    assert ell[ambient] == 0
    rho = next(h for h in range(1, ambient + 1) if ell[h] == 0)
    assert ambient % rho == 0
    zeros = [h for h in range(1, D + 1) if h % rho == 0]

    check = [0] * (D + 1)
    if D >= 1:
        check[1] = 1
    for h in range(1, D):
        check[h + 1] = (34 * check[h] - check[h - 1]) % p
    assert [h for h in range(1, D + 1) if check[h] == 0] == zeros
    return {
        "ambient_index": ambient,
        "rank": rho,
        "zero_rows": zeros,
        "ell_mod_p": check,
    }


def fraction_text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def localized_clique_gate(
    fibres: dict[tuple[int, int], list[int]], N: int, D: int, weighted: Fraction
) -> dict[str, object]:
    """Check every step of the localized inverse-square clique inequality."""

    localized_choose3 = 0
    localized_energy = Fraction(0)
    class_checks = 0
    for block in range(N // D + 1):
        lower = block * D
        upper = (block + 2) * D
        for positions in fibres.values():
            local = [position for position in positions if lower <= position < upper]
            if len(local) < 2:
                continue
            energy = Fraction(0)
            for i, left in enumerate(local):
                for right in local[i + 1 :]:
                    gap = right - left
                    if gap <= D:
                        energy += Fraction(1, gap * gap)
            triples = math.comb(len(local), 3)
            assert triples <= 11 * D * D * energy
            localized_choose3 += triples
            localized_energy += energy
            class_checks += 1

    assert localized_energy <= 2 * weighted
    return {
        "class_checks": class_checks,
        "localized_choose3": localized_choose3,
        "localized_energy": localized_energy,
    }


def apery_prime_gate(p: int) -> dict[str, object]:
    assert is_prime(p) and p >= 7
    N = p - 2
    D = floor_three_fifths(p)
    assert 2 <= D < N

    b, c, keys = apery_orbit(p)
    # Strong projective reflection on the complete regular germ.
    for n in range(1, N + 1):
        assert keys[n] == keys[p - 1 - n]

    polynomials = gap_polynomials(p, D)
    apparition = lc_apparition_data(p, D)
    ell = apparition["ell_mod_p"]
    assert isinstance(ell, list)

    degree_rows: list[int] = [0] * (D + 1)
    degree_drops: list[tuple[int, int, int]] = []
    for d in range(1, D + 1):
        polynomial = polynomials[d]
        assert polynomial is not None and polynomial != [0]
        nominal_degree = 3 * (d - 1)
        degree = len(polynomial) - 1
        degree_rows[d] = degree
        assert degree <= nominal_degree
        nominal_coefficient = (
            polynomial[nominal_degree]
            if nominal_degree < len(polynomial)
            else 0
        )
        assert nominal_coefficient == ell[d]
        if degree < nominal_degree:
            degree_drops.append((d, degree, nominal_degree))
    assert [row[0] for row in degree_drops] == apparition["zero_rows"]

    C = [0] * (D + 1)
    d_D = [0] * (N + 1)
    collision_gaps: list[list[int]] = [[] for _ in range(N + 1)]
    checked_pairs = 0
    determinant_checks = 0
    for r in range(1, N):
        max_d = min(D, N - r)
        if max_d == 0:
            continue
        previous = 1  # N_1(r)
        current = None
        denominator = 1
        for d in range(1, max_d + 1):
            if d == 1:
                value = previous
            elif d == 2:
                current = apery_P(r + 1) % p
                value = current
            else:
                assert current is not None
                following = (
                    apery_P(r + d - 1) * current
                    - (r + d - 1) ** 6 * previous
                ) % p
                previous, current = current, following
                value = current

            polynomial = polynomials[d]
            assert polynomial is not None
            assert value == poly_eval(polynomial, r, p)
            denominator = denominator * pow(r + d, 3, p) % p
            determinant = (b[r] * c[r + d] - b[r + d] * c[r]) % p
            assert determinant * denominator % p == value
            collision = keys[r] == keys[r + d]
            assert collision == (determinant == 0) == (value == 0)
            checked_pairs += 1
            determinant_checks += 1
            if collision:
                C[d] += 1
                d_D[r] += 1
                collision_gaps[r].append(d)

    # Root counts are also checked independently by dense polynomial evaluation.
    for d in range(1, D + 1):
        polynomial = polynomials[d]
        assert polynomial is not None
        full_field_root_count = sum(
            poly_eval(polynomial, r, p) == 0 for r in range(p)
        )
        dense_count = sum(
            poly_eval(polynomial, r, p) == 0 for r in range(1, N - d + 1)
        )
        assert dense_count == C[d]
        assert C[d] <= full_field_root_count <= degree_rows[d]
        assert C[d] <= 3 * (d - 1)

    # Every even row has its forced mirror root in the declared window.
    for d in range(2, D + 1, 2):
        mirror = (p - 1 - d) // 2
        assert 1 <= mirror <= N - d
        assert keys[mirror] == keys[mirror + d]
        polynomial = polynomials[d]
        assert polynomial is not None and poly_eval(polynomial, mirror, p) == 0

    # For p>13 none of the four omitted boundary values is an N_2 root.
    expected_C2 = 1 + 2 * (legendre_symbol(-51, p) == 1)
    assert C[2] == expected_C2

    S_D = sum(C[1:])
    assert S_D == sum(d_D)
    Q_D = sum(math.comb(multiplicity, 2) for multiplicity in d_D)
    max_multiplicity = max(d_D)
    assert sum(value * value for value in d_D) == S_D + 2 * Q_D
    assert S_D * S_D <= N * (S_D + 2 * Q_D)
    assert 2 * Q_D <= max(0, max_multiplicity - 1) * S_D
    # Banked restart scale M_D <= 4 D^(2/3), checked without floats.
    assert max_multiplicity**3 <= 64 * D * D

    fibres: dict[tuple[int, int], list[int]] = defaultdict(list)
    for r in range(1, N + 1):
        fibres[keys[r]].append(r)

    # Evaluate every continuant on F_p once, then exhaust the exact addition
    # law and its root-level restart corollary throughout the physical triangle.
    continuant_values = [[0] * p]
    for d in range(1, D + 1):
        polynomial = polynomials[d]
        assert polynomial is not None
        continuant_values.append(
            [poly_eval(polynomial, r, p) for r in range(p)]
        )
    addition_checks = 0
    restart_equivalence_checks = 0
    for r in range(1, N + 1):
        for a in range(1, min(D - 1, N - r) + 1):
            max_g = min(D - a, N - r - a)
            for g in range(1, max_g + 1):
                left = continuant_values[a + g][r]
                right = (
                    continuant_values[g][r + a]
                    * continuant_values[a + 1][r]
                    - pow(r + a + 1, 6, p)
                    * continuant_values[g - 1][r + a + 1]
                    * continuant_values[a][r]
                ) % p
                assert left == right
                addition_checks += 1
                if continuant_values[a][r] == 0:
                    assert continuant_values[a + 1][r] != 0
                    assert (left == 0) == (continuant_values[g][r + a] == 0)
                    restart_equivalence_checks += 1

    # Raw columns can acquire a zero ray after crossing the singular transfer
    # edge.  Check every such Apéry-triggered ray and that it starts strictly
    # beyond the nonwrapping physical triangle.
    polluted_rays: list[tuple[int, int]] = []
    for m in range(2, p):
        if b[m - 1] != 0:
            continue
        raw_base = (-m) % p
        raw_values = [0, 1]
        for h in range(1, p - 1):
            raw_values.append(
                (
                    apery_P(raw_base + h) * raw_values[h]
                    - (raw_base + h) ** 6 * raw_values[h - 1]
                )
                % p
            )
        assert all(raw_values[h] == 0 for h in range(m, p))
        polluted_rays.append((raw_base, m))
        if 1 <= raw_base <= N:
            assert N - raw_base == m - 2 < m

    # A projective fibre is the zero set of one homogeneous Apéry solution.
    # This directly checks that the zero-restart argument is target-independent.
    general_fibre_recurrence_checks = 0
    for positions in fibres.values():
        anchor = positions[0]
        solution = [
            (b[anchor] * c[n] - c[anchor] * b[n]) % p
            for n in range(N + 1)
        ]
        assert [n for n in range(1, N + 1) if solution[n] == 0] == positions
        for n in range(1, N):
            assert (
                (n + 1) ** 3 * solution[n + 1]
                - apery_P(n) * solution[n]
                + n**3 * solution[n - 1]
            ) % p == 0
            general_fibre_recurrence_checks += 1

    # Exhaust every cutoff H <= D.  This includes the projective-fibre window
    # theorem, primitive/renewal decomposition, moments, and all inequalities
    # later quoted in the report.
    cutoff_checks = 0
    fibre_window_checks = 0
    max_cutoff_row: dict[str, object] | None = None
    for H in range(1, D + 1):
        bound_H, optimizer_H = restart_exact_bound(H)
        assert bound_H**3 <= 64 * H * H
        for positions in fibres.values():
            left = 0
            max_occupancy = 0
            for right, position in enumerate(positions):
                while position - positions[left] > H:
                    left += 1
                max_occupancy = max(max_occupancy, right - left + 1)
            assert max_occupancy - 1 <= bound_H
            fibre_window_checks += 1

        multiplicities = [
            sum(gap <= H for gap in collision_gaps[r])
            for r in range(1, N + 1)
        ]
        S_H = sum(multiplicities)
        Q_H = sum(math.comb(value, 2) for value in multiplicities)
        M_H = max(multiplicities)
        assert S_H == sum(C[1 : H + 1])
        assert sum(value * value for value in multiplicities) == S_H + 2 * Q_H
        assert S_H * S_H <= N * (S_H + 2 * Q_H)
        assert 2 * Q_H <= max(0, M_H - 1) * S_H
        assert M_H <= bound_H

        first_return_endpoints: list[int] = []
        nonprimitive_count = 0
        split_multiplicity = 0
        split_histogram: Counter[int] = Counter()
        for r in range(1, N + 1):
            gaps = [gap for gap in collision_gaps[r] if gap <= H]
            if not gaps:
                continue
            first_return_endpoints.append(r + gaps[0])
            nonprimitive_count += len(gaps) - 1
            for index, gap in enumerate(gaps[1:], start=1):
                split_histogram[index] += 1
                for earlier_gap in gaps[:index]:
                    assert keys[r + earlier_gap] == keys[r + gap]
                    assert gap - earlier_gap in collision_gaps[r + earlier_gap]
                    split_multiplicity += 1
        primitive_count = len(first_return_endpoints)
        assert len(set(first_return_endpoints)) == primitive_count
        assert primitive_count <= N
        assert S_H == primitive_count + nonprimitive_count
        assert nonprimitive_count <= split_multiplicity == Q_H
        assert S_H <= N + Q_H

        weighted_H = sum(
            (Fraction(C[d], d * d) for d in range(1, H + 1)), Fraction(0)
        )
        assert Q_H <= 22 * H * H * weighted_H
        if H >= 2:
            assert weighted_H >= Fraction(1, 4)
        cutoff_checks += 1
        if H == D:
            max_cutoff_row = {
                "primitive_count": primitive_count,
                "nonprimitive_count": nonprimitive_count,
                "split_multiplicity": split_multiplicity,
                "split_histogram": dict(sorted(split_histogram.items())),
                "restart_exact_bound": bound_H,
                "restart_optimizer": optimizer_H,
            }
    assert max_cutoff_row is not None

    triple_count = 0
    for positions in fibres.values():
        for i, left in enumerate(positions):
            for j in range(i + 1, len(positions)):
                middle = positions[j]
                if middle - left >= D:
                    break
                for right in positions[j + 1 :]:
                    if right - left > D:
                        break
                    triple_count += 1
    assert triple_count == Q_D

    weighted = sum(
        (Fraction(C[d], d * d) for d in range(1, D + 1)), Fraction(0)
    )
    forced_weight = sum(
        (Fraction(1, d * d) for d in range(2, D + 1, 2)), Fraction(0)
    )
    assert weighted >= forced_weight >= Fraction(1, 4)
    inverse_square_upper = 22 * D * D * weighted
    assert Q_D <= inverse_square_upper

    localization = localized_clique_gate(fibres, N, D, weighted)
    assert Q_D <= localization["localized_choose3"]
    assert 11 * D * D * localization["localized_energy"] <= inverse_square_upper

    # Exact threshold at which the raw clique upper bound alone would certify Q_D <= N.
    certification_threshold = Fraction(N, 22 * D * D)
    assert weighted > certification_threshold
    Y = min(D, math.floor(math.log(p) ** 2))
    prefix_weight = sum(
        (Fraction(C[d], d * d) for d in range(1, Y + 1)), Fraction(0)
    )
    assert prefix_weight >= Fraction(1, 4) > certification_threshold
    chebotarev_benchmark = sum(
        (
            Fraction(1 + (d % 2 == 0), d * d)
            for d in range(1, D + 1)
        ),
        Fraction(0),
    )

    drop_rows_with_counts = [
        {
            "d": d,
            "degree": actual,
            "nominal_degree": nominal,
            "C_d": C[d],
        }
        for d, actual, nominal in degree_drops
    ]
    return {
        "p": p,
        "N": N,
        "D": D,
        "effective_L_squared": Fraction(D * D, N),
        "C": C[1:],
        "S_D": S_D,
        "Q_D": Q_D,
        "max_d_D": max_multiplicity,
        "primitive": max_cutoff_row,
        "polluted_rays": polluted_rays,
        "weighted_inverse_square_mass": weighted,
        "forced_mirror_weight": forced_weight,
        "nonforced_weight": weighted - forced_weight,
        "inverse_square_upper": inverse_square_upper,
        "certification_threshold_for_Q_le_N": certification_threshold,
        "threshold_failure_factor": weighted / certification_threshold,
        "polylog_prefix": {
            "Y": Y,
            "S_Y": sum(C[1 : Y + 1]),
            "weight": prefix_weight,
            "tail_weight": weighted - prefix_weight,
        },
        "fixed_d_chebotarev_benchmark_weight": chebotarev_benchmark,
        "leading_coefficient_apparition": {
            "ambient_index": apparition["ambient_index"],
            "rank": apparition["rank"],
            "degree_drop_rows": drop_rows_with_counts,
        },
        "gates": {
            "regular_pairs_checked": checked_pairs,
            "determinant_identity_checks": determinant_checks,
            "addition_identity_checks": addition_checks,
            "root_restart_checks": restart_equivalence_checks,
            "general_fibre_recurrence_checks": general_fibre_recurrence_checks,
            "all_cutoffs_checked": cutoff_checks,
            "fibre_windows_checked": fibre_window_checks,
            "localized_classes_checked": localization["class_checks"],
            "localized_choose3": localization["localized_choose3"],
            "localized_energy": localization["localized_energy"],
        },
    }


def countermodel_positions(q: int) -> list[int]:
    """Positions of the repeated color in the abstract countermodel."""

    assert is_prime(q) and q % 2 == 1
    return [2 * q * n + ((n % q) ** 2 % q) for n in range(q * q)]


def difference_counter_digest(R: Counter[int]) -> str:
    canonical = json.dumps(sorted(R.items()), separators=(",", ":"))
    return hashlib.sha256(canonical.encode("ascii")).hexdigest()


def next_prime_at_least(lower_bound: int) -> int:
    candidate = lower_bound if lower_bound % 2 else lower_bound + 1
    while not is_prime(candidate):
        candidate += 2
    return candidate


def countermodel_gate(q: int) -> dict[str, object]:
    positions = countermodel_positions(q)
    m = q * q
    assert len(positions) == m
    assert positions[0] == 0
    assert all(left < right for left, right in zip(positions, positions[1:]))
    H = positions[-1] - positions[0]
    assert H == 2 * q**3 - 2 * q + 1

    R: Counter[int] = Counter()
    for i, left in enumerate(positions):
        for right in positions[i + 1 :]:
            R[right - left] += 1
    assert sum(R.values()) == math.comb(m, 2)
    assert all(1 <= gap <= H for gap in R)
    # The requested gate is 3d; this construction in fact passes R_d <= d.
    assert all(multiplicity <= gap for gap, multiplicity in R.items())
    assert all(multiplicity <= 3 * gap for gap, multiplicity in R.items())

    # Directly verify that a difference d uniquely determines the index gap k.
    index_gap_for_difference: dict[int, int] = {}
    for i, left in enumerate(positions):
        for j in range(i + 1, m):
            difference = positions[j] - left
            index_gap = j - i
            old = index_gap_for_difference.setdefault(difference, index_gap)
            assert old == index_gap

    local_returns = [m - i - 1 for i in range(m)]
    max_local_multiplicity = max(local_returns)
    assert max_local_multiplicity == q * q - 1
    assert max_local_multiplicity**3 <= 64 * H * H

    S_H = sum(local_returns)
    Q_H = sum(math.comb(value, 2) for value in local_returns)
    assert S_H == math.comb(q * q, 2)
    assert Q_H == math.comb(q * q, 3)

    # In this one-color class the primitive edges are precisely consecutive
    # repeated-color positions.  Every nonprimitive edge has an intermediate
    # repeated-color point, and the total number of such renewal witnesses is
    # exactly Q_H (one for each ordered choice left < middle < right).
    primitive_edges = m - 1
    nonprimitive_edges = S_H - primitive_edges
    renewal_witness_mass = sum(
        max(0, j - i - 1) for i in range(m) for j in range(i + 1, m)
    )
    assert nonprimitive_edges == math.comb(m - 1, 2)
    assert nonprimitive_edges <= renewal_witness_mass == Q_H

    N = q**5
    assert positions[-1] < N
    L_squared = Fraction(H * H, N)
    Q_over_N = Fraction(Q_H, N)
    exact_ratio_formula = Fraction(q, 6) - Fraction(1, 2 * q) + Fraction(1, 3 * q**3)
    assert Q_over_N == exact_ratio_formula
    assert L_squared >= q
    if q >= 5:
        assert Q_over_N >= Fraction(q, 7)

    max_ratio_gap, max_ratio_multiplicity = max(
        R.items(), key=lambda item: (Fraction(item[1], item[0]), -item[0])
    )
    max_R_gap, max_R = max(R.items(), key=lambda item: (item[1], -item[0]))
    return {
        "q": q,
        "repeated_color_size": m,
        "H": H,
        "ambient_N": N,
        "S_H": S_H,
        "Q_H": Q_H,
        "primitive_edges": primitive_edges,
        "nonprimitive_edges": nonprimitive_edges,
        "renewal_witness_mass": renewal_witness_mass,
        "max_local_multiplicity": max_local_multiplicity,
        "nonzero_R_count": len(R),
        "max_R": max_R,
        "argmax_R": max_R_gap,
        "max_R_over_d": Fraction(max_ratio_multiplicity, max_ratio_gap),
        "R_sha256": difference_counter_digest(R),
        "L_squared": L_squared,
        "Q_over_N": Q_over_N,
    }


def symmetric_countermodel_gate(q: int) -> dict[str, object]:
    """Add the exact Apéry reflection pattern to the abstract word model.

    The repeated-color cluster is placed at the left edge and reflected to
    the right edge.  Every remaining reflection pair receives its own color.
    Thus there is exactly one isolated collision at every even gap, while the
    two clusters do not interact at gaps at most H.
    """

    base = countermodel_positions(q)
    repeated_size = q * q
    H = 2 * q**3 - 2 * q + 1
    p = next_prime_at_least(q**5)
    N = p - 2
    assert p < 2 * q**5  # finite instance of Bertrand's bound
    assert N > 3 * H + 1

    left_cluster = [1 + position for position in base]
    right_cluster = [N + 1 - position for position in reversed(left_cluster)]
    special = left_cluster + right_cluster
    special_set = set(special)
    assert len(special_set) == 2 * repeated_size
    assert [N + 1 - position for position in reversed(special)] == special
    assert right_cluster[0] - left_cluster[-1] > H

    base_R: Counter[int] = Counter()
    for i, left in enumerate(base):
        for right in base[i + 1 :]:
            base_R[right - left] += 1
    assert all(multiplicity <= gap for gap, multiplicity in base_R.items())

    special_R: Counter[int] = Counter()
    for i, left in enumerate(special):
        for right in special[i + 1 :]:
            gap = right - left
            if gap <= H:
                special_R[gap] += 1
    assert all(special_R[d] == 2 * base_R[d] for d in range(1, H + 1))

    C: Counter[int] = Counter()
    for d in range(1, H + 1):
        C[d] = special_R[d]
        if d % 2 == 0:
            mirror_left = (N + 1 - d) // 2
            mirror_right = N + 1 - mirror_left
            assert mirror_right - mirror_left == d
            assert mirror_left not in special_set and mirror_right not in special_set
            C[d] += 1
        assert C[d] <= 3 * (d - 1)

    S_H = sum(C.values())
    Q_H = 2 * math.comb(repeated_size, 3)
    max_local_multiplicity = repeated_size - 1
    primitive_edges = 2 * (repeated_size - 1) + H // 2
    nonprimitive_edges = 2 * math.comb(repeated_size - 1, 2)
    renewal_witness_mass = Q_H
    assert S_H == 2 * math.comb(repeated_size, 2) + H // 2
    assert S_H == primitive_edges + nonprimitive_edges
    assert nonprimitive_edges <= renewal_witness_mass
    assert max_local_multiplicity**3 <= 64 * H * H
    assert Q_H > 0

    L_squared = Fraction(H * H, N)
    Q_over_N = Fraction(Q_H, N)
    assert L_squared >= q
    if q >= 5:
        assert Q_over_N >= Fraction(q, 7)
    return {
        "q": q,
        "ambient_prime": p,
        "ambient_N": N,
        "H": H,
        "S_H": S_H,
        "Q_H": Q_H,
        "primitive_edges": primitive_edges,
        "nonprimitive_edges": nonprimitive_edges,
        "renewal_witness_mass": renewal_witness_mass,
        "max_local_multiplicity": max_local_multiplicity,
        "C_sha256": difference_counter_digest(C),
        "L_squared": L_squared,
        "Q_over_N": Q_over_N,
    }


def print_apery_result(result: dict[str, object]) -> None:
    print(
        f"p={result['p']} N={result['N']} D={result['D']} "
        f"S_D={result['S_D']} Q_D={result['Q_D']} max_d_D={result['max_d_D']}"
    )
    print(f"  C_d={result['C']}")
    primitive = result["primitive"]
    assert isinstance(primitive, dict)
    print(
        "  primitive/restart: "
        f"P={primitive['primitive_count']} "
        f"A={primitive['nonprimitive_count']} "
        f"split_mass={primitive['split_multiplicity']} "
        f"split_hist={primitive['split_histogram']} "
        f"restart_bound={primitive['restart_exact_bound']}"
        f"@K={primitive['restart_optimizer']}"
    )
    print(f"  cut-edge polluted raw rays={result['polluted_rays']}")
    print(
        "  inverse-square: "
        f"A_D={fraction_text(result['weighted_inverse_square_mass'])} "
        f"({float(result['weighted_inverse_square_mass']):.12f}), "
        f"mirror={fraction_text(result['forced_mirror_weight'])}, "
        f"nonforced={fraction_text(result['nonforced_weight'])}"
    )
    print(
        "  raw Q<=N certificate: "
        f"needed A_D<={fraction_text(result['certification_threshold_for_Q_le_N'])}; "
        f"failure factor={float(result['threshold_failure_factor']):.6f}"
    )
    prefix = result["polylog_prefix"]
    assert isinstance(prefix, dict)
    print(
        f"  polylog prefix Y={prefix['Y']}: S_Y={prefix['S_Y']} "
        f"A_Y={fraction_text(prefix['weight'])}, "
        f"tail={fraction_text(prefix['tail_weight'])}"
    )
    apparition = result["leading_coefficient_apparition"]
    assert isinstance(apparition, dict)
    print(
        "  lc apparition: "
        f"rho_p={apparition['rank']} in p-(Delta/p)={apparition['ambient_index']}; "
        f"degree drops={apparition['degree_drop_rows']}"
    )
    gates = result["gates"]
    assert isinstance(gates, dict)
    print(
        "  gates PASS: "
        f"pointwise={gates['regular_pairs_checked']}, "
        f"determinant={gates['determinant_identity_checks']}, "
        f"addition={gates['addition_identity_checks']}, "
        f"root_restart={gates['root_restart_checks']}, "
        f"fibre_recurrence={gates['general_fibre_recurrence_checks']}, "
        f"cutoffs={gates['all_cutoffs_checked']}, "
        f"fibre_windows={gates['fibre_windows_checked']}, "
        f"localized classes={gates['localized_classes_checked']}"
    )


def print_countermodel_result(result: dict[str, object]) -> None:
    print(
        f"q={result['q']} m={result['repeated_color_size']} H={result['H']} "
        f"N={result['ambient_N']} Q_H={result['Q_H']} "
        f"max_local={result['max_local_multiplicity']}"
    )
    print(
        "  primitive/renewal: "
        f"primitive={result['primitive_edges']} "
        f"nonprimitive={result['nonprimitive_edges']} "
        f"witness_mass={result['renewal_witness_mass']}"
    )
    print(
        "  differences PASS: "
        f"support={result['nonzero_R_count']} max_R={result['max_R']} "
        f"at d={result['argmax_R']}, "
        f"max(R_d/d)={float(result['max_R_over_d']):.6f}, "
        f"sha256={result['R_sha256']}"
    )
    print(
        f"  supercritical: L^2={fraction_text(result['L_squared'])} "
        f"({float(result['L_squared']):.6f}), "
        f"Q_H/N={fraction_text(result['Q_over_N'])} "
        f"({float(result['Q_over_N']):.6f})"
    )


def print_symmetric_countermodel_result(result: dict[str, object]) -> None:
    print(
        f"q={result['q']} p={result['ambient_prime']} N={result['ambient_N']} "
        f"H={result['H']} S_H={result['S_H']} Q_H={result['Q_H']} "
        f"max_local={result['max_local_multiplicity']}"
    )
    print(
        "  reflection/primitive: "
        f"primitive={result['primitive_edges']} "
        f"nonprimitive={result['nonprimitive_edges']} "
        f"witness_mass={result['renewal_witness_mass']} "
        f"C_sha256={result['C_sha256']}"
    )
    print(
        f"  supercritical: L^2={fraction_text(result['L_squared'])} "
        f"({float(result['L_squared']):.6f}), "
        f"Q_H/N={fraction_text(result['Q_over_N'])} "
        f"({float(result['Q_over_N']):.6f})"
    )


def main() -> None:
    print("=== Apéry last-wall exact gates ===")
    apery_results = [apery_prime_gate(p) for p in APERY_PRIMES]
    for result in apery_results:
        print_apery_result(result)

    print("\n=== Abstract one-color countermodel exact gates ===")
    countermodel_results = [countermodel_gate(q) for q in COUNTERMODEL_PRIMES]
    for result in countermodel_results:
        print_countermodel_result(result)

    print("\n=== Reflection-symmetric prime-length countermodel exact gates ===")
    symmetric_results = [
        symmetric_countermodel_gate(q) for q in COUNTERMODEL_PRIMES
    ]
    for result in symmetric_results:
        print_symmetric_countermodel_result(result)

    print("\nALL EXACT GATES PASS")
    print(
        "VERDICT: raw inverse-square small-d input is terminally insufficient: "
        "its required weighted mass is O(N/D^2), while forced even-gap mirror "
        "roots already contribute at least 1/4.  Fixed-d Chebotarev and "
        "leading-coefficient apparition do not alter that nonnegative mass."
    )
    print(
        "VERDICT: the abstract family has R_d<=d, restart-scale local "
        "multiplicity, and all renewal identities of a genuine word, yet "
        "D/sqrt(N)->infinity and Q_D/N->infinity."
    )
    print(
        "VERDICT: after reflection symmetrization at prime length, it also has "
        "the forced collision at every even gap and C_d<=3(d-1), while the "
        "same supercritical Q_D/N obstruction remains."
    )


if __name__ == "__main__":
    main()
