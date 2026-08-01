#!/usr/bin/env python3
"""Exact finite checks for ``CODEX_SPEC_strike2.md``.

The script has two independent parts.

* For the four requested Apery primes it constructs the projective orbit and
  the gap continuants modulo p, checks that both constructions give exactly
  the same nonwrapping collisions, and records first-return and dyadic-shell
  statistics at the two requested cutoffs.
* For q=5,7,11 it constructs the reflection-symmetric word from (4.10) of
  ``CODEX_LASTWALL_report.md`` literally, rather than just evaluating the
  closed formulas for its statistics.

All counts and finite-field operations are integer exact.  Decimal arithmetic
is used only to choose ceil(sqrt(p) log(p)); the result is checked far from an
integer boundary.  Displayed decimal ratios are derived from exact fractions.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from decimal import Decimal, getcontext
from fractions import Fraction
import hashlib
import math
from pathlib import Path
from typing import DefaultDict, Iterable, Sequence


APERY_PRIMES = (997, 1999, 4001, 7919)
WORD_QS = (5, 7, 11)
getcontext().prec = 80


def apery_P(n: int) -> int:
    return 34 * n**3 + 51 * n**2 + 27 * n + 5


def ceil_decimal(value: Decimal) -> int:
    integer = int(value)
    return integer if value == integer else integer + 1


def ceil_three_fifths(p: int) -> int:
    """Smallest D with D^5 >= p^3, using integer comparisons only."""

    target = p**3
    lower, upper = 0, p + 1
    while lower + 1 < upper:
        middle = (lower + upper) // 2
        if middle**5 >= target:
            upper = middle
        else:
            lower = middle
    assert (upper - 1) ** 5 < target <= upper**5
    return upper


def requested_cutoffs(p: int) -> tuple[tuple[str, int], ...]:
    real_cutoff = Decimal(p).sqrt() * Decimal(p).ln()
    fractional = real_cutoff % 1
    assert min(fractional, 1 - fractional) > Decimal("1e-50")
    return (
        ("ceil(sqrt(p) log(p))", ceil_decimal(real_cutoff)),
        ("ceil(p^(3/5))", ceil_three_fifths(p)),
    )


def inverse_table(p: int) -> list[int]:
    inverse = [0] * p
    inverse[1] = 1
    for value in range(2, p):
        inverse[value] = (p - (p // value) * inverse[p % value] % p) % p
    return inverse


def apery_orbit(p: int) -> tuple[list[int], list[int], list[int]]:
    """Return b, c and normalized projective keys on 0,...,p-2."""

    N = p - 2
    inverse = inverse_table(p)
    b = [0] * (N + 1)
    c = [0] * (N + 1)
    b[0], b[1] = 1, 5 % p
    c[0], c[1] = 0, 1
    for n in range(1, N):
        inv_cube = inverse[n + 1] ** 3 % p
        b[n + 1] = (
            (apery_P(n) * b[n] - n**3 * b[n - 1]) * inv_cube
        ) % p
        c[n + 1] = (
            (apery_P(n) * c[n] - n**3 * c[n - 1]) * inv_cube
        ) % p
    # [b:c] is represented by b/c, with p as the infinity sentinel.
    keys = [p if cn == 0 else bn * inverse[cn] % p for bn, cn in zip(b, c)]
    return b, c, keys


def dyadic_shells(D: int) -> list[tuple[int, int]]:
    """Disjoint shells (lo,hi] whose union is [1,D], largest first."""

    shells = []
    high = D
    while high:
        low = high // 2
        shells.append((low, high))
        high = low
    assert sorted(d for lo, hi in shells for d in range(lo + 1, hi + 1)) == list(
        range(1, D + 1)
    )
    return shells


def nearest_rank_quantiles(values: Sequence[int]) -> tuple[int, ...]:
    ordered = sorted(values)
    assert ordered
    result = [ordered[0]]
    for j in range(1, 11):
        result.append(ordered[(j * len(ordered) + 9) // 10 - 1])
    return tuple(result)


def equal_width_first_return_bins(first: Counter[int], D: int) -> list[tuple[int, int, int]]:
    result = []
    for j in range(10):
        low = j * D // 10
        high = (j + 1) * D // 10
        if high == low:
            continue
        result.append((low + 1, high, sum(first[d] for d in range(low + 1, high + 1))))
    assert sum(count for _, _, count in result) == sum(first.values())
    return result


def collision_lags(p: int, max_D: int) -> tuple[list[list[int]], int]:
    """Build all collisions twice and return the per-base lag lists.

    The first construction compares normalized orbit values.  The second
    propagates N_d(r) by the continuant recurrence simultaneously for all r.
    """

    N = p - 2
    _, _, keys = apery_orbit(p)
    for n in range(1, N + 1):
        assert keys[n] == keys[p - 1 - n]

    p_values = [apery_P(x) % p for x in range(p)]
    sixth = [pow(x, 6, p) for x in range(p)]
    previous = [0] * p
    current = [1] * p
    lags = [[] for _ in range(N + 1)]
    checks = 0
    for d in range(1, max_D + 1):
        for r in range(1, N - d + 1):
            orbit_collision = keys[r] == keys[r + d]
            continuant_collision = current[r] == 0
            assert orbit_collision == continuant_collision
            checks += 1
            if orbit_collision:
                lags[r].append(d)
        if d < max_D:
            following = [0] * p
            for r in range(p):
                x = (r + d) % p
                following[r] = (
                    p_values[x] * current[r] - sixth[x] * previous[r]
                ) % p
            previous, current = current, following
    return lags, checks


def analyze_cutoff(p: int, D: int, all_lags: Sequence[Sequence[int]]) -> dict[str, object]:
    N = p - 2
    lags = [[d for d in all_lags[r] if d <= D] for r in range(N + 1)]
    multiplicities = [len(lags[r]) for r in range(1, N + 1)]
    S = sum(multiplicities)
    P = sum(value > 0 for value in multiplicities)
    Q = sum(math.comb(value, 2) for value in multiplicities)
    first = Counter(values[0] for values in lags[1:] if values)
    assert sum(first.values()) == P
    assert P <= S

    # First-return endpoints are injective, and the return is to the next
    # occurrence of that projective value.
    _, _, keys = apery_orbit(p)
    endpoints: dict[int, int] = {}
    for r in range(1, N + 1):
        if not lags[r]:
            continue
        endpoint = r + lags[r][0]
        assert endpoint not in endpoints
        endpoints[endpoint] = r
        assert keys[r] == keys[endpoint]
        assert all(keys[r] != keys[middle] for middle in range(r + 1, endpoint))

    # Every intersection of two same-base zero sets is exactly a restarted
    # collision.  Summing these intersections is Q.
    cascade_checks = 0
    for r in range(1, N + 1):
        values = lags[r]
        for i, a in enumerate(values):
            for b in values[i + 1 :]:
                assert keys[r + a] == keys[r + b]
                assert b - a in all_lags[r + a]
                cascade_checks += 1
    assert cascade_checks == Q

    # Full finite inclusion-exclusion of the row zero sets.
    maximum = max(multiplicities, default=0)
    intersection_sums = {
        order: sum(math.comb(value, order) for value in multiplicities)
        for order in range(1, maximum + 1)
    }
    assert intersection_sums.get(1, 0) == S
    assert intersection_sums.get(2, 0) == Q
    assert sum((-1) ** (order + 1) * value for order, value in intersection_sums.items()) == P
    if maximum:
        # The optimal quadratic majorant of 1_{k>0} based on 1<=k<=M is
        # 1 <= k-2*binom(k,2)/M.  The triangle inequality is squared to
        # avoid irrational arithmetic.
        assert maximum * P <= maximum * S - 2 * Q
    assert 9 * Q * Q <= 2 * S**3

    shells = []
    shell_union_sum = 0
    for low, high in dyadic_shells(D):
        counts = [sum(low < d <= high for d in lags[r]) for r in range(1, N + 1)]
        row_sum = sum(counts)
        union = sum(value > 0 for value in counts)
        pair_intersections = sum(math.comb(value, 2) for value in counts)
        primitive = sum(first[d] for d in range(low + 1, high + 1))
        shell_union_sum += union
        shells.append(
            {
                "low": low,
                "high": high,
                "sum": row_sum,
                "union": union,
                "excess": row_sum - union,
                "pair_intersections": pair_intersections,
                "primitive": primitive,
            }
        )
    assert sum(int(shell["sum"]) for shell in shells) == S
    assert sum(int(shell["primitive"]) for shell in shells) == P
    assert shell_union_sum >= P

    ratio = Fraction(P * Q, N * N)
    first_values = list(first.elements())
    return {
        "p": p,
        "N": N,
        "D": D,
        "S": S,
        "P": P,
        "Q": Q,
        "ratio": ratio,
        "multiplicity_histogram": Counter(multiplicities),
        "first_quantiles": nearest_rank_quantiles(first_values),
        "first_bins": equal_width_first_return_bins(first, D),
        "intersection_sums": intersection_sums,
        "shells": shells,
        "shell_union_sum": shell_union_sum,
        "cross_shell_excess": shell_union_sum - P,
        "cascade_checks": cascade_checks,
    }


def is_prime(value: int) -> bool:
    if value < 2:
        return False
    if value % 2 == 0:
        return value == 2
    for divisor in range(3, math.isqrt(value) + 1, 2):
        if value % divisor == 0:
            return False
    return True


def first_prime_at_least(value: int) -> int:
    candidate = value if value % 2 else value + 1
    while not is_prime(candidate):
        candidate += 2
    return candidate


def bounded_fibre_statistics(
    fibres: Iterable[Sequence[int]], D: int
) -> tuple[int, int, int, list[int]]:
    bases: set[int] = set()
    outgoing: DefaultDict[int, int] = defaultdict(int)
    row_counts = [0] * (D + 1)
    S = 0
    for positions in fibres:
        ordered = sorted(positions)
        for index, left in enumerate(ordered):
            for right in ordered[index + 1 :]:
                gap = right - left
                if gap > D:
                    break
                bases.add(left)
                outgoing[left] += 1
                row_counts[gap] += 1
                S += 1
    Q = sum(math.comb(value, 2) for value in outgoing.values())
    return len(bases), Q, S, row_counts


def reflection_word_gate(q: int) -> dict[str, object]:
    """Construct (4.10) literally and audit all bounded differences."""

    assert is_prime(q) and q % 2 == 1
    p = first_prime_at_least(q**5)
    N = p - 2
    H = 2 * q**3 - 2 * q + 1
    left = [1 + 2 * q * n + ((n % q) ** 2 % q) for n in range(q * q)]
    assert left[-1] - left[0] == H
    special = set(left)
    special.update(N + 1 - position for position in left)
    assert len(special) == 2 * q * q
    assert max(left) + H < min(N + 1 - position for position in left)

    fibres: DefaultDict[int, list[int]] = defaultdict(list)
    colors = [None] * (N + 1)
    for position in range(1, N + 1):
        color = 0 if position in special else min(position, N + 1 - position)
        colors[position] = color
        fibres[color].append(position)
    assert all(colors[position] == colors[N + 1 - position] for position in range(1, N + 1))
    assert all(colors[position] != colors[position + 1] for position in range(1, N))

    P, Q, S, rows = bounded_fibre_statistics(fibres.values(), H)
    assert all(rows[d] <= 3 * (d - 1) for d in range(1, H + 1))
    expected_P = q**3 + 2 * q**2 - q - 2
    expected_Q = 2 * math.comb(q**2, 3)
    expected_S = 2 * math.comb(q**2, 2) + q**3 - q
    assert (P, Q, S) == (expected_P, expected_Q, expected_S)
    ratio = Fraction(P * Q, N * N)
    return {
        "q": q,
        "p": p,
        "N": N,
        "H": H,
        "P": P,
        "Q": Q,
        "S": S,
        "product": P * Q,
        "ratio": ratio,
        "nonzero_rows": sum(value > 0 for value in rows[1:]),
    }


def all_window_multiplicity_gate(fibres: Iterable[Sequence[int]]) -> int:
    """Check m_J(v)-1 <= 4 span(J)^(2/3) for every tight fibre window.

    Cubing removes floating point: (m-1)^3 <= 64 span^2.  It is enough to
    check windows whose endpoints are occurrences of the fibre value.
    """

    checks = 0
    for positions in fibres:
        ordered = sorted(positions)
        for left_index, left in enumerate(ordered):
            for right_index in range(left_index + 1, len(ordered)):
                span = ordered[right_index] - left
                excess = right_index - left_index
                assert excess**3 <= 64 * span**2
                checks += 1
    return checks


def padded_reflection_word_gate(q: int) -> dict[str, object]:
    """Pad the triangle-rich word by primitive-only four-point fibres.

    In the left half, a block of length 2d consists of d disjoint pairs at
    gap d.  Reflection creates a second copy.  Every pair has a private
    colour, so these blocks increase P but not Q.  Only large, distinct d are
    used, and a row is accepted only when its literal count remains at most
    3(d-1).  This is a finite hostile test of any purported word-level P-Q
    tradeoff.
    """

    base = reflection_word_gate(q)
    p = int(base["p"])
    N = int(base["N"])
    D = int(base["H"])
    left_special = [
        1 + 2 * q * n + ((n % q) ** 2 % q) for n in range(q * q)
    ]
    special = set(left_special)
    special.update(N + 1 - position for position in left_special)

    colors = [None] * (N + 1)
    for position in range(1, N + 1):
        colors[position] = (
            0 if position in special else min(position, N + 1 - position)
        )
    initial_fibres: DefaultDict[int, list[int]] = defaultdict(list)
    for position in range(1, N + 1):
        initial_fibres[int(colors[position])].append(position)
    _, _, _, initial_rows = bounded_fibre_statistics(initial_fibres.values(), D)

    cursor = max(left_special) + 1
    left_limit = (N - D) // 2
    next_color = N + 1
    # Match the exact Apery h=2 law C_2=1+2*1_{(-51/p)=1}.  The base word
    # already has its single forced mirror edge.  One reflected private pair
    # adds exactly two further gap-2 edges when required.
    legendre = pow((-51) % p, (p - 1) // 2, p)
    assert legendre in (1, p - 1)
    h2_patch = 0
    if legendre == 1:
        for position in (cursor, cursor + 2, N + 1 - cursor, N - 1 - cursor):
            assert position not in special
            colors[position] = next_color
        next_color += 1
        cursor += 3
        h2_patch = 2
    used_gaps = []
    padded_pairs = 0
    for d in range(D, 3, -1):
        if cursor + 2 * d - 1 > left_limit:
            continue
        if initial_rows[d] + 2 * d > 3 * (d - 1):
            continue
        for offset in range(d):
            left = cursor + offset
            right = left + d
            reflected_left = N + 1 - right
            reflected_right = N + 1 - left
            assert len({left, right, reflected_left, reflected_right}) == 4
            for position in (left, right, reflected_left, reflected_right):
                assert position not in special
                colors[position] = next_color
            next_color += 1
            padded_pairs += 2
        used_gaps.append(d)
        cursor += 2 * d

    assert used_gaps
    fibres: DefaultDict[int, list[int]] = defaultdict(list)
    for position in range(1, N + 1):
        fibres[int(colors[position])].append(position)
    assert all(colors[position] == colors[N + 1 - position] for position in range(1, N + 1))
    assert all(colors[position] != colors[position + 1] for position in range(1, N))
    P, Q, S, rows = bounded_fibre_statistics(fibres.values(), D)
    assert all(rows[d] <= 3 * (d - 1) for d in range(1, D + 1))
    assert rows[2] == 1 + (2 if legendre == 1 else 0)
    assert Q == int(base["Q"])
    assert P == int(base["P"]) + padded_pairs + h2_patch
    window_checks = all_window_multiplicity_gate(fibres.values())
    ratio = Fraction(P * Q, N * N)
    return {
        "q": q,
        "p": p,
        "N": N,
        "D": D,
        "L_squared": Fraction(D * D, N),
        "blocks": len(used_gaps),
        "gap_min": min(used_gaps),
        "gap_max": max(used_gaps),
        "padded_pairs": padded_pairs,
        "h2_patch": h2_patch,
        "P": P,
        "Q": Q,
        "S": S,
        "ratio": ratio,
        "window_checks": window_checks,
    }


def fraction_decimal(value: Fraction) -> str:
    return f"{Decimal(value.numerator) / Decimal(value.denominator):.12f}"


def print_apery(result: dict[str, object], label: str) -> None:
    ratio = result["ratio"]
    assert isinstance(ratio, Fraction)
    print(
        f"  {label}: D={result['D']} S={result['S']} P={result['P']} "
        f"Q={result['Q']} P*Q/N^2={ratio.numerator}/{ratio.denominator} "
        f"={fraction_decimal(ratio)}"
    )
    print(f"    multiplicities={dict(sorted(result['multiplicity_histogram'].items()))}")
    print(f"    first-return deciles={result['first_quantiles']}")
    print(f"    first-return equal-width bins={result['first_bins']}")
    print(f"    intersection sums={result['intersection_sums']}")
    for shell in result["shells"]:
        row_sum = int(shell["sum"])
        excess = int(shell["excess"])
        saving = Fraction(excess, row_sum) if row_sum else Fraction(0)
        print(
            f"    shell ({shell['low']},{shell['high']}]: "
            f"sum={row_sum} union={shell['union']} excess={excess} "
            f"saving={fraction_decimal(saving)} pairI={shell['pair_intersections']} "
            f"first={shell['primitive']}"
        )
    print(
        f"    sum(shell unions)={result['shell_union_sum']} "
        f"cross-shell excess={result['cross_shell_excess']} "
        f"cascade checks={result['cascade_checks']}"
    )


def main() -> None:
    for p in APERY_PRIMES:
        cutoffs = requested_cutoffs(p)
        all_lags, checks = collision_lags(p, max(D for _, D in cutoffs))
        print(f"APERY p={p} N={p - 2} orbit/continuant checks={checks}")
        for label, D in cutoffs:
            print_apery(analyze_cutoff(p, D, all_lags), label)

    for q in WORD_QS:
        result = reflection_word_gate(q)
        ratio = result["ratio"]
        assert isinstance(ratio, Fraction)
        print(
            f"WORD q={q} p={result['p']} N={result['N']} H={result['H']} "
            f"S={result['S']} P={result['P']} Q={result['Q']} "
            f"P*Q={result['product']} P*Q/N^2={ratio.numerator}/{ratio.denominator} "
            f"={fraction_decimal(ratio)} rows={result['nonzero_rows']}"
        )

    for q in WORD_QS:
        result = padded_reflection_word_gate(q)
        ratio = result["ratio"]
        L_squared = result["L_squared"]
        assert isinstance(ratio, Fraction) and isinstance(L_squared, Fraction)
        print(
            f"PADDED_WORD q={q} p={result['p']} N={result['N']} D={result['D']} "
            f"D^2/N={L_squared.numerator}/{L_squared.denominator}="
            f"{fraction_decimal(L_squared)} blocks={result['blocks']} "
            f"gaps=[{result['gap_min']},{result['gap_max']}] "
            f"paddedP={result['padded_pairs']} h2patch={result['h2_patch']} "
            f"S={result['S']} P={result['P']} "
            f"Q={result['Q']} P*Q/N^2={ratio.numerator}/{ratio.denominator}="
            f"{fraction_decimal(ratio)} windows={result['window_checks']}"
        )

    digest = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    print(f"SHA256 {digest}")
    print("ALL STRIKE2 GATES PASS")


if __name__ == "__main__":
    main()
