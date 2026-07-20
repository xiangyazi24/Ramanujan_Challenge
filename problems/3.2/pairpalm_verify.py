#!/usr/bin/env python3
"""Exact checks for the higher pair--Palm audit.

Only the Python standard library is used.  The finite calculations verify the
combinatorial identities, the Apéry transcription, and the counterexamples in
``pairpalm_result.tex``.  They are not substituted for the missing cross-prime
estimate (PP)_7.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
from math import comb, isqrt, lcm


def falling(value: int, order: int) -> int:
    """Return the falling factorial (value)_order, with its usual zeroes."""

    if order < 0:
        raise ValueError("negative falling-factorial order")
    if value < order:
        return 0
    answer = 1
    for offset in range(order):
        answer *= value - offset
    return answer


def primes_up_to(limit: int) -> list[int]:
    sieve = bytearray(b"\x01") * (limit + 1)
    if limit >= 0:
        sieve[0] = 0
    if limit >= 1:
        sieve[1] = 0
    for prime in range(2, isqrt(limit) + 1):
        if sieve[prime]:
            start = prime * prime
            sieve[start : limit + 1 : prime] = b"\x00" * (
                (limit - start) // prime + 1
            )
    return [number for number in range(2, limit + 1) if sieve[number]]


def apery_zeros_division(prime: int) -> tuple[int, ...]:
    """Compute b_0,...,b_{p-1} modulo p from the divided recurrence."""

    values = [1, 5 % prime]
    for index in range(1, prime - 1):
        coefficient = (
            34 * index**3 + 51 * index**2 + 27 * index + 5
        ) % prime
        numerator = (
            coefficient * values[index] - index**3 * values[index - 1]
        ) % prime
        denominator = pow(index + 1, 3, prime)
        values.append(numerator * pow(denominator, -1, prime) % prime)
    return tuple(index for index, value in enumerate(values) if value == 0)


def apery_zeros_cleared(prime: int) -> tuple[int, ...]:
    """Compute the same zeros using y_n=(n!)^3 b_n (no division)."""

    previous, current = 1 % prime, 5 % prime
    zeros: list[int] = []
    if previous == 0:
        zeros.append(0)
    if current == 0:
        zeros.append(1)
    for index in range(1, prime - 1):
        coefficient = (
            34 * index**3 + 51 * index**2 + 27 * index + 5
        ) % prime
        following = (
            coefficient * current
            - pow(index, 6, prime) * previous
        ) % prime
        previous, current = current, following
        if current == 0:
            zeros.append(index + 1)
    return tuple(zeros)


def apery_integers(length: int) -> list[int]:
    if length <= 0:
        return []
    values = [1]
    if length == 1:
        return values
    values.append(5)
    for index in range(1, length - 1):
        coefficient = 34 * index**3 + 51 * index**2 + 27 * index + 5
        numerator = (
            coefficient * values[index] - index**3 * values[index - 1]
        )
        denominator = (index + 1) ** 3
        assert numerator % denominator == 0
        values.append(numerator // denominator)
    return values


def build_actual_incidence(
    scale: int,
) -> tuple[list[list[int]], list[Fraction], list[int], list[tuple[int, ...]]]:
    """Build rows m<X^2 and columns p in (X,2X] from actual zero sets."""

    primes = [p for p in primes_up_to(2 * scale) if p > scale]
    zero_sets = [apery_zeros_cleared(prime) for prime in primes]
    rows: list[list[int]] = [[] for _ in range(scale * scale)]
    for column, (prime, zeros) in enumerate(zip(primes, zero_sets)):
        for residue in zeros:
            for integer in range(residue, scale * scale, prime):
                rows[integer].append(column)
    weights = [
        Fraction(len(zeros), prime)
        for prime, zeros in zip(primes, zero_sets)
    ]
    return rows, weights, primes, zero_sets


def factorial_moment(rows: list[list[int]], order: int) -> int:
    return sum(falling(len(row), order) for row in rows)


def pair_palm(
    rows: list[list[int]], weights: list[Fraction], order: int
) -> tuple[int, int, Fraction, Fraction, tuple[int, int, int]]:
    """Return F_k, sum T, positive excess, baseline sum, sign counts."""

    if order < 3:
        raise ValueError("pair-Palm extension requires k >= 3")
    extension_order = order - 2
    columns = len(weights)
    intersections = [[0] * columns for _ in range(columns)]
    extensions = [[0] * columns for _ in range(columns)]

    for row in rows:
        row_size = len(row)
        extension_count = falling(row_size - 2, extension_order)
        for first in row:
            for second in row:
                if first == second:
                    continue
                intersections[first][second] += 1
                extensions[first][second] += extension_count

    lam = sum(weights, Fraction(0))
    excess = Fraction(0)
    baseline_sum = Fraction(0)
    signs = Counter()
    for first in range(columns):
        for second in range(columns):
            if first == second:
                continue
            local_lam = lam - weights[first] - weights[second]
            baseline = (
                intersections[first][second]
                * local_lam**extension_order
            )
            difference = Fraction(extensions[first][second]) - baseline
            baseline_sum += baseline
            excess += max(difference, Fraction(0))
            signs[1 if difference > 0 else -1 if difference < 0 else 0] += 1

    direct = factorial_moment(rows, order)
    extension_total = sum(map(sum, extensions))
    return (
        direct,
        extension_total,
        excess,
        baseline_sum,
        (signs[1], signs[-1], signs[0]),
    )


def ordered_tuple_mass(zero_counts: list[int], order: int) -> int:
    """Sum prod Z(p_i) over ordered, pairwise-distinct prime labels."""

    elementary = [0] * (order + 1)
    elementary[0] = 1
    for count in zero_counts:
        for degree in range(order, 0, -1):
            elementary[degree] += count * elementary[degree - 1]
    return falling(order, order) * elementary[order]


def crt_factorial_moments(
    primes: list[int],
    zero_sets: list[tuple[int, ...]],
    limit: int,
    max_order: int = 8,
) -> dict[int, int]:
    """Independently enumerate canonical tuples by their two least primes."""

    active = [
        (prime, frozenset(zeros))
        for prime, zeros in zip(primes, zero_sets)
        if zeros
    ]
    canonical = Counter()
    for left_index, (prime, left_zeros) in enumerate(active):
        for right_index in range(left_index + 1, len(active)):
            other_prime, right_zeros = active[right_index]
            assert prime * other_prime > limit
            inverse = pow(prime, -1, other_prime)
            for left_residue in left_zeros:
                for right_residue in right_zeros:
                    representative = left_residue + prime * (
                        (right_residue - left_residue)
                        * inverse
                        % other_prime
                    )
                    if representative >= limit:
                        continue
                    later_hits = sum(
                        representative % later_prime in later_zeros
                        for later_prime, later_zeros in active[right_index + 1 :]
                    )
                    for order in range(2, min(max_order, later_hits + 2) + 1):
                        canonical[order] += comb(later_hits, order - 2)
    return {
        order: falling(order, order) * canonical[order]
        for order in range(2, max_order + 1)
    }


def anchored_star(
    scale: int,
) -> tuple[list[list[int]], list[Fraction], list[int], int]:
    """The reflection-invariant anchored family from hm3_result.tex."""

    bound = scale * scale
    anchor = 3 * bound // 4
    obstruction = (
        (anchor - 1)
        * anchor
        * (anchor + 1)
        * (anchor + 2)
        * (2 * anchor + 1)
    )
    primes = [
        prime
        for prime in primes_up_to(2 * scale)
        if prime > scale and obstruction % prime != 0
    ]
    rows: list[list[int]] = [[] for _ in range(bound)]
    weights: list[Fraction] = []
    for column, prime in enumerate(primes):
        residue = anchor % prime
        zeros = {residue, prime - 1 - residue}
        assert len(zeros) == 2
        weights.append(Fraction(2, prime))
        for zero in zeros:
            for integer in range(zero, bound, prime):
                rows[integer].append(column)
    return rows, weights, primes, anchor


def check_apery_transcription() -> str:
    expected_prefix = [
        1,
        5,
        73,
        1445,
        33001,
        819005,
        21460825,
        584307365,
        16367912425,
    ]
    assert apery_integers(len(expected_prefix)) == expected_prefix

    expected_zeros = {
        7: (),
        11: (5,),
        13: (),
        17: (3, 13),
        19: (8, 10),
        37: (17, 19),
        73: (2, 70),
        251: (114, 136),
    }
    for prime in primes_up_to(257):
        if prime < 5:
            continue
        divided = apery_zeros_division(prime)
        cleared = apery_zeros_cleared(prime)
        assert divided == cleared
        assert all((prime - 1 - residue) in cleared for residue in cleared)
        if prime in expected_zeros:
            assert cleared == expected_zeros[prime]
    return "Apéry recurrence, fixed zero sets, and reflection agree"


def check_anchor_identity_and_multiplicity() -> str:
    for row_size in range(0, 25):
        for order in range(3, 9):
            assert (
                falling(row_size, 2)
                * falling(row_size - 2, order - 2)
                == falling(row_size, order)
            )
            overlap = sum(
                comb(order - 2, shared) ** 2
                * falling(shared, shared)
                * falling(row_size, 2 * order - 2 - shared)
                for shared in range(order - 1)
            )
            assert falling(row_size, order) ** 2 == (
                falling(row_size, 2) * overlap
            )

    # An existential third-prime indicator is not the required multiplicity:
    # at K=4 it gives (4)_2=12, whereas the ordered third moment is 24.
    assert falling(4, 2) == 12
    assert falling(4, 3) == 24
    return "anchor and overlap identities hold; existential T fails at K=4"


def check_synthetic_pair_palm() -> str:
    rows = [
        list(range(8)),
        [8, 9],
        [0, 8],
        [1, 9],
        [2, 3, 8],
        [4, 5, 9],
    ] + [[] for _ in range(6)]
    weights = [Fraction(index + 1, 100) for index in range(10)]
    lam = sum(weights, Fraction(0))
    f2 = factorial_moment(rows, 2)
    expected = [348, 1680, 6720, 20160, 40320, 40320]

    for order, expected_moment in zip(range(3, 9), expected):
        direct, extension_total, excess, baseline, signs = pair_palm(
            rows, weights, order
        )
        assert direct == extension_total == expected_moment
        assert excess <= direct <= excess + baseline
        assert baseline <= lam ** (order - 2) * f2
        assert excess < direct < excess + lam ** (order - 2) * f2
        if order == 3:
            assert signs == (64, 6, 20)
        else:
            assert signs == (56, 14, 20)

        overlap_moment = sum(
            comb(order - 2, shared) ** 2
            * falling(shared, shared)
            * factorial_moment(rows, 2 * order - 2 - shared)
            for shared in range(order - 1)
        )
        assert direct**2 <= f2 * overlap_moment
    return "nondegenerate synthetic pair-Palm identities pass for k=3,...,8"


def check_actual_apery_block() -> str:
    scale = 128
    rows, weights, primes, zero_sets = build_actual_incidence(scale)
    lam = sum(weights, Fraction(0))
    histogram = Counter(map(len, rows))
    zero_counts = [len(zeros) for zeros in zero_sets]

    assert len(primes) == 23
    assert sum(count > 0 for count in zero_counts) == 13
    assert sum(zero_counts) == 30
    assert histogram == Counter({0: 13976, 1: 2242, 2: 156, 3: 10})
    assert factorial_moment(rows, 2) == 372
    assert factorial_moment(rows, 3) == 60
    assert all(factorial_moment(rows, order) == 0 for order in range(4, 9))
    assert lam == Fraction(
        85675682425804117942566386934,
        543404033598502288453922853793,
    )

    # The sharper unconditional second-moment constant is 4, since pq>X^2.
    bound = scale * scale
    f2 = factorial_moment(rows, 2)
    assert Fraction(f2) <= 4 * bound * lam**2

    finite_weights = []
    for prime, zeros in zip(primes, zero_sets):
        representatives = sum(
            (bound - 1 - residue) // prime + 1 for residue in zeros
        )
        finite_weights.append(Fraction(representatives, bound))
    mu = sum(finite_weights, Fraction(0))
    assert abs(mu - lam) <= Fraction(2, scale) * lam

    expected_tuple_bounds = [
        824,
        20592,
        464640,
        9377280,
        167270400,
        2597253120,
        34402959360,
    ]
    for order, expected_tuple_bound in zip(range(2, 9), expected_tuple_bounds):
        moment = factorial_moment(rows, order)
        tuple_bound = ordered_tuple_mass(zero_counts, order)
        assert tuple_bound == expected_tuple_bound
        assert moment <= tuple_bound <= sum(zero_counts) ** order
        assert sum(zero_counts) ** order <= (2 * scale * lam) ** order

    for order in range(3, 8):
        direct, extension_total, excess, baseline, _ = pair_palm(
            rows, weights, order
        )
        assert direct == extension_total
        assert excess <= direct <= excess + baseline
        assert baseline <= lam ** (order - 2) * f2
        assert direct <= excess + 4 * bound * lam**order

        finite = pair_palm(rows, finite_weights, order)
        assert finite[0] == finite[1] == direct
        assert finite[2] <= direct <= finite[2] + finite[3]
        assert finite[3] <= mu ** (order - 2) * f2
        assert direct <= finite[2] + 4 * (
            1 + Fraction(2, scale)
        ) ** (order - 2) * bound * lam**order
    return "actual X=128 block matches fixed moments and all Palm sandwiches"


def check_actual_apery_1024() -> str:
    """Reproduce the nonzero fourth-moment block used in the TeX audit."""

    scale = 1024
    bound = scale * scale
    rows, weights, primes, zero_sets = build_actual_incidence(scale)
    for prime, cleared in zip(primes, zero_sets):
        assert apery_zeros_division(prime) == cleared
        zero_set = set(cleared)
        assert all(prime - 1 - residue in zero_set for residue in cleared)
        assert all((residue + 1) % prime not in zero_set for residue in cleared)

    histogram = Counter(map(len, rows))
    assert histogram == Counter({0: 965700, 1: 79588, 2: 3210, 3: 77, 4: 1})
    moments = {order: factorial_moment(rows, order) for order in range(1, 9)}
    assert [moments[order] for order in range(1, 9)] == [
        86243,
        6894,
        486,
        24,
        0,
        0,
        0,
        0,
    ]
    existence_aggregate = sum(
        falling(len(row), 2) for row in rows if len(row) >= 3
    )
    assert existence_aggregate == 474 < moments[3] == 486

    crt_moments = crt_factorial_moments(primes, zero_sets, bound)
    assert all(crt_moments[order] == moments[order] for order in range(2, 9))

    lam = sum(weights, Fraction(0))
    assert f"{float(lam):.9f}" == "0.082244154"
    assert Fraction(moments[2]) <= 4 * bound * lam**2
    zero_counts = list(map(len, zero_sets))
    for order in range(2, 9):
        assert moments[order] <= ordered_tuple_mass(zero_counts, order)
        assert ordered_tuple_mass(zero_counts, order) <= (2 * scale * lam) ** order

    finite_weights = []
    for prime, zeros in zip(primes, zero_sets):
        representatives = sum(
            (bound - 1 - residue) // prime + 1 for residue in zeros
        )
        finite_weights.append(Fraction(representatives, bound))
    mu = sum(finite_weights, Fraction(0))
    assert abs(mu - lam) <= Fraction(2, scale) * lam

    periodic_four = pair_palm(rows, weights, 4)
    periodic_seven = pair_palm(rows, weights, 7)
    finite_four = pair_palm(rows, finite_weights, 4)
    assert periodic_four[0] == periodic_four[1] == moments[4]
    assert periodic_seven[0] == periodic_seven[1] == moments[7] == 0
    assert f"{float(periodic_four[2]):.9f}" == "23.448043852"
    assert f"{float(finite_four[2]):.9f}" == "23.448002897"
    for order in range(3, 9):
        direct, extension_total, excess, baseline, _ = pair_palm(
            rows, weights, order
        )
        assert direct == extension_total == moments[order]
        assert excess <= direct <= excess + baseline
        assert baseline <= lam ** (order - 2) * moments[2]

        finite = pair_palm(rows, finite_weights, order)
        assert finite[0] == finite[1] == direct
        assert finite[2] <= direct <= finite[2] + finite[3]
        assert finite[3] <= mu ** (order - 2) * moments[2]

    return (
        "actual X=1024 has F2..F8=(6894,486,24,0,0,0,0); "
        "CRT, multiplicity, and Palm checks agree"
    )


def check_cs_and_six_tuple_claim() -> str:
    f2 = falling(6, 2)
    f4 = falling(6, 4)
    f6 = falling(6, 6)
    assert (f2, f4, f6) == (30, 360, 720)
    assert f4**2 == 6 * f2 * f6
    assert falling(4, 4) > 0 == falling(4, 6)

    rows, weights, primes, anchor = anchored_star(64)
    lam = sum(weights, Fraction(0))
    bound = 64 * 64
    assert len(primes) == 12
    assert len(rows[anchor]) == 12
    assert factorial_moment(rows, 6) > bound * lam**6

    common_integer = 1
    for value in range(1, 65):
        common_integer = lcm(common_integer, value)
    assert all(common_integer % prime == 0 for prime in primes_up_to(64))
    assert common_integer.bit_length() < 128
    return (
        "Cauchy--Schwarz/F6 claims have exact counterexamples; "
        "height alone is non-distributional"
    )


def check_unconditional_envelope_and_star() -> str:
    rows, weights, primes, anchor = anchored_star(64)
    lam = sum(weights, Fraction(0))
    bound = 64 * 64
    active_columns = len(primes)
    f2 = factorial_moment(rows, 2)

    assert primes == [
        67,
        71,
        73,
        79,
        89,
        97,
        101,
        103,
        107,
        109,
        113,
        127,
    ]
    assert lam == Fraction(
        108674385686125190956832,
        412378374392817307962353,
    )
    assert f2 == 290
    assert len(rows[anchor]) == active_columns

    expected = [290, 1332, 11880, 95040, 665280, 3991680, 19958400]
    for order, expected_moment in zip(range(2, 9), expected):
        moment = factorial_moment(rows, order)
        assert moment == expected_moment
        assert moment >= falling(active_columns, order)
        assert moment <= falling(active_columns - 2, order - 2) * f2
        assert moment <= (2 * 64 * lam) ** order

    for order in range(3, 9):
        direct, _, excess, _, _ = pair_palm(rows, weights, order)
        assert excess <= direct <= excess + 4 * bound * lam**order
    return "anchored star attains the high-moment envelope through k=8"


def main() -> None:
    checks = [
        check_apery_transcription,
        check_anchor_identity_and_multiplicity,
        check_synthetic_pair_palm,
        check_actual_apery_block,
        check_actual_apery_1024,
        check_cs_and_six_tuple_claim,
        check_unconditional_envelope_and_star,
    ]
    for check in checks:
        print(f"PASS: {check()}")
    print(f"ALL {len(checks)} CHECKS PASS")


if __name__ == "__main__":
    main()
