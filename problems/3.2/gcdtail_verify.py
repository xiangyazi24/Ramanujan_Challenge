#!/usr/bin/env python3
"""Finite checks for ``gcdtail_result.tex``.

For every prime ``p <= 2000`` and every
``4 <= H <= floor(sqrt(p))``, this pure-Python script checks the exact
reflection, collision, column-pair, and deep-pair identities used in the
TeX audit.  It also checks three isolated examples used there: a
resultant-component mismatch in the scanned range, the same mismatch in
the deep range, and a deep polluted-boundary collision.

The finite scan is not used as evidence for an asymptotic estimate.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from hashlib import sha256
from math import comb, gcd, isqrt
import sys


PRIME_LIMIT = 2000
EXPECTED_DIGEST = "76b9e86cf7f8f007826efba856bf2728875eff0feddf10afd947cec050941a84"


def primes_up_to(limit: int) -> list[int]:
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[:2] = b"\x00\x00"
    for divisor in range(2, isqrt(limit) + 1):
        if not sieve[divisor]:
            continue
        start = divisor * divisor
        sieve[start : limit + 1 : divisor] = b"\x00" * (
            (limit - start) // divisor + 1
        )
    return [number for number, flag in enumerate(sieve) if flag]


def deep_cutoff(height: int) -> int:
    """Return floor(H^(3/4)) without floating-point arithmetic."""

    low, high = 0, height + 1
    target = height**3
    while high - low > 1:
        middle = (low + high) // 2
        if middle**4 <= target:
            low = middle
        else:
            high = middle
    assert low**4 <= target < (low + 1) ** 4
    return low


def apery_coefficient(argument: int) -> int:
    return (2 * argument + 1) * (
        17 * argument * argument + 17 * argument + 5
    )


def apery_values_mod(prime: int, max_index: int) -> list[int]:
    values = [1 % prime, 5 % prime]
    for index in range(1, max_index):
        numerator = (
            apery_coefficient(index) * values[index]
            - index**3 * values[index - 1]
        ) % prime
        denominator = (index + 1) ** 3 % prime
        values.append(numerator * pow(denominator, -1, prime) % prime)
    return values[: max_index + 1]


def root_data(
    prime: int, height: int
) -> tuple[list[list[int]], list[bytearray], list[list[int]]]:
    """Evaluate N_0,...,N_H pointwise over F_p."""

    roots = [[] for _ in range(height + 1)]
    masks = [bytearray(prime) for _ in range(height + 1)]
    levels_by_column = [[] for _ in range(prime)]
    if height < 2:
        return roots, masks, levels_by_column

    p_values = [apery_coefficient(value) % prime for value in range(prime)]
    sixth_powers = [pow(value, 6, prime) for value in range(prime)]
    previous = [0] * prime
    current = [1] * prime
    for index in range(1, height):
        following = [
            (
                p_values[(x + index) % prime] * current[x]
                - sixth_powers[(x + index) % prime] * previous[x]
            )
            % prime
            for x in range(prime)
        ]
        level = index + 1

        # N_h(-X-h-1)=(-1)^(h-1)N_h(X).
        sign = 1 if (level - 1) % 2 == 0 else -1
        assert all(
            following[(-x - level - 1) % prime]
            == sign * following[x] % prime
            for x in range(prime)
        )

        level_roots = [x for x, value in enumerate(following) if value == 0]
        if prime >= 7:
            assert len(level_roots) <= 3 * (level - 1)
        roots[level] = level_roots
        for x in level_roots:
            masks[level][x] = 1
            levels_by_column[x].append(level)
        previous, current = current, following
    return roots, masks, levels_by_column


def records_for_prime(prime: int) -> list[tuple[int, ...]]:
    max_height = isqrt(prime)
    roots, masks, levels_by_column = root_data(prime, max_height)
    apery = apery_values_mod(prime, max(1, max_height - 1))

    pair_witnesses: dict[tuple[int, int], set[int]] = {}
    witnesses: list[tuple[int, int, int]] = []
    for first_gap in range(2, max(2, max_height - 1)):
        for second_gap in range(2, max_height - first_gap + 1):
            direct = {
                x
                for x in roots[first_gap]
                if masks[second_gap][(x + first_gap) % prime]
            }

            # With y=x+d, reflection changes the shifted collision into
            # N_d(-y-1)=N_r(y)=0.
            dual = {
                y
                for y in roots[second_gap]
                if masks[first_gap][(-y - 1) % prime]
            }
            assert {(x + first_gap) % prime for x in direct} == dual
            assert len(direct) <= 3 * (min(first_gap, second_gap) - 1)

            pair_witnesses[first_gap, second_gap] = direct
            for x in direct:
                assert masks[first_gap + second_gap][x]
                witnesses.append((first_gap, second_gap, x))

    # The exact involution (d,r,x) -> (r,d,-x-d-r-1).
    for first_gap, second_gap, x in witnesses:
        reflected = (-x - first_gap - second_gap - 1) % prime
        assert reflected in pair_witnesses[second_gap, first_gap]
        assert (
            -reflected - first_gap - second_gap - 1
        ) % prime == x
    for first_gap, second_gap in pair_witnesses:
        assert len(pair_witnesses[first_gap, second_gap]) == len(
            pair_witnesses[second_gap, first_gap]
        )

    zero_pairs: list[tuple[int, int, int]] = []
    for x, levels in enumerate(levels_by_column):
        for index, lower in enumerate(levels):
            for upper in levels[index + 1 :]:
                zero_pairs.append((lower, upper, x))

    records = []
    for height in range(4, max_height + 1):
        cutoff = deep_cutoff(height)
        boundary = {(-endpoint) % prime for endpoint in range(2, height + 1)}
        polluted = {
            (-cut) % prime
            for cut in range(2, height + 1)
            if apery[cut - 1] == 0
        }
        assert all((-x - 1) % prime not in polluted for x in polluted)
        current_witnesses = [
            item for item in witnesses if item[0] + item[1] <= height
        ]
        current_zero_pairs = [item for item in zero_pairs if item[1] <= height]

        off_witnesses = [
            item for item in current_witnesses if item[2] not in boundary
        ]
        off_zero_pairs = [
            item for item in current_zero_pairs if item[2] not in boundary
        ]

        # Exact off-boundary identity and the global boundary sandwich.
        assert {
            (first_gap, first_gap + second_gap, x)
            for first_gap, second_gap, x in off_witnesses
        } == set(off_zero_pairs)
        assert len(current_witnesses) <= len(current_zero_pairs)

        deep = [
            item
            for item in current_witnesses
            if min(item[0], item[1]) > cutoff
        ]
        deep_off = [item for item in deep if item[2] not in boundary]
        deep_column_pairs = [
            (lower, upper, x)
            for lower, upper, x in off_zero_pairs
            if lower > cutoff and upper - lower > cutoff
        ]
        assert {
            (first_gap, first_gap + second_gap, x)
            for first_gap, second_gap, x in deep_off
        } == set(deep_column_pairs)

        pair_multiplicity = Counter(
            (first_gap, second_gap) for first_gap, second_gap, _ in deep
        )
        by_column = Counter(x for _, _, x in deep)
        by_first_gap = Counter(first_gap for first_gap, _, _ in deep)
        fiber: defaultdict[tuple[int, int], set[int]] = defaultdict(set)
        for first_gap, second_gap, x in deep:
            fiber[first_gap, x].add(second_gap)

        # The aligned three-polynomial gcd identity, evaluated through
        # its split roots: both sides count pairs of r-levels meeting the
        # same root of N_d.
        aligned_first = sum(len(values) for values in fiber.values())
        aligned_second = sum(comb(len(values), 2) for values in fiber.values())
        aligned_pairwise = 0
        for first_gap in range(cutoff + 1, height - 1):
            returns = list(range(cutoff + 1, height - first_gap + 1))
            for index, second_gap in enumerate(returns):
                for third_gap in returns[index + 1 :]:
                    aligned_pairwise += len(
                        pair_witnesses[first_gap, second_gap]
                        & pair_witnesses[first_gap, third_gap]
                    )
        assert aligned_first == len(deep)
        assert aligned_second == aligned_pairwise

        records.append(
            (
                prime,
                height,
                cutoff,
                len(current_witnesses),
                len(current_zero_pairs),
                len(deep),
                sum(x in boundary for _, _, x in deep),
                sum(first_gap == second_gap for first_gap, second_gap, _ in deep),
                len(pair_multiplicity),
                max(pair_multiplicity.values(), default=0),
                max(by_column.values(), default=0),
                max(by_first_gap.values(), default=0),
                max((len(values) for values in fiber.values()), default=0),
            )
        )
    return records


def trim(polynomial: list[int], prime: int) -> list[int]:
    polynomial = [coefficient % prime for coefficient in polynomial]
    while len(polynomial) > 1 and polynomial[-1] == 0:
        polynomial.pop()
    return polynomial


def poly_add(
    left: list[int], right: list[int], prime: int, sign: int = 1
) -> list[int]:
    result = [0] * max(len(left), len(right))
    for index, value in enumerate(left):
        result[index] += value
    for index, value in enumerate(right):
        result[index] += sign * value
    return trim(result, prime)


def poly_mul(left: list[int], right: list[int], prime: int) -> list[int]:
    result = [0] * (len(left) + len(right) - 1)
    for i, first in enumerate(left):
        for j, second in enumerate(right):
            result[i + j] += first * second
    return trim(result, prime)


def poly_shift(polynomial: list[int], shift: int, prime: int) -> list[int]:
    """Return f(X+shift) over F_p by Horner's rule."""

    result = [0]
    for coefficient in reversed(polynomial):
        result = poly_mul(result, [shift, 1], prime)
        result[0] = (result[0] + coefficient) % prime
    return trim(result, prime)


def poly_pow(polynomial: list[int], exponent: int, prime: int) -> list[int]:
    result = [1]
    base = polynomial
    while exponent:
        if exponent & 1:
            result = poly_mul(result, base, prime)
        base = poly_mul(base, base, prime)
        exponent >>= 1
    return result


def poly_gcd(left: list[int], right: list[int], prime: int) -> list[int]:
    left, right = trim(left, prime), trim(right, prime)
    while right != [0]:
        remainder = left[:]
        inverse_lead = pow(right[-1], -1, prime)
        while remainder != [0] and len(remainder) >= len(right):
            offset = len(remainder) - len(right)
            factor = remainder[-1] * inverse_lead % prime
            for index, value in enumerate(right):
                remainder[index + offset] -= factor * value
            remainder = trim(remainder, prime)
        left, right = right, remainder
    inverse_lead = pow(left[-1], -1, prime)
    return trim([inverse_lead * value for value in left], prime)


def poly_eval(polynomial: list[int], value: int, prime: int) -> int:
    result = 0
    for coefficient in reversed(polynomial):
        result = (result * value + coefficient) % prime
    return result


def gap_polynomials_mod(prime: int, max_height: int) -> list[list[int]]:
    polynomials = [[0], [1]]
    p_base = [5, 27, 51, 34]
    for height in range(1, max_height):
        following = poly_add(
            poly_mul(poly_shift(p_base, height, prime), polynomials[height], prime),
            poly_mul(poly_pow([height, 1], 6, prime), polynomials[height - 1], prime),
            prime,
            sign=-1,
        )
        polynomials.append(following)
    return polynomials


def center_values(even_block: int, max_index: int) -> list[int]:
    """Generate T_b^(a) from the exact integer center recurrence."""

    assert even_block >= 2 and even_block % 2 == 0
    values = [0, 1]
    for index in range(1, max_index):
        z = even_block + 2 * index - 1
        coefficient = 34 * z**3 + 102 * z**2 + 108 * z + 40
        values.append(coefficient * values[index] - z**6 * values[index - 1])
    return values


def check_resultant_component_obstructions() -> None:
    # A small mismatch inside the main p<=2000 scan.
    prime, first_gap = 257, 4
    polynomials = gap_polynomials_mod(prime, 12)
    ninth = poly_shift(polynomials[9], first_gap, prime)
    eleventh = poly_shift(polynomials[11], first_gap, prime)
    assert [polynomials[first_gap][-1], ninth[-1], eleventh[-1]] == [172, 91, 67]
    assert poly_gcd(polynomials[first_gap], ninth, prime) == [217, 1]
    assert poly_gcd(polynomials[first_gap], eleventh, prime) == [215, 1]
    assert poly_gcd(ninth, eleventh, prime) == [1]
    assert poly_eval(ninth, 40, prime) == 0
    assert poly_eval(eleventh, 42, prime) == 0
    assert [
        second_gap
        for second_gap in range(2, 13)
        if poly_eval(
            poly_shift(polynomials[second_gap], first_gap, prime), 40, prime
        )
        == 0
    ] == [9]
    assert [
        second_gap
        for second_gap in range(2, 13)
        if poly_eval(
            poly_shift(polynomials[second_gap], first_gap, prime), 42, prime
        )
        == 0
    ] == [11]

    # The same component mismatch wholly inside the requested deep region.
    prime, height, first_gap = 5683, 72, 33
    assert all(prime % divisor for divisor in range(2, isqrt(prime) + 1))
    assert height**2 < prime
    assert 29**4 > height**3
    polynomials = gap_polynomials_mod(prime, 39)
    twenty_ninth = poly_shift(polynomials[29], first_gap, prime)
    thirty_ninth = poly_shift(polynomials[39], first_gap, prime)
    assert [
        polynomials[first_gap][-1],
        twenty_ninth[-1],
        thirty_ninth[-1],
    ] == [4692, 1160, 3467]
    assert poly_gcd(polynomials[first_gap], twenty_ninth, prime) == [4912, 1]
    assert poly_gcd(polynomials[first_gap], thirty_ninth, prime) == [1497, 1]
    assert poly_gcd(twenty_ninth, thirty_ninth, prime) == [1]
    assert poly_eval(polynomials[first_gap], 771, prime) == 0
    assert poly_eval(twenty_ninth, 771, prime) == 0
    assert poly_eval(polynomials[first_gap], 4186, prime) == 0
    assert poly_eval(thirty_ninth, 4186, prime) == 0
    assert [
        second_gap
        for second_gap in range(deep_cutoff(height) + 1, height - first_gap + 1)
        if poly_eval(
            poly_shift(polynomials[second_gap], first_gap, prime), 771, prime
        )
        == 0
    ] == [29]
    assert [
        second_gap
        for second_gap in range(deep_cutoff(height) + 1, height - first_gap + 1)
        if poly_eval(
            poly_shift(polynomials[second_gap], first_gap, prime), 4186, prime
        )
        == 0
    ] == [39]

    # The center values are not a fixed-first-index divisibility chain.
    centers = center_values(2, 3)
    center_second = centers[2]
    center_third = centers[3]
    assert (center_second, center_third) == (2200, 16220375)
    assert gcd(center_second, center_third) == 25
    assert center_second % center_third and center_third % center_second


def check_deep_boundary_example() -> None:
    prime, height = 4283, 65
    assert all(prime % divisor for divisor in range(2, isqrt(prime) + 1))
    roots, masks, levels_by_column = root_data(prime, height)
    assert deep_cutoff(height) == 22
    assert apery_values_mod(prime, 18)[18] == 0
    assert levels_by_column[(-19) % prime] == list(range(19, height + 1))
    assert levels_by_column[(-39) % prime] == [24, 57]
    for first_gap, second_gap, witness in ((24, 33, -39), (33, 24, -19)):
        x = witness % prime
        assert min(first_gap, second_gap) > deep_cutoff(height)
        assert x in roots[first_gap]
        assert masks[second_gap][(x + first_gap) % prime]


def run() -> None:
    if not __debug__:
        raise RuntimeError("run without -O; assertions are verification checks")

    records: list[tuple[int, ...]] = []
    primes = primes_up_to(PRIME_LIMIT)
    for prime in primes:
        records.extend(records_for_prime(prime))

    payload = "\n".join(",".join(map(str, record)) for record in records)
    digest = sha256(payload.encode("ascii")).hexdigest()
    assert digest == EXPECTED_DIGEST
    assert len(primes) == 303
    assert len(records) == 7439
    assert tuple(
        sum(record[index] for record in records)
        for index in range(2, len(records[0]))
    ) == (60645, 2349, 3832, 57, 0, 3, 57, 30, 37, 30, 30)
    assert tuple(
        max(record[index] for record in records)
        for index in range(3, len(records[0]))
    ) == (16, 173, 2, 0, 1, 2, 1, 2, 1, 1)

    endpoint_records = [
        record for record in records if record[1] == isqrt(record[0])
    ]
    assert sum(record[5] for record in endpoint_records) == 12
    assert [record[0] for record in endpoint_records if record[5]] == [
        547,
        607,
        691,
        1109,
        1409,
        1721,
    ]

    check_resultant_component_obstructions()
    check_deep_boundary_example()

    print("GCDTAIL IDENTITIES PASS")
    print(f"primes<={PRIME_LIMIT}: {len(primes)}; height records: {len(records)}")
    print(f"canonical_sha256={digest}")
    print("deep records: total energy 57, maximum 2; endpoint total 12")
    print("RESULTANT COMPONENT OBSTRUCTIONS PASS")
    print("DEEP BOUNDARY EXAMPLE PASS")
    print("GCDTAIL_VERIFY PASS")


if __name__ == "__main__":
    try:
        run()
    except Exception:
        print("GCDTAIL_VERIFY FAIL", file=sys.stderr)
        raise
