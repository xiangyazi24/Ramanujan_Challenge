#!/usr/bin/env python3
"""Independently verify the first short projective Apéry phantom chains.

The prime-first C++ census groups indices by a projective solution state.  This
standard-library check takes its first ratio-above-one records as fixed input
and recomputes, independently:

* both recurrence solutions modulo ``p``;
* every raw gap value ``N_h(x)`` in the selected span;
* exclusion from the structural carrier ``U_s``;
* nonwrapping, noncentered, non-all-equal gaps; and
* failure of the start to lie in the distinguished Apéry zero fiber.

No resultant computation or output from the scanner is consumed at runtime.
"""

from __future__ import annotations

from collections import defaultdict


RECORDS = (
    (1297, 360, (5, 20, 10)),
    (1297, 901, (10, 20, 5)),
    (7411, 2655, (31, 40, 3)),
    (7411, 4681, (3, 40, 31)),
    (10427, 2286, (49, 12, 41)),
    (10427, 8038, (41, 12, 49)),
    (10993, 3061, (23, 63, 18)),
    (10993, 7827, (18, 63, 23)),
    (97553, 34550, (138, 58, 21)),
    (97553, 62785, (21, 58, 138)),
    (128047, 42375, (41, 86, 37)),
    (128047, 85507, (37, 86, 41)),
    (131591, 53366, (189, 81, 70)),
    (131591, 77884, (70, 81, 189)),
    (208759, 89452, (242, 63, 31)),
    (208759, 118970, (31, 63, 242)),
    (213023, 68949, (49, 55, 255)),
    (213023, 143714, (255, 55, 49)),
    (321341, 150730, (61, 31, 328)),
    (321341, 170190, (328, 31, 61)),
)


def apery_coefficient(index: int, prime: int) -> int:
    return (
        34 * index**3 + 51 * index**2 + 27 * index + 5
    ) % prime


def solution_pair(prime: int) -> tuple[list[int], list[int]]:
    """Return the distinguished and companion divided solutions modulo p."""

    apery = [1, 5 % prime]
    companion = [0, 1]
    for index in range(1, prime - 1):
        inverse_denominator = pow((index + 1) ** 3, -1, prime)
        coefficient = apery_coefficient(index, prime)
        apery.append(
            (
                coefficient * apery[index]
                - index**3 * apery[index - 1]
            )
            * inverse_denominator
            % prime
        )
        companion.append(
            (
                coefficient * companion[index]
                - index**3 * companion[index - 1]
            )
            * inverse_denominator
            % prime
        )
    return apery, companion


def return_offsets(prime: int, start: int, span: int) -> tuple[int, ...]:
    """Evaluate N_1(start),...,N_span(start) by the continuant recurrence."""

    previous, current = 0, 1
    offsets = [0]
    for height in range(1, span + 1):
        if current == 0:
            offsets.append(height)
        argument = (start + height) % prime
        following = (
            apery_coefficient(argument, prime) * current
            - pow(argument, 6, prime) * previous
        ) % prime
        previous, current = current, following
    return tuple(offsets)


def carrier_is_clean(prime: int, span: int, apery: list[int]) -> bool:
    """Check p does not divide product_{j<=s} j! b_j V_{j+1}."""

    assert span < prime  # Hence all factorial factors are units modulo p.
    previous, current = 0, 1  # V_0,V_1
    for index in range(1, span + 1):
        previous, current = current, (34 * current - previous) % prime
        if apery[index] == 0 or current == 0:  # b_index or V_(index+1)
            return False
    return True


def projective_state(
    prime: int, apery_value: int, companion_value: int
) -> int | None:
    assert apery_value or companion_value
    if companion_value == 0:
        return None
    return apery_value * pow(companion_value, -1, prime) % prime


def main() -> None:
    by_prime: dict[int, list[tuple[int, tuple[int, int, int]]]] = defaultdict(list)
    for prime, start, gaps in RECORDS:
        by_prime[prime].append((start, gaps))

    verified = 0
    maximum = (0, 1, None)
    for prime, prime_records in sorted(by_prime.items()):
        apery, companion = solution_pair(prime)
        for start, gaps in prime_records:
            a, b, c = gaps
            span = a + b + c
            selected = (0, a, a + b, span)
            assert start + span < prime
            assert not (a == b == c)
            assert carrier_is_clean(prime, span, apery)
            assert return_offsets(prime, start, span) == selected
            assert all(
                (2 * start + left + right + 1) % prime
                for left, right in zip(selected, selected[1:])
            )

            state = projective_state(
                prime, apery[start], companion[start]
            )
            assert state is not None and state != 0
            for offset in selected:
                assert projective_state(
                    prime,
                    apery[start + offset],
                    companion[start + offset],
                ) == state
            assert apery[start] != 0

            if prime * maximum[1] ** 2 > maximum[0] * span**2:
                maximum = (prime, span, (prime, start, gaps))
            verified += 1
            print(
                "PRIMITIVE_PROJECTIVE_PHANTOM PASS "
                f"p={prime} x={start} gaps={gaps} span={span} "
                f"p/span^2={prime}/{span**2}"
            )

        if len(prime_records) == 2:
            (left_start, left_gaps), (right_start, right_gaps) = prime_records
            assert right_gaps == tuple(reversed(left_gaps))
            assert right_start == prime - 1 - left_start - sum(left_gaps)

    print(
        "PRIMITIVE_PROJECTIVE_PRIME_VERIFY PASS "
        f"records={verified} max={maximum}"
    )


if __name__ == "__main__":
    main()
