#!/usr/bin/env python3
"""Verify the universal zero-to-zero companion transport law.

Let ``u`` and ``v`` be independent solutions of the divided Apéry
recurrence.  If ``u_r = u_{r+h} = 0`` in a nonwrapping interval, then

    v_{r+h} / v_r = -prod_{j=1}^{h-1}(r+j)^3
                      / ((r+h)^3 N_{h-1}(r)).

The right side is independent of the projective fibre cut out by ``u``.
Consequently the marked companion coordinate at an actual Apéry zero is a
universal return multiplier, not distinguished-orbit data.  This script
checks the formula on actual zeros and on the first short projective phantom
chains, checks that products along a four-return chain give only the direct
transport formula, and records a carrier-clean short phantom for which the
marked coordinate itself has a nonreflection collision.
"""

from __future__ import annotations


def middle(index: int, prime: int) -> int:
    return (
        34 * index**3 + 51 * index**2 + 27 * index + 5
    ) % prime


def solution_pair(prime: int) -> tuple[list[int], list[int]]:
    """Return b and c=a/6 modulo prime on the nonsingular base block."""

    apery = [1, 5 % prime]
    companion = [0, 1]
    for index in range(1, prime - 1):
        denominator_inverse = pow((index + 1) ** 3, -1, prime)
        coefficient = middle(index, prime)
        apery.append(
            (
                coefficient * apery[index]
                - index**3 * apery[index - 1]
            )
            * denominator_inverse
            % prime
        )
        companion.append(
            (
                coefficient * companion[index]
                - index**3 * companion[index - 1]
            )
            * denominator_inverse
            % prime
        )
    return apery, companion


def gap_values(prime: int, start: int, height: int) -> list[int]:
    """Return N_0(start),...,N_height(start) modulo prime."""

    values = [0, 1]
    for gap in range(1, height):
        argument = start + gap
        values.append(
            (
                middle(argument, prime) * values[gap]
                - pow(argument, 6, prime) * values[gap - 1]
            )
            % prime
        )
    return values[: height + 1]


def transport_multiplier(prime: int, start: int, gap: int) -> int:
    """Return the universal multiplier v_(start+gap)/v_start."""

    assert gap >= 2 and start + gap < prime
    continuants = gap_values(prime, start, gap)
    assert continuants[gap] == 0
    assert continuants[gap - 1] != 0
    numerator = 1
    for argument in range(start + 1, start + gap):
        numerator = numerator * pow(argument, 3, prime) % prime
    denominator = (
        pow(start + gap, 3, prime) * continuants[gap - 1]
    ) % prime
    return -numerator * pow(denominator, -1, prime) % prime


def verify_determinant_formula(
    prime: int,
    start: int,
    gap: int,
    apery: list[int],
    companion: list[int],
) -> None:
    """Verify the exact marked determinant formula modulo prime."""

    denominator_product = 1
    for argument in range(start + 2, start + gap + 1):
        denominator_product = (
            denominator_product * pow(argument, 3, prime)
        ) % prime
    endpoint_determinant = (
        apery[start + gap] * companion[start]
        - companion[start + gap] * apery[start]
    ) % prime
    expected = (
        -pow(start + 1, 3, prime)
        * denominator_product
        * endpoint_determinant
    ) % prime
    assert gap_values(prime, start, gap)[gap] == expected


def carrier_is_clean(prime: int, span: int, apery: list[int]) -> bool:
    """Check p does not divide prod_{j<=span} j! b_j V_(j+1)."""

    previous, current = 0, 1  # V_0,V_1
    for index in range(1, span + 1):
        previous, current = current, (34 * current - previous) % prime
        if apery[index] == 0 or current == 0:
            return False
    return True


def verify_marked_rows(
    prime: int,
    indices: tuple[int, ...],
    apery: list[int],
    companion: list[int],
) -> None:
    """Check R_r=c_r L and transverse q-derivatives at Apéry zeros.

    Relative to the standard initial coordinates (u_0,u_1), every solution
    is u = u_0 b + (u_1-5u_0)c.  Hence the evaluation row at index r is
    R_r=(b_r-5c_r,c_r).  At an Apéry zero it equals c_r(-5,1).
    """

    for index in indices:
        assert apery[index] == 0
        assert companion[index] != 0
        evaluation_row = (
            (apery[index] - 5 * companion[index]) % prime,
            companion[index],
        )
        assert evaluation_row == (
            (-5 * companion[index]) % prime,
            companion[index],
        )
        # d/dq (b_r-qc_r) at q=0 is nonzero.
        assert (-companion[index]) % prime != 0


def verify_chain(
    prime: int,
    indices: tuple[int, ...],
    expected_state: int,
) -> tuple[int, ...]:
    """Verify one projective fibre and all transport/product identities."""

    apery, companion = solution_pair(prime)
    state = apery[indices[0]] * pow(companion[indices[0]], -1, prime) % prime
    assert state == expected_state
    labels = []
    for index in indices:
        assert companion[index] != 0
        assert (apery[index] - state * companion[index]) % prime == 0
        labels.append(companion[index])

    product = 1
    for left, right, left_label, right_label in zip(
        indices, indices[1:], labels, labels[1:]
    ):
        verify_determinant_formula(
            prime, left, right - left, apery, companion
        )
        multiplier = transport_multiplier(prime, left, right - left)
        assert right_label == multiplier * left_label % prime
        product = product * multiplier % prime

    verify_determinant_formula(
        prime, indices[0], indices[-1] - indices[0], apery, companion
    )
    direct = transport_multiplier(prime, indices[0], indices[-1] - indices[0])
    assert product == direct
    assert labels[-1] == direct * labels[0] % prime
    return tuple(labels)


def main() -> None:
    # The first actual off-center consecutive four-zero window.  Its large
    # span is irrelevant here: it checks the same transport identity on the
    # distinguished state q=0.
    actual = (99, 868, 1011, 1294)
    actual_labels = verify_chain(3727, actual, 0)
    actual_apery, actual_companion = solution_pair(3727)
    verify_marked_rows(
        3727, actual, actual_apery, actual_companion
    )
    print(
        "ACTUAL_TRANSPORT PASS "
        f"p=3727 indices={actual} labels={actual_labels}"
    )
    print("MARKED_ROW_TRANSVERSALITY PASS p=3727")

    # Short formal four-return chains outside the distinguished fibre.
    phantom_records = (
        (1297, (360, 365, 385, 395), 454),
        (7411, (4681, 4684, 4724, 4755), 2717),
        (128047, (42375, 42416, 42502, 42539), 90334),
    )
    for prime, indices, state in phantom_records:
        labels = verify_chain(prime, indices, state)
        assert state != 0
        reflected = tuple(prime - 1 - index for index in reversed(indices))
        reflected_labels = verify_chain(prime, reflected, state)
        assert reflected_labels == tuple(reversed(labels))
        print(
            "PHANTOM_TRANSPORT PASS "
            f"p={prime} state={state} indices={indices} labels={labels}"
        )

    # The marked coordinate is not even injective modulo reflection on all
    # clean short projective fibres.  At p=709 the complete solution vector,
    # hence both its projective state and its companion label, repeats after
    # 18 steps.  The gap is below sqrt(p) and the standard carrier U_18 is a
    # unit.  Its reflected copy is the second pair below.
    prime = 709
    apery, companion = solution_pair(prime)
    indices = (282, 300, 408, 426)
    assert 18**2 < prime
    assert carrier_is_clean(prime, 18, apery)
    assert (apery[282], companion[282]) == (apery[300], companion[300])
    assert (apery[408], companion[408]) == (apery[426], companion[426])
    assert 282 + 300 != prime - 1
    state = apery[282] * pow(companion[282], -1, prime) % prime
    assert state == 67
    assert verify_chain(prime, (282, 300), state) == (
        companion[282], companion[300]
    )
    print(
        "NONREFLECTION_LABEL_COLLISION PASS "
        f"p={prime} state={state} indices={indices} gap=18 carrier_clean=1"
    )
    print("DISTINGUISHED_LABEL_VERIFY PASS")


if __name__ == "__main__":
    main()
