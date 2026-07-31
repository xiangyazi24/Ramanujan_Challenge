#!/usr/bin/env sage
"""Exact maximal-margin boundary-increment audit.

For the diagonal-margin vectors

    U_m = G_{D-m+1,N+m-2}(Y,W),
    V_m = G_{D,N+m-2}(Y,W),

the last step adds one point on each side.  This script evaluates the
primitive high-difference increments

    U_m-U_{m-1},  V_m-V_{m-1}

before their Pascal binomial multipliers.  It also checks the constant-term
identity

    Delta^L C_M(D)
      = sum_{kappa in P(Z)} CT[
          Lambda^M X^(-D*kappa) (X^(-kappa)-1)^L
        ]

by expanding the right side back into coefficient calls.  The latter check
is deliberately independent of the shell finite-difference computation.
"""

from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

# Loading the gauge file constructs and certifies the moment operator.  Its
# public objects are used exactly as in the multi-margin audit.
load(str(HERE / "q32_doubled_period_gauge_audit.sage"))
from q32_cartier_packet_audit import coefficient, polytope_points

common_origin_operator = P.lclm(apery_operator)
assert common_origin_operator.order() == 5
print(
    "COMMON_ORIGIN_OPERATOR",
    common_origin_operator.order(),
    tuple(
        (
            K(common_origin_operator[shift]).numerator().degree(),
            K(common_origin_operator[shift]).denominator().degree(),
        )
        for shift in range(common_origin_operator.order() + 1)
    ),
)


CASES = (
    (200, 128, 63),
    (272, 180, 63),
    (300, 180, 57),
    (321, 168, 53),
)


def finite_difference(values, start, order):
    return ZZ(
        sum(
            (-1) ** (order - shift)
            * binomial(order, shift)
            * values[start + shift]
            for shift in range(order + 1)
        )
    )


def coefficient_difference(moment, start, order):
    """Expand the constant-term formula coefficient by coefficient."""

    total = ZZ(0)
    for kappa in polytope_points(moment // start):
        for shift in range(order + 1):
            node = start + shift
            total += (
                (-1) ** (order - shift)
                * binomial(order, shift)
                * coefficient(
                    moment,
                    node * kappa[0],
                    node * kappa[1],
                    node * kappa[2],
                )
            )
    return total


def primitive_content(values):
    answer = ZZ(0)
    for value in values:
        answer = gcd(answer, ZZ(value))
    return abs(answer)


def endpoint(values_Y, values_W, d, length, D):
    left_Y = newton_carrier(values_Y, d, length)
    left_W = newton_carrier(values_W, d, length)
    right_Y = newton_carrier(values_Y, D, length)
    right_W = newton_carrier(values_W, D, length)
    raw = left_Y * right_W - left_W * right_Y
    count = D - d
    weights = [
        binomial(d + shift + length + 1, length)
        for shift in range(count)
    ]
    content = gcd_many(weights)
    assert raw % content == 0
    return raw // content


def euclidean_profile(left, right, limit=10):
    a, b = abs(ZZ(left)), abs(ZZ(right))
    out = []
    while b and len(out) < limit:
        quotient, remainder = a.quo_rem(b)
        out.append(
            (
                quotient
                if abs(quotient).nbits() <= 40
                else f"<{abs(quotient).nbits()} bits>",
                abs(remainder).nbits(),
            )
        )
        a, b = b, remainder
    return tuple(out)


def boundary_record(index, D, N):
    moment = index - 1
    left_margin = D - moment // 2
    right_margin = moment - D - N + 2
    margin = min(left_margin, right_margin)
    assert margin >= 2
    length = N + margin - 2
    left = D - margin + 1
    right = D + length

    nodes = list(range(left, right + 1))
    required = list(range(left - 1, right + 2))
    Y = shell_batch(moment, nodes)
    moment_shells = {
        shift: shell_batch(index + shift, required)
        for shift in range(origin_operator.order() + 1)
    }
    Z = {
        node: sum(
            ZZ(integer_coefficients[shift](n=index))
            * (
                moment_shells[shift][node - 1]
                + moment_shells[shift][node + 1]
            )
            for shift in range(origin_operator.order() + 1)
        )
        for node in nodes
    }
    mu = ZZ(multiplier(n=index))
    W = {node: Z[node] - mu * Y[node] for node in nodes}

    left_Y = finite_difference(Y, left, length)
    left_W = finite_difference(W, left, length)
    right_Y = finite_difference(Y, D, length)
    right_W = finite_difference(W, D, length)

    current = endpoint(Y, W, left, length, D)
    previous = endpoint(Y, W, left + 1, length - 1, D)
    common = gcd(abs(current), abs(previous))

    # The direct coefficient expansion is affordable for Y at the active
    # right boundary and provides an independent exact CT audit.
    right_ct = None
    if right_margin <= left_margin:
        assert right == moment
        right_ct = coefficient_difference(moment, D, length)
        assert right_ct == right_Y

    record = {
        "n": index,
        "core": (D, N),
        "margin": margin,
        "active": (
            "both"
            if left_margin == right_margin
            else ("left" if left_margin < right_margin else "right")
        ),
        "window": (left, right),
        "difference_bits": (
            abs(left_Y).nbits(),
            abs(left_W).nbits(),
            abs(right_Y).nbits(),
            abs(right_W).nbits(),
        ),
        "left_content": primitive_content((left_Y, left_W)),
        "right_content": primitive_content((right_Y, right_W)),
        "cross_bits": abs(left_Y * right_W - left_W * right_Y).nbits(),
        "endpoint_bits": (abs(previous).nbits(), abs(current).nbits()),
        "endpoint_gcd_bits": common.nbits(),
        "euclidean_head": euclidean_profile(current, previous),
        "right_ct_checked": right_ct is not None,
    }
    print("BOUNDARY_INCREMENT", record)
    return record


for case in CASES:
    boundary_record(*case)

print("Q32_BOUNDARY_INCREMENT_AUDIT=PASS")
