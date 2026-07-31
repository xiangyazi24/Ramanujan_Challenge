#!/usr/bin/env sage
"""Compare the corrected quotient residual with the W terminal packet.

This is an exact finite experiment.  It asks whether the real
moment-raised source left by the order-66 quotient operator is already
a short constant-coefficient convolution of the origin-cancelled W
coordinate at fixed M.
"""

from fractions import Fraction
from math import comb
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

# Provides integer_coefficients, multiplier, and origin_operator.
load(str(HERE / "q32_doubled_period_gauge_audit.sage"))
from q32_cartier_packet_audit import coefficient, polytope_points


def convolution(left, right):
    out = [QQ.zero()] * (len(left) + len(right) - 1)
    for i, aa in enumerate(left):
        for j, bb in enumerate(right):
            out[i + j] += aa * bb
    return out


quadratics = (
    (1, -6, 1),
    (-1, 2, 1),
    (2, -4, 1),
    (QQ(-1) / 2, 0, 1),
    (QQ(7) / 4, -3, 1),
    (1, 2, 1),
    (-1, -2, 1),
    (QQ(1) / 2, -2, 1),
    (QQ(1) / 4, -1, 1),
    (-2, 0, 1),
    (QQ(-1) / 4, -1, 1),
)
q = [QQ.one()]
for factor in quadratics:
    q = convolution(q, [QQ(entry) for entry in factor])
stable_operator = convolution(convolution(q, q), q)
assert len(stable_operator) == 67


def fixed_first_cell_shell(moment, endpoint):
    return ZZ(
        sum(
            coefficient(
                moment,
                endpoint * point[0],
                endpoint * point[1],
                endpoint * point[2],
            )
            for point in polytope_points(1)
        )
    )


def packet_values(moment, endpoint, maximum_order):
    raw = [
        fixed_first_cell_shell(moment, endpoint - residue)
        for residue in range(maximum_order + 1)
    ]
    return [
        ZZ(
            sum(
                (-1)**residue
                * binomial(order, residue)
                * raw[residue]
                for residue in range(order + 1)
            )
        )
        for order in range(maximum_order + 1)
    ]


def solve_and_test_convolution(residuals, source, order):
    rows = []
    rhs = []
    for start in range(order + 1):
        rows.append(
            [source[start + shift] for shift in range(order + 1)]
        )
        rhs.append(residuals[start])
    matrix_value = matrix(QQ, rows)
    if matrix_value.rank() < order + 1:
        return None
    coefficients_value = matrix_value.solve_right(vector(QQ, rhs))
    good = all(
        residuals[start]
        == sum(
            coefficients_value[shift] * source[start + shift]
            for shift in range(order + 1)
        )
        for start in range(len(residuals))
    )
    return good, coefficients_value


def polynomial_shift_span(
    residuals, source, shift_order, degree, prime=1000003
):
    field = GF(prime)
    columns = [
        (shift, power)
        for shift in range(shift_order + 1)
        for power in range(degree + 1)
    ]
    matrix_value = matrix(
        field,
        [
            [
                field(start)**power
                * field(source[start + shift])
                for shift, power in columns
            ]
            for start in range(len(residuals))
        ],
    )
    rhs = vector(field, [field(value) for value in residuals])
    return matrix_value.rank() == matrix_value.augment(rhs).rank()


for moment in (80, 90):
    global_index = moment + 1
    max_start = 45
    max_w_shift = 6
    maximum_order = 66 + max_start + max_w_shift
    y_values = packet_values(moment, moment, maximum_order)

    shifted_pairs = []
    for shift in range(origin_operator.order() + 1):
        shifted_moment = global_index + shift
        left = packet_values(
            shifted_moment, moment - 1, maximum_order
        )
        right = packet_values(
            shifted_moment, moment + 1, maximum_order
        )
        shifted_pairs.append(
            [left[index] + right[index]
             for index in range(maximum_order + 1)]
        )

    w_values = [
        sum(
            ZZ(integer_coefficients[shift](n=global_index))
            * shifted_pairs[shift][order]
            for shift in range(origin_operator.order() + 1)
        )
        - ZZ(multiplier(n=global_index)) * y_values[order]
        for order in range(maximum_order + 1)
    ]
    residuals = [
        sum(
            stable_operator[shift] * y_values[start + shift]
            for shift in range(67)
        )
        for start in range(max_start + 1)
    ]
    print("MOMENT", moment)
    print("RESIDUAL_BITS",
          [abs(value.numerator()).nbits() for value in residuals])
    print("W_BITS", [abs(value).nbits() for value in w_values[:11]])
    print("SCALAR_RATIOS",
          [residuals[start] * w_values[0]
           == residuals[0] * w_values[start]
           for start in range(1, 5)])
    for order in range(1, max_w_shift + 1):
        result = solve_and_test_convolution(
            residuals, w_values, order
        )
        print(
            "W_CONVOLUTION",
            order,
            None if result is None else result[0],
        )
    print(
        "POLYNOMIAL_W_SPANS",
        [
            (shift_order, degree)
            for shift_order in range(max_w_shift + 1)
            for degree in range(5)
            if polynomial_shift_span(
                residuals, w_values, shift_order, degree
            )
        ],
    )
