#!/usr/bin/env sage-python
"""Temporary exact search for the Wilson-to-challenge rectangular gauge."""

from sage.all import QQ, PolynomialRing, binomial, factorial, matrix, vector


def generalized_binomial(value, count):
    answer = QQ.one()
    for index in range(count):
        answer *= value - index
        answer /= index + 1
    return answer


def wilson_pair(index):
    if index == 0:
        return QQ.one(), QQ.zero()
    z = 4 * index
    a = QQ(z - 1) / 2
    denominator = sum(
        binomial(index, j) * generalized_binomial(a, j)
        * generalized_binomial(a + j, j)
        for j in range(index + 1)
    )
    q_polynomial = QQ(z) / 4 * sum(
        binomial(index, j) * sum(
            generalized_binomial(a + j, j - k)
            * generalized_binomial(a - k, j - k)
            * QQ((-1) ** (k - 1)) / (k**2 * binomial(j, k)**2)
            for k in range(1, j + 1)
        )
        for j in range(1, index + 1)
    )
    partial = sum(QQ((-1) ** k) / (2 * k + 1)**2
                  for k in range(2 * index))
    numerator = denominator * partial + q_polynomial / (2 * z)
    return denominator, numerator


def challenge_matrix(n):
    return matrix(QQ, 3, 3, [
        (-2*n-5)*(n+3)**2*(136*n**4+1424*n**3+5548*n**2+9551*n+6141),
        384*n**6+6384*n**5+44168*n**4+162698*n**3+336377*n**2+369933*n+169011,
        -(480*n**4+4980*n**3+19210*n**2+32690*n+20730),
        (n+2)**2*(n+3)**2*(4*n+10)*(48*n**3+386*n**2+1017*n+879),
        (n+2)**2*(-272*n**5-3848*n**4-21732*n**3-61184*n**2-85761*n-47808),
        (n+2)**2*(320*n**3+2540*n**2+6610*n+5640),
        (-4*n-10)*(n+2)**2*(n+3)**2*(32*n**4+302*n**3+1037*n**2+1530*n+813),
        (n+2)**2*(192*n**6+2984*n**5+19116*n**4+64452*n**3+120256*n**2+117279*n+46476),
        (n+2)**2*(-16*n**5-408*n**4-2912*n**3-8884*n**2-12254*n-6240),
    ])


count = 36
wilson = [wilson_pair(index) for index in range(count + 1)]
target_numerator = vector(QQ, [30921, -32972, 8240])
target_denominator = vector(QQ, [33750, -36000, 9000])
gauge = QQ.one()
values = [[[], []] for _ in range(3)]
projection_values = {(left, right): [[], []]
                     for left in range(3) for right in range(left + 1, 3)}
signs = [1, -1, 1]
for index in range(count):
    if index > 0:
        u, v = wilson[index]
        u_next, v_next = wilson[index + 1]
        wronskian = u * v_next - u_next * v
        for column in range(3):
            q = signs[column] * target_denominator[column] / gauge
            p = signs[column] * target_numerator[column] / gauge
            values[column][0].append((index, (q * v_next - u_next * p) / wronskian))
            values[column][1].append((index, (u * p - q * v) / wronskian))
        positive_q = [signs[column] * target_denominator[column] / gauge
                      for column in range(3)]
        positive_p = [signs[column] * target_numerator[column] / gauge
                      for column in range(3)]
        for left, right in projection_values:
            determinant = (positive_q[left] * positive_p[right]
                           - positive_q[right] * positive_p[left])
            c_left = (u * positive_p[right] - v * positive_q[right]) / determinant
            c_right = (v * positive_q[left] - u * positive_p[left]) / determinant
            projection_values[(left, right)][0].append((index, c_left))
            projection_values[(left, right)][1].append((index, c_right))
    target_numerator *= challenge_matrix(index)
    target_denominator *= challenge_matrix(index)
    gauge *= -2*(index+2)**2*(index+3)**2*(2*index+5)*(2*index+7)**2

# Compare challenge envelopes at stage N with the r=2 Wilson approximant at n=N+2.
target_numerator = vector(QQ, [30921, -32972, 8240])
target_denominator = vector(QQ, [33750, -36000, 9000])
for index in range(8):
    u, v = wilson[index + 2]
    ratios = [target_numerator[column] / target_denominator[column]
              for column in range(3)]
    print("sandwich", index, [ratio - v / u for ratio in ratios], flush=True)
    target_numerator *= challenge_matrix(index)
    target_denominator *= challenge_matrix(index)


ring = PolynomialRing(QQ, "n")
n = ring.gen()


def rational_fit(data, maximum_total=32):
    for total in range(maximum_total + 1):
        for numerator_degree in range(total + 1):
            denominator_degree = total - numerator_degree
            unknowns = total + 2
            rows = []
            for point, value in data[:unknowns + 2]:
                rows.append(
                    [QQ(point)**power for power in range(numerator_degree + 1)]
                    + [-value * QQ(point)**power
                       for power in range(denominator_degree + 1)]
                )
            kernel = matrix(QQ, rows).right_kernel_matrix()
            if kernel.nrows() != 1:
                continue
            coefficients = list(kernel[0])
            numerator = sum(coefficients[power] * n**power
                            for power in range(numerator_degree + 1))
            denominator = sum(
                coefficients[numerator_degree + 1 + power] * n**power
                for power in range(denominator_degree + 1)
            )
            if denominator == 0:
                continue
            if all(numerator(point) == value * denominator(point)
                   for point, value in data):
                return numerator.factor() / denominator.factor(), (
                    numerator_degree, denominator_degree)
    return None


for column in range(3):
    for coordinate, name in enumerate(("a", "b")):
        print("fitting", column, name, flush=True)
        print(rational_fit(values[column][coordinate]), flush=True)

for pair, pair_values in projection_values.items():
    for coordinate, name in enumerate(("left", "right")):
        print("projecting", pair, name, flush=True)
        print(rational_fit(pair_values[coordinate]), flush=True)
