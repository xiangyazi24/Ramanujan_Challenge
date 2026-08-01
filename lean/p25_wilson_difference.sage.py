#!/usr/bin/env sage-python
"""Temporary exact search for Wilson/challenge cross-difference products."""

from sage.all import QQ, PolynomialRing, binomial, matrix, vector


def generalized_binomial(value, count):
    result = QQ.one()
    for index in range(count):
        result *= value - index
        result /= index + 1
    return result


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


count = 42
target_numerator = vector(QQ, [30921, -32972, 8240])
target_denominator = vector(QQ, [33750, -36000, 9000])
signs = [1, -1, 1]
differences = [[], [], []]
gauge = QQ.one()
for index in range(count):
    u, v = wilson_pair(index + 2)
    for column in range(3):
        p = signs[column] * target_numerator[column]
        q = signs[column] * target_denominator[column]
        differences[column].append((p * u - q * v) / gauge)
    target_numerator *= challenge_matrix(index)
    target_denominator *= challenge_matrix(index)
    gauge *= -2*(index+2)**2*(index+3)**2*(2*index+5)*(2*index+7)**2

R = PolynomialRing(QQ, "n")
n = R.gen()


def rational_fit(data, maximum_total=32):
    points = list(enumerate(data))
    for total in range(maximum_total + 1):
        for numerator_degree in range(total + 1):
            denominator_degree = total - numerator_degree
            unknowns = total + 2
            rows = []
            for point, value in points[:unknowns + 4]:
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
            if denominator and all(numerator(point) == value * denominator(point)
                                   for point, value in points):
                return numerator.factor() / denominator.factor(), (
                    numerator_degree, denominator_degree)
    return None


for column in range(3):
    ratios = [differences[column][index + 1] / differences[column][index]
              for index in range(count - 1)]
    print("column", column, "signs", [value.sign() for value in differences[column][:8]])
    print("first", differences[column][0].factor())
    print("ratios", [value.factor() for value in ratios[:5]])
    print("ratio fit", rational_fit(ratios), flush=True)
