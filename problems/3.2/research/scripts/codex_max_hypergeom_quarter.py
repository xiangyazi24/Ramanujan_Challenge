#!/usr/bin/env python3
"""Exact checks for the tau/sigma pullbacks and their coefficient formulae.

Everything in this file is over Q.  It verifies, through the requested finite
range, the algebraic pullback from the Apéry series to the Franel series and
the Lagrange-inversion formula used in CODEX_JACOBSTHAL_DEEP.md.  It also
gives an exact obstruction to the most restrictive "one Pochhammer term"
interpretation of the half-integer recurrences.
"""

from fractions import Fraction as Q
from math import comb


ORDER = 40


def add(a, b, n=ORDER):
    out = [Q(0)] * n
    for i in range(n):
        out[i] = (a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)
    return out


def mul(a, b, n=ORDER):
    out = [Q(0)] * n
    for i, ai in enumerate(a[:n]):
        for j, bj in enumerate(b[: n - i]):
            out[i + j] += ai * bj
    return out


def power(a, exponent, n=ORDER):
    out = [Q(1)] + [Q(0)] * (n - 1)
    base = a[:n] + [Q(0)] * max(0, n - len(a))
    while exponent:
        if exponent & 1:
            out = mul(out, base, n)
        base = mul(base, base, n)
        exponent //= 2
    return out


def inverse(a, n=ORDER):
    assert a and a[0] != 0
    out = [1 / a[0]] + [Q(0)] * (n - 1)
    for k in range(1, n):
        out[k] = -sum(a[i] * out[k - i] for i in range(1, min(k, len(a) - 1) + 1)) / a[0]
    return out


def compose(a, xseries, n=ORDER):
    out = [Q(0)] * n
    xp = [Q(1)] + [Q(0)] * (n - 1)
    for ai in a[:n]:
        for j in range(n):
            out[j] += ai * xp[j]
        xp = mul(xp, xseries, n)
    return out


def binomial(alpha, n):
    out = Q(1)
    for j in range(n):
        out *= (alpha - j) / (j + 1)
    return out


def pochhammer(alpha, n):
    out = Q(1)
    for j in range(n):
        out *= alpha + j
    return out


def apery(n):
    b = [1]
    for m in range(n - 1):
        previous = b[m - 1] if m else 0
        numerator = (2 * m + 1) * (17 * m * m + 17 * m + 5) * b[m] - m**3 * previous
        quotient, remainder = divmod(numerator, (m + 1) ** 3)
        assert remainder == 0
        b.append(quotient)
    return [Q(v) for v in b]


def half_integer_sequences(n):
    tau = [Q(1), Q(5, 2)]
    sigma = [Q(1), Q(39, 2)]
    for m in range(1, n - 1):
        tau.append(
            (
                2 * (68 * m * m + 34 * m + 5) * tau[m]
                - (2 * m - 1) ** 2 * tau[m - 1]
            )
            / (4 * (m + 1) ** 2)
        )
        sigma.append(
            (
                2 * (68 * m * m + 102 * m + 39) * sigma[m]
                - (2 * m + 1) ** 2 * sigma[m - 1]
            )
            / (4 * (m + 1) ** 2)
        )
    return tau[:n], sigma[:n]


def franel(n):
    return [Q(sum(comb(m, k) ** 3 for k in range(m + 1))) for m in range(n)]


def solve_linear(matrix, rhs):
    a = [list(row) + [value] for row, value in zip(matrix, rhs)]
    size = len(a)
    for col in range(size):
        pivot = next(row for row in range(col, size) if a[row][col])
        a[col], a[pivot] = a[pivot], a[col]
        factor = a[col][col]
        a[col] = [value / factor for value in a[col]]
        for row in range(size):
            if row != col and a[row][col]:
                factor = a[row][col]
                a[row] = [a[row][j] - factor * a[col][j] for j in range(size + 1)]
    return [a[i][-1] for i in range(size)]


def pochhammer_obstruction(name, sequence):
    ratios = [sequence[n + 1] / sequence[n] for n in range(5)]
    # Put L=lambda, S=lambda(A+B), U=lambda AB.  The equation
    # r_n(n+1)(n+C)=L n^2+S n+U is linear in L,S,U,C.
    matrix = []
    rhs = []
    for n, ratio in enumerate(ratios[:4]):
        matrix.append([Q(n * n), Q(n), Q(1), -ratio * (n + 1)])
        rhs.append(ratio * n * (n + 1))
    lam, scaled_sum, scaled_product, cpar = solve_linear(matrix, rhs)
    predicted = (lam * 16 + scaled_sum * 4 + scaled_product) / ((4 + cpar) * 5)
    assert predicted != ratios[4]
    print(
        f"{name}: one-Pochhammer fit fails at n=4; "
        f"predicted={predicted}, actual={ratios[4]}, difference={ratios[4] - predicted}"
    )


def main():
    b = apery(ORDER)
    tau, sigma = half_integer_sequences(ORDER)

    assert mul(tau, tau) == b
    q_sigma_sq = mul([Q(1), Q(-34), Q(1)], mul(sigma, sigma))
    assert q_sigma_sq == b

    # t=x(1-8x)/(1+x), phi=x/t=(1+x)/(1-8x).
    one_plus_x = [Q(1), Q(1)]
    one_minus_8x = [Q(1), Q(-8)]
    t_of_x = mul([Q(0), Q(1), Q(-8)], inverse(one_plus_x))
    phi = mul(one_plus_x, inverse(one_minus_8x))
    h = franel(ORDER)
    sqrt_one_plus_x = [binomial(Q(1, 2), n) for n in range(ORDER)]
    one_plus_x_3_2 = [binomial(Q(3, 2), n) for n in range(ORDER)]
    tau_pullback = mul(sqrt_one_plus_x, h)
    sigma_pullback = mul(
        mul(one_plus_x_3_2, h),
        inverse([Q(1), Q(-16), Q(-8)]),
    )
    assert compose(tau, t_of_x) == tau_pullback
    assert compose(sigma, t_of_x) == sigma_pullback

    # q(t(x))=((1-16x-8x^2)/(1+x))^2, checked as a polynomial identity.
    numerator_t = [Q(0), Q(1), Q(-8)]
    numerator_q = add(
        add(mul(one_plus_x, one_plus_x), [-34 * v for v in mul(numerator_t, one_plus_x)]),
        mul(numerator_t, numerator_t),
    )
    assert numerator_q[:5] == mul([Q(1), Q(-16), Q(-8)], [Q(1), Q(-16), Q(-8)])[:5]

    # Franel's classical 2F1 pullback.
    z = mul(
        [Q(0), Q(0), Q(27)],
        power(inverse([Q(1), Q(-2)]), 3),
    )
    hyper = []
    factorial = 1
    for n in range(ORDER):
        if n:
            factorial *= n
        hyper.append(
            pochhammer(Q(1, 3), n)
            * pochhammer(Q(2, 3), n)
            / (factorial * factorial)
        )
    assert mul(inverse([Q(1), Q(-2)]), compose(hyper, z)) == h

    # Lagrange-Bueremann: [t^n]G(x(t))=(1/n)[x^(n-1)]G'(x)phi(x)^n.
    for name, expected, g in (("tau", tau, tau_pullback), ("sigma", sigma, sigma_pullback)):
        derivative = [Q(k + 1) * g[k + 1] for k in range(ORDER - 1)]
        recovered = [Q(1)]
        for n in range(1, ORDER):
            recovered.append(mul(derivative, power(phi, n), ORDER)[n - 1] / n)
        assert recovered == expected
        print(f"{name}: recurrence, pullback, and Lagrange formula agree through n={ORDER - 1}")

    pochhammer_obstruction("tau", tau)
    pochhammer_obstruction("sigma", sigma)
    print("Franel 2F1 pullback and q-pullback identities: VERIFIED")


if __name__ == "__main__":
    main()
