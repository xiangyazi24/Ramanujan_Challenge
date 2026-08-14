#!/usr/bin/env python3
"""Exact verifier for chatgpt_q32_turan_crossray.

This standard-library script checks the finite identities used in the Q8376
cross-ray report.  It does not claim to verify the asymptotic impossibility
statement, whose recurrence-theoretic input is the already banked
G_diff = SL_2 / no-Riccati theorem.

Checks:
  * the absolute endpoint has a nonorigin ray of coefficient at least 2^M;
  * that ray occurs in no other shell grid d with M/2 < d <= M;
  * Newton evaluation rows have augmentation one, while differences have zero;
  * the integral cross-ray quotient Z^r/ker(<-,s>) is rank one;
  * the Apéry recurrence and canonical homogeneous Casoratian;
  * the exact modular/Eichler counterexample (p,n)=(31,8):
        b_8 = 0, kappa_8 = 13, Xi_8 = 28 (mod 31),
    and the corresponding Green seam is nonzero.

Run from the repository root:
    python3 problems/3.2/research/scripts/chatgpt_q32_turan_crossray_verify.py
"""

from __future__ import annotations

from collections import defaultdict
from fractions import Fraction
from functools import reduce
from math import comb, gcd
from typing import Dict, Iterable, List, Sequence, Tuple


Exponent = Tuple[int, int, int]
Laurent = Dict[Exponent, int]


def multiply(left: Laurent, right: Laurent) -> Laurent:
    out: Dict[Exponent, int] = defaultdict(int)
    for u, a in left.items():
        for v, b in right.items():
            out[(u[0] + v[0], u[1] + v[1], u[2] + v[2])] += a * b
    return dict(out)


def lambda_polynomial() -> Laurent:
    one = {(0, 0, 0): 1}
    x = {(0, 0, 0): 1, (1, 0, 0): 1}
    y = {(0, 0, 0): 1, (0, 1, 0): 1}
    z = {(0, 0, 0): 1, (0, 0, 1): 1}
    bracket = multiply(y, z)
    bracket[(1, 1, 1)] = bracket.get((1, 1, 1), 0) + 1
    numerator = one
    for factor in (x, y, z, bracket):
        numerator = multiply(numerator, factor)
    return {
        (u[0] - 1, u[1] - 1, u[2] - 1): coefficient
        for u, coefficient in numerator.items()
    }


def power(base: Laurent, exponent: int) -> Laurent:
    out: Laurent = {(0, 0, 0): 1}
    for _ in range(exponent):
        out = multiply(out, base)
    return out


def endpoint_ray_checks() -> int:
    lam = lambda_polynomial()
    ray = (-1, 0, -1)
    assert lam[ray] == 2
    assert ray != (0, 0, 0)

    checks = 0
    for moment in range(1, 9):
        packet = power(lam, moment)
        endpoint = tuple(moment * coordinate for coordinate in ray)
        # Repeating the coefficient-2 monomial M times is one positive
        # contribution to this coefficient.
        assert packet[endpoint] >= 2**moment
        checks += 1

    # The point M*ray belongs to a d-grid iff d divides M.  A proper divisor
    # of M cannot lie strictly above M/2.
    for moment in range(2, 81):
        endpoint = tuple(moment * coordinate for coordinate in ray)
        for d in range(moment // 2 + 1, moment + 1):
            on_grid = all(coordinate % d == 0 for coordinate in endpoint)
            assert on_grid == (d == moment)
            checks += 1
    return checks


def newton_weights(start: int, order: int) -> List[int]:
    return [
        (-1) ** i
        * comb(start + i, i)
        * comb(start + order + 1, order - i)
        for i in range(order + 1)
    ]


def augmentation_checks() -> int:
    checks = 0
    for start in range(1, 15):
        for order in range(0, 12):
            weights = newton_weights(start, order)
            assert sum(weights) == 1
            checks += 1
            if order:
                difference = [
                    (-1) ** (order - i) * comb(order, i)
                    for i in range(order + 1)
                ]
                assert sum(difference) == 0
                checks += 1
    return checks


def dot(left: Sequence[int], right: Sequence[int]) -> int:
    return sum(a * b for a, b in zip(left, right))


def bezout_vector(values: Sequence[int]) -> Tuple[int, List[int]]:
    """Return positive gcd g and a vector a with a dot values = g."""

    coefficients = [0] * len(values)
    current = 0
    for index, value in enumerate(values):
        old = current
        g = gcd(old, value)
        if g == 0:
            continue
        # Extended Euclid for old*x + value*y = g.
        a0, b0 = abs(old), abs(value)
        x0, x1, y0, y1 = 1, 0, 0, 1
        while b0:
            q, a0, b0 = a0 // b0, b0, a0 % b0
            x0, x1 = x1, x0 - q * x1
            y0, y1 = y1, y0 - q * y1
        x = x0 if old >= 0 else -x0
        y = y0 if value >= 0 else -y0
        coefficients = [x * coefficient for coefficient in coefficients]
        coefficients[index] += y
        current = g
    if current < 0:
        current = -current
        coefficients = [-coefficient for coefficient in coefficients]
    assert dot(coefficients, values) == current
    return current, coefficients


def crossray_quotient_checks() -> int:
    samples = [
        ([1, 1, 1, 1], [7, -4, 9, 2]),
        ([2, 3, 5], [11, -8, 6]),
        ([6, 10, 14, 22], [-3, 5, 7, -2]),
        ([-3, 6, 9], [8, 1, -5]),
    ]
    checks = 0
    for s, c in samples:
        g, anchor = bezout_vector(s)
        assert g == reduce(gcd, (abs(value) for value in s), 0)
        assert dot(anchor, s) == g
        assert dot(c, s) % g == 0
        quotient = dot(c, s) // g
        zero_ray = [
            c_i - quotient * anchor_i
            for c_i, anchor_i in zip(c, anchor)
        ]
        assert dot(zero_ray, s) == 0
        reconstructed = [
            zero + quotient * anchor_i
            for zero, anchor_i in zip(zero_ray, anchor)
        ]
        assert reconstructed == c
        checks += 1
    return checks


def apery_values(limit: int) -> List[int]:
    if limit == 0:
        return [1]
    values = [1, 5]
    for n in range(1, limit):
        p_n = 34 * n**3 + 51 * n**2 + 27 * n + 5
        numerator = p_n * values[n] - n**3 * values[n - 1]
        denominator = (n + 1) ** 3
        assert numerator % denominator == 0
        values.append(numerator // denominator)
    return values[: limit + 1]


def homogeneous_companion(limit: int) -> List[Fraction]:
    if limit == 0:
        return [Fraction(0)]
    values = [Fraction(0), Fraction(1)]
    for n in range(1, limit):
        p_n = 34 * n**3 + 51 * n**2 + 27 * n + 5
        values.append(
            (p_n * values[n] - n**3 * values[n - 1]) / (n + 1) ** 3
        )
    return values[: limit + 1]


def series_mul(left: Sequence[Fraction], right: Sequence[Fraction], limit: int) -> List[Fraction]:
    out = [Fraction(0)] * (limit + 1)
    for i, a in enumerate(left):
        if i > limit:
            break
        for j, b in enumerate(right):
            if i + j > limit:
                break
            out[i + j] += a * b
    return out


def series_inv(values: Sequence[Fraction], limit: int) -> List[Fraction]:
    assert values[0] != 0
    out = [Fraction(1, 1) / values[0]]
    for n in range(1, limit + 1):
        known = sum(values[k] * out[n - k] for k in range(1, n + 1))
        out.append(-known / values[0])
    return out


def inverse_sqrt_discriminant(limit: int) -> List[Fraction]:
    """Q(t)=(1-34t+t^2)^(-1/2), solved from D*Q^2=1."""

    q = [Fraction(1)]

    def coefficient(candidate: Sequence[Fraction], n: int) -> Fraction:
        square = series_mul(candidate, candidate, n)
        answer = square[n]
        if n >= 1:
            answer -= 34 * square[n - 1]
        if n >= 2:
            answer += square[n - 2]
        return answer

    for n in range(1, limit + 1):
        candidate = q + [Fraction(0)]
        known = coefficient(candidate, n)
        # The unknown q_n occurs as 2*q_n in Q^2.
        q.append(-known / 2)
        assert coefficient(q, n) == 0
    return q


def eichler_data(limit: int) -> Tuple[List[Fraction], List[Fraction]]:
    b_int = apery_values(limit)
    b = [Fraction(value) for value in b_int]
    q = inverse_sqrt_discriminant(limit)
    b_square = series_mul(b, b, limit)
    g = series_mul(q, series_inv(b_square, limit), limit)

    kappa = [Fraction(0), Fraction(-36)]
    for n in range(2, limit + 1):
        p_prev = 34 * (n - 1) ** 3 + 51 * (n - 1) ** 2 + 27 * (n - 1) + 5
        kappa.append(
            (
                p_prev * kappa[n - 1]
                - (n - 1) ** 3 * kappa[n - 2]
                - 5 * g[n]
            )
            / n**3
        )
    return g, kappa


def mod_fraction(value: Fraction, prime: int) -> int:
    denominator = value.denominator % prime
    assert denominator
    return value.numerator % prime * pow(denominator, -1, prime) % prime


def recurrence_and_modular_checks() -> int:
    limit = 12
    b = apery_values(limit)
    u = homogeneous_companion(limit)

    for n in range(1, limit + 1):
        casoratian = n**3 * (Fraction(b[n - 1]) * u[n] - Fraction(b[n]) * u[n - 1])
        assert casoratian == 1

    g, kappa = eichler_data(limit)
    assert [int(g[i]) for i in range(5)] == [1, 7, 192, 5520, 165168]

    prime, n = 31, 8
    assert b[n] % prime == 0
    assert b[n + 1] % prime != 0
    assert mod_fraction(u[n], prime) != 0
    assert mod_fraction(kappa[n], prime) == 13

    xi = Fraction(n**3) * (
        Fraction(b[n - 1]) * kappa[n] - Fraction(b[n]) * kappa[n - 1]
    )
    assert mod_fraction(xi, prime) == 28

    green_seam = Fraction(b[n]) * kappa[n + 1] - Fraction(b[n + 1]) * kappa[n]
    assert mod_fraction(green_seam, prime) == (
        -b[n + 1] * mod_fraction(kappa[n], prime)
    ) % prime
    assert mod_fraction(green_seam, prime) != 0
    return limit + 6


def main() -> None:
    result = {
        "endpoint_ray_checks": endpoint_ray_checks(),
        "augmentation_checks": augmentation_checks(),
        "crossray_quotient_checks": crossray_quotient_checks(),
        "recurrence_and_modular_checks": recurrence_and_modular_checks(),
    }
    print("CHATGPT_Q32_TURAN_CROSSRAY_VERIFY=PASS", result)


if __name__ == "__main__":
    main()
