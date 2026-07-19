#!/usr/bin/env python3
"""Exact checks for the range nonvanishing theorem.

Only the Python standard library is used.  Finite computations here check the
symbolic identities used in nv_theorem.tex; they are not substituted for the
proof.  In particular, the script checks both the theorem's natural range
``k < p`` and the exact degeneration at ``(p, k) = (7, 21)``.
"""

from __future__ import annotations

from fractions import Fraction
from math import isqrt
import random
import sys
from typing import Union


Number = Union[int, Fraction]
Poly = tuple[Number, ...]  # coefficients in increasing order
ZERO: Poly = (0,)
ONE: Poly = (1,)
P_BASE: Poly = (5, 27, 51, 34)


def trim(coeffs: list[Number] | tuple[Number, ...]) -> Poly:
    out = list(coeffs)
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return tuple(out)


def poly_add(left: Poly, right: Poly) -> Poly:
    out: list[Number] = [0] * max(len(left), len(right))
    for index, value in enumerate(left):
        out[index] += value
    for index, value in enumerate(right):
        out[index] += value
    return trim(out)


def poly_neg(poly: Poly) -> Poly:
    return trim([-value for value in poly])


def poly_sub(left: Poly, right: Poly) -> Poly:
    return poly_add(left, poly_neg(right))


def poly_mul(left: Poly, right: Poly) -> Poly:
    if left == ZERO or right == ZERO:
        return ZERO
    out: list[Number] = [0] * (len(left) + len(right) - 1)
    for i, a_value in enumerate(left):
        for j, b_value in enumerate(right):
            out[i + j] += a_value * b_value
    return trim(out)


def poly_pow(poly: Poly, exponent: int) -> Poly:
    result = ONE
    base = poly
    power = exponent
    while power:
        if power & 1:
            result = poly_mul(result, base)
        base = poly_mul(base, base)
        power >>= 1
    return result


def linear(shift: Number) -> Poly:
    return (shift, 1)


def poly_shift(poly: Poly, shift: Number) -> Poly:
    """Return poly(x + shift), exactly."""

    result = ZERO
    for coefficient in reversed(poly):
        result = poly_add(poly_mul(result, linear(shift)), (coefficient,))
    return result


def poly_mod(poly: Poly, modulus: int) -> tuple[int, ...]:
    out: list[int] = []
    for coefficient in poly:
        if isinstance(coefficient, Fraction):
            value = (
                coefficient.numerator
                * pow(coefficient.denominator, -1, modulus)
            ) % modulus
        else:
            value = coefficient % modulus
        out.append(value)
    return trim(out)  # type: ignore[return-value]


def poly_add_mod(left: Poly, right: Poly, modulus: int) -> tuple[int, ...]:
    return poly_mod(poly_add(left, right), modulus)


def poly_mul_mod(left: Poly, right: Poly, modulus: int) -> tuple[int, ...]:
    return poly_mod(poly_mul(left, right), modulus)


def poly_pow_mod(poly: Poly, exponent: int, modulus: int) -> tuple[int, ...]:
    result: tuple[int, ...] = (1,)
    base = poly_mod(poly, modulus)
    power = exponent
    while power:
        if power & 1:
            result = poly_mul_mod(result, base, modulus)
        base = poly_mul_mod(base, base, modulus)
        power >>= 1
    return result


def poly_shift_mod(poly: Poly, shift: int, modulus: int) -> tuple[int, ...]:
    result: tuple[int, ...] = (0,)
    shifted_linear = (shift % modulus, 1)
    for coefficient in reversed(poly):
        result = poly_mul_mod(result, shifted_linear, modulus)
        result = poly_add_mod(result, (int(coefficient) % modulus,), modulus)
    return result


def coefficient(poly: Poly, exponent: int) -> Number:
    if exponent < 0 or exponent >= len(poly):
        return 0
    return poly[exponent]


def build_polynomials(max_m: int) -> tuple[list[Poly], list[Poly], list[Poly]]:
    """Build N_m, Pi_m, and B_m over Z[x]."""

    n_polys = [ZERO for _ in range(max_m + 1)]
    pi_polys = [ONE for _ in range(max_m + 1)]
    b_polys = [ZERO for _ in range(max_m + 1)]
    n_polys[0], n_polys[1] = ZERO, ONE
    b_polys[0], b_polys[1] = ONE, ZERO

    for m in range(1, max_m + 1):
        pi_polys[m] = poly_mul(
            pi_polys[m - 1], poly_pow(linear(m), 3)
        )

    for m in range(1, max_m):
        p_shifted = poly_shift(P_BASE, m)
        sixth = poly_pow(linear(m), 6)
        n_polys[m + 1] = poly_sub(
            poly_mul(p_shifted, n_polys[m]),
            poly_mul(sixth, n_polys[m - 1]),
        )
        b_polys[m + 1] = poly_sub(
            poly_mul(p_shifted, b_polys[m]),
            poly_mul(sixth, b_polys[m - 1]),
        )
    return n_polys, pi_polys, b_polys


def delta_poly(
    h: int, k: int, n_polys: list[Poly], pi_polys: list[Poly]
) -> Poly:
    d = k - h
    return poly_add(
        poly_sub(
            poly_mul(n_polys[h], poly_shift(pi_polys[d], h)),
            n_polys[k],
        ),
        poly_mul(pi_polys[h], poly_shift(n_polys[d], h)),
    )


def ell_values(max_m: int) -> list[int]:
    values = [0, 1]
    for _ in range(1, max_m):
        values.append(34 * values[-1] - values[-2])
    return values


def u_closed(m: int, ell: list[int]) -> Fraction:
    if m == 0:
        return Fraction(0)
    return Fraction(
        -(m - 1) * (32 * m**2 - 64 * m - 11) * ell[m]
        + 5 * m * ell[m - 1],
        256,
    )


def v_closed(m: int, ell: list[int]) -> Fraction:
    if m == 0:
        return Fraction(0)
    return Fraction(
        (m - 1)
        * (
            5120 * m**5
            - 31744 * m**4
            + 62016 * m**3
            - 19264 * m**2
            - 37024 * m
            - 2095
        )
        * ell[m],
        655360,
    ) - Fraction(
        m * (320 * m**3 - 1600 * m**2 - 320 * m + 621) * ell[m - 1],
        131072,
    )


def is_prime(number: int) -> bool:
    if number < 2:
        return False
    if number % 2 == 0:
        return number == 2
    divisor = 3
    while divisor <= isqrt(number):
        if number % divisor == 0:
            return False
        divisor += 2
    return True


def primes_up_to(limit: int) -> list[int]:
    return [number for number in range(2, limit + 1) if is_prime(number)]


def check_centered_coefficients() -> str:
    max_m = 20
    n_polys, pi_polys, _ = build_polynomials(max_m)
    ell = ell_values(max_m)

    for m in range(1, max_m + 1):
        centered_n = poly_shift(n_polys[m], -Fraction(m + 1, 2))
        assert coefficient(centered_n, 3 * m - 3) == ell[m]
        assert coefficient(centered_n, 3 * m - 5) == u_closed(m, ell)
        assert coefficient(centered_n, 3 * m - 7) == v_closed(m, ell)

        centered_pi = poly_shift(pi_polys[m], -Fraction(m + 1, 2))
        assert coefficient(centered_pi, 3 * m) == 1
        assert coefficient(centered_pi, 3 * m - 1) == 0
        assert coefficient(centered_pi, 3 * m - 2) == Fraction(
            -m * (m**2 - 1), 8
        )

    delta_checks = 0
    for k in range(2, 13):
        for h in range(1, k):
            d = k - h
            degree = 3 * (k - 1)
            delta = delta_poly(h, k, n_polys, pi_polys)
            centered = poly_shift(delta, -Fraction(k + 1, 2))
            expected_l = ell[h] + ell[d] - ell[k]
            expected_c1 = Fraction(3, 2) * (d * ell[h] - h * ell[d])
            expected_c2 = (
                u_closed(h, ell)
                + u_closed(d, ell)
                - u_closed(k, ell)
                + Fraction(
                    d * (-d**2 + 1 + 12 * d - 3 * h * k) * ell[h]
                    + h * (-h**2 + 1 + 12 * h - 3 * d * k) * ell[d],
                    8,
                )
            )
            assert coefficient(centered, degree) == expected_l
            assert coefficient(centered, degree - 1) == expected_c1
            assert coefficient(centered, degree - 2) == expected_c2
            delta_checks += 1

    return (
        f"u_m,v_m and centered Pi_m for m<=20; "
        f"{delta_checks} exact L,C1,C2 expansions"
    )


def check_closed_recurrences_and_elimination() -> str:
    ell = ell_values(101)
    for m in range(1, 100):
        expected_u = (
            34 * u_closed(m, ell)
            - u_closed(m - 1, ell)
            - Fraction(3, 4) * (17 * m**2 - 17 * m - 2) * ell[m]
            + Fraction(3, 4) * m * (m - 2) * ell[m - 1]
        )
        assert u_closed(m + 1, ell) == expected_u

        expected_v = (
            34 * v_closed(m, ell)
            - v_closed(m - 1, ell)
            - Fraction(3, 4) * (17 * m**2 - 17 * m - 36) * u_closed(m, ell)
            + Fraction(3, 4) * (m**2 - 2 * m - 4) * u_closed(m - 1, ell)
            + Fraction(
                3 * (m - 1) * (17 * m**3 + 51 * m**2 - 90 * m - 24),
                64,
            )
            * ell[m]
            - Fraction(3 * m * (m - 2) * (m**2 - 6), 16) * ell[m - 1]
        )
        assert v_closed(m + 1, ell) == expected_v

    for h in range(1, 50):
        for d in range(1, 50):
            assert ell[h + d] == (
                34 * ell[h] * ell[d]
                - ell[h] * ell[d - 1]
                - ell[h - 1] * ell[d]
            )
            assert ell[h + d - 1] == (
                ell[h] * ell[d] - ell[h - 1] * ell[d - 1]
            )
    for m in range(1, 101):
        assert ell[m - 1] ** 2 - 34 * ell[m] * ell[m - 1] + ell[m] ** 2 == 1

    rng = random.Random(32021)
    for _ in range(2000):
        h = rng.randrange(1, 100)
        d = rng.randrange(1, 100)
        k = h + d
        u_value = rng.randrange(-1000, 1001)
        b_value = rng.randrange(-1000, 1001)
        h_times_d_value = d * (34 * u_value - b_value) - k
        left = (
            h_times_d_value**2
            - 34 * d * u_value * h_times_d_value
            + d**2 * u_value**2
            - h**2
            - d**2 * (b_value**2 - 34 * u_value * b_value + u_value**2 - 1)
        )
        right = 2 * d * k * (b_value - 17 * u_value + 1)
        assert left == right
        assert (
            (17 * u_value - 1) ** 2
            - 34 * u_value * (17 * u_value - 1)
            + u_value**2
            - 1
        ) == -288 * u_value**2

    return "u_m,v_m recurrences; addition, Cassini, and elimination identities"


def u_mod(m: int, ell: list[int], p: int) -> int:
    numerator = (
        -(m - 1) * (32 * m**2 - 64 * m - 11) * ell[m]
        + 5 * m * ell[m - 1]
    )
    return numerator * pow(256, -1, p) % p


def check_modular_classification() -> str:
    pairs = 0
    leading_candidates = 0
    dangerous = 0
    for p in primes_up_to(251):
        if p < 7:
            continue
        ell = [0, 1]
        for _ in range(1, p):
            ell.append((34 * ell[-1] - ell[-2]) % p)
        inv_two = pow(2, -1, p)
        inv_eight = pow(8, -1, p)
        inv_2048 = pow(2048, -1, p)
        for k in range(2, p):
            for h in range(1, k):
                pairs += 1
                d = k - h
                a_value = ell[h]
                c_value = ell[d]
                leading = (a_value + c_value - ell[k]) % p
                c1 = 3 * inv_two * (d * a_value - h * c_value) % p
                if leading != 0 or c1 != 0:
                    continue
                leading_candidates += 1
                assert a_value == 0 and c_value == 0
                b_value = ell[h - 1]
                d_value = ell[d - 1]
                assert b_value in (1, p - 1)
                assert d_value in (1, p - 1)
                c2 = (
                    u_mod(h, ell, p)
                    + u_mod(d, ell, p)
                    - u_mod(k, ell, p)
                    + inv_eight
                    * (
                        d * (-d**2 + 1 + 12 * d - 3 * h * k) * a_value
                        + h * (-h**2 + 1 + 12 * h - 3 * d * k) * c_value
                    )
                ) % p
                if c2 != 0:
                    continue
                dangerous += 1
                assert b_value == p - 1 and d_value == p - 1
                c4 = -75 * h * d * k * inv_2048 % p
                assert c4 != 0

    return (
        f"{pairs} pairs with p<=251; {leading_candidates} leading "
        f"candidates, {dangerous} C2-dangerous cases, all killed by C4"
    )


def check_stress_examples_and_random_deltas() -> str:
    n_polys, pi_polys, _ = build_polynomials(50)
    examples = [
        (9369319, 19, 38, 1454304),
        (6771937, 23, 46, 1989688),
        (45245801, 25, 50, 13210276),
    ]
    for p, h, k, documented_coefficient in examples:
        assert is_prime(p)
        delta = delta_poly(h, k, n_polys, pi_polys)
        centered = poly_shift_mod(delta, (-(k + 1) * pow(2, -1, p)) % p, p)
        degree = 3 * (k - 1)
        assert len(centered) - 1 == degree - 4
        expected = -75 * h * (k - h) * k * pow(2048, -1, p) % p
        assert centered[-1] == expected == documented_coefficient

    rng = random.Random(7321)
    tested = 0
    primes = [p for p in primes_up_to(997) if p >= 7]
    for _ in range(1000):
        p = rng.choice(primes)
        k = rng.randrange(2, min(p, 51))
        h = rng.randrange(1, k)
        delta = delta_poly(h, k, n_polys, pi_polys)
        assert poly_mod(delta, p) != (0,)
        tested += 1
    return f"3 large-prime drop-four stress cases and {tested} random full Deltas"


Matrix = tuple[tuple[Poly, Poly], tuple[Poly, Poly]]


def matrix_mul_mod(left: Matrix, right: Matrix, p: int) -> Matrix:
    entries: list[list[Poly]] = [[ZERO, ZERO], [ZERO, ZERO]]
    for i in range(2):
        for j in range(2):
            value = ZERO
            for r in range(2):
                value = poly_add_mod(
                    value,
                    poly_mul_mod(left[i][r], right[r][j], p),
                    p,
                )
            entries[i][j] = value
    return (
        (entries[0][0], entries[0][1]),
        (entries[1][0], entries[1][1]),
    )


def matrix_pow_mod(matrix: Matrix, exponent: int, p: int) -> Matrix:
    result: Matrix = ((ONE, ZERO), (ZERO, ONE))
    base = matrix
    power = exponent
    while power:
        if power & 1:
            result = matrix_mul_mod(result, base, p)
        base = matrix_mul_mod(base, base, p)
        power >>= 1
    return result


def check_exact_counterexample() -> str:
    p = 7
    transfer: Matrix = ((ONE, ZERO), (ZERO, ONE))
    for j in range(1, 8):
        matrix: Matrix = (
            (
                poly_mod(poly_shift(P_BASE, j), p),
                poly_mod(poly_neg(poly_pow(linear(j), 6)), p),
            ),
            (ONE, ZERO),
        )
        transfer = matrix_mul_mod(matrix, transfer, p)

    q = trim([0, -1, 0, 0, 0, 0, 0, 1])
    q3 = poly_pow_mod(q, 3, p)
    q6 = poly_pow_mod(q, 6, p)
    q9 = poly_pow_mod(q, 9, p)
    trace = poly_add_mod(transfer[0][0], transfer[1][1], p)
    determinant = poly_sub(
        poly_mul_mod(transfer[0][0], transfer[1][1], p),
        poly_mul_mod(transfer[0][1], transfer[1][0], p),
    )
    assert poly_mod(trace, p) == poly_mod(poly_neg(q3), p)
    assert poly_mod(determinant, p) == q6
    assert matrix_pow_mod(transfer, 3, p) == ((q9, ZERO), (ZERO, q9))

    n_polys, pi_polys, b_polys = build_polynomials(22)
    assert poly_mod(n_polys[21], p) == (0,)
    assert poly_mod(b_polys[21], p) == q9
    assert poly_mod(pi_polys[21], p) == q9
    for h in range(1, 21):
        assert poly_mod(delta_poly(h, 21, n_polys, pi_polys), p) == (0,)
    return "T^3=(x^7-x)^9 I and Delta_(h,21)=0 mod 7 for all 1<=h<21"


def main() -> int:
    if not __debug__:
        print("FAIL: run without python -O so that exact assertions remain active")
        return 1
    checks = [
        ("a", check_centered_coefficients),
        ("b", check_closed_recurrences_and_elimination),
        ("c", check_modular_classification),
        ("d", check_stress_examples_and_random_deltas),
        ("e", check_exact_counterexample),
    ]
    failures = 0
    for label, check in checks:
        try:
            detail = check()
        except Exception as error:
            failures += 1
            print(f"FAIL ({label}) {type(error).__name__}: {error}")
        else:
            print(f"PASS ({label}) {detail}")

    if failures:
        print(f"OVERALL FAIL: {failures} check(s) failed")
        return 1
    print("OVERALL PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
