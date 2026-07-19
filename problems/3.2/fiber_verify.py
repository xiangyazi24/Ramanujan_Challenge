#!/usr/bin/env python3
"""Exact checks for the bordered Apéry fiber certificate.

This script uses only the Python standard library.  It verifies every finite
calculation requested in CODEX_SPEC_fiber.md, but it does not promote those
calculations to a proof of the missing range-content theorem; see the STALL
REPORT in fiber_bound.tex.
"""

from __future__ import annotations

from fractions import Fraction
from math import ceil, factorial, isqrt
import sys


Poly = tuple[int, ...]  # coefficients in increasing order
ZERO: Poly = (0,)
ONE: Poly = (1,)
X: Poly = (0, 1)
P_BASE: Poly = (5, 27, 51, 34)


def trim(coeffs: list[int] | tuple[int, ...]) -> Poly:
    out = list(coeffs)
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return tuple(out)


def poly_add(f: Poly, g: Poly) -> Poly:
    out = [0] * max(len(f), len(g))
    for i, value in enumerate(f):
        out[i] += value
    for i, value in enumerate(g):
        out[i] += value
    return trim(out)


def poly_neg(f: Poly) -> Poly:
    return trim([-value for value in f])


def poly_sub(f: Poly, g: Poly) -> Poly:
    return poly_add(f, poly_neg(g))


def poly_scale(c: int, f: Poly) -> Poly:
    return trim([c * value for value in f])


def poly_mul(f: Poly, g: Poly) -> Poly:
    if f == ZERO or g == ZERO:
        return ZERO
    out = [0] * (len(f) + len(g) - 1)
    for i, left in enumerate(f):
        for j, right in enumerate(g):
            out[i + j] += left * right
    return trim(out)


def poly_pow(f: Poly, exponent: int) -> Poly:
    result = ONE
    base = f
    n = exponent
    while n:
        if n & 1:
            result = poly_mul(result, base)
        base = poly_mul(base, base)
        n >>= 1
    return result


def linear(shift: int) -> Poly:
    """The polynomial x + shift."""

    return (shift, 1)


def poly_shift(f: Poly, shift: int) -> Poly:
    """Return f(x + shift), exactly over Z."""

    result = ZERO
    for coefficient in reversed(f):
        result = poly_add(poly_mul(result, linear(shift)), (coefficient,))
    return result


def poly_eval(f: Poly, value: int) -> int:
    result = 0
    for coefficient in reversed(f):
        result = result * value + coefficient
    return result


def poly_eval_mod(f: Poly, value: int, modulus: int) -> int:
    result = 0
    for coefficient in reversed(f):
        result = (result * value + coefficient) % modulus
    return result


def poly_eval_fraction(f: Poly, value: Fraction) -> Fraction:
    result = Fraction(0)
    for coefficient in reversed(f):
        result = result * value + coefficient
    return result


def poly_mod(f: Poly, modulus: int) -> Poly:
    return trim([coefficient % modulus for coefficient in f])


def divide_x_minus_root_mod(f: Poly, root: int, modulus: int) -> tuple[Poly, int]:
    """Synthetic division by x-root over F_modulus."""

    f = poly_mod(f, modulus)
    if f == ZERO:
        return ZERO, 0
    degree = len(f) - 1
    if degree == 0:
        return ZERO, f[0]
    quotient = [0] * degree
    quotient[-1] = f[-1]
    for i in range(degree - 1, 0, -1):
        quotient[i - 1] = (f[i] + root * quotient[i]) % modulus
    remainder = (f[0] + root * quotient[0]) % modulus
    return trim(quotient), remainder


def build_polynomials(max_m: int) -> tuple[list[Poly], list[Poly], list[Poly]]:
    """Build N_m, Pi_m, and the recurrence solution B_m over Z[x]."""

    n_polys = [ZERO for _ in range(max_m + 1)]
    pi_polys = [ONE for _ in range(max_m + 1)]
    b_polys = [ZERO for _ in range(max_m + 1)]
    n_polys[0] = ZERO
    n_polys[1] = ONE
    b_polys[0] = ONE
    b_polys[1] = ZERO

    for m in range(1, max_m + 1):
        pi_polys[m] = poly_mul(pi_polys[m - 1], poly_pow(linear(m), 3))

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


def b_closed(m: int, n_polys: list[Poly]) -> Poly:
    assert m >= 1
    return poly_neg(
        poly_mul(poly_pow(linear(1), 6), poly_shift(n_polys[m - 1], 1))
    )


def delta_poly(
    h: int, k: int, n_polys: list[Poly], pi_polys: list[Poly]
) -> Poly:
    d = k - h
    shifted_pi = poly_shift(pi_polys[d], h)
    shifted_n = poly_shift(n_polys[d], h)
    return poly_add(
        poly_sub(poly_mul(n_polys[h], shifted_pi), n_polys[k]),
        poly_mul(pi_polys[h], shifted_n),
    )


def raw_d_poly(
    h: int,
    k: int,
    n_polys: list[Poly],
    pi_polys: list[Poly],
    b_polys: list[Poly],
) -> Poly:
    return_row_h = poly_sub(pi_polys[h], b_polys[h])
    return_row_k = poly_sub(pi_polys[k], b_polys[k])
    return poly_sub(
        poly_mul(n_polys[h], return_row_k),
        poly_mul(n_polys[k], return_row_h),
    )


def apery_integers(length: int) -> list[int]:
    if length <= 0:
        return []
    if length == 1:
        return [1]
    values = [1, 5]
    for n in range(1, length - 1):
        p_n = 34 * n**3 + 51 * n**2 + 27 * n + 5
        numerator = p_n * values[n] - n**3 * values[n - 1]
        denominator = (n + 1) ** 3
        assert numerator % denominator == 0
        values.append(numerator // denominator)
    return values


def apery_mod_prime(p: int) -> list[int]:
    values = [1, 5 % p]
    for n in range(1, p - 1):
        p_n = (34 * n**3 + 51 * n**2 + 27 * n + 5) % p
        numerator = (p_n * values[n] - n**3 * values[n - 1]) % p
        denominator = pow(n + 1, 3, p)
        values.append(numerator * pow(denominator, -1, p) % p)
    return values


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    divisor = 3
    while divisor <= isqrt(n):
        if n % divisor == 0:
            return False
        divisor += 2
    return True


def primes_between(low: int, high: int) -> list[int]:
    return [n for n in range(low, high + 1) if is_prime(n)]


def fibers(values: list[int]) -> dict[int, list[int]]:
    result: dict[int, list[int]] = {}
    for index, value in enumerate(values):
        result.setdefault(value, []).append(index)
    return result


def pi_value_mod(t: int, m: int, p: int) -> int:
    result = 1
    for j in range(1, m + 1):
        result = result * pow(t + j, 3, p) % p
    return result


def forced_midpoint_constants(h: int, k: int) -> list[int]:
    """Constants c for the documented factors 2x+c of Delta_{h,k}."""

    constants: list[int] = []
    if h % 2 == 0:
        constants.append(h + 1)
    if k % 2 == 0:
        constants.append(k + 1)
    if (k - h) % 2 == 0:
        constants.append(h + k + 1)
    if (h, k) == (1, 2):
        # Here the sole midpoint factor is cubic, not merely linear.
        constants.extend([3, 3])
    return constants


def check_solution_formula() -> str:
    n_polys, pi_polys, b_polys = build_polynomials(8)

    # The coefficient pair of (Y_1,Y_0) obtained directly from the recurrence.
    coeff_y1 = [ZERO for _ in range(9)]
    coeff_y0 = [ZERO for _ in range(9)]
    coeff_y1[0], coeff_y1[1] = ZERO, ONE
    coeff_y0[0], coeff_y0[1] = ONE, ZERO
    for m in range(1, 8):
        p_shifted = poly_shift(P_BASE, m)
        sixth = poly_pow(linear(m), 6)
        coeff_y1[m + 1] = poly_sub(
            poly_mul(p_shifted, coeff_y1[m]),
            poly_mul(sixth, coeff_y1[m - 1]),
        )
        coeff_y0[m + 1] = poly_sub(
            poly_mul(p_shifted, coeff_y0[m]),
            poly_mul(sixth, coeff_y0[m - 1]),
        )

    for m in range(9):
        assert coeff_y1[m] == n_polys[m]
        assert coeff_y0[m] == b_polys[m]
    for m in range(1, 9):
        assert b_polys[m] == b_closed(m, n_polys)

    ell = [0, 1]
    for _ in range(1, 8):
        ell.append(34 * ell[-1] - ell[-2])
    for m in range(1, 9):
        assert len(n_polys[m]) - 1 == 3 * (m - 1)
        assert n_polys[m][-1] == ell[m]

    # Also verify the exact bordered factorization in small symbolic cases.
    factor_checks = 0
    for k in range(2, 9):
        for h in range(1, k):
            delta = delta_poly(h, k, n_polys, pi_polys)
            raw = raw_d_poly(h, k, n_polys, pi_polys, b_polys)
            assert raw == poly_mul(pi_polys[h], delta)
            assert len(delta) - 1 == 3 * (k - 1)
            assert delta[-1] == ell[h] + ell[k - h] - ell[k]
            for constant in forced_midpoint_constants(h, k):
                assert poly_eval_fraction(delta, Fraction(-constant, 2)) == 0
            factor_checks += 1

    return (
        "solution coefficients and B_m=-(x+1)^6 N_{m-1}(x+1) for m<=8; "
        f"{factor_checks} bordered factorizations"
    )


def check_triple_determinants() -> str:
    n_polys, pi_polys, b_polys = build_polynomials(25)
    tested_primes = primes_between(50, 300)
    triple_count = 0

    for p in tested_primes:
        value_fibers = fibers(apery_mod_prime(p))
        for a, positions in value_fibers.items():
            if a == 0:
                continue
            size = len(positions)
            for left in range(size - 2):
                t = positions[left]
                for middle in range(left + 1, size - 1):
                    h = positions[middle] - t
                    if h >= 25:
                        break
                    for right in range(middle + 1, size):
                        k = positions[right] - t
                        if k > 25:
                            break
                        n_h = poly_eval_mod(n_polys[h], t, p)
                        n_k = poly_eval_mod(n_polys[k], t, p)
                        b_h = poly_eval_mod(b_polys[h], t, p)
                        b_k = poly_eval_mod(b_polys[k], t, p)
                        pi_h = pi_value_mod(t, h, p)
                        pi_k = pi_value_mod(t, k, p)
                        determinant = (
                            n_h * (pi_k - b_k) - n_k * (pi_h - b_h)
                        ) % p
                        assert determinant == 0, (p, a, t, h, k, determinant)
                        triple_count += 1

    assert len(tested_primes) == 47
    return (
        f"all {len(tested_primes)} primes in [50,300], every a!=0, "
        f"{triple_count} in-fiber triples with k<=25"
    )


def check_finite_nonvanishing() -> str:
    n_polys, pi_polys, _ = build_polynomials(25)
    deltas = {
        (h, k): delta_poly(h, k, n_polys, pi_polys)
        for k in range(2, 26)
        for h in range(1, k)
    }
    tested_primes = primes_between(50, 300)
    reductions = 0

    for p in tested_primes:
        inv_two = pow(2, -1, p)
        for (h, k), delta in deltas.items():
            quotient = poly_mod(delta, p)
            for constant in forced_midpoint_constants(h, k):
                root = (-constant * inv_two) % p
                quotient, remainder = divide_x_minus_root_mod(quotient, root, p)
                assert remainder == 0, (p, h, k, constant)
            assert quotient != ZERO, (p, h, k, "zero polynomial")
            reductions += 1

    return (
        f"{reductions} reduced Delta_(h,k), h<k<=25; raw Pi_h and all "
        "documented midpoint factors removed"
    )


def check_endpoint_formulas() -> str:
    n_polys, pi_polys, _ = build_polynomials(25)
    apery = apery_integers(25)
    evaluations = 0

    for k in range(2, 26):
        for h in range(1, k):
            delta = delta_poly(h, k, n_polys, pi_polys)
            raw = poly_mul(pi_polys[h], delta)
            for r in range(1, k + 1):
                factorial_factor = (factorial(r - 1) * factorial(k - r)) ** 3
                if r <= h:
                    expected = (
                        (-1) ** (r - 1)
                        * apery[r - 1]
                        * (apery[h - r] - apery[k - r])
                        * factorial_factor
                    )
                    assert poly_eval(raw, -r) == 0
                else:
                    expected = (
                        (-1) ** (r - 1)
                        * apery[k - r]
                        * (apery[r - h - 1] - apery[r - 1])
                        * factorial_factor
                    )
                assert poly_eval(delta, -r) == expected, (h, k, r)
                evaluations += 1

    # A concrete witness that the proposed two-evaluation argument is invalid.
    p = 131
    values = apery_mod_prime(p)
    assert values[10] == values[53] == 15
    assert values[11] == values[54] == 15
    r = 3
    delta_minus_three = (
        values[r - 1]
        * (values[12 - r] - values[55 - r])
        * pow(factorial(r - 1) * factorial(55 - r), 3, p)
    ) % p
    assert delta_minus_three == 2

    return (
        f"{evaluations} exact endpoint evaluations for h<k<=25; "
        "p=131,(h,k)=(12,55) two-point obstruction confirmed"
    )


def conditional_exact_bound(p: int) -> tuple[int, int]:
    h_cutoff = ceil((2 * p / 3) ** 0.25)
    short = h_cutoff * (h_cutoff - 1) * (2 * h_cutoff - 1) // 2
    bound = 2 * ((p - 1) // (h_cutoff + 1)) + short + 2
    return h_cutoff, bound


def check_empirical_fibers() -> str:
    tested_primes = primes_between(7, 2000)
    max_nonzero = 0
    max_nonzero_primes: list[int] = []
    max_zero = 0
    max_all = 0
    worst_ratio = 0.0

    for p in tested_primes:
        counts = {a: len(indices) for a, indices in fibers(apery_mod_prime(p)).items()}
        nonzero = max((count for a, count in counts.items() if a != 0), default=0)
        zero = counts.get(0, 0)
        overall = max(counts.values())
        _, bound = conditional_exact_bound(p)
        assert nonzero <= bound
        # The zero fiber has a separate proved O(p^(2/3)) bound; this comparison
        # is only the empirical curve requested by the specification.
        assert overall <= bound
        ratio = nonzero / bound
        worst_ratio = max(worst_ratio, ratio)

        if nonzero > max_nonzero:
            max_nonzero = nonzero
            max_nonzero_primes = [p]
        elif nonzero == max_nonzero:
            max_nonzero_primes.append(p)
        max_zero = max(max_zero, zero)
        max_all = max(max_all, overall)

    assert len(tested_primes) == 300
    assert max_nonzero == 12
    assert 1231 in max_nonzero_primes and 1933 in max_nonzero_primes
    assert max_zero == 6
    return (
        f"{len(tested_primes)} primes p<=2000; max_(a!=0)={max_nonzero} "
        f"at {max_nonzero_primes}, max_zero={max_zero}, max_all={max_all}, "
        f"largest observed/conditional-curve ratio={worst_ratio:.6f}"
    )


def main() -> int:
    checks = [
        ("a", check_solution_formula),
        ("b", check_triple_determinants),
        ("c", check_finite_nonvanishing),
        ("d", check_endpoint_formulas),
        ("e", check_empirical_fibers),
    ]
    failures = 0
    for label, check in checks:
        try:
            detail = check()
        except Exception as error:  # Print a useful failure and keep the exit code nonzero.
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
