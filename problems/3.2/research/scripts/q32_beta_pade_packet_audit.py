#!/usr/bin/env python3
"""Exact audit of the beta--Padé packet identity for Apéry shells.

For

    Q(T) = sum_i (-1)^i binom(D-1+i,i) binom(D+N,N-i) T^i

the incomplete-beta identity is

    T^D Q(T) + (1-T)^(N+1) R(T) = 1.

After multiplying by T^(-1), inserting the coefficient array of the
Apéry Laurent polynomial, and summing over the fixed cube
||kappa||_infty <= floor(M/(D-1)), this becomes

    A_{D,N} + H_{M,D,N} = S_{M,K},

where A is the Newton shell carrier and

    S_{M,K} = sum_{||kappa||_infty <= K} c_M(-kappa).

In particular S is a near-origin coefficient packet, not the shell
C_M(D-1).  The script also audits the exact local obstruction to
cancelling a packet shared by several carriers: on a prime node common
to all their intervals, a coefficient combination with sum zero is
coefficientwise divisible by that prime, independently of the shell
values.  Thus its guaranteed prime factor is presentation content.

This is an identity/regression audit.  It proves no asymptotic radical
bound.
"""

from __future__ import annotations

from functools import lru_cache
from itertools import product
from math import comb


Exponent = tuple[int, int, int]
Polynomial = dict[int, int]


def comb_zero(n: int, k: int) -> int:
    return comb(n, k) if 0 <= k <= n else 0


@lru_cache(maxsize=None)
def coefficient(m: int, u: int, v: int, w: int) -> int:
    """[x^u y^v z^w] Lambda^m, independently using formula (49.1)."""
    return sum(
        comb(m, t)
        * comb_zero(m, t - u)
        * comb_zero(2 * m - t, m - v)
        * comb_zero(2 * m - t, m - w)
        for t in range(m + 1)
    )


def scale(vector: Exponent, scalar: int) -> Exponent:
    return tuple(scalar * entry for entry in vector)  # type: ignore[return-value]


def cube(radius: int):
    return product(range(-radius, radius + 1), repeat=3)


def shell(m: int, spacing: int) -> int:
    assert spacing >= 1
    radius = m // spacing
    return sum(coefficient(m, *scale(kappa, spacing)) for kappa in cube(radius))


def packet(m: int, radius: int, step: int = -1) -> int:
    """The truncated packet sum_kappa c_M(step*kappa)."""
    return sum(coefficient(m, *scale(kappa, step)) for kappa in cube(radius))


def add_poly(left: Polynomial, right: Polynomial) -> Polynomial:
    out = dict(left)
    for exponent, value in right.items():
        out[exponent] = out.get(exponent, 0) + value
        if out[exponent] == 0:
            del out[exponent]
    return out


def scale_poly(poly: Polynomial, scalar: int) -> Polynomial:
    return {
        exponent: scalar * value
        for exponent, value in poly.items()
        if scalar * value
    }


def shift_poly(poly: Polynomial, shift: int) -> Polynomial:
    return {exponent + shift: value for exponent, value in poly.items()}


def multiply_poly(left: Polynomial, right: Polynomial) -> Polynomial:
    out: Polynomial = {}
    for i, a in left.items():
        for j, b in right.items():
            out[i + j] = out.get(i + j, 0) + a * b
    return {exponent: value for exponent, value in out.items() if value}


def one_minus_t_power(degree: int) -> Polynomial:
    return {
        exponent: (-1) ** exponent * comb(degree, exponent)
        for exponent in range(degree + 1)
    }


def q_polynomial(d_start: int, length: int) -> Polynomial:
    return {
        i: (
            (-1) ** i
            * comb(d_start - 1 + i, i)
            * comb(d_start + length, length - i)
        )
        for i in range(length + 1)
    }


def r_polynomial(d_start: int, length: int) -> Polynomial:
    out: Polynomial = {}
    for r in range(d_start):
        term = shift_poly(one_minus_t_power(d_start - 1 - r), r)
        out = add_poly(
            out,
            scale_poly(term, comb(d_start + length, r)),
        )
    return out


def pair_with_cube(
    m: int,
    radius: int,
    polynomial: Polynomial,
) -> int:
    return sum(
        value * coefficient(m, *scale(kappa, exponent))
        for kappa in cube(radius)
        for exponent, value in polynomial.items()
    )


def carrier(m: int, d_start: int, length: int) -> int:
    q_poly = q_polynomial(d_start, length)
    return sum(
        weight * shell(m, d_start - 1 + i)
        for i, weight in q_poly.items()
    )


def is_prime(value: int) -> bool:
    if value < 2:
        return False
    divisor = 2
    while divisor * divisor <= value:
        if value % divisor == 0:
            return value == divisor
        divisor += 1
    return True


def carrier_row(
    d_start: int,
    length: int,
    lower: int,
    upper: int,
) -> list[int]:
    row = [0] * (upper - lower + 1)
    for i, weight in q_polynomial(d_start, length).items():
        node = d_start - 1 + i
        row[node - lower] = weight
    return row


def audit_beta_packet() -> tuple[int, int, int]:
    polynomial_checks = 0
    packet_checks = 0
    shifted_checks = 0

    for m in range(4, 10):
        for d_start in range(3, m + 2):
            for length in range(1, min(4, d_start - 1) + 1):
                q_poly = q_polynomial(d_start, length)
                r_poly = r_polynomial(d_start, length)
                bezout = add_poly(
                    shift_poly(q_poly, d_start),
                    multiply_poly(
                        one_minus_t_power(length + 1),
                        r_poly,
                    ),
                )
                assert bezout == {0: 1}
                assert q_poly[0] == comb(d_start + length, length)
                assert sum(q_poly.values()) == 1
                polynomial_checks += 1

                radius = m // (d_start - 1)
                value = carrier(m, d_start, length)
                kernel_value = pair_with_cube(
                    m,
                    radius,
                    shift_poly(q_poly, d_start - 1),
                )
                assert value == kernel_value

                h_poly = shift_poly(
                    multiply_poly(
                        one_minus_t_power(length + 1),
                        r_poly,
                    ),
                    -1,
                )
                high_difference = pair_with_cube(m, radius, h_poly)
                near_origin = packet(m, radius)
                assert value + high_difference == near_origin

                # This catches the indexing error that identifies the
                # near-origin cube with the (D-1)-spaced shell.
                if m == 4 and d_start == 3 and length == 2:
                    assert near_origin == 1_826_539
                    assert shell(m, d_start - 1) == 320_000
                    assert near_origin != shell(m, d_start - 1)
                packet_checks += 1

                # The full T^{-h} family is an exact Laurent identity.
                for h in range(4):
                    left = pair_with_cube(
                        m,
                        radius,
                        shift_poly(q_poly, d_start - h),
                    )
                    companion = pair_with_cube(
                        m,
                        radius,
                        shift_poly(
                            multiply_poly(
                                one_minus_t_power(length + 1),
                                r_poly,
                            ),
                            -h,
                        ),
                    )
                    right = packet(m, radius, step=-h)
                    assert left + companion == right
                    shifted_checks += 1

                for q in range(d_start + 1, d_start + length + 1):
                    if not is_prime(q):
                        continue
                    j = q - d_start
                    for i in range(length + 1):
                        expected = 1 if i == j else 0
                        assert (q_poly.get(i, 0) - expected) % q == 0
                    assert value % q == shell(m, q - 1) % q

    return polynomial_checks, packet_checks, shifted_checks


def audit_same_radius_cancellation() -> tuple[int, int]:
    common_node_checks = 0
    universal_content_checks = 0

    for m in range(8, 30):
        candidates: list[tuple[int, int, int]] = []
        for d_start in range(3, m + 1):
            for length in range(1, min(5, d_start - 1) + 1):
                radius = m // (d_start - 1)
                candidates.append((d_start, length, radius))

        for index, (d1, n1, radius1) in enumerate(candidates):
            for d2, n2, radius2 in candidates[index + 1 :]:
                if radius1 != radius2:
                    continue
                lower = min(d1 - 1, d2 - 1)
                upper = max(d1 + n1 - 1, d2 + n2 - 1)
                row1 = carrier_row(d1, n1, lower, upper)
                row2 = carrier_row(d2, n2, lower, upper)

                interval_left = max(d1 + 1, d2 + 1)
                interval_right = min(d1 + n1, d2 + n2)
                for q in range(interval_left, interval_right + 1):
                    if not is_prime(q):
                        continue
                    node_index = q - 1 - lower
                    for row in (row1, row2):
                        for position, value in enumerate(row):
                            expected = 1 if position == node_index else 0
                            assert (value - expected) % q == 0
                    common_node_checks += 1

                    # Packet cancellation has coefficient sum 1-1=0.
                    # Its q-factor is therefore present for every input
                    # sequence, not only when the marked shell vanishes.
                    difference = [
                        left - right for left, right in zip(row1, row2)
                    ]
                    assert all(value % q == 0 for value in difference)
                    universal_content_checks += 1

    return common_node_checks, universal_content_checks


def main() -> None:
    polynomial, packet_count, shifted = audit_beta_packet()
    common, universal = audit_same_radius_cancellation()
    print("Q32_BETA_PADE_PACKET_AUDIT=PASS")
    print("POLYNOMIAL_BEZOUT_CHECKS", polynomial)
    print("EXACT_PACKET_CHECKS", packet_count)
    print("SHIFTED_PACKET_CHECKS", shifted)
    print("COMMON_NODE_CHECKS", common)
    print("UNIVERSAL_DIFFERENCE_CONTENT_CHECKS", universal)


if __name__ == "__main__":
    main()
