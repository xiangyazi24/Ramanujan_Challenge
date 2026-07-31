#!/usr/bin/env python3
"""Audit the terminal Bernstein transform and its constant-tail obstruction."""

from functools import lru_cache
from math import comb


def C(n, k):
    return comb(n, k) if 0 <= k <= n else 0


# Representatives and multiplicities under interchange of y and z.
RAY_CLASSES = (
    ((-1, -1, -1), 1),
    ((-1, -1, 0), 2),
    ((-1, -1, 1), 2),
    ((-1, 0, 0), 1),
    ((-1, 0, 1), 2),
    ((-1, 1, 1), 1),
    ((0, -1, -1), 1),
    ((0, -1, 0), 2),
    ((0, -1, 1), 2),
    ((0, 0, 1), 2),
    ((0, 1, 1), 1),
    ((1, 0, 0), 1),
    ((1, 0, 1), 2),
    ((1, 1, 1), 1),
)
assert sum(multiplicity for _, multiplicity in RAY_CLASSES) == 21


@lru_cache(maxsize=None)
def apery(moment):
    return sum(
        C(moment, index) ** 2
        * C(moment + index, index) ** 2
        for index in range(moment + 1)
    )


@lru_cache(maxsize=None)
def ray_value(moment, residue, point):
    """Return c_M((M-residue) point) from the one-fold formula."""

    u, v, w = point
    node = moment - residue
    return sum(
        C(moment, index)
        * C(moment, index - node * u)
        * C(2 * moment - index, moment - node * v)
        * C(2 * moment - index, moment - node * w)
        for index in range(moment + 1)
    )


@lru_cache(maxsize=None)
def ray_correction(moment, residue):
    return sum(
        multiplicity * ray_value(moment, residue, point)
        for point, multiplicity in RAY_CLASSES
    )


@lru_cache(maxsize=None)
def shell_fast(moment, node):
    quotient = moment // node
    out = 0
    for index in range(moment + 1):
        base = moment - index
        x_packet = sum(
            C(moment, base + node * u)
            for u in range(-quotient, quotient + 1)
        )
        yz_packet = sum(
            C(2 * moment - index, base + node * v)
            for v in range(-quotient, quotient + 1)
        )
        out += C(moment, index) * x_packet * yz_packet**2
    return out


def terminal_carrier(moment, order):
    node = moment - order
    return sum(
        (-1) ** index
        * C(node + index, index)
        * C(moment + 1, order - index)
        * shell_fast(moment, node + index)
        for index in range(order + 1)
    )


def terminal_transform(moment, order):
    return apery(moment) + sum(
        (-1) ** (order - residue)
        * C(moment - residue, order - residue)
        * C(moment + 1, residue)
        * ray_correction(moment, residue)
        for residue in range(order + 1)
    )


def boundary_from_rays(moment, order):
    out = apery(moment) if order == 0 else 0
    return out + sum(
        multiplicity
        * (-1) ** residue
        * C(order, residue)
        * ray_value(moment, residue, point)
        for point, multiplicity in RAY_CLASSES
        for residue in range(order + 1)
    )


def add_bernstein_term(poly, coefficient, residue, power):
    degree = residue + power
    if len(poly) <= degree:
        poly.extend([0] * (degree + 1 - len(poly)))
    for index in range(power + 1):
        poly[residue + index] += (
            coefficient * (-1) ** index * C(power, index)
        )


def normalized(poly):
    out = list(poly) or [0]
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def multiply_by_one_minus_z(poly):
    out = [0] * (len(poly) + 1)
    for index, value in enumerate(poly):
        out[index] += value
        out[index + 1] -= value
    return normalized(out)


def E_polynomial(moment):
    out = [0] * (moment + 1)
    for residue in range(moment + 1):
        add_bernstein_term(
            out,
            C(moment + 1, residue)
            * ray_correction(moment, residue),
            residue,
            moment - residue,
        )
    return out


def N_polynomial(moment):
    out = [0] * (moment + 2)
    for point, multiplicity in RAY_CLASSES:
        for residue in range(moment + 2):
            add_bernstein_term(
                out,
                multiplicity
                * C(moment + 1, residue)
                * ray_value(moment, residue, point),
                residue,
                moment + 1 - residue,
            )
    return out


def audit(max_moment=16):
    checks = 0
    for moment in range(2, max_moment + 1):
        last = (moment - 1) // 2
        for residue in range(last + 1):
            assert (
                shell_fast(moment, moment - residue)
                == apery(moment) + ray_correction(moment, residue)
            )
            checks += 1

        for order in range(last + 1):
            assert terminal_carrier(moment, order) == terminal_transform(
                moment, order
            )
            checks += 1
            if order:
                assert (
                    terminal_carrier(moment, order)
                    - terminal_carrier(moment, order - 1)
                    == (-1) ** order
                    * C(moment + 1, order)
                    * boundary_from_rays(moment, order)
                )
                checks += 1

        polynomial_e = E_polynomial(moment)
        polynomial_n = N_polynomial(moment)
        tail = sum(
            multiplicity * ray_value(moment, moment + 1, point)
            for point, multiplicity in RAY_CLASSES
        )
        assert sum(polynomial_n) == tail
        right = list(polynomial_n)
        right[moment + 1] -= tail
        assert multiply_by_one_minus_z(polynomial_e) == normalized(right)
        checks += 1

        running = 0
        for order in range(moment + 5):
            if order < len(polynomial_n):
                running += polynomial_n[order]
            expected = (
                polynomial_e[order] if order <= moment else tail
            )
            assert running == expected
            checks += 1

        assert tail >= 4 * moment * 5 ** (moment - 1)
        checks += 1
    return checks


if __name__ == "__main__":
    print("TERMINAL_BERNSTEIN_CHECKS", audit())
    print("Q32_TERMINAL_BERNSTEIN_AUDIT=PASS")
