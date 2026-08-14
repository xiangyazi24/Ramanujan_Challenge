#!/usr/bin/env python3
"""Exact stdlib verifier for the Q8375 intermediate-cube arithmetic audit.

This script checks only algebraic identities used in the report:

* direct/reflected fixed-quotient affine prime-cube formulas;
* barycentric weights and R/p_i == V_i (mod p_i);
* the exact 8-node Vandermonde determinant det = -8 B Delta;
* full Boolean interpolation dimension;
* CRT freedom for local quotient/boundary residues;
* solvability of a triangular first-divided local lift;
* aggregate exponent bookkeeping.

Synthetic residues below are witnesses for algebraic surjectivity only.  They
are not claimed to be actual Apery rows or arithmetic counterexamples.
"""

from __future__ import annotations

from fractions import Fraction
from math import gcd, prod


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2
    return True


def inv_mod(a: int, m: int) -> int:
    return pow(a % m, -1, m)


def crt(residues: list[int], moduli: list[int]) -> int:
    assert len(residues) == len(moduli) and residues
    for i, left in enumerate(moduli):
        for right in moduli[i + 1 :]:
            assert gcd(left, right) == 1
    modulus = prod(moduli)
    answer = 0
    for residue, local_modulus in zip(residues, moduli):
        complement = modulus // local_modulus
        answer += (
            residue
            * complement
            * inv_mod(complement, local_modulus)
        )
    return answer % modulus


def masks() -> range:
    return range(8)


def subset_sum(mask: int, steps: tuple[int, int, int]) -> int:
    return sum(steps[j] for j in range(3) if mask & (1 << j))


def quotient_cube(
    a: int,
    eps: int,
    m: int,
    u: int,
    steps: tuple[int, int, int],
) -> tuple[int, ...]:
    assert eps in (0, 1)
    c = a + eps
    tau = 2 * eps - 1
    shifted_m = m + eps
    assert c > 0
    assert all(step % c == 0 for step in steps)
    assert (shifted_m + tau * u) % c == 0
    base_prime = (shifted_m + tau * u) // c
    scaled_steps = tuple(step // c for step in steps)

    values = []
    for mask in masks():
        h = u + subset_sum(mask, steps)
        numerator = shifted_m + tau * h
        assert numerator % c == 0
        p = numerator // c
        expected = base_prime + tau * subset_sum(mask, scaled_steps)
        assert p == expected
        # The examples are chosen inside one ordinary quotient cell.
        assert m // p == a
        values.append(p)
    return tuple(values)


def barycentric_v(nodes: tuple[int, ...], i: int) -> int:
    return prod(nodes[j] - nodes[i] for j in range(len(nodes)) if j != i)


def vandermonde(nodes: tuple[int, ...]) -> int:
    return prod(
        nodes[j] - nodes[i]
        for i in range(len(nodes))
        for j in range(i + 1, len(nodes))
    )


def det_bareiss(matrix: list[list[int]]) -> int:
    """Fraction-free determinant with exact divisions."""

    n = len(matrix)
    assert n and all(len(row) == n for row in matrix)
    a = [row[:] for row in matrix]
    sign = 1
    previous = 1

    for k in range(n - 1):
        if a[k][k] == 0:
            swap = next((r for r in range(k + 1, n) if a[r][k]), None)
            if swap is None:
                return 0
            a[k], a[swap] = a[swap], a[k]
            sign *= -1
        pivot = a[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                numerator = a[i][j] * pivot - a[i][k] * a[k][j]
                assert numerator % previous == 0
                a[i][j] = numerator // previous
            a[i][k] = 0
        previous = pivot
    return sign * a[-1][-1]


def interpolation_leading_coefficient(
    nodes: tuple[int, ...], values: tuple[int, ...]
) -> Fraction:
    assert len(nodes) == len(values)
    total = Fraction(0)
    for i, (node, value) in enumerate(zip(nodes, values)):
        derivative = prod(
            node - nodes[j] for j in range(len(nodes)) if j != i
        )
        total += Fraction(value, derivative)
    return total


def submasks(mask: int) -> list[int]:
    answer = []
    current = mask
    while True:
        answer.append(current)
        if current == 0:
            return answer
        current = (current - 1) & mask


def boolean_coefficients(values: tuple[int, ...]) -> tuple[int, ...]:
    """Möbius coefficients of the unique multilinear polynomial on {0,1}^3."""

    assert len(values) == 8
    coefficients = []
    for mask in masks():
        total = 0
        for submask in submasks(mask):
            parity = (mask.bit_count() - submask.bit_count()) & 1
            total += (-1 if parity else 1) * values[submask]
        coefficients.append(total)
    return tuple(coefficients)


def boolean_reconstruct(coefficients: tuple[int, ...], mask: int) -> int:
    return sum(coefficients[submask] for submask in submasks(mask))


def exponent_sub(
    left: tuple[int, int, int], right: tuple[int, int, int]
) -> tuple[int, int, int]:
    return tuple(a - b for a, b in zip(left, right))


def exponent_scale(
    value: tuple[int, int, int], factor: int
) -> tuple[int, int, int]:
    return tuple(factor * entry for entry in value)


def main() -> None:
    # Unified quotient geometry.  These are exact integer examples of a fixed
    # ordinary quotient cell; primality is not needed for this formula check.
    direct = quotient_cube(
        a=2,
        eps=0,
        m=404,
        u=4,
        steps=(12, 24, 60),
    )
    reflected = quotient_cube(
        a=2,
        eps=1,
        m=595,
        u=4,
        steps=(18, 36, 90),
    )
    assert direct == (200, 194, 188, 182, 170, 164, 158, 152)
    assert reflected == (200, 206, 212, 218, 230, 236, 242, 248)

    # An exact affine 3-cube all of whose nodes are prime.  This is used only
    # for the CRT/barycentric algebra, not as an asserted legal Apery cube.
    prime_nodes = tuple(11 + subset_sum(mask, (6, 12, 30)) for mask in masks())
    assert prime_nodes == (11, 17, 23, 29, 41, 47, 53, 59)
    assert all(is_prime(p) for p in prime_nodes)
    assert len(set(prime_nodes)) == 8

    R = prod(prime_nodes)
    weights = tuple(barycentric_v(prime_nodes, i) for i in masks())
    for p, weight in zip(prime_nodes, weights):
        assert weight % p != 0
        assert (R // p - weight) % p == 0

    # Canonical interpolation determinant: the span-controlled Vandermonde
    # factor survives, but the free global quotient B remains as -8 B.
    B = 37
    y_values = tuple(B * weight for weight in weights)
    leading = interpolation_leading_coefficient(prime_nodes, y_values)
    assert leading == -8 * B

    matrix = [
        [1, p, p**2, p**3, p**4, p**5, p**6, y]
        for p, y in zip(prime_nodes, y_values)
    ]
    delta = vandermonde(prime_nodes)
    assert det_bareiss(matrix) == -8 * B * delta

    H_over_c = 6 + 12 + 30
    assert all(abs(weight) <= H_over_c**7 for weight in weights)
    assert abs(delta) <= H_over_c**28

    # Ordinary Boolean interpolation has all eight degrees of freedom.
    arbitrary_cube_values = (3, 1, 4, 1, 5, 9, 2, 6)
    coefficients = boolean_coefficients(arbitrary_cube_values)
    reconstructed = tuple(
        boolean_reconstruct(coefficients, mask) for mask in masks()
    )
    assert reconstructed == arbitrary_cube_values

    # CRT makes one global B mod R equivalent to arbitrary independent B_i.
    local_B = tuple((i * i + 3 * i + 7) % p for i, p in enumerate(prime_nodes))
    global_B = crt(list(local_B), list(prime_nodes))
    assert all(global_B % p == residue for p, residue in zip(prime_nodes, local_B))

    # The same is true for a single integral boundary state modulo R^2.
    square_moduli = [p * p for p in prime_nodes]
    local_x = [
        (17 * i + 11) % modulus
        for i, modulus in enumerate(square_moduli)
    ]
    local_y = [
        (29 * i * i + 5) % modulus
        for i, modulus in enumerate(square_moduli)
    ]
    global_x = crt(local_x, square_moduli)
    global_y = crt(local_y, square_moduli)
    assert all(global_x % modulus == value for modulus, value in zip(square_moduli, local_x))
    assert all(global_y % modulus == value for modulus, value in zip(square_moduli, local_y))

    # Model the proved first-divided structure B_i V_i = beta z_i + Gamma_i.
    # beta=5 is the direct quotient-one coefficient.  Gamma_i are arbitrary
    # local constants here: the point is exact solvability for every B_i.
    beta = 5
    divided_digits = []
    for i, (p, weight, b_residue) in enumerate(
        zip(prime_nodes, weights, local_B)
    ):
        gamma = (i**3 + 2 * i + 1) % p
        z = (
            (b_residue * (weight % p) - gamma)
            * inv_mod(beta, p)
        ) % p
        assert (beta * z + gamma - b_residue * weight) % p == 0
        divided_digits.append(z)
    assert len(divided_digits) == 8

    # Adding p_i t_i to a lifted nodal value leaves arbitrary independent
    # cofactor directions; there is no common-field interpolation constraint.
    lift_multipliers = tuple(2 * i + 1 for i in masks())
    lifted_values = tuple(
        B * weight + p * t
        for p, weight, t in zip(prime_nodes, weights, lift_multipliers)
    )
    assert all(
        (lifted - B * weight) % p == 0
        for lifted, weight, p in zip(lifted_values, weights, prime_nodes)
    )
    assert any(lifted != B * weight for lifted, weight in zip(lifted_values, weights))

    # Exponent bookkeeping in the order (power of L, power of X, power of log X).
    forced = (15, -14, -14)
    cells = (-1, 1, 2)  # X log^2 X / L
    pointwise = exponent_sub(forced, cells)
    second_moment = exponent_sub(exponent_scale(forced, 2), cells)
    assert pointwise == (16, -15, -16)
    assert second_moment == (31, -29, -30)

    print("PASS: direct/reflected fixed-cell affine cube formulas")
    print("PASS: barycentric R/p_i congruences and span bounds")
    print("PASS: det = -8*B*Vandermonde exactly")
    print("PASS: full Boolean interpolation has eight degrees of freedom")
    print("PASS: CRT boundary/B residues are componentwise free")
    print("PASS: triangular first-divided lifts solve for local digits")
    print("PASS: aggregate exponent bookkeeping")
    print("Q8375 INTERMEDIATE CUBE ARITHMETIC VERIFIER: PASS")


if __name__ == "__main__":
    main()
