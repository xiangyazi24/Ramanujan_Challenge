#!/usr/bin/env python3
"""Verify the carry-free Gamma/Jacobi reconstruction of the rank-two branches.

For the relevant branch degree D, all binomial factors in the
Lagrange--Buermann formula lie in a carry-free range.  They can therefore be
written both as Morita-Gamma ratios modulo p and as reductions of Jacobi sums.
The number of summands still grows with the coefficient index; this script
does not assert an O(1)-term Gross--Koblitz formula.

The two requested tests are exhaustive:

* p=13: sigma_0,...,sigma_5 and b_0,...,b_12;
* p=29: tau_0,...,tau_14 and b_0,...,b_28.
"""

from __future__ import annotations

from fractions import Fraction as Q
from functools import lru_cache
from math import comb


PRIMES = (13, 29)


def legendre(value: int, prime: int) -> int:
    residue = pow(value % prime, (prime - 1) // 2, prime)
    return -1 if residue == prime - 1 else residue


def fraction_residue(value: Q, prime: int) -> int:
    return value.numerator % prime * pow(value.denominator % prime, -1, prime) % prime


def gamma_p_mod(value: Q, prime: int) -> int:
    """Morita Gamma_p(value) modulo p, using local constancy modulo p."""

    residue = fraction_residue(value, prime)
    if residue == 0:
        return 1
    return (-1 if residue & 1 else 1) * factorial_mod(residue - 1, prime) % prime


def factorial_mod(number: int, prime: int) -> int:
    result = 1
    for factor in range(2, number + 1):
        result = result * factor % prime
    return result


def gamma_binomial(alpha: Q, lower: int, prime: int) -> int:
    """Return binom(alpha,lower) via Gamma_p in the no-crossing range."""

    assert 0 <= lower < prime
    assert all(fraction_residue(alpha - offset, prime) for offset in range(lower))
    numerator = -gamma_p_mod(alpha + 1, prime)
    denominator = gamma_p_mod(alpha - lower + 1, prime) * gamma_p_mod(
        Q(lower + 1), prime
    )
    return numerator * pow(denominator % prime, -1, prime) % prime


def jacobi_binomial(top: int, lower: int, prime: int) -> int:
    """Return binom(top,lower) as a residual Teichmueller Jacobi sum.

    Jbar_p(lower,top) = sum_x x^(p-1-lower) (1-x)^top, so
    binom(top,lower) = (-1)^(lower+1) Jbar_p(lower,top) modulo p.
    """

    assert 0 <= lower <= top < prime - 1
    if top == 0:
        # Avoid the convention-dependent J(epsilon,epsilon) endpoint.
        return 1
    jacobi = sum(
        pow(x, prime - 1 - lower, prime) * pow((1 - x) % prime, top, prime)
        for x in range(prime)
    ) % prime
    return ((-1 if lower % 2 == 0 else 1) * jacobi) % prime


def gross_koblitz_binomial(top: int, lower: int, prime: int) -> int:
    """The unit Gross--Koblitz Gamma quotient for an interior binomial.

    The endpoint Jacobi sums have a trivial character, so they are handled
    separately.  For 0 < lower < top < p-1, Gross--Koblitz applied to
    J(omega^(-lower), omega^top) gives the displayed three-Gamma quotient.
    """

    assert 0 <= lower <= top < prime - 1
    if lower in (0, top):
        return 1
    sign = -1 if lower & 1 else 1
    numerator = (
        sign
        * gamma_p_mod(Q(lower, prime - 1), prime)
        * gamma_p_mod(Q(prime - 1 - top, prime - 1), prime)
    )
    denominator = gamma_p_mod(Q(prime - 1 - top + lower, prime - 1), prime)
    return numerator * pow(denominator, -1, prime) % prime


class CarryFreeAtoms:
    def __init__(self, prime: int):
        self.prime = prime
        self.checked: set[tuple[Q, int]] = set()

    @lru_cache(maxsize=None)
    def choose(self, alpha: Q, lower: int) -> int:
        top = fraction_residue(alpha, self.prime)
        assert lower <= top < self.prime - 1
        gamma_value = gamma_binomial(alpha, lower, self.prime)
        jacobi_value = jacobi_binomial(top, lower, self.prime)
        gross_koblitz_value = gross_koblitz_binomial(top, lower, self.prime)
        direct = 1
        for offset in range(lower):
            direct = direct * fraction_residue(alpha - offset, self.prime) % self.prime
            direct = direct * pow(offset + 1, -1, self.prime) % self.prime
        assert gamma_value == jacobi_value == gross_koblitz_value == direct
        self.checked.add((alpha, lower))
        return gamma_value


def direct_branch(branch: str, degree: int, prime: int) -> list[int]:
    if branch == "tau":
        values = [1, 5 * pow(2, -1, prime) % prime]
    else:
        values = [1, 39 * pow(2, -1, prime) % prime]
    for index in range(1, degree):
        if branch == "tau":
            numerator = (
                2 * (68 * index * index + 34 * index + 5) * values[index]
                - (2 * index - 1) ** 2 * values[index - 1]
            )
        else:
            numerator = (
                2 * (68 * index * index + 102 * index + 39) * values[index]
                - (2 * index + 1) ** 2 * values[index - 1]
            )
        denominator = 4 * (index + 1) ** 2
        values.append(numerator * pow(denominator % prime, -1, prime) % prime)
    return values[: degree + 1]


def lagrange_branch(branch: str, degree: int, prime: int, atoms: CarryFreeAtoms) -> list[int]:
    choose = atoms.choose

    franel = []
    for index in range(degree + 1):
        franel.append(
            sum(pow(choose(Q(index), lower), 3, prime) for lower in range(index + 1))
            % prime
        )

    g_values = []
    for total in range(degree + 1):
        if branch == "tau":
            value = sum(
                franel[index] * choose(Q(1, 2), total - index)
                for index in range(total + 1)
            )
        else:
            value = 0
            for index in range(total + 1):
                for half_index in range(total - index + 1):
                    remainder = total - index - half_index
                    for quadratic_index in range(remainder // 2 + 1):
                        linear_index = remainder - 2 * quadratic_index
                        value += (
                            franel[index]
                            * choose(Q(3, 2), half_index)
                            * choose(Q(linear_index + quadratic_index), linear_index)
                            * pow(16, linear_index, prime)
                            * pow(8, quadratic_index, prime)
                        )
        g_values.append(value % prime)

    def phi_coefficient(index: int, power_index: int) -> int:
        return sum(
            choose(Q(index), numerator_index)
            * choose(
                Q(index + power_index - numerator_index - 1),
                power_index - numerator_index,
            )
            * pow(8, power_index - numerator_index, prime)
            for numerator_index in range(min(index, power_index) + 1)
        ) % prime

    values = [1]
    for index in range(1, degree + 1):
        value = sum(
            inner * g_values[inner] * phi_coefficient(index, index - inner)
            for inner in range(1, index + 1)
        )
        values.append(value * pow(index, -1, prime) % prime)
    return values


def apery_mod(index: int, prime: int) -> int:
    return sum(
        comb(index, lower) ** 2 * comb(index + lower, lower) ** 2
        for lower in range(index + 1)
    ) % prime


def convolution(values: list[int], index: int, prime: int) -> int:
    degree = len(values) - 1
    return sum(
        values[left] * values[index - left]
        for left in range(max(0, index - degree), min(degree, index) + 1)
    ) % prime


def check_prime(prime: int) -> None:
    character = legendre(-6, prime)
    branch = "tau" if character == 1 else "sigma"
    degree = (prime - 1) // 2 if branch == "tau" else (prime - 3) // 2
    atoms = CarryFreeAtoms(prime)

    direct = direct_branch(branch, degree, prime)
    reconstructed = lagrange_branch(branch, degree, prime, atoms)
    assert reconstructed == direct

    branch_square = [convolution(reconstructed, index, prime) for index in range(prime)]
    recovered_apery = []
    for index in range(prime):
        value = branch_square[index]
        if branch == "sigma":
            if index >= 1:
                value -= 34 * branch_square[index - 1]
            if index >= 2:
                value += branch_square[index - 2]
        recovered_apery.append(value % prime)
    expected_apery = [apery_mod(index, prime) for index in range(prime)]
    assert recovered_apery == expected_apery

    print(
        f"p={prime}: chi={character:+d}, branch={branch}, D={degree}, "
        f"all {degree + 1} branch coefficients and all {prime} b_r verified"
    )
    print(f"  branch residues: {reconstructed}")
    print(f"  distinct carry-free Gamma/Jacobi atoms checked: {len(atoms.checked)}")


def main() -> None:
    for prime in PRIMES:
        check_prime(prime)
    print("quadratic Jacobi convolution agrees with the direct Apéry sum: VERIFIED")


if __name__ == "__main__":
    main()
