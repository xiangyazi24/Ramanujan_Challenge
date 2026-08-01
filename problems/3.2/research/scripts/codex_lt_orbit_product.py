#!/usr/bin/env python3
"""Exact checks behind the Galois-orbit and Parseval arguments.

The cyclotomic computations use the genuine integer Frobenius traces of the
Franel elliptic pencil

    E_u: y^2 + (1-2u)xy + u^2 y = x^3.

Only M is needed to test the Galois mechanism.  The same computation applies
to M-T whenever the pointwise amplitudes defining T are rational integers.
"""

from math import gcd

import sympy as sp


def legendre(value: int, prime: int) -> int:
    value %= prime
    if value == 0:
        return 0
    return 1 if pow(value, (prime - 1) // 2, prime) == 1 else -1


def elliptic_trace(parameter: int, prime: int) -> int:
    """Return a_p for the smooth fiber E_parameter over F_p."""
    a1 = (1 - 2 * parameter) % prime
    a3 = parameter * parameter % prime
    points = 1
    for x_value in range(prime):
        linear = (a1 * x_value + a3) % prime
        discriminant = (linear * linear + 4 * x_value**3) % prime
        points += 1 + legendre(discriminant, prime)
    return prime + 1 - points


def primitive_root(prime: int) -> int:
    return int(sp.primitive_root(prime))


def fiber_amplitudes(prime: int) -> tuple[dict[int, int], list[int]]:
    """Return B_t=sum_{phi(x)=t} a_{p,x}^2 and the individual a_x^2."""
    singular = {0, prime - 1, pow(8, -1, prime)}
    grouped = {value: 0 for value in range(1, prime)}
    individual = []
    for parameter in range(prime):
        if parameter in singular:
            continue
        image = (
            parameter
            * (1 - 8 * parameter)
            * pow(1 + parameter, -1, prime)
        ) % prime
        if image == 0:
            continue
        square_trace = elliptic_trace(parameter, prime) ** 2
        grouped[image] += square_trace
        individual.append(square_trace)
    return grouped, individual


def canonical(poly: sp.Expr, variable: sp.Symbol, cyclotomic: sp.Poly) -> sp.Poly:
    return sp.rem(sp.Poly(sp.expand(poly), variable, domain=sp.ZZ), cyclotomic)


def mellin_polynomials(prime: int) -> tuple[list[sp.Poly], dict[int, int], int, sp.Symbol, sp.Poly]:
    order = prime - 1
    generator = primitive_root(prime)
    logarithm = {pow(generator, exponent, prime): exponent for exponent in range(order)}
    grouped, _ = fiber_amplitudes(prime)
    zeta = sp.Symbol("z")
    cyclotomic = sp.Poly(sp.cyclotomic_poly(order, zeta), zeta, domain=sp.ZZ)
    values = []
    for r in range(order):
        expression = sum(
            amplitude * zeta ** ((-r * logarithm[t]) % order)
            for t, amplitude in grouped.items()
        )
        values.append(canonical(expression, zeta, cyclotomic))
    return values, grouped, generator, zeta, cyclotomic


def orbit(index: int, order: int) -> list[int]:
    return sorted({(unit * index) % order for unit in range(order) if gcd(unit, order) == 1})


def check_prime(prime: int) -> None:
    order = prime - 1
    values, grouped, generator, zeta, cyclotomic = mellin_polynomials(prime)

    # sigma_a(zeta)=zeta^a sends M_r to M_{ar}; rational trace amplitudes stay fixed.
    for unit in range(order):
        if gcd(unit, order) != 1:
            continue
        for r in range(order):
            conjugate = canonical(values[r].as_expr().subs(zeta, zeta**unit), zeta, cyclotomic)
            assert conjugate == values[(unit * r) % order]

    # The orbit product is fixed by the full Galois group and hence is an integer.
    for r in range(1, order):
        product = sp.Poly(1, zeta, domain=sp.ZZ)
        for s in orbit(r, order):
            product = canonical(product.as_expr() * values[s].as_expr(), zeta, cyclotomic)
        assert product.degree() <= 0

        # At the chosen split prime, zeta -> generator.  k hits at this same
        # prime give at least k factors of p in the rational orbit product.
        residues = [int(values[s].eval(generator)) % prime for s in orbit(r, order)]
        hits = sum(residue == 0 for residue in residues)
        integer_product = int(product.nth(0))
        if integer_product:
            quotient = integer_product
            valuation = 0
            while quotient % prime == 0:
                valuation += 1
                quotient //= prime
            assert valuation >= hits

    # Exact cyclotomic Parseval, including the r=0 character.
    parseval_left = sp.Poly(0, zeta, domain=sp.ZZ)
    for r in range(order):
        term = values[r].as_expr() * values[(-r) % order].as_expr()
        parseval_left = canonical(parseval_left.as_expr() + term, zeta, cyclotomic)
    grouped_right = order * sum(amplitude * amplitude for amplitude in grouped.values())
    assert parseval_left == sp.Poly(grouped_right, zeta, domain=sp.ZZ)

    # The ungrouped formula drops cross terms inside the quadratic phi-fibers.
    _, individual = fiber_amplitudes(prime)
    ungrouped_right = order * sum(amplitude * amplitude for amplitude in individual)
    assert grouped_right != ungrouped_right

    print(
        f"VERIFIED p={prime}: exact cyclotomic Galois action, integral orbit "
        "products, and p-adic hit transfer"
    )
    print(
        f"VERIFIED p={prime}: grouped Parseval identity; ungrouped RHS differs "
        f"by {grouped_right - ungrouped_right}"
    )


def main() -> None:
    for prime in (13, 17):
        check_prime(prime)


if __name__ == "__main__":
    main()
