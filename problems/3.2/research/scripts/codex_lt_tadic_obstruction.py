#!/usr/bin/env python3
"""Finite guards for the two elementary obstructions to the T-adic route."""

from math import comb


def check_no_zp_interpolation() -> None:
    """The integer function n -> zeta^n is not Z_p-continuous for zeta != 1."""
    for prime in (5, 7, 13, 17, 29):
        tame_order = prime - 1
        # In exponent notation choose zeta itself, a primitive tame root.
        for exponent in range(1, 9):
            assert pow(prime, exponent, tame_order) == 1
        # Thus p^j -> 0 in Z_p, while zeta^(p^j)=zeta != 1=zeta^0.
    print(
        "VERIFIED tame exponent obstruction: p^j=1 mod (p-1), so "
        "zeta^(p^j)=zeta != 1 although p^j tends p-adically to 0"
    )


def check_legendre_hasse_degree() -> None:
    """The Legendre Hasse polynomial has degree (p-1)/2, not O(1)."""
    for prime in (5, 7, 11, 13, 17, 19, 23, 29, 31):
        half = (prime - 1) // 2
        coefficients = [comb(half, index) ** 2 % prime for index in range(half + 1)]
        assert coefficients[-1] == 1
        assert len(coefficients) - 1 == half
    print(
        "VERIFIED Legendre Hasse polynomial degree=(p-1)/2 for "
        "p=5,7,11,13,17,19,23,29,31"
    )


def main() -> None:
    check_no_zp_interpolation()
    check_legendre_hasse_degree()


if __name__ == "__main__":
    main()
