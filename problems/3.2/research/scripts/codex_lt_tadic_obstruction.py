#!/usr/bin/env python3
"""Finite guards for elementary obstructions to the T-adic route.

These checks do not replace the cited T-adic literature.  They guard the
finite-etale tame parameter space, the failure of a Z_p interpolation of its
indices, and the extra hypotheses needed to identify trace divisibility with
non-ordinarity.
"""

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


def check_split_tame_parameter_space() -> None:
    """Check that U^(p-1)-1 is the product of p-1 distinct linear factors."""
    for prime in (5, 7, 13, 17, 29):
        roots = list(range(1, prime))
        assert len({root % prime for root in roots}) == prime - 1
        for value in range(prime):
            left = (pow(value, prime - 1, prime) - 1) % prime
            right = 1
            for root in roots:
                right = right * (value - root) % prime
            assert left == right

        # Delta functions are arbitrary regular functions on the split
        # zero-dimensional fiber.  This is the finite-field shadow of
        # Z_p[U]/(U^(p-1)-1) being a product of p-1 copies of Z_p.
        for root in roots:
            evaluations = [
                (1 - pow(value - root, prime - 1, prime)) % prime
                for value in roots
            ]
            assert evaluations == [int(value == root) for value in roots]
    print(
        "VERIFIED tame character special fiber splits into p-1 independent points "
        "at p=5,7,13,17,29"
    )


def check_trace_ordinarity_caveats() -> None:
    """Give exact rank-two counterexamples to the trace-zero slogan."""
    for prime in (5, 7, 13, 17, 29):
        # Two slope-zero eigenvalues may cancel although the object is
        # ordinary for Hodge slopes {0,0}.
        unit_roots = (1, -1)
        assert sum(unit_roots) % prime == 0
        assert all(root % prime != 0 for root in unit_roots)

        # If the Hodge slopes are {1,2}, the ordinary eigenvalues p,p^2
        # already have trace divisible by p.
        positive_roots = (prime, prime**2)
        assert sum(positive_roots) % prime == 0
        assert positive_roots[0] * positive_roots[1] == prime**3

        # With the additional unique-unit-root normalization {0,3}, an
        # ordinary model has unit trace.  This is the narrow setting in
        # which trace divisibility can detect loss of the unit root.
        unique_unit_roots = (1, prime**3)
        assert sum(unique_unit_roots) % prime == 1
        assert unique_unit_roots[0] * unique_unit_roots[1] == prime**3
    print(
        "VERIFIED trace mod p is not a general ordinarity test; the "
        "unique-unit-root Hodge/determinant hypotheses are load-bearing"
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
    check_split_tame_parameter_space()
    check_trace_ordinarity_caveats()
    check_legendre_hasse_degree()


if __name__ == "__main__":
    main()
