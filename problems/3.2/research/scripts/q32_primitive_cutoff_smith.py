#!/usr/bin/env python3
"""Audit primitive all-cutoff Smith divisors.

The raw all-cutoff matrix has rows

    h_J(c) = T_(n,J)(c) - T_(n,J-1)(c)
           = K_(n,J)(c) g_J(c).

Most of every maximal minor is a universal factorial/binomial factor.  Since
g_J is primitive (its top coefficient is +/-F_0=+/-1), Gauss's lemma says
that content(h_J)=content(K_(n,J)).  Divide each row by this exact universal
content and take the top determinantal divisor of the resulting primitive
row matrix.

This is a stronger normalization than the consecutive-minor audit in
q32_all_cutoff_determinant.py.  It produces surprisingly small numbers at
first, and every aligned q=1 bad prime in the tested range survives.
However it is not selective: good q=1 candidates and even primes above n
also survive.  The explicit counterexamples printed below close the hope
that row-content saturation alone isolates the target radical.

SymPy is used only to compute exact Smith normal forms in this finite audit.
"""

from __future__ import annotations

from math import gcd, prod

from sympy import Matrix, ZZ, factorint
from sympy.matrices.normalforms import smith_normal_form

from q32_all_cutoff_determinant import difference_rows
from q32_newton import apery_numbers
from q32_strehl_gcd import franel_numbers, primes_up_to


LIMIT = 42


def primitive_rows(n: int, franel: list[int]) -> list[list[int]]:
    result: list[list[int]] = []
    for row in difference_rows(n, franel):
        row_content = 0
        for entry in row:
            row_content = gcd(row_content, entry)
        assert row_content
        result.append([entry // row_content for entry in row])
    return result


def top_determinantal_divisor(rows: list[list[int]]) -> int:
    smith = smith_normal_form(Matrix(rows), domain=ZZ)
    rank = len(rows)
    invariants = [abs(int(smith[index, index])) for index in range(rank)]
    assert all(invariants)
    return prod(invariants)


def main() -> None:
    apery = apery_numbers(LIMIT)
    franel = franel_numbers(LIMIT)
    primes = primes_up_to(2 * LIMIT)
    counterexamples: dict[int, list[int]] = {}

    for n in range(6, LIMIT + 1):
        divisor = top_determinantal_divisor(primitive_rows(n, franel))
        target = 1
        for prime in primes:
            if not n / 2 < prime <= n:
                continue
            folded = min(n - prime, 2 * prime - 1 - n)
            if apery[folded] % prime == 0:
                target *= prime
        assert divisor % target == 0, (n, divisor, target)

        contaminating = []
        for prime in factorint(divisor):
            if prime <= n / 2:
                continue
            aligned_bad = False
            if prime <= n:
                folded = min(n - prime, 2 * prime - 1 - n)
                aligned_bad = apery[folded] % prime == 0
            if not aligned_bad:
                contaminating.append(prime)
        if contaminating:
            counterexamples[n] = contaminating

    # Good candidate contamination: at n=23 neither 17 nor 19 is aligned
    # with an Apéry zero, yet both divide the primitive determinantal divisor.
    assert 17 in counterexamples[23]
    assert 19 in counterexamples[23]

    # Above-n contamination also occurs, so a top-support cutoff alone does
    # not repair the construction.
    assert 47 in counterexamples[38]
    assert 53 in counterexamples[42]

    print(
        "primitive all-cutoff Smith audit verified through "
        f"n={LIMIT}; contamination examples: "
        f"n=23->{counterexamples[23]}, "
        f"n=38->{counterexamples[38]}, "
        f"n=42->{counterexamples[42]}"
    )


if __name__ == "__main__":
    main()
