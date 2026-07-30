#!/usr/bin/env python3
"""Exact audit of the block Newton--Pascal target carrier.

For an integer sequence Y define

    G_{d,L}(Y) =
        sum_{i=0}^L (-1)^i binom(d+i,i)
        binom(d+L+1,L-i) Y_{d+i}.

If D > N and q is prime in (D,D+N], then

    G_{D-1,N}(Y) == Y_{q-1} (mod q)

and q occurs exactly once in binom(D+N,N).  Consequently the part of
gcd(G_{D-1,N}(Y), binom(D+N,N)) supported on that prime interval is
exactly the target radical.

The second half independently computes the fixed-moment Apéry shell
C_M(d) from its coefficient formula and checks that its prime targets
agree with the corresponding Apéry residues inside fixed-quotient
cells.  This script proves no asymptotic bound; it audits only the exact
interface to the remaining primitive-content problem.
"""

from __future__ import annotations

from functools import lru_cache
from math import comb, gcd


def prime_flags(limit: int) -> bytearray:
    flags = bytearray(b"\x01") * (limit + 1)
    if limit >= 0:
        flags[0] = 0
    if limit >= 1:
        flags[1] = 0
    p = 2
    while p * p <= limit:
        if flags[p]:
            flags[p * p : limit + 1 : p] = b"\x00" * (
                (limit - p * p) // p + 1
            )
        p += 1
    return flags


def apery_values(limit: int) -> list[int]:
    if limit == 0:
        return [1]
    values = [1, 5]
    for n in range(1, limit):
        numerator = (
            (34 * n**3 + 51 * n**2 + 27 * n + 5) * values[n]
            - n**3 * values[n - 1]
        )
        denominator = (n + 1) ** 3
        value, remainder = divmod(numerator, denominator)
        assert remainder == 0
        values.append(value)
    return values[: limit + 1]


def newton_carrier(values: list[int], d: int, length: int) -> int:
    assert d >= 0 and length >= 0
    assert d + length < len(values)
    return sum(
        (-1) ** i
        * comb(d + i, i)
        * comb(d + length + 1, length - i)
        * values[d + i]
        for i in range(length + 1)
    )


def interval_radical_identity(
    values: list[int],
    d_start: int,
    length: int,
    is_prime: bytearray,
) -> tuple[int, int]:
    """Return the target radical and the interval-prime modulus."""
    assert d_start > length >= 1
    carrier = newton_carrier(values, d_start - 1, length)
    pascal = comb(d_start + length, length)
    target_radical = 1
    interval_modulus = 1

    for q in range(d_start + 1, d_start + length + 1):
        if not is_prime[q]:
            continue
        assert q > length
        assert carrier % q == values[q - 1] % q
        assert pascal % q == 0
        assert pascal % (q * q) != 0
        interval_modulus *= q
        if values[q - 1] % q == 0:
            target_radical *= q

    assert gcd(carrier, interval_modulus) == target_radical
    assert gcd(gcd(carrier, pascal), interval_modulus) == target_radical
    return target_radical, interval_modulus


def comb_zero(n: int, k: int) -> int:
    return comb(n, k) if 0 <= k <= n else 0


@lru_cache(maxsize=None)
def apery_shell_coefficient(m: int, u: int, v: int, w: int) -> int:
    """Coefficient [x^u y^v z^w] Lambda(x,y,z)^m from (49.1)."""
    return sum(
        comb(m, t)
        * comb_zero(m, t - u)
        * comb_zero(2 * m - t, m - v)
        * comb_zero(2 * m - t, m - w)
        for t in range(m + 1)
    )


@lru_cache(maxsize=None)
def apery_shell(m: int, d: int) -> int:
    """C_M(d): sum of coefficients on the lattice d Z^3."""
    assert d >= 1
    radius = m // d
    return sum(
        apery_shell_coefficient(m, d * i, d * j, d * k)
        for i in range(-radius, radius + 1)
        for j in range(-radius, radius + 1)
        for k in range(-radius, radius + 1)
    )


def generic_sequence_audit() -> int:
    limit = 160
    is_prime = prime_flags(limit)
    apery = apery_values(limit)
    sequences = [
        [3 * n**4 - 7 * n**2 + 11 for n in range(limit + 1)],
        [(-5) ** n + 2 * n + 9 for n in range(limit + 1)],
        apery,
    ]
    checks = 0
    for values in sequences:
        for d_start in range(3, 80):
            for length in range(1, min(d_start - 1, 12) + 1):
                interval_radical_identity(values, d_start, length, is_prime)
                checks += 1
    return checks


def fixed_moment_shell_audit() -> tuple[int, int, int]:
    max_m = 24
    is_prime = prime_flags(max_m + 3)
    apery = apery_values(max_m)
    shell_lucas_checks = 0
    block_checks = 0
    nontrivial_targets = 0

    for m in range(6, max_m + 1):
        # The constant shell is the Apéry number itself.
        assert apery_shell(m, m + 1) == apery[m]

        values = [0] + [apery_shell(m, d) for d in range(1, m + 1)]
        for d_start in range(3, m + 1):
            max_length = min(4, d_start - 1, m - d_start + 1)
            for length in range(1, max_length + 1):
                nodes = range(d_start, d_start + length)
                quotients = {m // d for d in nodes}
                if len(quotients) != 1:
                    continue
                quotient = next(iter(quotients))
                if quotient == 0:
                    continue
                residues = [m - quotient * d for d in nodes]
                if any(not (1 <= r <= d - 1) for r, d in zip(residues, nodes)):
                    continue

                radical, _ = interval_radical_identity(
                    values, d_start, length, is_prime
                )
                direct_radical = 1
                for d, residue in zip(nodes, residues):
                    q = d + 1
                    if not is_prime[q]:
                        continue
                    shell_zero = values[d] % q == 0
                    apery_zero = apery[residue] % q == 0
                    assert shell_zero == apery_zero
                    shell_lucas_checks += 1
                    if apery_zero:
                        direct_radical *= q
                        nontrivial_targets += 1
                assert radical == direct_radical
                block_checks += 1

    return shell_lucas_checks, block_checks, nontrivial_targets


def main() -> None:
    generic_checks = generic_sequence_audit()
    shell_lucas, shell_blocks, targets = fixed_moment_shell_audit()
    print("Q32_BLOCK_PASCAL_GCD_AUDIT=PASS")
    print("GENERIC_BLOCK_CHECKS", generic_checks)
    print("SHELL_LUCAS_CHECKS", shell_lucas)
    print("FIXED_MOMENT_BLOCK_CHECKS", shell_blocks)
    print("NONTRIVIAL_SHELL_TARGETS", targets)


if __name__ == "__main__":
    main()
