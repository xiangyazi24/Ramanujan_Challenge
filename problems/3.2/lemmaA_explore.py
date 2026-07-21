#!/usr/bin/env python3
"""Exact computations for Lemma A (Jacobi skeleton and decision gate).

The script uses only Python's standard library.  Characters are represented by
their exponents modulo h = p - 1, and their values at zero are always zero,
including for the trivial character.  The marked-coordinate convention at
j = 0 is handled separately, as required by oracleC_result.tex.
"""

from __future__ import annotations

import argparse
from collections import Counter
from math import comb


def primes_up_to(limit: int) -> list[int]:
    """Return all primes <= limit."""
    if limit < 2:
        return []
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[0:2] = b"\x00\x00"
    for q in range(2, int(limit**0.5) + 1):
        if sieve[q]:
            start = q * q
            sieve[start : limit + 1 : q] = b"\x00" * (
                (limit - start) // q + 1
            )
    return [q for q in range(2, limit + 1) if sieve[q]]


def apery_residues(p: int) -> list[int]:
    """Return b_j mod p for 0 <= j <= p-2, using the Apéry recurrence."""
    h = p - 1
    values = [0] * h
    values[0] = 1
    if h == 1:
        return values
    values[1] = 5 % p
    for n in range(1, h - 1):
        coefficient = (34 * n**3 + 51 * n**2 + 27 * n + 5) % p
        rhs = (coefficient * values[n] - n**3 * values[n - 1]) % p
        denominator = pow(n + 1, 3, p)
        values[n + 1] = rhs * pow(denominator, -1, p) % p
    return values


def apery_binomial_mod(j: int, p: int) -> int:
    """Compute b_j from its defining binomial sum, reduced modulo p."""
    return sum(comb(j, k) ** 2 * comb(j + k, k) ** 2 for k in range(j + 1)) % p


def character_mod(p: int, exponent: int, value: int) -> int:
    """Reduction of omega^exponent(value), with every character zero at 0."""
    value %= p
    if value == 0:
        return 0
    return pow(value, exponent % (p - 1), p)


def jacobi_mod_direct(p: int, a: int, b: int) -> int:
    """Direct finite-field reduction of J(omega^a, omega^b)."""
    return sum(
        character_mod(p, a, t) * character_mod(p, b, 1 - t)
        for t in range(p)
    ) % p


def jacobi_mod(p: int, a: int, b: int) -> int:
    """O(1) closed formula for J(omega^a, omega^b) modulo the prime above p."""
    h = p - 1
    a %= h
    b %= h
    if a == 0:
        return (p - 2) if b == 0 else (p - 1)
    if b == 0:
        return p - 1
    r = h - a
    if r > b:
        return 0
    sign = -1 if r % 2 else 1
    return (-sign * comb(b, r)) % p


def skeleton_mod(p: int, j: int) -> int:
    """Reduction of the exact four-Jacobi skeleton for c_{p,j}^{tor}.

    For 1 <= j <= p-2 this is

      -1/(p-1) sum_k J(-k,j)^2 J(-k,j+k)^2.

    At j=0 the zero-extended trivial-character skeleton counts only the open
    set U.  The exact endpoint correction #U-(p-1)^3=-4p^2+14p-13 converts it
    to the prescribed raw zeroth power, including zeros of Lambda.
    """
    h = p - 1
    if not 0 <= j <= p - 2:
        raise ValueError("j must lie in [0,p-2]")
    total = 0
    for k in range(h):
        first = jacobi_mod(p, -k, j)
        second = jacobi_mod(p, -k, j + k)
        total = (total + first * first * second * second) % p
    result = -pow(h, -1, p) * total
    if j == 0:
        result += -4 * p * p + 14 * p - 13
    return result % p


def lambda_value_mod(p: int, x: int, y: int, z: int) -> int:
    """The marked Laurent polynomial Lambda evaluated in F_p^3."""
    numerator = (
        (1 + x)
        * (1 + y)
        * (1 + z)
        * ((1 + y) * (1 + z) + x * y * z)
    )
    return numerator * pow((x * y * z) % p, -1, p) % p


def torus_coordinate_mod(p: int, j: int) -> int:
    """Direct marked torus sum, with the prescribed zero^0 = 1 convention."""
    total = 0
    for x in range(1, p):
        for y in range(1, p):
            for z in range(1, p):
                total += pow(lambda_value_mod(p, x, y, z), j, p)
    return (-total) % p


def digit_fraction_numerator(exponent: int, h: int) -> int:
    """Numerator in d(a)=[-a]_h/h from Gross--Koblitz."""
    return (-exponent) % h


def jacobi_valuation(p: int, a: int, b: int) -> int:
    """v_p(J(omega^a,omega^b)), normalized by v_p(p)=1.

    The final correction is the inverse-character degeneracy.  Omitting it
    gives the wrong answer on a+b=0.
    """
    h = p - 1
    a0 = a % h
    b0 = b % h
    numerator = (
        digit_fraction_numerator(a0, h)
        + digit_fraction_numerator(b0, h)
        - digit_fraction_numerator(a0 + b0, h)
    )
    if numerator % h:
        raise AssertionError("fractional-part expression is not integral")
    correction = int(a0 != 0 and b0 != 0 and (a0 + b0) % h == 0)
    value = numerator // h - correction
    if value not in (0, 1):
        raise AssertionError(f"unexpected Jacobi valuation {value}")
    return value


def term_valuation(p: int, j: int, k: int) -> int:
    """Valuation of T_{j,k}=J(-k,j)^2 J(-k,j+k)^2."""
    if not 1 <= j <= p - 2:
        raise ValueError("the Kummer skeleton requires 1 <= j <= p-2")
    return 2 * jacobi_valuation(p, -k, j) + 2 * jacobi_valuation(p, -k, j + k)


def term_valuation_simplified(p: int, j: int, k: int) -> int:
    """The explicit carry formula 2[k>j]+2[k>p-1-j]."""
    h = p - 1
    k %= h
    return 2 * int(k > j) + 2 * int(k > h - j)


def unit_indices(p: int, j: int) -> list[int]:
    """Indices of valuation-zero Jacobi summands."""
    if not 1 <= j <= p - 2:
        raise ValueError("U(j) is defined only on the nontrivial-character range")
    return list(range(min(j, p - 1 - j) + 1))


def unit_count(p: int, j: int) -> int:
    """U(j)=1+min(j,p-1-j)."""
    if not 1 <= j <= p - 2:
        raise ValueError("U(j) is defined only on the nontrivial-character range")
    return 1 + min(j, p - 1 - j)


def eta_product_coefficients(limit: int) -> list[int]:
    """Coefficients of eta(2z)^4 eta(4z)^4 through q^limit."""
    if limit < 0:
        return []
    # eta(2z)^4 eta(4z)^4 = q prod_n (1-q^(2n))^4(1-q^(4n))^4.
    product_limit = max(0, limit - 1)
    product = [0] * (product_limit + 1)
    product[0] = 1

    def multiply_fourth_power(step: int) -> None:
        old = product[:]
        for degree in range(product_limit + 1):
            value = old[degree]
            for power, coefficient in ((1, -4), (2, 6), (3, -4), (4, 1)):
                source = degree - power * step
                if source >= 0:
                    value += coefficient * old[source]
            product[degree] = value

    for n in range(1, product_limit // 2 + 1):
        multiply_fourth_power(2 * n)
    for n in range(1, product_limit // 4 + 1):
        multiply_fourth_power(4 * n)

    coefficients = [0] * (limit + 1)
    for degree in range(1, limit + 1):
        coefficients[degree] = product[degree - 1]
    return coefficients


def decision_gate(limit: int) -> dict[str, object]:
    """Compute the empirical A3 classification through the given prime bound."""
    primes = [p for p in primes_up_to(limit) if p >= 5]
    distribution: Counter[int] = Counter()
    total_zeros = 0
    collision_zeros = 0
    nonordinary: list[tuple[int, int]] = []
    zero_free = 0
    max_zeros = 0
    records: list[tuple[int, int, int]] = []
    for p in primes:
        row = apery_residues(p)
        zeros = [j for j, value in enumerate(row) if value == 0]
        midpoint = (p - 1) // 2
        central = int(midpoint in zeros)
        if central:
            nonordinary.append((p, midpoint))
        collision_zeros += len(zeros) - central
        total_zeros += len(zeros)
        distribution[len(zeros)] += 1
        zero_free += int(not zeros)
        max_zeros = max(max_zeros, len(zeros))
        records.append((p, len(zeros), p - 2))
    return {
        "limit": limit,
        "prime_count": len(primes),
        "total_zeros": total_zeros,
        "collision_zeros": collision_zeros,
        "nonordinary": nonordinary,
        "zero_free": zero_free,
        "max_zeros": max_zeros,
        "distribution": dict(sorted(distribution.items())),
        "records": records,
    }


def format_distribution(distribution: dict[int, int]) -> str:
    return ", ".join(f"Z={z}:{count}" for z, count in distribution.items())


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--limit", type=int, default=5000, help="A3 prime bound")
    parser.add_argument(
        "--skeleton-limit",
        type=int,
        default=200,
        help="verify the A1 skeleton through this prime bound",
    )
    args = parser.parse_args()

    checked_pairs = 0
    for p in primes_up_to(args.skeleton_limit):
        if p < 5:
            continue
        row = apery_residues(p)
        for j, value in enumerate(row):
            if skeleton_mod(p, j) != value:
                raise AssertionError(f"skeleton mismatch at p={p}, j={j}")
            checked_pairs += 1
    print(
        f"A1 PASS: skeleton=b_j for all {checked_pairs} pairs "
        f"with 5<=p<={args.skeleton_limit}."
    )

    report = decision_gate(args.limit)
    print(
        f"A3 p<={args.limit}: primes={report['prime_count']}, "
        f"sum_Z={report['total_zeros']}, zero_free={report['zero_free']}, "
        f"max_Z={report['max_zeros']}"
    )
    print(f"distribution: {format_distribution(report['distribution'])}")
    print(
        f"classification: U=0:0, U>=2 collisions={report['collision_zeros']}, "
        f"nonordinary={len(report['nonordinary'])} {report['nonordinary']}"
    )
    print("container: |{1<=j<=p-2: U(j)!=1}|=p-2 for every prime (S2).")

    records = report["records"]
    if records:
        print("sample p,Z(p),container:")
        sample_indices = sorted(
            {0, len(records) // 4, len(records) // 2, 3 * len(records) // 4, len(records) - 1}
        )
        for index in sample_indices:
            p, zeros, container = records[index]
            print(f"  {p:5d} {zeros:2d} {container:5d}")


if __name__ == "__main__":
    main()
