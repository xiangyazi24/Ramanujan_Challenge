#!/usr/bin/env python3
"""Audit the mod-p quarter law, reversal signs, and complete zero sets.

The proof in CODEX_JACOBSTHAL_DEEP.md is uniform in p.  This program checks
its exact identities and records the additional (currently conjectural) zero
patterns through all primes below LIMIT.
"""

from collections import Counter, defaultdict
from math import isqrt


LIMIT = 3000


def primes_below(n):
    sieve = bytearray(b"\x01") * n
    sieve[:2] = b"\x00\x00"
    for p in range(2, isqrt(n - 1) + 1):
        if sieve[p]:
            sieve[p * p : n : p] = b"\x00" * (((n - 1 - p * p) // p) + 1)
    return [p for p in range(2, n) if sieve[p]]


def legendre(a, p):
    value = pow(a % p, (p - 1) // 2, p)
    return -1 if value == p - 1 else value


def apery_truncation_mod(p):
    """Return b_0,...,b_(p-1) modulo p from the Apéry recurrence."""
    values = [1, 5]
    for n in range(1, p - 1):
        numerator = (
            (2 * n + 1) * (17 * n * n + 17 * n + 5) * values[n]
            - n**3 * values[n - 1]
        ) % p
        values.append(numerator * pow((n + 1) ** 3, -1, p) % p)
    return values[:p]


def branch_sequence(kind, degree, p):
    inv2 = pow(2, -1, p)
    if kind == "tau":
        values = [1, 5 * inv2 % p]
        for n in range(1, degree):
            numerator = (
                2 * (68 * n * n + 34 * n + 5) * values[n]
                - (2 * n - 1) ** 2 * values[n - 1]
            ) % p
            values.append(numerator * pow(4 * (n + 1) ** 2, -1, p) % p)
    else:
        values = [1, 39 * inv2 % p]
        for n in range(1, degree):
            numerator = (
                2 * (68 * n * n + 102 * n + 39) * values[n]
                - (2 * n + 1) ** 2 * values[n - 1]
            ) % p
            values.append(numerator * pow(4 * (n + 1) ** 2, -1, p) % p)
    return values[: degree + 1]


def representation(p, kind):
    for y in range(isqrt(p) + 1):
        remainder = p - (6 * y * y if kind == "principal" else 3 * y * y)
        if remainder < 0:
            break
        if kind == "principal":
            x = isqrt(remainder)
            if x * x == remainder:
                return x, y
        elif remainder % 2 == 0:
            x = isqrt(remainder // 2)
            if 2 * x * x == remainder:
                return x, y
    return None


def main():
    table = defaultdict(lambda: [0, 0])
    size_distributions = defaultdict(Counter)
    extra_tau = []
    p23_examples = []
    p13_examples = []
    raw_exceptions = []
    ordinary_double_roots = 0

    for p in primes_below(LIMIT):
        residue = p % 24
        if p < 5 or residue not in (1, 5, 7, 11, 13, 17, 19, 23):
            continue
        kind = "tau" if legendre(-6, p) == 1 else "sigma"
        degree = (p - 1) // 2 if kind == "tau" else (p - 3) // 2
        values = branch_sequence(kind, degree, p)
        epsilon = legendre(-2, p)
        assert all(values[degree - j] == epsilon * values[j] % p for j in range(degree + 1))

        # This directly audits the Sun + ordinary-point input in the proof:
        # A_p(1) has order exactly two precisely in the anti-reciprocal cases.
        apery_values = apery_truncation_mod(p)
        at_one = sum(apery_values) % p
        first = sum(n * value for n, value in enumerate(apery_values)) % p
        second = sum(n * (n - 1) * value for n, value in enumerate(apery_values)) % p
        assert (at_one == 0) == (epsilon == -1)
        assert (2 * first - (p - 1) * at_one) % p == 0
        if epsilon == -1:
            assert first == 0 and second != 0
            ordinary_double_roots += 1

        quarter = (p - 1) // 4 if kind == "tau" else (p - 3) // 4
        zero = values[quarter] == 0
        expected = residue in (5, 23)
        table[residue][0 if zero else 1] += 1
        assert zero == expected

        zeros = [j for j, value in enumerate(values) if value == 0]
        size_distributions[(residue, kind)][len(zeros)] += 1
        assert zeros == sorted(degree - j for j in zeros)

        if residue == 5 and len(zeros) > 1:
            extra_tau.append((p, representation(p, "nonprincipal"), quarter, zeros))
        if residue == 23:
            forced = [(p - 7) // 8, (p - 3) // 4, (3 * p - 5) // 8]
            assert all(j in zeros for j in forced)
            assert representation(p, "principal") is None
            assert representation(p, "nonprincipal") is None
            if len(p23_examples) < 4:
                p23_examples.append((p, forced, zeros))
        if residue == 13:
            forced = [(p - 5) // 8, (3 * p - 7) // 8]
            assert all(j in zeros for j in forced)
            if len(p13_examples) < 4:
                p13_examples.append((p, forced, zeros))

        # Audit a stronger raw-series interpretation, not used in the theorem.
        # It is false: p=71 already gives an additional tau floor-quarter zero.
        tau_index = (p - 1) // 4
        sigma_index = (p - 3) // 4
        tau_raw = branch_sequence("tau", tau_index, p)[tau_index]
        sigma_raw = branch_sequence("sigma", sigma_index, p)[sigma_index]
        if (tau_raw == 0) != (residue == 5) or (sigma_raw == 0) != (residue == 23):
            raw_exceptions.append((p, tau_raw, sigma_raw))

        if residue in (1, 7):
            assert representation(p, "principal") is not None
            assert representation(p, "nonprincipal") is None
        elif residue in (5, 11):
            assert representation(p, "principal") is None
            assert representation(p, "nonprincipal") is not None

    # The two non-vanishing-centre recurrence coefficients in the proof.
    for p in primes_below(LIMIT):
        if p > 3 and p % 4 == 1:
            j = (p - 1) // 4
            assert (4 * (j + 1) ** 2 + (2 * j - 1) ** 2) % p == 9 * pow(2, -1, p) % p
        if p > 3 and p % 4 == 3:
            j = (p - 3) // 4
            assert (4 * (j + 1) ** 2 + (2 * j + 1) ** 2) % p == pow(2, -1, p)

    print(f"Relevant-branch quarter law and reversal checked for every prime p < {LIMIT}")
    print(f"A_p(1) ordinary double-root cases checked: {ordinary_double_roots}")
    for residue in (1, 5, 7, 11, 13, 17, 19, 23):
        print(f"  p mod 24 = {residue:2d}: zero/nonzero = {table[residue]}")
    print("Complete relevant-branch zero-count distributions:")
    for key in sorted(size_distributions):
        print(f"  {key}: {dict(sorted(size_distributions[key].items()))}")
    print("p = 5 mod 24 cases with an extra symmetric zero pair:")
    for row in extra_tau:
        print(" ", row)
    print("First p = 23 mod 24 forced-eighth examples (forced, complete):")
    for row in p23_examples:
        print(" ", row)
    print("First p = 13 mod 24 forced-eighth examples (forced, complete):")
    for row in p13_examples:
        print(" ", row)
    print("Counterexamples to the stronger raw-both-series floor interpretation:")
    print(" ", raw_exceptions[:20])
    print("Discriminant -24 form classification and p=23 non-representability: VERIFIED")


if __name__ == "__main__":
    main()
