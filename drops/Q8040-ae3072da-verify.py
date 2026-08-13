#!/usr/bin/env python3
"""Verifier for the exponent ledger and separation model in Q8040.

This is an exact finite/combinatorial verifier.  It does not claim to verify
any unproved statement about the actual Apéry zero sets.
"""

from fractions import Fraction
from math import isqrt, log


def falling(n: int, k: int) -> int:
    out = 1
    for j in range(k):
        out *= n - j
    return out


def primes_between(lo: int, hi: int) -> list[int]:
    sieve = bytearray(b"\x01") * (hi + 1)
    sieve[0:2] = b"\x00\x00"
    for q in range(2, isqrt(hi) + 1):
        if sieve[q]:
            start = q * q
            sieve[start : hi + 1 : q] = b"\x00" * (((hi - start) // q) + 1)
    return [p for p in range(lo + 1, hi + 1) if sieve[p]]


def fold(p: int, r: int) -> int:
    return min(r, p - 1 - r)


def isolated_folds(p: int, row: set[int], H: int) -> set[int]:
    ys = sorted({fold(p, r) for r in row})
    return {
        s
        for s in ys
        if all(t == s or abs(t - s) > H for t in ys)
    }


def exponent_ledger() -> None:
    # With a carrier cost H^2 and isolated density X/(H log X), an order-t
    # moment gives exponent theta_t = 2/3 + 4/(3t).
    def theta(t: int) -> Fraction:
        return Fraction(2, 3) + Fraction(4, 3 * t)

    assert theta(3) == Fraction(10, 9)
    assert theta(4) == Fraction(1, 1)
    assert theta(5) == Fraction(14, 15)

    alpha = Fraction(7, 15)
    assert 2 * alpha == Fraction(14, 15)                 # carrier load H^2
    assert 1 - alpha == Fraction(8, 15)                  # isolated density
    assert Fraction(2, 5) + (1 - alpha) == Fraction(14, 15)

    eta = Fraction(1, 30)
    threshold = Fraction(14, 15) + eta
    assert threshold == Fraction(29, 30)
    assert threshold < 1

    rhs_power = 2 + 5 * Fraction(8, 15)
    one_high_output_power = 5 * threshold
    assert rhs_power == Fraction(14, 3)
    assert one_high_output_power == Fraction(29, 6)
    assert one_high_output_power > rhs_power


def separation_model() -> None:
    # Every active row is reflected and has a short nonreflection neighbour:
    # Z_p^* = {2,4,p-5,p-3}.  At m=2 all active rows hit, but the isolated
    # subbank is empty.  The active count is capped by H^2, matching the
    # deterministic carrier conclusion, while full HM7 fails by ~X^5.
    X = 10_000
    H = max(2, int((X ** (7 / 15)) / (log(X) ** (1 / 3))))
    bank = primes_between(X, 2 * X)
    R = min(H * H, len(bank))
    assert R >= 7
    active = bank[:R]

    m0 = 2
    for p in active:
        row = {2, 4, p - 5, p - 3}
        assert len(row) == 4
        assert {p - 1 - r for r in row} == row
        ordered = sorted(row)
        assert all(ordered[j + 1] - ordered[j] != 1 for j in range(3))
        assert m0 % p in row
        assert isolated_folds(p, row, H) == set()

    K_near = R
    K_iso = 0
    assert K_near <= H * H
    assert K_iso == 0

    lam = sum((Fraction(4, p) for p in active), Fraction(0, 1))
    F2_at_anchor = falling(R, 2)
    F7_at_anchor = falling(R, 7)

    # The existing HM2 scale is respected.
    hm2_rhs = 5 * X * X * lam * lam
    assert Fraction(F2_at_anchor, 1) <= hm2_rhs

    # The isolated high-load fifth-moment hypothesis is vacuous.
    isolated_tail_fifth = 0
    assert isolated_tail_fifth == 0

    # Full HM7 at the independence scale is violated.
    hm7_rhs = X * X * lam**7
    assert Fraction(F7_at_anchor, 1) > hm7_rhs

    # A robust symbolic upper comparison: lambda <= 4R/X, so the HM7 RHS is
    # at most 4^7 R^7 / X^5, whereas (R)_7 is comparable to R^7.
    assert lam <= Fraction(4 * R, X)
    assert Fraction(F7_at_anchor, 1) > Fraction((4**7) * (R**7), X**5)

    print("X", X)
    print("H", H)
    print("active_rows", R)
    print("lambda", float(lam))
    print("HM2_ratio", float(Fraction(F2_at_anchor, 1) / hm2_rhs))
    print("HM7_failure_ratio", float(Fraction(F7_at_anchor, 1) / hm7_rhs))


def main() -> None:
    exponent_ledger()
    separation_model()
    print("Q8040 verifier: PASS")


if __name__ == "__main__":
    main()
