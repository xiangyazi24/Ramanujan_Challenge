#!/usr/bin/env python3
"""Self-contained deterministic verifier for Q7707.

No third-party packages.  This checks finite instances of the exact algebraic
identities proved in chatgpt_q7707_all_index_mellin_attack.md:

  * Morita-Gamma factorial identity with the sign cancellation;
  * the 0/2 valuation break and slope-zero sum;
  * reflection b_r == b_{p-1-r} (mod p);
  * the explicit mod-p^2 reflection defect;
  * the all-quotient reindexing n=qp+r, m=n-q;
  * the n=11 counterexample to "q=1 suffices";
  * finite checks of the gap-polynomial reflection symmetry and the automatic
    central factor for even reflection gaps.

The mathematical audit proves the identities symbolically; this script is only
an independent finite regression checker.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from math import comb


def primes_upto(n: int) -> list[int]:
    if n < 2:
        return []
    sieve = bytearray(b"\x01") * (n + 1)
    sieve[0:2] = b"\x00\x00"
    for d in range(2, int(n**0.5) + 1):
        if sieve[d]:
            sieve[d * d : n + 1 : d] = b"\x00" * (((n - d * d) // d) + 1)
    return [i for i in range(2, n + 1) if sieve[i]]


def apery(n: int) -> int:
    return sum((comb(n, k) * comb(n + k, k)) ** 2 for k in range(n + 1))


def apery_term(r: int, k: int) -> int:
    return (comb(r, k) * comb(r + k, k)) ** 2


def gamma_p_integer(p: int, s: int) -> int:
    """Morita convention Gamma_p(s)=(-1)^s prod_{1<=j<s,p\nmid j} j."""
    ans = -1 if (s & 1) else 1
    for j in range(1, s):
        if j % p:
            ans *= j
    return ans


def morita_u(p: int, r: int, k: int) -> Fraction:
    return Fraction(
        gamma_p_integer(p, r + k + 1),
        gamma_p_integer(p, k + 1) ** 2 * gamma_p_integer(p, r - k + 1),
    )


def frac_mod(x: Fraction, modulus: int) -> int:
    den = x.denominator % modulus
    return (x.numerator % modulus) * pow(den, -1, modulus) % modulus


def harmonic(t: int) -> Fraction:
    return sum((Fraction(1, j) for j in range(1, t + 1)), Fraction(0, 1))


def check_morita_and_reflection(max_p: int) -> None:
    for p in primes_upto(max_p):
        if p < 5:
            continue
        b = [apery(r) for r in range(p)]
        for r in range(p):
            a = min(r, p - 1 - r)
            theta = 0
            for k in range(r + 1):
                A = comb(r, k) * comb(r + k, k)
                carry = int(r + k >= p)
                U = morita_u(p, r, k)
                assert Fraction(A, p**carry) == U, ("Morita", p, r, k)

                T = A * A
                if k <= a:
                    assert T % p != 0, ("slope-zero", p, r, k)
                    theta = (theta + frac_mod(U * U, p)) % p
                else:
                    assert T % (p * p) == 0, ("slope-two", p, r, k)

            assert b[r] % p == theta, ("theta", p, r)
            assert b[r] % p == b[a] % p, ("reflection", p, r, a)
            if r in (0, p - 1):
                assert b[r] % p != 0, ("endpoint", p, r)


def check_second_digit(max_p: int) -> None:
    for p in primes_upto(max_p):
        if p < 5:
            continue
        modulus = p * p
        for a in range((p - 2) // 2 + 1):
            r = p - 1 - a
            D = Fraction(0, 1)
            for k in range(a + 1):
                T = apery_term(a, k)
                D += T * (harmonic(a + k) - harmonic(a - k))
            lhs = (
                apery(r)
                - apery(a)
                + 2 * p * frac_mod(D, modulus)
            ) % modulus
            assert lhs == 0, ("p2-reflection", p, a, r, D)


def check_all_q(max_p: int, max_n: int) -> None:
    primes = [p for p in primes_upto(max_p) if p >= 5]
    b_cache = {r: apery(r) for r in range(max_p)}
    for n in range(1, max_n + 1):
        for p in primes:
            if p > n:
                break
            q, r = divmod(n, p)
            assert q >= 1
            m = n - q
            assert m == q * (p - 1) + r
            assert (m - r) % (p - 1) == 0

            a = min(r, p - 1 - r)
            assert b_cache[r] % p == b_cache[a] % p

            if 2 * r <= p - 1:
                assert a == r
                assert (n - a) % p == 0
                assert p * q == n - a
            else:
                assert a == p - 1 - r
                assert (n + 1 + a) % p == 0
                assert p * (q + 1) == n + 1 + a

            if r in (0, p - 1):
                assert b_cache[r] % p != 0

    # Explicit counterexample to the false implication "q=1 controls all q".
    n = 11
    p = 5
    q, r = divmod(n, p)
    assert (q, r, apery(r) % p) == (2, 1, 0)
    top = []
    for ell in primes_upto(n):
        if ell >= 5 and n // ell == 1:
            rr = n % ell
            top.append((ell, rr, apery(rr) % ell))
    assert top == [(7, 4, 3), (11, 0, 1)], top


# Polynomial helpers, coefficients in ascending order.
def trim(a: list[int]) -> list[int]:
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return a


def padd(a: list[int], b: list[int], sign: int = 1) -> list[int]:
    out = [0] * max(len(a), len(b))
    for i, x in enumerate(a):
        out[i] += x
    for i, x in enumerate(b):
        out[i] += sign * x
    return trim(out)


def pmul(a: list[int], b: list[int]) -> list[int]:
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return trim(out)


def ppow_linear(c: int, e: int) -> list[int]:
    # (x+c)^e
    return [comb(e, j) * c ** (e - j) for j in range(e + 1)]


def P_shift(c: int) -> list[int]:
    # P(x+c)=34(x+c)^3+51(x+c)^2+27(x+c)+5
    out = [5]
    out = padd(out, [27 * x for x in ppow_linear(c, 1)])
    out = padd(out, [51 * x for x in ppow_linear(c, 2)])
    out = padd(out, [34 * x for x in ppow_linear(c, 3)])
    return trim(out)


def peval(a: list[int], x: Fraction | int) -> Fraction:
    ans = Fraction(0, 1)
    xx = Fraction(x, 1) if isinstance(x, int) else x
    for c in reversed(a):
        ans = ans * xx + c
    return ans


def gap_polynomials(max_h: int) -> dict[int, list[int]]:
    # N_1=1, N_2=P(x+1),
    # N_{h+1}=P(x+h)N_h-(x+h)^6 N_{h-1}.
    N: dict[int, list[int]] = {1: [1]}
    if max_h >= 2:
        N[2] = P_shift(1)
    for h in range(2, max_h):
        N[h + 1] = padd(
            pmul(P_shift(h), N[h]),
            pmul(ppow_linear(h, 6), N[h - 1]),
            sign=-1,
        )
    return N


def check_gap_center(max_h: int, max_p: int) -> None:
    N = gap_polynomials(max_h)
    for h, poly in N.items():
        # The symbolic proof is in the audit.  At even h, oddness around the
        # center forces 2x+h+1 as a factor; center evaluation is an exact check.
        if h % 2 == 0:
            center = Fraction(-(h + 1), 2)
            assert peval(poly, center) == 0, ("gap-center", h)

    # For reflection gaps h=p-1-2a (necessarily even), N_h(a)=0 mod p
    # automatically because 2a+h+1=p.  Check every instance available in N.
    for p in primes_upto(max_p):
        if p < 5:
            continue
        for a in range((p - 1) // 2 + 1):
            h = p - 1 - 2 * a
            if 2 <= h <= max_h:
                assert h % 2 == 0
                assert 2 * a + h + 1 == p
                assert peval(N[h], a).denominator == 1
                assert int(peval(N[h], a)) % p == 0, ("gap-auto", p, a, h)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-p", type=int, default=43)
    ap.add_argument("--max-n", type=int, default=250)
    ap.add_argument("--max-h", type=int, default=14)
    args = ap.parse_args()

    if args.max_p < 11 or args.max_n < 11 or args.max_h < 2:
        raise SystemExit("use --max-p >= 11, --max-n >= 11, --max-h >= 2")

    check_morita_and_reflection(args.max_p)
    check_second_digit(args.max_p)
    check_all_q(args.max_p, args.max_n)
    check_gap_center(args.max_h, args.max_p)
    print("Q7707_VERIFY PASS")


if __name__ == "__main__":
    main()
