#!/usr/bin/env sage
"""
Exact characteristic-zero verifier for Q7693.

It constructs the Apéry numbers b_r from the stated recurrence, constructs

g(t) = 1 / (F(t)^2 * sqrt(1 - 34*t + t^2))

in QQ[[t]], checks that the required coefficients g_r are integral, and then
forms

    Xi_0 = -1,
    Xi_r = -1 - 5 * sum_{m=1}^r g_m b_{m-1}.

For each r it factors gcd(b_r, Xi_r) with Sage's proof flag enabled.  Thus a
reported prime p is an exact characteristic-zero common prime, not a modular
or floating-point candidate.  The scan records repeated prime labels instead
of silently deduplicating them.

Typical use:

    sage q7693_transverse_support.sage --max-r 10000 --expect-known-10000

To search for a violation of p <= 5 r:

    sage q7693_transverse_support.sage --max-r 20000 --C 5 --offset 0

The optional offset checks p <= C*r + offset.  No finite run proves an
all-index support theorem.
"""

from sage.all import (
    ZZ,
    QQ,
    PowerSeriesRing,
    factor,
    gcd,
    proof,
)
import argparse
from collections import defaultdict


def apery_numbers(max_r):
    """Return b_0,...,b_max_r exactly in ZZ."""
    if max_r == 0:
        return [ZZ(1)]
    b = [ZZ(1), ZZ(5)]
    for n in range(1, max_r):
        nZ = ZZ(n)
        num = (
            (34 * nZ**3 + 51 * nZ**2 + 27 * nZ + 5) * b[n]
            - nZ**3 * b[n - 1]
        )
        den = ZZ(n + 1) ** 3
        q, rem = num.quo_rem(den)
        assert rem == 0, (n, num, den)
        b.append(q)
    return b


def g_coefficients(b, max_r):
    """Return g_0,...,g_max_r from the exact power-series definition."""
    prec = max_r + 2
    R = PowerSeriesRing(QQ, "t", default_prec=prec)
    t = R.gen()
    F = R([QQ(x) for x in b], prec=prec)
    D = R(1 - 34 * t + t**2, prec=prec)
    sqrtD = D.sqrt(prec=prec)
    G = (F * F * sqrtD) ** (-1)

    out = []
    for n in range(max_r + 1):
        c = QQ(G[n])
        assert c.denominator() == 1, (n, c)
        out.append(ZZ(c))
    assert out[0] == 1
    return out


def xi_numbers(b, g, max_r):
    """Return Xi_0,...,Xi_max_r exactly in ZZ."""
    xi = [ZZ(-1)]
    running = ZZ(0)
    for r in range(1, max_r + 1):
        running += g[r] * b[r - 1]
        xi.append(ZZ(-1) - 5 * running)
    return xi


def scan(max_r, C, offset, expect_known_10000=False):
    proof.all(True)
    b = apery_numbers(max_r)
    g = g_coefficients(b, max_r)
    xi = xi_numbers(b, g, max_r)

    # Endpoint checks used in the audit.
    if max_r >= 2:
        assert b[0] == 1 and b[1] == 5 and b[2] == 73
        assert g[0] == 1 and g[1] == 7 and g[2] == 192
        assert xi[0] == -1 and xi[1] == -36 and xi[2] == -4836
    # Since Xi_r = -1 - 5*(integer), p=5 never occurs.
    for r in range(max_r + 1):
        assert xi[r] % 5 == 4

    common = []
    violations = []
    labels = defaultdict(list)

    for r in range(max_r + 1):
        G = gcd(abs(b[r]), abs(xi[r]))
        if G <= 1:
            continue
        for p, exponent in factor(G):
            p = ZZ(p)
            exponent = ZZ(exponent)
            assert p.is_prime(proof=True)
            assert b[r] % p == 0
            assert xi[r] % p == 0
            if p > r:
                row = (int(p), int(r), int(exponent))
                common.append(row)
                labels[int(p)].append(int(r))
                if p > C * r + offset:
                    violations.append(row)
                    print(
                        "VIOLATION:",
                        f"p={p}, r={r}, v_p(gcd)={exponent}, "
                        f"p-Cr={p-C*r}"
                    )

    common.sort(key=lambda x: (x[1], x[0]))
    repeated = {p: rs for p, rs in labels.items() if len(rs) > 1}

    print(f"max_r={max_r}; C={C}; offset={offset}")
    print("common pairs p>r (p,r,v_p(gcd)):")
    for row in common:
        print("  ", row)
    print("repeated high-prime labels:", repeated)
    print("violations:", violations)

    if expect_known_10000:
        assert max_r == 10000, "use exactly --max-r 10000 with this regression"
        got = {(p, r) for p, r, _ in common}
        expected = {(17, 13), (2237, 492)}
        assert got == expected, (got, expected)
        assert repeated == {}, repeated
        assert violations == [], violations
        print("KNOWN-10000 REGRESSION VERIFIED")

    return common, repeated, violations


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-r", type=int, default=1000)
    parser.add_argument("--C", type=int, default=5)
    parser.add_argument("--offset", type=int, default=0)
    parser.add_argument("--expect-known-10000", action="store_true")
    args = parser.parse_args()
    if args.max_r < 0:
        raise ValueError("--max-r must be nonnegative")
    if args.C <= 0:
        raise ValueError("--C must be positive")
    scan(args.max_r, ZZ(args.C), ZZ(args.offset), args.expect_known_10000)


if __name__ == "__main__":
    main()
