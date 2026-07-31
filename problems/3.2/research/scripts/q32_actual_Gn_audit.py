#!/usr/bin/env python3
"""The headline quantity itself: G_n = gcd(d_n a_n, d_n b_n), d_n = lcm(1..n)^3.

Also splits it as G_n = (d_n/D_n) * gcd(A_n, b_n) with a_n = A_n/D_n in lowest terms,
so that the two factors of the reduction can be compared empirically:
  - the OVER-CLEARANCE factor d_n/D_n (is the classical denominator envelope wasteful?)
  - the INTRINSIC numerator gcd gcd(A_n, b_n).
Prints log(.)/n for each, so the o(n) question is directly visible.
"""
import sys
from fractions import Fraction
from math import gcd, log, lcm


def run(N):
    a0, a1 = Fraction(0), Fraction(6)
    b0, b1 = Fraction(1), Fraction(5)
    A = {0: a0, 1: a1}
    B = {0: b0, 1: b1}
    for n in range(1, N):
        c = Fraction(34 * n ** 3 + 51 * n ** 2 + 27 * n + 5, (n + 1) ** 3)
        d = Fraction(n ** 3, (n + 1) ** 3)
        A[n + 1] = c * A[n] - d * A[n - 1]
        B[n + 1] = c * B[n] - d * B[n - 1]
    L = 1
    print("  n   log(G_n)/n   log(over-clear)/n   log gcd(A,b)/n   G_n bits", flush=True)
    for n in range(1, N + 1):
        L = lcm(L, n)
        dn = L ** 3
        an, bn = A[n], B[n]
        assert bn.denominator == 1
        num = dn * an
        assert num.denominator == 1
        num = int(num)
        G = gcd(abs(num), int(bn) * dn)
        Dn = an.denominator
        over = dn // Dn
        intr = gcd(abs(int(an.numerator)), int(bn))
        if n % 10 == 0 or n <= 5:
            print(f"{n:4d}  {log(G)/n if G>1 else 0:10.5f}   {log(over)/n if over>1 else 0:15.5f}"
                  f"   {log(intr)/n if intr>1 else 0:13.5f}   {G.bit_length():8d}", flush=True)


if __name__ == "__main__":
    run(int(sys.argv[1]) if len(sys.argv) > 1 else 200)
