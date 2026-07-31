#!/usr/bin/env python3
"""The exact denominator D_n of Apery's a_n: find the law.

d_n = lcm(1..n)^3, v_p(d_n) = 3*floor(log_p n).  Define the DEFECT
    e_p(n) = v_p(d_n) - v_p(D_n) >= 0,
so that G_n = d_n/D_n = prod_p p^{e_p(n)} (given gcd(A_n,b_n)=1, which we also recheck).

Question: what is the rule for e_p(n)?  Prints, for each n, the list of (p, e_p) with e_p > 0,
plus whether p | b_n, plus the range bucket of p (top window, middle, small).
"""
import sys
from fractions import Fraction
from math import gcd, lcm, log


def primes_upto(n):
    s = bytearray([1]) * (n + 1)
    s[0:2] = b'\x00\x00'
    for i in range(2, int(n ** .5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(n + 1) if s[i]]


def vp(x, p):
    e = 0
    while x % p == 0:
        x //= p
        e += 1
    return e


def run(N):
    A = {0: Fraction(0), 1: Fraction(6)}
    B = {0: Fraction(1), 1: Fraction(5)}
    for n in range(1, N):
        c = Fraction(34 * n ** 3 + 51 * n ** 2 + 27 * n + 5, (n + 1) ** 3)
        d = Fraction(n ** 3, (n + 1) ** 3)
        A[n + 1] = c * A[n] - d * A[n - 1]
        B[n + 1] = c * B[n] - d * B[n - 1]
    ps = primes_upto(N)
    print("n | defects (p^e, tag)   [tag: T=top window & p|b_n, t=top window not dividing, M=middle, S=small]", flush=True)
    for n in range(2, N + 1):
        an, bn = A[n], int(B[n])
        Dn = an.denominator
        out = []
        for p in ps:
            if p > n:
                break
            vd = 3 * int(log(n) / log(p) + 1e-9)
            e = vd - vp(Dn, p)
            if e:
                if p > n / 2:
                    tag = 'T' if bn % p == 0 else 't'
                elif p * p > n:
                    tag = 'M'
                else:
                    tag = 'S'
                out.append((p, e, tag))
        # sanity: gcd(A_n,b_n)
        g = gcd(abs(an.numerator), bn)
        if out or g > 1:
            print(f"{n:4d} | {out}  gcd(A,b)={g}", flush=True)


if __name__ == "__main__":
    run(int(sys.argv[1]) if len(sys.argv) > 1 else 120)
