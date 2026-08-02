#!/usr/bin/env python3
"""Left-pole residue terms for the P2.5 Meijer-G recessive solution."""

import sympy as s

n, k = s.symbols("n k", integer=True, nonnegative=True)


def bracket(nn, kk):
    return s.simplify(
        -s.digamma(nn+kk+s.Rational(7, 2))
        -s.digamma(nn+kk+s.Rational(5, 2))
        -s.digamma(nn+kk+3)-s.digamma(nn+kk+2)
        +s.digamma(2*nn+kk+s.Rational(9, 2))
        +s.digamma(2*nn+kk+5)+2*s.digamma(kk+1))


for nn in range(5):
    print("n", nn)
    for kk in range(8):
        b = s.expand_func(bracket(nn, kk)).simplify()
        print(kk, b, float(b.evalf()), "minus", float((-b).evalf()))
