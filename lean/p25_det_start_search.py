#!/usr/bin/env python3
"""Temporary determinant search for a direct Catalan 3F2 trajectory."""

import itertools
from fractions import Fraction

import sympy as sp

n = sp.symbols("n")
target = (n+1)*(n+sp.Rational(3, 2))**2 / (
    (n+3)*(n+sp.Rational(7, 2))**2)


def determinant(a, c, d, e):
    x = a-2*n
    z = c+2*n

    def f(u):
        return (u-d+1)*(u-e+1)/(2*u**2)

    return sp.cancel(f(z)*f(z+1)/(f(x-1)*f(x-2)))


def value_at(a, c, d, e, index):
    x = a-2*index
    z = c+2*index

    def f(u):
        return (u-d+1)*(u-e+1)/(2*u*u)

    return f(z)*f(z+1)/(f(x-1)*f(x-2))


sample_indices = (0, 1, 2, 3)
target_values = [Fraction(2*k+2, 2*k+6)*Fraction(2*k+3, 2*k+7)**2
                 for k in sample_indices]
target_floats = list(map(float, target_values))

for residues in itertools.product((Fraction(0), Fraction(1, 2)), repeat=4):
  for offsets in itertools.product(range(-4, 5), repeat=4):
    aa, cc, dd, ee = offsets
    a = float(residues[0]+aa)
    c = float(residues[1]+cc)
    d = float(residues[2]+dd)
    e = float(residues[3]+ee)
    try:
        values = [value_at(a, c, d, e, k) for k in sample_indices]
        if any(abs(value-wanted) > 1e-10
               for value, wanted in zip(values, target_floats)):
            continue
    except ZeroDivisionError:
        continue
    a = residues[0]+aa
    c = residues[1]+cc
    d = residues[2]+dd
    e = residues[3]+ee
    if [value_at(a, c, d, e, k) for k in sample_indices] != target_values:
        continue
    a, c, d, e = map(sp.Rational, (a, c, d, e))
    candidate = determinant(a, c, d, e)
    if sp.cancel(candidate-target) == 0:
        print(residues, offsets, flush=True)
