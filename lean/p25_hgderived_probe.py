#!/usr/bin/env python3
"""Exact probe of the derived 3F2 trajectory matching P2.5."""

import sympy as sp

from ramanujantools import Position
from ramanujantools.cmf import known_cmfs


n = sp.Symbol("n")
x0, x1, x2, y0, y1 = sp.symbols("x0 x1 x2 y0 y1")

cmf = known_cmfs.hypergeometric_derived_3F2().subs({sp.Symbol("z"): -1})
base = Position(
    {
        x0: sp.Rational(1, 2),
        x1: sp.Rational(1, 2),
        x2: 1,
        y0: sp.Rational(1, 2),
        y1: sp.Rational(1, 2),
    }
)
trajectory = {x0: 0, x1: 0, x2: 0, y0: 2, y1: -2}
T = cmf.trajectory_matrix(trajectory, base, n)
print("T=")
print(sp.sstr(T))
print("det", sp.factor(T.det()))
print("char infinity")
d = max(sp.degree(sp.together(entry).as_numer_denom()[0], n) - sp.degree(sp.together(entry).as_numer_denom()[1], n) for entry in T)
Tinfty = T.applyfunc(lambda entry: sp.limit(entry / n**d, n, sp.oo))
print("degree", d)
print(Tinfty)
print(sp.factor(Tinfty.charpoly().as_expr()))
