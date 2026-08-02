#!/usr/bin/env python3.12
import sys

sys.path.insert(0, "/private/tmp/p25-py312-deps")
sys.path.insert(0, "/private/tmp/p25-researchtools")

import sympy as sp
from sympy import Rational
from sympy.abc import n
from ramanujantools import Matrix, Position
from ramanujantools.cmf import pFq

x = sp.symbols("x:4")
y = sp.symbols("y:3")
start = Position({x[0]: 4, x[1]: 4, x[2]: -Rational(1, 2), x[3]: 0,
                  y[0]: Rational(7, 2), y[1]: Rational(7, 2), y[2]: 4})
trajectory = Position({x[0]: 11, x[1]: 9, x[2]: -2, x[3]: -11,
                       y[0]: -1, y[1]: 0, y[2]: 2})

raw = pFq(4, 3, 1).trajectory_matrix(trajectory, start, n)
source = raw.inv().T

entries = [
    (-2*n-5)*(n+3)**2*(136*n**4+1424*n**3+5548*n**2+9551*n+6141),
    384*n**6+6384*n**5+44168*n**4+162698*n**3+336377*n**2+369933*n+169011,
    -(480*n**4+4980*n**3+19210*n**2+32690*n+20730),
    (n+2)**2*(n+3)**2*(4*n+10)*(48*n**3+386*n**2+1017*n+879),
    (n+2)**2*(-272*n**5-3848*n**4-21732*n**3-61184*n**2-85761*n-47808),
    (n+2)**2*(320*n**3+2540*n**2+6610*n+5640),
    (-4*n-10)*(n+2)**2*(n+3)**2*(32*n**4+302*n**3+1037*n**2+1530*n+813),
    (n+2)**2*(192*n**6+2984*n**5+19116*n**4+64452*n**3+120256*n**2+117279*n+46476),
    (n+2)**2*(-16*n**5-408*n**4-2912*n**3-8884*n**2-12254*n-6240),
]
delta = -2*(n+2)**2*(n+3)**2*(2*n+5)*(2*n+7)**2
target = Matrix(3, 3, entries) / delta

for label, value in (("source", source), ("target", target)):
    print(label, "companionizing", flush=True)
    companion = value.as_companion(n)
    gauge = value.companion_coboundary_matrix(n)
    print(label, "rank", companion.rows, "companion", flush=True)
    for row in range(companion.rows):
        for column in range(companion.cols):
            entry = sp.cancel(companion[row, column])
            numerator, denominator = entry.as_numer_denom()
            print(row, column, "deg", sp.degree(numerator, n),
                  sp.degree(denominator, n), sp.factor(entry), flush=True)
    print(label, "gauge shape", gauge.shape, flush=True)
    for row in range(gauge.rows):
        for column in range(gauge.cols):
            entry = sp.cancel(gauge[row, column])
            numerator, denominator = entry.as_numer_denom()
            print("g", row, column, "deg", sp.degree(numerator, n),
                  sp.degree(denominator, n), flush=True)
