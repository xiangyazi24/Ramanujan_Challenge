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
start = Position({
    x[0]: 4,
    x[1]: 4,
    x[2]: -Rational(1, 2),
    x[3]: 0,
    y[0]: Rational(7, 2),
    y[1]: Rational(7, 2),
    y[2]: 4,
})
trajectory = Position({
    x[0]: 1,
    x[1]: 1,
    x[2]: -1,
    x[3]: -1,
    y[0]: 0,
    y[1]: 0,
    y[2]: 0,
})

raw_matrix = pFq(4, 3, 1).trajectory_matrix(trajectory, start, n)
matrix = sp.Matrix(raw_matrix.tolist())
print("det", sp.factor(matrix.det()), flush=True)
for row in range(3):
    for column in range(3):
        value = sp.factor(matrix[row, column])
        numerator, denominator = sp.cancel(value).as_numer_denom()
        print(row, column, "degrees", sp.degree(numerator, n),
              sp.degree(denominator, n), flush=True)
        print(value, flush=True)

challenge_entries = [
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
challenge_delta = -2*(n+2)**2*(n+3)**2*(2*n+5)*(2*n+7)**2
challenge = sp.Matrix(3, 3, challenge_entries) / challenge_delta

lam = sp.symbols("lam")
raw_balanced = sp.Matrix(3, 3, lambda i, j:
                         sp.limit(matrix[i, j] * n**(i-j), n, sp.oo))
dual_balanced = raw_balanced.inv().T
target_balanced = sp.Matrix(3, 3, lambda i, j:
                            sp.limit(challenge[i, j] * n**(j-i), n, sp.oo))
for label, value in (("raw_limit", raw_balanced),
                     ("dual_limit", dual_balanced),
                     ("target_limit", target_balanced)):
    print(label, value, flush=True)
    print(label + "_charpoly", sp.factor(value.charpoly(lam).as_expr()), flush=True)

source_companion = raw_matrix.inv().T.as_companion(n)
target_companion = (Matrix(3, 3, challenge_entries) /
                    challenge_delta).as_companion(n)
for row in range(3):
    for column in range(3):
        print("companion_difference", row, column,
              sp.factor(source_companion[row, column] -
                        target_companion[row, column]), flush=True)

for depth in (1, 2, 3, 5, 10):
    walk = sp.Matrix(pFq(4, 3, 1).trajectory_matrix(
        trajectory, start, n).walk({n: 1}, depth, {n: 0}).tolist()).inv().T
    column = walk[:, 0]
    estimate = sp.cancel((4915*column[1] + 754*column[2]) /
                         (6930*column[1] + 1260*column[2]))
    print("depth", depth, "column", [sp.N(v/column[0], 30) for v in column],
          "estimate", sp.N(estimate, 50), flush=True)
