#!/usr/bin/env python3
"""Exact probe for a parity-folded Wilson/challenge gauge."""

from fractions import Fraction as F
import ast
import sympy as sp


source = open("p25_cross_fit_codex.py").read()
tree = ast.parse(source)
keep = []
for node in tree.body:
    if isinstance(node, (ast.Import, ast.ImportFrom)) or (
        isinstance(node, ast.FunctionDef) and node.name == "challenge_matrix"
    ):
        keep.append(node)
namespace = {}
exec(compile(ast.Module(keep, type_ignores=[]), "probe", "exec"), namespace)
challenge_matrix = namespace["challenge_matrix"]


def wilson_coefficients(m):
    a4 = 12265 + 29296*m + 26176*m**2 + 10368*m**3 + 1536*m**4
    c4 = 313 + 1904*m + 4288*m**2 + 4224*m**3 + 1536*m**4
    b10 = (
        111992515 + 1144683736*m + 5147619352*m**2
        + 13412393984*m**3 + 22433518592*m**4
        + 25185342464*m**5 + 19235018752*m**6
        + 9876373504*m**7 + 3265527808*m**8
        + 628359168*m**9 + 53477376*m**10
    )
    a = 4*(m+1)**2*(4*m+1)**2*(4*m+3)**2*a4
    c = 4*(m+2)**2*(4*m+5)**2*(4*m+7)**2*c4
    return F(a), F(b10), F(c)


u = [F(1), F(19, 4)]
v = [F(0), F(313, 72)]
for m in range(70):
    a, b, c = wilson_coefficients(m)
    u.append((b*u[-1] - a*u[-2]) / c)
    v.append((b*v[-1] - a*v[-2]) / c)

p = [F(30921), F(-32972), F(8240)]
q = [F(33750), F(-36000), F(9000)]
gauge = F(1)
challenge = []
for N in range(80):
    challenge.append(([q[j] / gauge for j in range(3)],
                      [p[j] / gauge for j in range(3)]))
    matrix = challenge_matrix(N)
    p = [sum(p[i] * matrix[i][j] for i in range(3)) for j in range(3)]
    q = [sum(q[i] * matrix[i][j] for i in range(3)) for j in range(3)]
    gauge *= F(-2*(N+2)**2*(N+3)**2*(2*N+5)*(2*N+7)**2)


x = sp.symbols("r")


def rational_fit(values, train=22):
    points = [(sp.Integer(i), sp.Rational(y.numerator, y.denominator))
              for i, y in enumerate(values)]
    for numerator_degree in range(train):
        candidate = sp.cancel(
            sp.rational_interpolate(points[:train], numerator_degree, x)
        )
        if all(sp.cancel(candidate.subs(x, i) - y) == 0
               for i, y in points[train:]):
            return sp.factor(candidate)
    return None


for phase in (0, 1):
    for offset in range(7):
        coefficients = [[] for _ in range(6)]
        for r in range(32):
            N = 2*r + phase
            m = r + offset
            determinant = u[m]*v[m+1] - u[m+1]*v[m]
            qq, pp = challenge[N]
            for j in range(3):
                coefficients[2*j].append(
                    (qq[j]*v[m+1] - u[m+1]*pp[j]) / determinant
                )
                coefficients[2*j+1].append(
                    (u[m]*pp[j] - qq[j]*v[m]) / determinant
                )
        base = coefficients[0]
        quotient_fits = [
            rational_fit([coefficients[j][r] / base[r] for r in range(32)])
            for j in range(1, 6)
        ]
        if any(value is not None for value in quotient_fits):
            print("phase", phase, "offset", offset,
                  "mask", [value is not None for value in quotient_fits],
                  flush=True)
            for j, value in enumerate(quotient_fits, start=1):
                if value is not None:
                    print(j, value, flush=True)
