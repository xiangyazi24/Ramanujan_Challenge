#!/usr/bin/env python3
"""Exact coefficients of G, log(2), and 1 in the Meijer recessive solution."""

from fractions import Fraction as Q
import sympy as s


def transition(n):
    n = Q(n)
    return [
        [4*(n+2)*(17*n**3+111*n**2+240*n+171)/((n+1)*(n+3)*(2*n+3)*(2*n+5)),
         (n+2)*(24*n**2+101*n+102)/((n+1)*(2*n+3)),
         (n+2)*(2*n+5)*(16*n**2+81*n+90)/(2*(n+1)*(2*n+3))],
        [(96*n**4+780*n**3+2384*n**2+3273*n+1723)/((n+1)*(n+2)*(n+3)*(2*n+3)*(2*n+5)),
         (68*n**3+398*n**2+778*n+523)/(2*(n+1)*(n+2)*(2*n+3)),
         (96*n**4+884*n**3+2970*n**2+4360*n+2403)/(4*(n+1)*(n+2)*(2*n+3))],
        [-5*(24*n**2+117*n+143)/((n+1)*(n+2)*(n+3)*(2*n+3)*(2*n+5)),
         -5*(16*n+41)/(2*(n+1)*(n+2)*(2*n+3)),
         (8*n**3-44*n**2-478*n-801)/(4*(n+1)*(n+2)*(2*n+3))],
    ]


def mm(a, b):
    return [[sum(a[i][k]*b[k][j] for k in range(3)) for j in range(3)]
             for i in range(len(a))]


initial_f = [
    [Q(150), Q(-128), Q(-146, 3)],
    [Q(24745, 4), Q(-14624, 3), Q(-823511, 360)],
    [Q(7225281, 32), Q(-886784, 5), Q(-2818419551, 33600)],
]

# Recover Y_0 coefficient rows from F-vector = Y_0 U(0).
t0 = transition(0)
t1 = transition(1)
e = [[Q(1)], [Q(0)], [Q(0)]]
u0 = [[Q(i == j) for j in range(3)] for i in range(3)]
col1 = [t0[i][0] for i in range(3)]
t01 = [[sum(t0[i][k]*t1[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
col2 = [t01[i][0] for i in range(3)]
u0 = [[Q(i == 0), col1[i], col2[i]] for i in range(3)]


def inverse(a):
    x = [row[:] + [Q(i == j) for j in range(3)] for i, row in enumerate(a)]
    for j in range(3):
        p = next(i for i in range(j, 3) if x[i][j])
        x[j], x[p] = x[p], x[j]
        v = x[j][j]
        x[j] = [z/v for z in x[j]]
        for i in range(3):
            if i != j:
                v = x[i][j]
                x[i] = [x[i][k]-v*x[j][k] for k in range(6)]
    return [row[3:] for row in x]


# Rows are F_n, columns constants; Y coefficients = U^-T? Each constant's
# F-vector is a row, so transpose the data for propagation.
f_by_constant = [[initial_f[r][c] for r in range(3)] for c in range(3)]
y = mm(f_by_constant, inverse(u0))
values = [[] for _ in range(3)]
for n in range(31):
    for c in range(3):
        values[c].append(y[c][0])  # F_n coefficient = Y_n first coordinate
    y = mm(y, transition(n))

x = s.symbols("n")
for c, label in enumerate(("G", "L2", "one")):
    print(label)
    for n in range(8):
        value = values[c][n]
        print(n, s.factor(s.Rational(value.numerator, value.denominator)))
    ratios = [values[c][n+1]/values[c][n] for n in range(30)
              if values[c][n]]
    points = [(s.Integer(n), s.Rational(v.numerator, v.denominator))
              for n, v in enumerate(ratios)]
    for degree in range(12):
        candidate = s.factor(s.rational_interpolate(points[:20], degree, X=x))
        if all(s.cancel(candidate.subs(x, n)-v) == 0 for n, v in points[20:]):
            print("hypergeometric ratio", candidate)
            break
