#!/usr/bin/env python3
"""Exact probe for the Euler-transformed Catalan squeeze."""

from fractions import Fraction as F
import sympy as sp


def positive_matrix(n):
    return [
        [(2*n+5)*(n+3)**2*(136*n**4+1424*n**3+5548*n**2+9551*n+6141),
         384*n**6+6384*n**5+44168*n**4+162698*n**3+336377*n**2+369933*n+169011,
         480*n**4+4980*n**3+19210*n**2+32690*n+20730],
        [(n+2)**2*(n+3)**2*(4*n+10)*(48*n**3+386*n**2+1017*n+879),
         (n+2)**2*(272*n**5+3848*n**4+21732*n**3+61184*n**2+85761*n+47808),
         (n+2)**2*(320*n**3+2540*n**2+6610*n+5640)],
        [(4*n+10)*(n+2)**2*(n+3)**2*(32*n**4+302*n**3+1037*n**2+1530*n+813),
         (n+2)**2*(192*n**6+2984*n**5+19116*n**4+64452*n**3+120256*n**2+117279*n+46476),
         (n+2)**2*(16*n**5+408*n**4+2912*n**3+8884*n**2+12254*n+6240)],
    ]


def gauge(n):
    return 2*(n+2)**2*(n+3)**2*(2*n+5)*(2*n+7)**2


def advance_euler(a, odd, partial, count=6):
    # a_k, O_k, E_k on input; add c_k,...,c_{k+count-1}.
    for _ in range(count):
        partial += a*odd
        k = advance_euler.index
        a *= F(k+1, 2*k+3)
        odd += F(1, 2*k+3)
        advance_euler.index += 1
    return a, odd, partial


def values(count=36, offset=12, stride=6):
    p = [F(30921), F(32972), F(8240)]
    q = [F(33750), F(36000), F(9000)]
    a, odd, partial = F(1, 2), F(1), F(0)
    advance_euler.index = 0
    a, odd, partial = advance_euler(a, odd, partial, offset)
    output = [[] for _ in range(3)]
    aux = []
    for n in range(count):
        errors = [q[j]*partial-p[j] for j in range(3)]
        for j in range(3):
            output[j].append(errors[j])
        aux.append((a, odd, partial, p[:], q[:]))
        matrix, d = positive_matrix(n), F(gauge(n))
        p = [sum(p[i]*matrix[i][j] for i in range(3))/d for j in range(3)]
        q = [sum(q[i]*matrix[i][j] for i in range(3))/d for j in range(3)]
        a, odd, partial = advance_euler(a, odd, partial, stride)
    return output, aux


errors, aux = values()
for n in range(12):
    a, odd, partial, p, q = aux[n]
    print("n", n, "sign", [v > 0 for v in (errors[j][n] for j in range(3))],
          "scaled", [float(errors[j][n]/(q[j]*a)) for j in range(3)])

x = sp.symbols("n")
for j in range(3):
    ratios = [errors[j][n+1]/errors[j][n] for n in range(len(errors[j])-1)]
    points = [(sp.Integer(n), sp.Rational(v.numerator, v.denominator))
              for n, v in enumerate(ratios)]
    print("column", j, "ratio limits", [float(v) for v in ratios[-3:]])
    found = False
    for numdeg in range(18):
        candidate = sp.factor(sp.rational_interpolate(points[:26], numdeg, X=x))
        if all(sp.cancel(candidate.subs(x, n)-v) == 0 for n, v in points[26:]):
            print("hypergeometric", j, sp.factor(candidate))
            found = True
            break
    if not found:
        print("not hypergeometric", j)
