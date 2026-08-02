#!/usr/bin/env python3
"""Inspect signs of exact local order-six recurrences for Lima gaps."""

from fractions import Fraction as F


def challenge(n):
    e = [
        (2*n+5)*(n+3)**2*(136*n**4+1424*n**3+5548*n**2+9551*n+6141),
        384*n**6+6384*n**5+44168*n**4+162698*n**3+336377*n**2+369933*n+169011,
        480*n**4+4980*n**3+19210*n**2+32690*n+20730,
        (n+2)**2*(n+3)**2*(4*n+10)*(48*n**3+386*n**2+1017*n+879),
        (n+2)**2*(272*n**5+3848*n**4+21732*n**3+61184*n**2+85761*n+47808),
        (n+2)**2*(320*n**3+2540*n**2+6610*n+5640),
        (4*n+10)*(n+2)**2*(n+3)**2*(32*n**4+302*n**3+1037*n**2+1530*n+813),
        (n+2)**2*(192*n**6+2984*n**5+19116*n**4+64452*n**3+120256*n**2+117279*n+46476),
        (n+2)**2*(16*n**5+408*n**4+2912*n**3+8884*n**2+12254*n+6240),
    ]
    return [[F(e[3*i+j]) for j in range(3)] for i in range(3)]


def lima_ratio(m):
    return -F((m+1)**3*(3*m+5), (2*m+3)**3*(3*m+2))


def transition(n, kind):
    m = 2*n+4
    a = lima_ratio(m)
    s = a*lima_ratio(m+1)
    forcing = 1+a if kind == "lower" else -(a+s)
    c = challenge(n)
    return [
        [(c[i][j] if r == 0 and q == 0 else
          forcing*c[i][j] if r == 1 and q == 0 else
          s*c[i][j] if r == 1 and q == 1 else F(0))
         for q in range(2) for j in range(3)]
        for r in range(2) for i in range(3)
    ]


def matmul(a, b):
    return [[sum(a[i][k]*b[k][j] for k in range(len(b)))
             for j in range(len(b[0]))] for i in range(len(a))]


def solve(a, b):
    aug = [row[:] + [b[i]] for i, row in enumerate(a)]
    for col in range(len(aug)):
        pivot = next(r for r in range(col, len(aug)) if aug[r][col])
        aug[col], aug[pivot] = aug[pivot], aug[col]
        scale = aug[col][col]
        aug[col] = [v/scale for v in aug[col]]
        for r in range(len(aug)):
            if r != col and aug[r][col]:
                scale = aug[r][col]
                aug[r] = [aug[r][j]-scale*aug[col][j]
                          for j in range(len(aug)+1)]
    return [aug[i][-1] for i in range(len(aug))]


def coefficients(n, kind, coordinate):
    # Match Sage's Krylov columns: product(n)...product(n+k-1) e_coordinate.
    eye = [[F(int(i == j)) for j in range(6)] for i in range(6)]
    seed = [[F(int(i == coordinate))] for i in range(6)]
    columns = [seed]
    product = eye
    for shift in range(6):
        product = matmul(product, transition(n+shift, kind))
        columns.append(matmul(product, seed))
    basis = [[columns[j][i][0] for j in range(6)] for i in range(6)]
    return solve(basis, [columns[6][i][0] for i in range(6)])


for kind, coordinate in (("lower", 2), ("upper", 0)):
    print(kind)
    for n in list(range(12)) + [20, 50, 100]:
        cs = coefficients(n, kind, coordinate)
        print(n, "".join("+" if c > 0 else "-" if c < 0 else "0" for c in cs),
              [float(c) for c in cs])
