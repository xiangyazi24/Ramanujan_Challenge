#!/usr/bin/env python3
"""Exact experiments for a moving cone around the neutral P2.5 solution."""

from fractions import Fraction as F
import sympy as s
import mpmath as mp


n = s.symbols("n", integer=True, nonnegative=True)
R = (1806156 + 3929280*n + 3483853*n**2 + 1615610*n**3
     + 414064*n**4 + 55680*n**5 + 3072*n**6)
S = s.expand(R.subs(n, n + 1))
T = (78407415225 + 226466421477*n + 286416594222*n**2
     + 208287850700*n**3 + 96028512072*n**4
     + 29119642544*n**5 + 5810223840*n**6
     + 735843584*n**7 + 53692416*n**8 + 1720320*n**9)
U = (195670909710 + 603914277213*n + 825643834707*n**2
     + 659435701854*n**3 + 341146546318*n**4
     + 119560001580*n**5 + 28770039448*n**6
     + 4696686576*n**7 + 498082432*n**8
     + 30999552*n**9 + 860160*n**10)

p0 = -4*(n+1)*(n+3)**2*(n+4)*(2*n+3)*(2*n+5)*(2*n+7)*S
p1 = (n+2)**2*(n+4)*(2*n+7)*T
p2 = -2*(n+2)*(n+3)*(2*n+9)*U
p3 = 4*(n+2)*(n+3)*(n+4)**2*(2*n+7)*(2*n+11)**2*R

c = s.factor(-(p0+p1+p2+p3)/p3)
minus_d = s.factor(3+(p1+2*p2)/p3)
e = s.factor(-2-p2/p3)


def matrix(k):
    return [
        [(2*k+5)*(k+3)**2*(136*k**4+1424*k**3+5548*k**2+9551*k+6141),
         384*k**6+6384*k**5+44168*k**4+162698*k**3+336377*k**2+369933*k+169011,
         480*k**4+4980*k**3+19210*k**2+32690*k+20730],
        [(k+2)**2*(k+3)**2*(4*k+10)*(48*k**3+386*k**2+1017*k+879),
         (k+2)**2*(272*k**5+3848*k**4+21732*k**3+61184*k**2+85761*k+47808),
         (k+2)**2*(320*k**3+2540*k**2+6610*k+5640)],
        [(4*k+10)*(k+2)**2*(k+3)**2*(32*k**4+302*k**3+1037*k**2+1530*k+813),
         (k+2)**2*(192*k**6+2984*k**5+19116*k**4+64452*k**3+120256*k**2+117279*k+46476),
         (k+2)**2*(16*k**5+408*k**4+2912*k**3+8884*k**2+12254*k+6240)],
    ]


def divisor(k):
    return 2*(k+2)**2*(k+3)**2*(2*k+5)*(2*k+7)**2


def row_step(row, k):
    m = matrix(k)
    d0 = F(divisor(k))
    return [sum(row[i]*m[i][j] for i in range(3))/d0 for j in range(3)]


def numeric_values(count=120):
    mp.mp.dps = 160
    g = mp.catalan
    p = list(map(F, (30921, 32972, 8240)))
    q = list(map(F, (33750, 36000, 9000)))
    values = []
    for k in range(count + 3):
        pp = mp.mpf(p[0].numerator)/p[0].denominator
        qq = mp.mpf(q[0].numerator)/q[0].denominator
        values.append(pp-g*qq)
        p = row_step(p, k)
        q = row_step(q, k)
    return values


if __name__ == "__main__":
    print("c =", c)
    print("d =", s.factor(-minus_d))
    print("e =", e)
    vals = numeric_values()
    for k in list(range(16)) + [20, 30, 50, 80, 100, 120]:
        z0 = vals[k]
        z1 = vals[k]-vals[k+1]
        z2 = vals[k]-2*vals[k+1]+vals[k+2]
        z3 = vals[k]-3*vals[k+1]+3*vals[k+2]-vals[k+3]
        print(k, "a", mp.nstr(z0, 12),
              "u=n*z1/z0", mp.nstr((k+1)*z1/z0, 12),
              "v=n2*z2/z0", mp.nstr((k+1)**2*z2/z0, 12),
              "w=n3*z3/z0", mp.nstr((k+1)**3*z3/z0, 12))
