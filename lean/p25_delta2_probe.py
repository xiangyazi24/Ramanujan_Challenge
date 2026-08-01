#!/usr/bin/env python3
"""Inspect second differences and asymptotic slopes of Delannoy coefficients."""

from fractions import Fraction as F
from math import comb
import mpmath as mp


def matrix(n):
    return [
        [(-2*n-5)*(n+3)**2*(136*n**4+1424*n**3+5548*n**2+9551*n+6141),
         384*n**6+6384*n**5+44168*n**4+162698*n**3+336377*n**2+369933*n+169011,
         -(480*n**4+4980*n**3+19210*n**2+32690*n+20730)],
        [(n+2)**2*(n+3)**2*(4*n+10)*(48*n**3+386*n**2+1017*n+879),
         (n+2)**2*(-272*n**5-3848*n**4-21732*n**3-61184*n**2-85761*n-47808),
         (n+2)**2*(320*n**3+2540*n**2+6610*n+5640)],
        [(-4*n-10)*(n+2)**2*(n+3)**2*(32*n**4+302*n**3+1037*n**2+1530*n+813),
         (n+2)**2*(192*n**6+2984*n**5+19116*n**4+64452*n**3+120256*n**2+117279*n+46476),
         (n+2)**2*(-16*n**5-408*n**4-2912*n**3-8884*n**2-12254*n-6240)],
    ]


def delta(n):
    return -2*(n+2)**2*(n+3)**2*(2*n+5)*(2*n+7)**2


def basis(n, k):
    return F(2**k*comb(2*k, k)*comb(n, k)*comb(n+k, k))


count = 90
rows = [[F(1), F(0), F(0)], [F(0), F(1), F(0)], [F(0), F(0), F(1)]]
q0 = [F(33750), F(-36000), F(9000)]
p0 = [F(30921), F(-32972), F(8240)]
qvalues = []
pvalues = []
for n in range(count+1):
    qvalues.append(sum(q0[i]*rows[i][0] for i in range(3)))
    pvalues.append(sum(p0[i]*rows[i][0] for i in range(3)))
    if n == count:
        break
    raw = matrix(n)
    d = F(delta(n))
    transition = [[F(raw[i][j], d) for j in range(3)] for i in range(3)]
    rows = [[sum(row[i]*transition[i][j] for i in range(3)) for j in range(3)]
            for row in rows]


def invert(values):
    answer = []
    for n, value in enumerate(values):
        residue = value-sum(answer[k]*basis(n, k) for k in range(n))
        answer.append(residue/basis(n, n))
    return answer


f = invert(qvalues)
g = invert(pvalues)
d2f = [f[k+2]-2*f[k+1]+f[k] for k in range(count-1)]
d2g = [g[k+2]-2*g[k+1]+g[k] for k in range(count-1)]

print("first f", f[:5])
print("first g", g[:5])
for k in range(12):
    print("d2", k, d2f[k], d2g[k], "ratios",
          d2f[k+1]/d2f[k] if d2f[k] else None,
          d2g[k+1]/d2g[k] if d2g[k] else None)

mp.mp.dps = 100
def real(value):
    return mp.mpf(value.numerator)/value.denominator

for end in (10, 20, 40, 60, 80):
    af = f[1]-f[0]+sum(d2f[:end])
    ag = g[1]-g[0]+sum(d2g[:end])
    print("slope", end, mp.nstr(real(af), 60), mp.nstr(real(ag), 60),
          mp.nstr(real(ag)/real(af), 60),
          "Gerr", mp.nstr(real(ag)/real(af)-mp.catalan, 8))

af = real(f[-1]-f[-2])
ag = real(g[-1]-g[-2])
print("last differences", mp.nstr(af, 100), mp.nstr(ag, 100))
print("ratio", mp.nstr(ag/af, 100))
print("Af constants")
for c, name in [(1, "1"), (mp.pi, "pi"), (mp.pi**2, "pi2"),
                (mp.sqrt(2), "sqrt2"), (mp.log(2), "log2"),
                (mp.ellipk(mp.mpf('0.5')), "Khalf")]:
    print(name, mp.nstr(af/c, 70))
