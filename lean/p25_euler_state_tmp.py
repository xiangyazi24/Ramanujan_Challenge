#!/usr/bin/env python3
from fractions import Fraction as F


def matrix(n):
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


def divisor(n):
    return 2*(n+2)**2*(n+3)**2*(2*n+5)*(2*n+7)**2


p = list(map(F, (30921, 32972, 8240)))
q = list(map(F, (33750, 36000, 9000)))
a, odd, partial = F(1, 2), F(1), F(0)
k = 0
for _ in range(12):
    partial += a*odd
    a *= F(k+1, 2*k+3)
    odd += F(1, 2*k+3)
    k += 1

for n in range(101):
    c = a*odd
    x = [q[j]*partial-p[j] for j in range(3)]
    if n < 12 or n in (20, 30, 40, 60, 80, 100):
        print(n, "odd", float(odd), "x/t", [float(x[j]/(q[j]*c)) for j in range(3)],
              "qproj", [float(q[j]/q[0]) for j in range(3)])
    m, d = matrix(n), F(divisor(n))
    p = [sum(p[i]*m[i][j] for i in range(3))/d for j in range(3)]
    q = [sum(q[i]*m[i][j] for i in range(3))/d for j in range(3)]
    for _ in range(6):
        partial += a*odd
        a *= F(k+1, 2*k+3)
        odd += F(1, 2*k+3)
        k += 1
