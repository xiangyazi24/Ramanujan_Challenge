#!/usr/bin/env python3
"""Exact bracket experiments with the Lima--Guillera Catalan partial sums."""

from fractions import Fraction as F
from math import comb


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


def term(m):
    return F((-1)**m*(3*m+2)*8**m,
             2*(2*m+1)**3*comb(2*m, m)**3)


partials = [F(0)]
for m in range(250):
    partials.append(partials[-1]+term(m))


def trajectories(count):
    p = [F(30921), F(32972), F(8240)]
    q = [F(33750), F(36000), F(9000)]
    answer = []
    for n in range(count):
        answer.append((p[:], q[:]))
        matrix = positive_matrix(n)
        p = [sum(p[i]*matrix[i][j] for i in range(3)) for j in range(3)]
        q = [sum(q[i]*matrix[i][j] for i in range(3)) for j in range(3)]
    return answer


data = trajectories(70)
for lower_shift in range(0, 16, 2):
    lower = [q[2]*partials[2*n+lower_shift]-p[2] for n, (p, q) in enumerate(data)]
    print("lower shift", lower_shift, "all", all(x >= 0 for x in lower),
          "first bad", next((n for n, x in enumerate(lower) if x < 0), None),
          "min float", float(min(lower)))

for upper_shift in range(1, 16, 2):
    upper = [p[0]-q[0]*partials[2*n+upper_shift] for n, (p, q) in enumerate(data)]
    print("upper shift", upper_shift, "all", all(x >= 0 for x in upper),
          "first bad", next((n for n, x in enumerate(upper) if x < 0), None),
          "min float", float(min(upper)))

lower_shift, upper_shift = 4, 5
print("\nscaled states for shifts", lower_shift, upper_shift)
for n, (p, q) in enumerate(data[:20]):
    lo = [q[j]*partials[2*n+lower_shift]-p[j] for j in range(3)]
    up = [p[j]-q[j]*partials[2*n+upper_shift] for j in range(3)]
    scale = abs(term(2*n+lower_shift))
    print(n,
          "lo/qterm", [float(x/(q[j]*scale)) for j, x in enumerate(lo)],
          "up/qterm", [float(x/(q[j]*scale)) for j, x in enumerate(up)])

print("\nendpoint movements / Lima interval movements")
for n in range(20):
    p, q = data[n]
    pn, qn = data[n+1]
    lowmove = p[2]/q[2]
    lowmove = pn[2]/qn[2]-lowmove
    highmove = p[0]/q[0]-pn[0]/qn[0]
    m = 2*n+4
    lower_increment = partials[m+2]-partials[m]
    upper_decrement = partials[m+1]-partials[m+3]
    r0, r2 = p[0]/q[0], p[2]/q[2]
    theta_lower = (partials[m]-r2)/(r0-r2)
    theta_upper = (partials[m+1]-r2)/(r0-r2)
    print(n, float(lowmove/lower_increment), float(highmove/upper_decrement),
          "theta", float(theta_lower), float(theta_upper))
