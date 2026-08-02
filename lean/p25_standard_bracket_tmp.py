#!/usr/bin/env python3
"""Search exact moving standard-Catalan partial-sum brackets."""

from fractions import Fraction as F


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


count = 120
partial = [F(0)]
for k in range(12*count + 100):
    partial.append(partial[-1] + F((-1)**k, (2*k+1)**2))

p = [F(30921), F(32972), F(8240)]
q = [F(33750), F(36000), F(9000)]
data = []
for n in range(count):
    data.append(([p[j]/q[j] for j in range(3)], p[:], q[:]))
    matrix = positive_matrix(n)
    p = [sum(p[i]*matrix[i][j] for i in range(3)) for j in range(3)]
    q = [sum(q[i]*matrix[i][j] for i in range(3)) for j in range(3)]

for slope in range(0, 13):
    for offset in range(0, 31):
        indices = [slope*n+offset for n in range(count)]
        if any(index % 2 for index in indices):
            continue
        lower = [ratios[2] <= partial[index]
                 for (ratios, _, _), index in zip(data, indices)]
        if all(lower):
            print("LOWER", slope, offset)
for slope in range(0, 13):
    for offset in range(0, 31):
        indices = [slope*n+offset for n in range(count)]
        if any(index % 2 == 0 for index in indices):
            continue
        upper = [partial[index] <= ratios[0]
                 for (ratios, _, _), index in zip(data, indices)]
        if all(upper):
            print("UPPER", slope, offset)

print("sample")
for n in range(12):
    ratios, _, _ = data[n]
    candidates = []
    for index in range(2*n, 6*n+31):
        if index % 2 == 0 and ratios[2] <= partial[index]:
            lower_gap = partial[index]-ratios[2]
            break
    for index2 in range(2*n+1, 6*n+32):
        if index2 % 2 == 1 and partial[index2] <= ratios[0]:
            upper_gap = ratios[0]-partial[index2]
            break
    print(n, index, index2, float(lower_gap), float(upper_gap))
