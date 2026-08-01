"""Temporary numerical probe for a projective cone of the Catalan error row."""

import mpmath as mp

mp.mp.dps = 1200
G = mp.catalan
error = [G * 33750 - 30921, G * 36000 - 32972, G * 9000 - 8240]
points = []


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


for n in range(101):
    a, b, c = -error[0], error[1], error[2]
    x = (n + 1)**2 * b / a
    y = (n + 1)**2 * c / a
    scale = (n + mp.mpf(5) / 2) / (n + mp.mpf(3) / 4)
    x *= scale
    y *= scale
    points.append((x, y, n))
    if n < 15 or n in [20, 30, 50, 100]:
        print(n, mp.nstr(x, 20), mp.nstr(y, 20), mp.nstr(y / x, 14),
              mp.nstr(error[0] / max(abs(v) for v in error), 8))
    if n == 100:
        break
    matrix = positive_matrix(n)
    error = [sum(error[i] * matrix[i][j] for i in range(3))
             for j in range(3)]


def cross(origin, left, right):
    return ((left[0] - origin[0]) * (right[1] - origin[1])
            - (left[1] - origin[1]) * (right[0] - origin[0]))


hull = []
for point in sorted(points[1:], key=lambda item: (item[0], item[1])):
    while len(hull) >= 2 and cross(hull[-2], hull[-1], point) <= 0:
        hull.pop()
    hull.append(point)
upper = []
for point in reversed(sorted(points[1:], key=lambda item: (item[0], item[1]))):
    while len(upper) >= 2 and cross(upper[-2], upper[-1], point) <= 0:
        upper.pop()
    upper.append(point)
print("hull indices", [point[2] for point in hull[:-1] + upper[:-1]])
