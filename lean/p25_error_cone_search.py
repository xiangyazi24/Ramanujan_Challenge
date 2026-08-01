#!/usr/bin/env python3
"""Temporary numerical search for shrinking error-cone rectangles."""


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


def bounds(n, constants):
    a, b, c, d = constants
    return (1.25-a/(n+1), 1.25-b/(n+1),
            2-c/(n+1)**2, 2-d/(n+1)**2)


def test(constants, start=2, stop=10000):
    worst = [float("inf")]*5
    where = [None]*5
    for n in range(start, stop+1):
        matrix = positive_matrix(n)
        h = (n+1)**2*(4*n+10)/(4*n+3)
        hn = (n+2)**2*(4*n+14)/(4*n+7)
        lx, ux, ly, uy = bounds(n, constants)
        lxn, uxn, lyn, uyn = bounds(n+1, constants)
        for x in (lx, ux):
            for y in (ly, uy):
                t = [h*matrix[0][j]-x*matrix[1][j]-y*matrix[2][j]
                     for j in range(3)]
                tests = [t[0], hn*(-t[1])-lxn*t[0],
                         uxn*t[0]-hn*(-t[1]),
                         hn*(-t[2])-lyn*t[0],
                         uyn*t[0]-hn*(-t[2])]
                scale = max(abs(value) for value in tests) or 1
                for i, value in enumerate(tests):
                    normalized = value/scale
                    if normalized < worst[i]:
                        worst[i] = normalized
                        where[i] = (n, x, y, value)
                    if value < 0:
                        return False, worst, where
    return True, worst, where


if __name__ == "__main__":
    candidates = [
        (1/4, 1/16, 1, 1/4),
        (1/3, 1/32, 3/2, 1/8),
        (1/2, 1/64, 2, 1/16),
    ]
    for candidate in candidates:
        print(candidate, test(candidate))

    import mpmath as mp
    mp.mp.dps = 100
    error = [mp.catalan*33750-30921, mp.catalan*36000-32972,
             mp.catalan*9000-8240]
    points = []
    for n in range(201):
        h = (n+1)**2*mp.mpf(4*n+10)/(4*n+3)
        x = h*error[1]/(-error[0])
        y = h*error[2]/(-error[0])
        u = (n+1)*(mp.mpf(5)/4-x)
        v = (n+1)**2*(2-y)
        points.append((float(u), float(v), n))
        if n in list(range(12))+[20,50,100,200]:
            print("uv", n, mp.nstr(u, 15), mp.nstr(v, 15))
        matrix = positive_matrix(n)
        error = [sum(error[i]*matrix[i][j] for i in range(3))
                 for j in range(3)]

    def cross(origin, left, right):
        return ((left[0]-origin[0])*(right[1]-origin[1])
                -(left[1]-origin[1])*(right[0]-origin[0]))

    ordered = sorted(points[1:])
    lower = []
    for point in ordered:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], point) <= 0:
            lower.pop()
        lower.append(point)
    upper = []
    for point in reversed(ordered):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], point) <= 0:
            upper.pop()
        upper.append(point)
    print("uv hull", lower[:-1]+upper[:-1])

    def hull(points2):
        ordered2 = sorted(set(points2))
        if len(ordered2) <= 1:
            return ordered2
        lo = []
        for point in ordered2:
            while len(lo) >= 2 and cross(lo[-2], lo[-1], point) <= 0:
                lo.pop()
            lo.append(point)
        hi = []
        for point in reversed(ordered2):
            while len(hi) >= 2 and cross(hi[-2], hi[-1], point) <= 0:
                hi.pop()
            hi.append(point)
        return lo[:-1]+hi[:-1]

    def map_uv(n, u, v):
        matrix = positive_matrix(n)
        h = (n+1)**2*(4*n+10)/(4*n+3)
        hn = (n+2)**2*(4*n+14)/(4*n+7)
        x = 1.25-u/(n+1)
        y = 2-v/(n+1)**2
        t = [h*matrix[0][j]-x*matrix[1][j]-y*matrix[2][j]
             for j in range(3)]
        if t[0] <= 0:
            return None
        xp = hn*(-t[1])/t[0]
        yp = hn*(-t[2])/t[0]
        return ((n+2)*(1.25-xp), (n+2)**2*(2-yp))

    polygon = [(0.07, 0.35), (0.16, 0.35), (0.16, 0.55),
               (0.14, 0.85), (0.07, 0.85)]
    for iteration in range(10):
        images = []
        failed = False
        for n in list(range(1, 100))+[100, 200, 500, 1000, 10000]:
            for u, v in polygon:
                mapped = map_uv(n, u, v)
                if mapped is None:
                    failed = True
                else:
                    images.append(mapped)
        new_polygon = hull(polygon+images)
        print("closure", iteration, "failed", failed, "vertices", new_polygon)
        if len(new_polygon) == len(polygon) and all(
                abs(a-b) < 1e-10 for p, q in zip(new_polygon, polygon)
                for a, b in zip(p, q)):
            break
        polygon = new_polygon
