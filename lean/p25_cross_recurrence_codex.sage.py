#!/usr/bin/env sage-python
"""Guess scalar recurrences for exact challenge/Wilson cross differences."""

import sys
import time

from sage.all import QQ, PolynomialRing, matrix, vector

sys.path.insert(0, "/Users/huangx/Library/SageMath-10-9/lib/python3.14/site-packages")
from ore_algebra import OreAlgebra, guess


def challenge_matrix(n):
    entries = [
        (-2*n-5)*(n+3)**2*(136*n**4+1424*n**3+5548*n**2+9551*n+6141),
        384*n**6+6384*n**5+44168*n**4+162698*n**3+336377*n**2+369933*n+169011,
        -(480*n**4+4980*n**3+19210*n**2+32690*n+20730),
        (n+2)**2*(n+3)**2*(4*n+10)*(48*n**3+386*n**2+1017*n+879),
        (n+2)**2*(-272*n**5-3848*n**4-21732*n**3-61184*n**2-85761*n-47808),
        (n+2)**2*(320*n**3+2540*n**2+6610*n+5640),
        (-4*n-10)*(n+2)**2*(n+3)**2*(32*n**4+302*n**3+1037*n**2+1530*n+813),
        (n+2)**2*(192*n**6+2984*n**5+19116*n**4+64452*n**3+120256*n**2+117279*n+46476),
        (n+2)**2*(-16*n**5-408*n**4-2912*n**3-8884*n**2-12254*n-6240),
    ]
    delta = -2*(n+2)**2*(n+3)**2*(2*n+5)*(2*n+7)**2
    return matrix(QQ, 3, 3, entries) / delta


def wilson_coefficients(m):
    a4 = 12265 + 29296*m + 26176*m**2 + 10368*m**3 + 1536*m**4
    c4 = 313 + 1904*m + 4288*m**2 + 4224*m**3 + 1536*m**4
    b10 = (
        111992515 + 1144683736*m + 5147619352*m**2
        + 13412393984*m**3 + 22433518592*m**4
        + 25185342464*m**5 + 19235018752*m**6
        + 9876373504*m**7 + 3265527808*m**8
        + 628359168*m**9 + 53477376*m**10
    )
    a = 4*(m+1)**2*(4*m+1)**2*(4*m+3)**2*a4
    c = 4*(m+2)**2*(4*m+5)**2*(4*m+7)**2*c4
    return QQ(a), QQ(b10), QQ(c)


def extend_wilson(length):
    # Pairs at z = 4m. These exact initial values come from Wilson's finite sums.
    u = [QQ(1), QQ(19)/4]
    v = [QQ(0), QQ(313)/72]
    for m in range(length - 2):
        a, b, c = wilson_coefficients(m)
        u.append((b*u[-1] - a*u[-2]) / c)
        v.append((b*v[-1] - a*v[-2]) / c)
    return u, v


COUNT = int(sys.argv[1]) if len(sys.argv) > 1 else 260
t0 = time.time()
u, v = extend_wilson(COUNT + 2)
assert u[2] == QQ(5545)/64
assert v[2] == QQ(22398503)/282240
print("Wilson terms", len(u), "seconds", time.time() - t0, flush=True)

p = vector(QQ, [30921, -32972, 8240])
q = vector(QQ, [33750, -36000, 9000])
signs = [1, -1, 1]
cross = [[] for _ in range(3)]
for N in range(COUNT):
    for j in range(3):
        cross[j].append(signs[j] * (p[j]*u[N+2] - q[j]*v[N+2]))
    transition = challenge_matrix(N)
    p *= transition
    q *= transition
print("cross terms", COUNT, "seconds", time.time() - t0, flush=True)
print("initial signs", [[x.sign() for x in seq[:12]] for seq in cross], flush=True)
for column, sequence in enumerate(cross):
    print("column", column, "first", [float(value) for value in sequence[:6]], flush=True)
    print("column", column, "tail ratios",
          [float(sequence[index+1] / sequence[index])
           for index in range(max(0, COUNT-8), COUNT-1)], flush=True)

if len(sys.argv) > 2 and sys.argv[2] == "noguess":
    sys.exit(0)

R = PolynomialRing(QQ, "n")
n = R.gen()
A = OreAlgebra(R, "Sn")


def verify(operator, sequence, start=0):
    for index in range(start, len(sequence) - operator.order()):
        if sum(operator[j](index) * sequence[index+j]
               for j in range(operator.order()+1)) != 0:
            return index
    return None


for column, sequence in enumerate(cross):
    print("guess column", column, flush=True)
    operator = guess(sequence[:COUNT*3//4], A)
    print("order", operator.order(), "degree", operator.degree(), flush=True)
    print("holdout bad", verify(operator, sequence, COUNT*3//4), flush=True)
    for shift in range(operator.order()+1):
        print("c", shift, operator[shift].factor(), flush=True)
    print("operator", operator, flush=True)
    try:
        print("factor", operator.factor(), flush=True)
    except Exception as error:
        print("factor error", repr(error), flush=True)
