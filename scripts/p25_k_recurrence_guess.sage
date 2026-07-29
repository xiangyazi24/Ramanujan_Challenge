#!/usr/bin/env sage
"""
P2.5: ore_algebra.guess() on the Delannoy-basis coefficients f(k).
Much stronger structured guesser than the naive (r,d) grid.
"""
import sys, time
from ore_algebra import OreAlgebra, guess

KMAX = 400

def M_entries(n):
    m11 = (-2*n-5)*(n+3)**2*(136*n**4+1424*n**3+5548*n**2+9551*n+6141)
    m12 = 384*n**6+6384*n**5+44168*n**4+162698*n**3+336377*n**2+369933*n+169011
    m13 = -(480*n**4+4980*n**3+19210*n**2+32690*n+20730)
    m21 = (n+2)**2*(n+3)**2*(4*n+10)*(48*n**3+386*n**2+1017*n+879)
    m22 = (n+2)**2*(-272*n**5-3848*n**4-21732*n**3-61184*n**2-85761*n-47808)
    m23 = (n+2)**2*(320*n**3+2540*n**2+6610*n+5640)
    m31 = (-4*n-10)*(n+2)**2*(n+3)**2*(32*n**4+302*n**3+1037*n**2+1530*n+813)
    m32 = (n+2)**2*(192*n**6+2984*n**5+19116*n**4+64452*n**3+120256*n**2+117279*n+46476)
    m33 = (n+2)**2*(-16*n**5-408*n**4-2912*n**3-8884*n**2-12254*n-6240)
    return [[m11,m12,m13],[m21,m22,m23],[m31,m32,m33]]

def delta_H(n):
    return -2*(n+2)**2*(n+3)**2*(2*n+5)*(2*n+7)**2

t0 = time.time()
r1 = vector(QQ, [1, 0, 0])
Qe1 = [QQ(1)]
for N in range(KMAX):
    M = matrix(QQ, M_entries(N)) / delta_H(N)
    r1 = r1 * M
    Qe1.append(r1[0])
print("trajectory done (%.1fs)" % (time.time()-t0)); sys.stdout.flush()

def Bsummand(N, k):
    return 2**k * binomial(2*k, k) * binomial(N, k) * binomial(N+k, k)

f = []
for K in range(KMAX+1):
    rhs = Qe1[K]
    for k in range(K):
        rhs -= f[k] * Bsummand(K, k)
    f.append(rhs / Bsummand(K, K))
assert f[1] == 5749/3136
print("inversion done (%.1fs)" % (time.time()-t0)); sys.stdout.flush()

R = QQ['k']
A = OreAlgebra(R, 'Sk')
for cut in (150, 250, KMAX+1):
    data = f[:cut]
    try:
        t1 = time.time()
        L = guess(data, A)
        print("FOUND with %d terms (%.1fs):" % (cut, time.time()-t1))
        print("  order =", L.order(), " degree =", L.degree())
        print(L)
        lc = L.leading_coefficient()
        print("  leading coeff factored:", lc.factor())
        # Poincare polynomial
        coeffs = [L[i] for i in range(L.order()+1)]
        D = max(c.degree() for c in coeffs)
        xi = polygen(QQ, 'xi')
        poin = sum(c[D] * xi**i for i, c in enumerate(coeffs))
        print("  Poincare polynomial:", poin, "=", poin.factor())
        print("  roots:", poin.roots(CC, multiplicities=False))
        break
    except ValueError as e:
        print("no operator with %d terms (%.1fs): %s"
              % (cut, time.time()-t1, e))
        sys.stdout.flush()
print("total %.1fs" % (time.time()-t0))
