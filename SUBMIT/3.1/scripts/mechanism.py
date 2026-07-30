#!/usr/bin/env python3
"""Verify the structural mechanism forcing D(V) + D(W) = 0 at complex embeddings.

Claim.  At every complex embedding of the alpha field:
   (i)   |a| = 1, hence X = a^4 satisfies conj(X) = 1/X;
   (ii)  u is real;
   (iii) T and U are real, so D(T) = D(U) = 0;
   (iv)  therefore conj(V) = conj(u/X) = u*X  and  W = 1/(1 - u X) = 1/(1 - conj(V));
   (v)   Bloch-Wigner satisfies D(1/(1-z)) = D(z) and D(conj z) = -D(z),
         so D(W) = D(conj V) = -D(V).
   Hence sum_j D(z_j) = 0 EXACTLY, not numerically.
"""
import mpmath as mp
mp.mp.dps = 60

coeffs = [1, -3, 4, -5, 6, -7, 7, -7, 6, -5, 4, -3, 1]
roots = mp.polyroots(coeffs, maxsteps=500, extraprec=300)

def data(a):
    M = a**2; L = a; X = M*M
    u = (L + X**3) / (X*(L + X))
    r = -(1 + mp.sqrt(1 + 4*u**2)) / (2*u)
    return dict(X=X, u=u, T=1-r**2, U=u, V=u/X, W=1/(1-u*X))

def D(z):
    z = mp.mpc(z)
    if abs(z) < mp.mpf('1e-40') or abs(z-1) < mp.mpf('1e-40'): return mp.mpf(0)
    return mp.im(mp.polylog(2, z)) + mp.arg(1-z)*mp.log(abs(z))

tol = mp.mpf('1e-40')
comp = [a for a in roots if abs(mp.im(a)) > tol]
print("complex embeddings:", len(comp))
print()
print("%-3s %-12s %-12s %-12s %-14s %s" %
      ("#", "| |a|-1 |", "Im(u)", "Im(T)", "|W - 1/(1-Vbar)|", "D(V)+D(W)"))
allok = True
for i, a in enumerate(comp):
    d = data(a)
    X, u, V, W, T = d['X'], d['u'], d['V'], d['W'], d['T']
    c1 = abs(abs(a) - 1)
    c2 = abs(mp.im(u))
    c3 = abs(mp.im(T))
    pred = 1/(1 - mp.conj(V))
    c4 = abs(W - pred)
    s  = D(V) + D(W)
    ok = c1 < tol and c2 < tol and c3 < tol and c4 < tol
    allok = allok and ok
    print("%-3d %-12s %-12s %-12s %-14s %s" %
          (i, mp.nstr(c1,4), mp.nstr(c2,4), mp.nstr(c3,4), mp.nstr(c4,4), mp.nstr(s,6)))

print()
print("all four structural conditions hold at every complex embedding:", allok)
print()
print("check conj(X) == 1/X on the unit circle:")
for a in comp[:3]:
    X = data(a)['X']
    print("   |conj(X) - 1/X| =", mp.nstr(abs(mp.conj(X) - 1/X), 6))
print()
print("check the Bloch-Wigner identity D(1/(1-z)) == D(z) numerically:")
for z in [mp.mpc(0.3, 0.7), mp.mpc(-1.2, 0.4), mp.mpc(2.5, -1.1)]:
    print("   z =", mp.nstr(z,8), "  D(z) =", mp.nstr(D(z),10),
          "  D(1/(1-z)) =", mp.nstr(D(1/(1-z)),10))
