#!/usr/bin/env python3
"""Identify the exact algebraic relation between the shapes V and W that forces
D(V) + D(W) = 0, and check T, U are real at every embedding.

Bloch-Wigner six-fold symmetry:
    D(z) = D(1 - 1/z) = D(1/(1-z))            (sign preserved)
    D(1/z) = D(1-z) = D(z/(z-1)) = -D(z)      (sign reversed)
So D(V) + D(W) = 0 identically if W is any of  1/V,  1-V,  V/(V-1).
"""
import mpmath as mp
mp.mp.dps = 60

coeffs = [1, -3, 4, -5, 6, -7, 7, -7, 6, -5, 4, -3, 1]
roots = mp.polyroots(coeffs, maxsteps=500, extraprec=300)

def shapes(a):
    M = a**2; L = a; X = M*M
    u = (L + X**3) / (X*(L + X))
    r = -(1 + mp.sqrt(1 + 4*u**2)) / (2*u)
    return {'T': 1 - r**2, 'U': u, 'V': u/X, 'W': 1/(1 - u*X), 'X': X, 'u': u}

cands = {
    "1/V"        : lambda v: 1/v,
    "1-V"        : lambda v: 1-v,
    "V/(V-1)"    : lambda v: v/(v-1),
    "1-1/V"      : lambda v: 1-1/v,
    "1/(1-V)"    : lambda v: 1/(1-v),
    "conj(V)"    : lambda v: mp.conj(v),
}

print("Testing W = phi(V) at every embedding of the alpha field:")
print()
hits = {k: 0 for k in cands}
n = 0
for a in roots:
    S = shapes(a)
    V, W = S['V'], S['W']
    n += 1
    for k, f in cands.items():
        try:
            if abs(f(V) - W) < mp.mpf('1e-40'):
                hits[k] += 1
        except Exception:
            pass

for k, c in hits.items():
    mark = "  <== HOLDS AT ALL EMBEDDINGS" if c == n else ""
    print("  W = %-10s : %2d/%d%s" % (k, c, n, mark))

print()
print("Is T real at every embedding?  U real?")
tr = ur = 0
for a in roots:
    S = shapes(a)
    if abs(mp.im(S['T'])) < mp.mpf('1e-40'): tr += 1
    if abs(mp.im(S['U'])) < mp.mpf('1e-40'): ur += 1
print("   T real at %d/%d embeddings" % (tr, len(roots)))
print("   U real at %d/%d embeddings" % (ur, len(roots)))

print()
print("Algebraic identity check:  W = 1/(1-uX),  V = u/X.")
print("If  W = V/(V-1)  then  1/(1-uX) = (u/X)/((u/X)-1) = u/(u-X),")
print("i.e.  u - X = u(1 - uX) = u - u^2 X,  i.e.  -X = -u^2 X,  i.e.  u^2 = 1.")
print("If  W = 1/V     then  1/(1-uX) = X/u,  i.e.  u = X - u X^2.")
print("If  W = 1-V     then  1/(1-uX) = 1 - u/X = (X-u)/X,")
print("i.e.  X = (X-u)(1-uX) = X - uX^2 - u + u^2 X,")
print("i.e.  0 = -uX^2 - u + u^2 X,  i.e.  u X^2 + u = u^2 X,  i.e.  X^2 + 1 = u X.")
print()
for a in roots[:3]:
    S = shapes(a)
    u, X = S['u'], S['X']
    print("  a =", mp.nstr(a, 10))
    print("     u^2 - 1        =", mp.nstr(u**2 - 1, 8))
    print("     u - X + u X^2  =", mp.nstr(u - X + u*X**2, 8))
    print("     X^2 + 1 - u X  =", mp.nstr(X**2 + 1 - u*X, 8))
