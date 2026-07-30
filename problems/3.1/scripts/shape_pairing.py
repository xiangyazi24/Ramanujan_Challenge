#!/usr/bin/env python3
"""Why does the Bloch-Wigner sum vanish at every complex embedding?

D(z) is real-analytic with  D(zbar) = -D(z)  and  D(z) = 0 for real z.
So the sum over the four shapes vanishes identically if, at each embedding,
the shapes are either real or occur in complex-conjugate pairs.

If that is the mechanism it is an EXACT structural fact, not a numerical
coincidence, and it proves xi_alpha is torsion.
"""
import mpmath as mp
mp.mp.dps = 60

coeffs = [1, -3, 4, -5, 6, -7, 7, -7, 6, -5, 4, -3, 1]
roots = mp.polyroots(coeffs, maxsteps=500, extraprec=300)

def shapes(a):
    M = a**2; L = a; X = M*M
    u = (L + X**3) / (X*(L + X))
    r = -(1 + mp.sqrt(1 + 4*u**2)) / (2*u)
    tau = 1 - r**2
    return {'T': tau, 'U': u, 'V': u/X, 'W': 1/(1 - u*X)}

def D(z):
    z = mp.mpc(z)
    if abs(z) < mp.mpf('1e-40') or abs(z-1) < mp.mpf('1e-40'):
        return mp.mpf(0)
    return mp.im(mp.polylog(2, z)) + mp.arg(1-z)*mp.log(abs(z))

def close(p, q, tol=mp.mpf('1e-30')):
    return abs(p-q) < tol

print("For each embedding: are the four shapes real, or conjugate-paired?")
print()
for a in roots:
    S = shapes(a)
    names = list(S)
    vals = [S[n] for n in names]
    imags = [abs(mp.im(v)) for v in vals]
    nreal = sum(1 for t in imags if t < mp.mpf('1e-30'))
    # look for conjugate pairs among the non-real ones
    nonreal = [(n, v) for n, v in S.items() if abs(mp.im(v)) >= mp.mpf('1e-30')]
    pairs = []
    used = set()
    for i in range(len(nonreal)):
        if i in used: continue
        for j in range(i+1, len(nonreal)):
            if j in used: continue
            if close(nonreal[i][1], mp.conj(nonreal[j][1])):
                pairs.append((nonreal[i][0], nonreal[j][0])); used.add(i); used.add(j)
                break
    tot = sum(D(v) for v in vals)
    unpaired = [nonreal[k][0] for k in range(len(nonreal)) if k not in used]
    print("a = %-34s  real:%d  conj-pairs:%s  unpaired:%s   sum D = %s"
          % (mp.nstr(a, 12), nreal,
             pairs if pairs else "-", unpaired if unpaired else "-",
             mp.nstr(tot, 6)))

print()
print("individual D values at one complex embedding, to see the cancellation:")
a = [r for r in roots if abs(mp.im(r)) > 1e-30][0]
S = shapes(a)
for n, v in S.items():
    print("   %s = %-40s   D = %s" % (n, mp.nstr(v, 14), mp.nstr(D(v), 14)))
print("   sum =", mp.nstr(sum(D(v) for v in S.values()), 14))
