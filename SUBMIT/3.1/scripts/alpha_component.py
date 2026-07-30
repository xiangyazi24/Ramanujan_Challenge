#!/usr/bin/env python3
"""Decide whether the alpha component carries the discrete faithful character.

If some complex embedding of the alpha field gives Bloch-Wigner value equal to
+/- Vol(S^3_{-1/2}(7_2)), then xi_alpha is NON-torsion and the torsion route is
dead.  If instead every complex embedding gives 0, xi_alpha is plausibly torsion
and the route is alive.

This is exactly the check the audit said had to be done before interpreting the
10^-50 table, and which was never done.
"""
import snappy
import mpmath as mp

mp.mp.dps = 60

# --- the hyperbolic volume of the filled manifold -------------------------
K = snappy.Manifold('7_2')
K.dehn_fill((-1, 2))          # slope -1/2
print("filled manifold :", K)
print("volume          :", K.volume())
print("solution type   :", K.solution_type())
try:
    print("is hyperbolic   :", K.solution_type() == 'all tetrahedra positively oriented')
except Exception:
    pass
VOL = mp.mpf(repr(K.volume()))
print()

# --- the alpha field and its embeddings ----------------------------------
# f_alpha(a) = a^12 - 3a^11 + 4a^10 - 5a^9 + 6a^8 - 7a^7 + 7a^6 - 7a^5
#              + 6a^4 - 5a^3 + 4a^2 - 3a + 1
coeffs = [1, -3, 4, -5, 6, -7, 7, -7, 6, -5, 4, -3, 1]   # descending
roots = mp.polyroots(coeffs, maxsteps=500, extraprec=200)
print("alpha field: 12 roots")
reals = [r for r in roots if abs(mp.im(r)) < mp.mpf('1e-40')]
comps = [r for r in roots if abs(mp.im(r)) >= mp.mpf('1e-40')]
print("  real roots     :", [mp.nstr(mp.re(r), 12) for r in reals])
print("  complex roots  :", len(comps), "( =", len(comps)//2, "conjugate pairs )")
print()

# --- shapes and Bloch-Wigner at each embedding ---------------------------
def shapes(a):
    M = a**2
    L = a
    X = M*M
    u = (L + X**3) / (X*(L + X))
    r = -(1 + mp.sqrt(1 + 4*u**2)) / (2*u)
    tau = 1 - r**2
    return [tau, u, u/X, 1/(1 - u*X)]

def bloch_wigner(z):
    z = mp.mpc(z)
    if abs(z) < mp.mpf('1e-40') or abs(z - 1) < mp.mpf('1e-40'):
        return mp.mpf(0)
    return mp.im(mp.polylog(2, z)) + mp.arg(1 - z) * mp.log(abs(z))

print("Bloch-Wigner sum  D_sigma = sum_j D(z_j)  at each embedding:")
print("%-28s %-24s %s" % ("root a", "D_sigma", "|D_sigma| vs Vol"))
vals = []
seen = set()
for r in roots:
    key = (round(float(mp.re(r)), 25), round(abs(float(mp.im(r))), 25))
    if key in seen:
        continue
    seen.add(key)
    try:
        S = shapes(r)
        D = sum(bloch_wigner(z) for z in S)
    except Exception as e:
        print("  (failed at", mp.nstr(r, 10), ":", e, ")")
        continue
    tag = ""
    if abs(abs(D) - VOL) < mp.mpf('1e-6'):
        tag = "   <== EQUALS THE HYPERBOLIC VOLUME"
    elif abs(D) < mp.mpf('1e-30'):
        tag = "   (vanishes)"
    print("%-28s %-24s %s" % (mp.nstr(r, 14), mp.nstr(D, 12), tag))
    vals.append(D)

print()
print("Vol(S^3_{-1/2}(7_2)) =", mp.nstr(VOL, 20))
print()
nonzero = [v for v in vals if abs(v) > mp.mpf('1e-25')]
if nonzero:
    print(">>> Some embedding has NONZERO Bloch-Wigner value.")
    print(">>> xi_alpha is NON-TORSION; the torsion route is dead.")
    for v in nonzero:
        print("      D =", mp.nstr(v, 20), "   D/Vol =", mp.nstr(v/VOL, 12))
else:
    print(">>> EVERY embedding gives Bloch-Wigner 0 to the tested precision.")
    print(">>> xi_alpha is plausibly TORSION -- the torsion route is ALIVE,")
    print(">>> and what is needed is an EXACT proof of these vanishings.")
