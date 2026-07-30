#!/usr/bin/env python3
"""SUPERSEDED -- kept only as a record of an argument that turned out to be wrong.

This script supported the claim that the Dehn-filling correction makes no
contribution because lambda_core is purely imaginary, hence lambda_core^2 is
real.  That reasoning is INVALID: Neumann's correction is LINEAR in lambda_core,
not quadratic (the quadratic lambda^2/4 belongs to the Neumann-Zagier deformation
asymptotics, not the exact surgery formula).  For lambda_core = i*theta the
correction is -theta*pi/2, which is real and not in general a rational multiple
of pi^2.

The claim it was meant to support is now established a different way, and the
correct script is  flattening_ambiguity.py :  at both endpoints all four shapes
are REAL, so ANY flattening change moves Re[Delta R]/pi^2 by a half-integer.  The
denominator bound therefore only doubles, and no flattening has to be exhibited.

Original docstring follows.
"""

"""Does the Dehn-filling correction contribute to Bloch-Wigner at conjugate
embeddings?

Neumann Thm 14.5: for a filled cusp the closed class is the relative four-shape
class plus a solid-torus term whose extended-Bloch image is governed by the core
holonomy lambda_core.  Its Bloch-Wigner (volume) contribution is what matters
for the torsion test.

Decisive empirical test: for a REAL representation the total volume is 0.  At a
Galois conjugate embedding sigma, the total volume is vol(sigma rho).  If the
four-shape sum already gives 0 at every embedding (proved structurally), then
the filling term contributes 0 exactly when vol(sigma rho) = 0 for all sigma.

We test this by computing, at each embedding, the peripheral data:
  M = a^2, L = a  (alpha chart), core holonomy for slope (-1,2).
The correction's volume contribution vanishes iff lambda_core is REAL.
"""
import mpmath as mp
mp.mp.dps = 60

coeffs = [1, -3, 4, -5, 6, -7, 7, -7, 6, -5, 4, -3, 1]
roots = mp.polyroots(coeffs, maxsteps=800, extraprec=400)

print("slope (p,q) = (-1,2):  filling relation  L^2 M^{-1} = 1")
print("dual pair (gamma,delta) with p*delta - q*gamma = 1:  (0,-1)")
print("core holonomy  lambda_core = gamma*log M + delta*log L = -log L")
print()
print("%-40s %-24s %-24s" % ("embedding a", "log L = log a", "Im(lambda_core)"))
allreal = True
for a in roots:
    L = a
    logL = mp.log(mp.mpc(L))
    lam = -logL
    im = mp.im(lam)
    # for a real positive L the log is real; on the unit circle log a is purely imaginary
    tag = ""
    if abs(im) < mp.mpf('1e-40'):
        tag = "  real"
    else:
        tag = "  PURELY IMAGINARY" if abs(mp.re(lam)) < mp.mpf('1e-40') else "  complex"
        allreal = False
    print("%-40s %-24s %-24s%s" % (mp.nstr(a, 12), mp.nstr(logL, 12), mp.nstr(im, 10), tag))

print()
print("On the unit circle a = e^{i theta}, so log a = i theta is PURELY IMAGINARY,")
print("hence lambda_core is purely imaginary at every complex embedding.")
print()
print("Neumann's correction term has the shape  (1/2)(lambda_core)^2 / (2 pi i)  or")
print("similar; a purely imaginary lambda_core makes lambda_core^2 REAL NEGATIVE,")
print("so the correction is REAL -- it contributes to Chern-Simons, NOT to volume.")
print()
print("Consistency check: the total volume of a conjugate representation.")
print("If the four-shape Bloch-Wigner sum is 0 (proved) and the filling term is")
print("real, then vol(sigma rho) = 0 at every embedding, i.e. the closed class is")
print("torsion.  Test: is lambda_core^2 real at every complex embedding?")
for a in roots:
    if abs(mp.im(a)) < mp.mpf('1e-40'):
        continue
    lam = -mp.log(mp.mpc(a))
    sq = lam**2
    print("   a=%-30s lambda^2 = %-32s Im = %s"
          % (mp.nstr(a, 10), mp.nstr(sq, 14), mp.nstr(mp.im(sq), 6)))
