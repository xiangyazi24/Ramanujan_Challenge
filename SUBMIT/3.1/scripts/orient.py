#!/usr/bin/env python3
"""Determine the tetrahedron orientation signs for the 7_2 ideal triangulation,
and test whether the four-shape regulator sum is signed or unsigned.

The project chart uses four shapes (T,U,V,W).  If any tetrahedron enters with
epsilon = -1 the Bloch-Wigner cancellation D(V) + D(W) = 0 becomes
D(V) - D(W) = 2 D(V) != 0 and the torsion argument fails.
"""
import snappy
import mpmath as mp
mp.mp.dps = 50

K = snappy.Manifold('7_2')
print("7_2 triangulation:")
print("  num_tetrahedra :", K.num_tetrahedra())
print("  solution type  :", K.solution_type())
print("  volume         :", K.volume())
print()
print("geometric shapes (complete structure):")
for i, z in enumerate(K.tetrahedra_shapes('rect')):
    zc = complex(z)
    print("   tet %d: z = %-34s  Im = %+.6f  %s"
          % (i, f"{zc.real:.12f}{zc.imag:+.12f}j", zc.imag,
             "positively oriented" if zc.imag > 0 else "NEGATIVELY oriented"))
print()

print("gluing equations (edge + cusp), as exponent vectors:")
eqs = K.gluing_equations('rect')
for row in eqs:
    print("   ", row)
print()

# retriangulate / check the 4-tetrahedron chart
print("filled -1/2:")
K2 = snappy.Manifold('7_2')
K2.dehn_fill((-1, 2))
print("  volume:", K2.volume())
print("  shapes:")
for i, z in enumerate(K2.tetrahedra_shapes('rect')):
    zc = complex(z)
    print("   tet %d: z = %-34s  Im = %+.6f" % (i, f"{zc.real:.12f}{zc.imag:+.12f}j", zc.imag))

print()
print("--- the project chart shapes at the alpha endpoint (all REAL) ---")
a = mp.mpf('0.59098942867025644049')
M = a**2; L = a; X = M*M
u = (L + X**3)/(X*(L+X))
r = -(1+mp.sqrt(1+4*u**2))/(2*u)
T = 1-r**2; U = u; V = u/X; W = 1/(1-u*X)
for n, v in [('T',T),('U',U),('V',V),('W',W)]:
    print("   %s = %-24s   sign of (v), (1-v), (1-1/v):  %+d %+d %+d"
          % (n, mp.nstr(v, 18), mp.sign(v), mp.sign(1-v), mp.sign(1-1/v)))
print()
print("For a REAL shape z the flat tetrahedron has a well-defined sign given by")
print("the cyclic order of (0, 1, z, infinity) on the real line, equivalently by")
print("which of the three cross-ratio branches z, 1/(1-z), 1-1/z lies in (1,inf).")
for n, v in [('T',T),('U',U),('V',V),('W',W)]:
    br = [v, 1/(1-v), 1-1/v]
    which = [i for i, x in enumerate(br) if x > 1]
    print("   %s: branches %s -> in (1,inf): %s"
          % (n, [mp.nstr(x, 8) for x in br], which))
