#!/usr/bin/env python3
"""The flattening ambiguity contributes a HALF-INTEGER multiple of pi^2, and
therefore cannot spoil the denominator bound.

Neumann's extended Rogers dilogarithm of a flattened simplex is

    Rhat([z;p,q]) = Li_2(z) + (1/2) log(z) log(1-z)
                    + (pi i/2)(q log z + p log(1-z)) - pi^2/6.

The (p,q)-dependent part is  (pi i/2)(q log z + p log(1-z)).  Writing the
principal branches as

    log z     = A + i theta_1,      log(1-z) = B + i theta_2,

with A, B real and theta_k in (-pi, pi], the (p,q)-part equals

    (pi i/2)(qA + pB) - (pi/2)(q theta_1 + p theta_2),

so its REAL part is  -(pi/2)(q theta_1 + p theta_2).

At BOTH endpoints of our arc all four shapes are REAL, and then theta_1, theta_2
are each either 0 or pi:

    z < 0      =>  theta_1 = pi,  theta_2 = 0     (1 - z > 1)
    z > 1      =>  theta_1 = 0,   theta_2 = pi    (1 - z < 0)

Hence each tetrahedron contributes an INTEGER multiple of pi^2/2 to Re Rhat, and
so does the sum over the four tetrahedra, and so does the difference between the
two endpoints.

CONSEQUENCE.  Changing the flattening changes Re[Delta R]/pi^2 by a half-integer.
So if the torsion argument bounds the denominator by Q, the flattening ambiguity
at worst doubles it: the denominator divides 2Q = 4080.  Since the reconstruction
tolerates any bound up to 1.75e150, this is irrelevant to the conclusion -- but it
means the flattening choice does not have to be exhibited at all.

This script verifies the claim numerically at both endpoints, for a range of
flattening integers.

Run:  python3 flattening_ambiguity.py
"""

import mpmath as mp

mp.mp.dps = 60

ALPHA = ["-0.157875280712128451", "6.8157987008302631",
         "55.8724744252461376", "5.93292171293941164"]
BETA = ["-0.258287827210", "4.342962293300",
        "26.241958221000", "3.555514130300"]


def Rhat(z, p, q):
    """Neumann's extended Rogers dilogarithm of the flattened simplex [z;p,q]."""
    z = mp.mpc(z)
    lz = mp.log(z)
    l1z = mp.log(1 - z)
    return (mp.polylog(2, z) + mp.mpf('0.5') * lz * l1z
            + (mp.pi * 1j / 2) * (q * lz + p * l1z) - mp.pi ** 2 / 6)


def branch_data(z):
    """(theta_1, theta_2) for the principal logs of z and 1-z."""
    z = mp.mpc(z)
    return mp.im(mp.log(z)), mp.im(mp.log(1 - z))


def check(name, shapes):
    print("=" * 70)
    print(name)
    print("=" * 70)
    print("  shape                   theta_1/pi   theta_2/pi")
    for z in shapes:
        t1, t2 = branch_data(z)
        print("  %-22s  %+.6f    %+.6f"
              % (z[:20], float(t1 / mp.pi), float(t2 / mp.pi)))

    base = sum(Rhat(z, 0, 0) for z in shapes)
    print()
    print("  Re sum Rhat with all (p,q) = (0,0):  %.20f" % float(mp.re(base)))
    print()
    print("  effect of flattening changes on Re(sum), in units of pi^2:")
    worst = mp.mpf(0)
    for j in range(len(shapes)):
        for (dp, dq) in [(1, 0), (0, 1), (-1, 0), (0, -1), (2, 3), (-1, 2)]:
            pq = [[0, 0] for _ in shapes]
            pq[j] = [dp, dq]
            tot = sum(Rhat(z, pq[i][0], pq[i][1]) for i, z in enumerate(shapes))
            delta = (mp.re(tot) - mp.re(base)) / mp.pi ** 2
            twice = delta * 2
            dev = abs(twice - mp.nint(twice))
            worst = max(worst, dev)
            print("     tet %d, (dp,dq) = (%+d,%+d):  Delta/pi^2 = %+8.4f"
                  "   2*Delta/pi^2 = %+6.1f   integer? %s"
                  % (j, dp, dq, float(delta), float(twice), dev < mp.mpf('1e-30')))
    print()
    print("  worst deviation of 2*Delta/pi^2 from an integer: %.3e" % float(worst))
    print("  => every flattening change moves Re(sum)/pi^2 by a HALF-INTEGER:",
          worst < mp.mpf('1e-30'))
    print()
    return worst < mp.mpf('1e-30')


print(__doc__)
ok_a = check("ALPHA endpoint", ALPHA)
ok_b = check("BETA endpoint", BETA)

print("=" * 70)
print("CONCLUSION")
print("=" * 70)
print("""
At both endpoints every flattening change shifts Re[sum Rhat]/pi^2 by a
half-integer.  Hence it shifts Re[Delta R]/pi^2 by a half-integer too.

So the flattening choice never has to be exhibited: whatever it is, the
denominator of Re[Delta R]/pi^2 divides 2Q = 4080 rather than Q = 2040, and the
301-digit reconstruction tolerates any bound below 1.75e150.

This replaces the earlier claim that the Dehn-filling core correction is
'absorbed into the flattening data' -- a claim that would have needed the filled
flattening to be constructed explicitly.
""")
print("alpha:", ok_a, "  beta:", ok_b)
