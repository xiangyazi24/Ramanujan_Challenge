"""The non-alternating twin: Layer C analogue, and why it is not a sign flip.

Two things are checked here, both of which decide how the last open certificate
input (`quadraticEulerTerm24`) has to be formalized.

1. The Layer C analogue holds:

       sum quadraticEulerTerm24 = int_0^1 (-log x)/x * Q(x) dx
                                = quadraticEulerValue24

   The already-proved ALTERNATING case (with Q(-x)) is run alongside as a
   control, so a setup error would show up in both rather than being mistaken
   for a result.

2. The integrand is NOT bounded at x = 1.  Q(x) = 2 J(x)/(1-x) has a pole there.
   The alternating chain evaluates Q(-x) for x in (0,1) -- argument in (-1,0) --
   and never meets it.  The non-alternating chain runs straight into it.  Only
   half the pole cancels: -log x ~ (1-x) kills the 1/(1-x), but J itself grows
   like log^2(1-x), so the integrand behaves like 2 log^2(1-x).

   Consequence: an integrable log^2 singularity at x = 1 that the alternating
   integrand does not have, so the twin needs endpoint integrability work the
   alternating chain never required.  "Mirror the alternating chain" is the
   wrong brief for this task.

Definitions are transcribed from the Lean source, not from memory:
  quadAltMclosed x = log(1-x)^2/2 + log(1+x)^2 + 2 log(1-x) log(1+x)
                     + 2 log2 log(1+x) + pi^2/6 - log^2 2 - 2 Li2((1+x)/2)
  quadAltJclosed x = quadAltMclosed x + Li2(x^2)
  quadAltQclosed x = 2 * quadAltJclosed x / (1 - x)
"""
from mpmath import mp, mpf, log, pi, polylog, quad, exp, nstr, zeta

mp.dps = 30

Z2 = pi ** 2 / 6
L = log(2)
Z3 = zeta(3)
A4 = polylog(4, mpf(1) / 2)


def Mclosed(x):
    return (log(1 - x) ** 2 / 2 + log(1 + x) ** 2
            + 2 * log(1 - x) * log(1 + x) + 2 * L * log(1 + x)
            + Z2 - L ** 2 - 2 * polylog(2, (1 + x) / 2))


def Jclosed(x):
    return Mclosed(x) + polylog(2, x ** 2)


def Qclosed(x):
    return 2 * Jclosed(x) / (1 - x)


def coeff_integral(sign):
    """int_0^1 (-log x)/x * Q(sign*x) dx, with the x=0 end substituted away."""
    f = lambda x: -log(x) / x * Qclosed(sign * x)
    return quad(lambda v: f(exp(-v)) * exp(-v), [0, 1, 8, 40])


NONALT = (20 * A4 + mpf(5) / 6 * L ** 4 + 7 * L ** 2 * Z2
          - mpf(59) / 10 * Z2 ** 2)                       # quadraticEulerValue24
ALT = (-22 * A4 - mpf(11) / 12 * L ** 4 - mpf(13) / 2 * L ** 2 * Z2
       - mpf(7) / 4 * L * Z3 + mpf(67) / 10 * Z2 ** 2)    # alternating, PROVED

print("1. Layer C analogue")
for name, sign, target in (("non-alternating  Q(+x)", 1, NONALT),
                           ("alternating      Q(-x)  [control, proved]", -1, ALT)):
    got = coeff_integral(sign)
    print(f"   {name:<44} {nstr(got, 18):>22}")
    print(f"   {'declared value':<44} {nstr(target, 18):>22}   diff {nstr(abs(got - target), 3)}")

print("\n2. Behaviour of the non-alternating integrand at x = 1")
print("   growth is log^2(1-x); the ratio below is the coefficient, and it")
print("   creeps toward 1 rather than sitting at a round constant")
for e in (mpf(10) ** -6, mpf(10) ** -10, mpf(10) ** -14, mpf(10) ** -20):
    x = 1 - e
    val = -log(x) / x * Qclosed(x)
    ref = log(e) ** 2
    print(f"   x = 1 - 1e{int(log(e)/log(10)):<4} integrand {nstr(val, 12):>16}"
          f"   log^2(1-x) = {nstr(ref, 12):>16}   ratio {nstr(val/ref, 6)}")
print("\n   -> integrable log^2 singularity at x = 1; the alternating integrand"
      "\n      has none, since Q(-x) never approaches the pole.")
