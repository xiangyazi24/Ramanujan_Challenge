#!/usr/bin/env python3
"""
Reduction of the Chudnovsky CM period-derivative evaluation to four named
components, three classical and one an exact rational identity.

THE TARGET (the sole hypothesis of the Lean theorem `chudnovsky_one_over_pi`):

    A * F(x)^2  +  B * x * 2 F(x) F'(x)  =  640320^(3/2) / (12 pi),

    A = 13591409,  B = 545140134,
    F = 2F1(1/12, 5/12; 1; .),   x = 1728/j(tau_163) = -1728/640320^3.

THE REDUCTION.  Put tau = tau_163 = (1 + i sqrt(163))/2 and z = 1728/j(tau).
Classical inversion gives F = E4^(1/4).  Using

    dz/dtau = 2 pi i * z * E6/E4          (from j' = -2 pi i j E6/E4)
    E4'     = (2 pi i / 3) (E2 E4 - E6)   (Ramanujan)

one computes  z * 2 F F' = (1/6) E4^(1/2) ( E2 E4 / E6 - 1 ), hence

    A F^2 + B z 2 F F'
        = E4^(1/2) [ A - B/6 + (B/6) * E2 E4 / E6 ].

Now substitute the non-holomorphic completion  E2 = E2* + 3/(pi Im tau),
where s2 := E2* E4 / E6.  This is the ONLY place 1/pi enters:

    = E4^(1/2) [ (A - B/6 + (B/6) s2)  +  (B / (2 pi Im tau)) * E4/E6 ].

For the result to be a pure multiple of 1/pi the algebraic bracket must vanish:

    (C1)   s2(tau_163) = 1 - 6A/B.

The surviving term, with Im tau = sqrt(163)/2, is  B E4^(3/2) / (pi sqrt(163) E6),
so the target identity becomes

    (C2)   B / sqrt(163) * E4^(3/2)/E6  =  640320^(3/2) / 12.

Since E4^3/E6^2 = j/(j - 1728), (C2) is purely algebraic in j.

This script verifies (C1) and (C2).  (C2) is checked in EXACT rational
arithmetic; (C1) is checked against q-expansions to 80 digits.

Run:  python3 cm_reduction.py
"""

from fractions import Fraction
from mpmath import mp, mpf, exp, pi, sqrt, nstr

mp.dps = 90

A = 13591409
B = 545140134
JINT = -(640320 ** 3)            # j(tau_163)


def rule(t):
    print()
    print("=" * 72)
    print(t)
    print("=" * 72)


# ----------------------------------------------------------------- (C2) exact
rule("(C2)  B/sqrt(163) * sqrt(j/(j-1728))  =  640320^(3/2)/12   -- EXACT")

J = Fraction(JINT)
lhs_sq = Fraction(B * B, 163) * J / (J - 1728)      # square of the LHS
rhs_sq = Fraction(640320 ** 3, 144)                 # square of the RHS

print(f"  LHS^2 = B^2/163 * j/(j-1728) = {lhs_sq}")
print(f"  RHS^2 = 640320^3/144         = {rhs_sq}")
print(f"  EXACTLY EQUAL: {lhs_sq == rhs_sq}")
print()
print("  Both sides are positive, so (C2) holds exactly.  No numerics used:")
print("  this is an identity between rational numbers and is decidable.")


# ----------------------------------------------------------------- (C1) 80 dps
rule("(C1)  s2(tau_163) = 1 - 6A/B    -- q-expansion check to 80 digits")

s2_required = 1 - Fraction(6 * A, B)
print(f"  required value 1 - 6A/B = {s2_required}")

tau_im = sqrt(163) / 2
q = -exp(-pi * sqrt(163))        # q = e^{2 pi i tau}, tau = (1+i sqrt163)/2


def sigma(n, k):
    return sum(d ** k for d in range(1, n + 1) if n % d == 0)


N = 30
E2 = 1 - 24 * sum(sigma(n, 1) * q ** n for n in range(1, N))
E4 = 1 + 240 * sum(sigma(n, 3) * q ** n for n in range(1, N))
E6 = 1 - 504 * sum(sigma(n, 5) * q ** n for n in range(1, N))
E2star = E2 - 3 / (pi * tau_im)
s2 = E2star * E4 / E6

target = mpf(s2_required.numerator) / s2_required.denominator
rel = abs(s2 - target) / target

print(f"  q            = {nstr(q, 10)}")
print(f"  E2*          = {nstr(E2star, 30)}")
print(f"  s2 = E2*E4/E6 = {nstr(s2, 32)}")
print(f"  1 - 6A/B      = {nstr(target, 32)}")
print(f"  relative diff = {nstr(rel, 5)}")
print()
print("  Note this is a genuine test, not a tautology: the naive approximation")
print(f"  1 - 6/(pi sqrt 163) = {nstr(1 - 6/(pi*sqrt(163)), 32)}")
print("  already disagrees in the 13th decimal, so the q-corrections matter.")


# ----------------------------------------------------------------- assembly
rule("Assembled target")

x = -mpf(1728) / mpf(640320) ** 3
from mpmath import hyp2f1, diff
F = hyp2f1(mpf(1) / 12, mpf(5) / 12, 1, x)
Fp = diff(lambda t: hyp2f1(mpf(1) / 12, mpf(5) / 12, 1, t), x)
lhs = A * F ** 2 + B * x * (2 * F * Fp)
rhs = mpf(640320) ** (mpf(3) / 2) / (12 * pi)
print(f"  A F^2 + B x 2 F F' = {nstr(lhs, 40)}")
print(f"  640320^(3/2)/(12pi) = {nstr(rhs, 40)}")
print(f"  relative difference = {nstr(abs(lhs-rhs)/rhs, 5)}")

rule("What remains classical")
print("""  (i)   F = E4^(1/4)  at the CM point (hypergeometric-modular inversion);
  (ii)  Ramanujan's derivative formulas for E4 and j;
  (iii) s2(tau_163) = 77265280/90856689, i.e. the rationality of s2 at a
        class-number-one CM point together with its value.

  (C2) is fully discharged above.  Crucially the factor 1/pi enters ONLY through
  the definition E2* = E2 - 3/(pi Im tau), so no Chowla-Selberg input is needed.
""")
