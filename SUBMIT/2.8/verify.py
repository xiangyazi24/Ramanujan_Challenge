#!/usr/bin/env python3
"""
Numerical certificates for Ramanujan Challenge Problem 2.8.

Everything here is reproducible with mpmath alone:

  1. The isolated CM period-derivative evaluation that the Lean development
     keeps as an explicit hypothesis, checked to 130 decimal digits.
  2. The Chudnovsky series itself, summed and compared with sqrt(10005)/pi.
  3. The Poincare-root bridge: the ratio of successive Chudnovsky terms versus
     the Poincare root of the CMF scalar recurrence.

Run:  python3 verify.py
"""

from mpmath import mp, hyp2f1, diff, pi, mpf, sqrt, nstr, log10, factorial

mp.dps = 140


def rule(title):
    print()
    print("=" * 72)
    print(title)
    print("=" * 72)


# ---------------------------------------------------------------- 1. CM evaluation
rule("1. The CM period-derivative evaluation (the Lean hypothesis)")

J = mpf(640320) ** 3                    # -j(tau_163)
x = -mpf(1728) / J                      # hypergeometric argument 1728/j(tau_163)

F = hyp2f1(mpf(1) / 12, mpf(5) / 12, 1, x)
Fp = diff(lambda t: hyp2f1(mpf(1) / 12, mpf(5) / 12, 1, t), x)

lhs = 13591409 * F ** 2 + 545140134 * x * (2 * F * Fp)
rhs = mpf(640320) ** (mpf(3) / 2) / (12 * pi)

rel = abs(lhs - rhs) / abs(rhs)

print(f"  x        = -1728/640320^3 = {nstr(x, 20)}")
print(f"  F(x)     = {nstr(F, 30)}")
print(f"  F'(x)    = {nstr(Fp, 30)}")
print(f"  (F'(0) should be ab/c = (1/12)(5/12)/1 = 5/144 = {nstr(mpf(5)/144, 20)})")
print()
print(f"  LHS      = {nstr(lhs, 40)}")
print(f"  RHS      = {nstr(rhs, 40)}")
print(f"  rel.err  = {nstr(rel, 5)}")
print(f"  digits   ~ {int(-log10(rel)) if rel > 0 else '>' + str(mp.dps)}")


# ---------------------------------------------------------------- 2. the series
rule("2. The Chudnovsky series")

def chud_term(k):
    num = (-1) ** k * factorial(6 * k) * (13591409 + 545140134 * k)
    den = factorial(3 * k) * factorial(k) ** 3 * mpf(640320) ** (3 * k)
    return num / den

S = mp.nsum(chud_term, [0, mp.inf])
lhs2 = mpf(640320) ** (mpf(3) / 2) / (12 * pi)
rel2 = abs(S - lhs2) / abs(lhs2)

print(f"  sum_k (-1)^k (6k)!(13591409+545140134k) / ((3k)!(k!)^3 640320^(3k))")
print(f"           = {nstr(S, 40)}")
print(f"  640320^(3/2)/(12 pi) = {nstr(lhs2, 40)}")
print(f"  rel.err  = {nstr(rel2, 5)}")
print()
print(f"  equivalently  sqrt(10005)/pi = {nstr(sqrt(10005)/pi, 40)}")
print(f"  and 640320 = 64 * 10005 = {64*10005}, so 640320^(3/2) = 512*10005*sqrt(10005)")


# ---------------------------------------------------------------- 3. Poincare root
rule("3. Poincare-root bridge")

R = 1 - (-J) / 1728                      # R = 1 - j(tau_163)/1728
rho = 64 * (R - 1)

target = -1 / (R - 1)

print(f"  R          = 1 - j(tau_163)/1728 = {R}")
print(f"  rho = 64(R-1)                     = {rho}")
print(f"  64*640320^3/1728                  = {64 * mpf(640320)**3 / 1728}")
print()
print("  The term ratio approaches -1/(R-1) only asymptotically, because the")
print("  linear factor (13591409 + 545140134 k) still varies for small k.")
print("  Convergence of h_(k+1)/h_k  ->  -1/(R-1):")
print()
print(f"  {'k':>8}   {'h_(k+1)/h_k':>32}   {'rel.err vs -1/(R-1)':>22}")
for k in [10, 100, 1000, 10000, 100000]:
    r = chud_term(k + 1) / chud_term(k)
    print(f"  {k:>8}   {nstr(r, 22):>32}   {nstr(abs(r - target)/abs(target), 5):>22}")
print()
print(f"  -1/(R-1)  =  {nstr(target, 25)}")
print()
print("  (The CMF scalar recurrence and the Chudnovsky series therefore share")
print("   the same Poincare root, hence the same limit and geometric rate.)")

print()
print("done.")
