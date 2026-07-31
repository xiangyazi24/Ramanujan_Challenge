#!/usr/bin/env python3
"""
Ramanujan Challenge, Problem 2.1 -- reproducible verification.

The point of this script is FAITHFULNESS: it reimplements, in Python, exactly
the convergent recursion that the Lean development uses (`cfP`, `cfQ`), feeds
it the challenge's own coefficients, and checks that the resulting convergents
really do tend to 6/(3-pi).  It then does the same for Cohen's Entry 5.3.22 and
checks the sign-flip identity convergent by convergent.

Checks:
  (1) the index-shift identities  a_n = -alpha(n+1),  b_n = beta(n)   (exact)
  (2) Cohen Entry 5.3.22 converges to pi                              (numeric)
  (3) the challenge PCF converges to 6/(3-pi)                         (numeric)
  (4) the sign-flip identity  Ptilde_k/Qtilde_k = -(P_k/Q_k)          (exact)
  (5) the displayed values of Cohen's entry: 42, 396, 1047, 38400     (exact)

Usage:  python3 verify.py     (mpmath optional but recommended)
"""

from fractions import Fraction as F
import sys

# ---------------------------------------------------------------- coefficients
# transcribed directly from the problem statement
def a21(n):
    return -220 * n**3 - 484 * n**2 - 301 * n - 42

def b21(n):
    return 4 * n**2 * (2*n + 1)**2 * (5*n - 4) * (5*n + 6)

# transcribed directly from Cohen, arXiv:2607.06581, Entry 5.3.22:
#   [()->Pi,[3, 220*n^3-176*n^2-7*n+5], [6, 4*n^2*(2*n+1)^2*(5*n-4)*(5*n+6)]]
def alphaC(n):
    return 220 * n**3 - 176 * n**2 - 7 * n + 5

def betaC(n):
    return 4 * n**2 * (2*n + 1)**2 * (5*n - 4) * (5*n + 6)

# Cohen's tail T = alpha(1) + beta(1)/(alpha(2) + beta(2)/(alpha(3) + ...))
def cohenC(n):
    return alphaC(n + 1)

def cohenD(n):
    return betaC(n)

# ------------------------------------------- the SAME recursion as the Lean file
# cfP c d 0 = 1, cfP c d 1 = c 0, cfP c d (n+2) = c(n+1)*cfP(n+1) + d(n+1)*cfP(n)
# cfQ c d 0 = 0, cfQ c d 1 = 1,   cfQ c d (n+2) = c(n+1)*cfQ(n+1) + d(n+1)*cfQ(n)

def cf_convergents(c, d, K, exact=True):
    one = F(1) if exact else 1.0
    P = [one, F(c(0)) if exact else float(c(0))]
    Q = [F(0) if exact else 0.0, one]
    for n in range(0, K):
        cn = F(c(n + 1)) if exact else float(c(n + 1))
        dn = F(d(n + 1)) if exact else float(d(n + 1))
        P.append(cn * P[n + 1] + dn * P[n])
        Q.append(cn * Q[n + 1] + dn * Q[n])
    return P, Q

# ------------------------------------------------------------------- check (1)

def check_shift(N=200):
    ok_a = all(a21(n) == -alphaC(n + 1) for n in range(N))
    ok_b = all(b21(n) == betaC(n) for n in range(N))
    print("(1) index-shift identities, n = 0 .. %d" % (N - 1))
    print("    a_n = -alpha(n+1) : %s" % ok_a)
    print("    b_n = beta(n)     : %s" % ok_b)
    return ok_a and ok_b

# ------------------------------------------------------------------- check (5)

def check_displayed():
    vals = {"alpha(1)": (alphaC(1), 42), "beta(1)": (betaC(1), 396),
            "alpha(2)": (alphaC(2), 1047), "beta(2)": (betaC(2), 38400),
            "alpha(3)": (alphaC(3), 4340), "a_0": (a21(0), -42),
            "a_1": (a21(1), -1047), "b_1": (b21(1), 396)}
    ok = True
    print("(5) values against Cohen's displayed continued fraction")
    for k, (got, want) in vals.items():
        good = (got == want)
        ok = ok and good
        print("    %-9s = %-8d expected %-8d %s" % (k, got, want, "ok" if good else "MISMATCH"))
    return ok

# ------------------------------------------------------------------- check (4)

def check_signflip(K=40):
    """The challenge's convergents are exactly minus Cohen's, term by term."""
    Pc, Qc = cf_convergents(cohenC, cohenD, K)
    Pa, Qa = cf_convergents(a21, b21, K)
    bad = []
    for k in range(K + 2):
        # Ptilde_k = (-1)^k P_k  and  Qtilde_k = -(-1)^k Q_k
        if Pa[k] != (-1)**k * Pc[k]:
            bad.append(("P", k))
        if Qa[k] != -((-1)**k) * Qc[k]:
            bad.append(("Q", k))
    print("(4) sign-flip identity on convergents, k = 0 .. %d" % (K + 1))
    print("    Ptilde_k = (-1)^k P_k and Qtilde_k = -(-1)^k Q_k : %s" % (not bad))
    if bad:
        print("      first mismatches:", bad[:5])
    return not bad

# --------------------------------------------------------------- checks (2),(3)

def check_limits(K=60):
    print("(2)/(3) convergence")
    try:
        import mpmath as mp
        mp.mp.dps = 160
        Pc, Qc = cf_convergents(cohenC, cohenD, K)
        Pa, Qa = cf_convergents(a21, b21, K)
        pi = mp.pi
        target_cohen = 3 + 6 / (mp.mpf(Pc[K + 1].numerator) / mp.mpf(Pc[K + 1].denominator)
                                / (mp.mpf(Qc[K + 1].numerator) / mp.mpf(Qc[K + 1].denominator))
                                ) if False else None
        print("    Cohen 5.3.22:  3 + 6/T_k  ->  pi")
        for k in (10, 20, 40, K + 1):
            T = mp.mpf(int(Pc[k])) / mp.mpf(int(Qc[k]))
            v = 3 + 6 / T
            print("      k=%3d  3+6/T = %s   err vs pi = %s"
                  % (k, mp.nstr(v, 30), mp.nstr(v - pi, 6)))
        tgt = 6 / (3 - pi)
        print("    challenge PCF  ->  6/(3-pi) = %s" % mp.nstr(tgt, 30))
        for k in (10, 20, 40, K + 1):
            v = mp.mpf(int(Pa[k])) / mp.mpf(int(Qa[k]))
            print("      k=%3d  P/Q   = %s   err = %s"
                  % (k, mp.nstr(v, 30), mp.nstr(v - tgt, 6)))
        errC = abs(3 + 6 / (mp.mpf(int(Pc[K + 1])) / mp.mpf(int(Qc[K + 1]))) - pi)
        errA = abs(mp.mpf(int(Pa[K + 1])) / mp.mpf(int(Qa[K + 1])) - tgt)
        phi = (1 + mp.sqrt(5)) / 2
        print("    predicted rate phi^-10 per term = %s digits/term"
              % mp.nstr(10 * mp.log10(phi), 6))
        return errC < mp.mpf(10) ** (-100) and errA < mp.mpf(10) ** (-100)
    except ImportError:
        print("    mpmath not installed; skipping the numerical limit checks")
        return None

# ------------------------------------------------------------------------ main

def main():
    print("Ramanujan Challenge, Problem 2.1 -- verification")
    print("=" * 68)
    r1 = check_shift(); print()
    r5 = check_displayed(); print()
    r4 = check_signflip(); print()
    r23 = check_limits(); print()
    print("=" * 68)
    results = [("index shift", r1), ("displayed values", r5),
               ("sign flip", r4), ("limits", r23)]
    for name, r in results:
        print("  %-18s %s" % (name, "PASS" if r else ("SKIP" if r is None else "FAIL")))
    if any(r is False for _, r in results):
        sys.exit(1)
    print("\nAll executed checks passed.")


if __name__ == "__main__":
    main()
