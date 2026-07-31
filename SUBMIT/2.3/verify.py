#!/usr/bin/env python3
"""
Ramanujan Challenge, Problem 2.3 -- reproducible verification.

Three independent checks, all in exact arithmetic except the last:

  (1) Forward-solve the challenge recurrence over Q from the challenge's own
      initial values, and confirm it agrees with the claimed closed forms
          q_n = A_{n+2} D_{n+3},   p_n = 4 B_{n+2} D_{n+3} + A_{n+2} (n+3)!
      This is the faithfulness check: nothing but the problem statement is
      used on the left-hand side.

  (2) The tensor theorem with FOUR FREE INITIAL VALUES: for symbolic
      X_{-1}, X_0, Y_0, Y_1, the product X_{n+2} Y_{n+3} annihilates the
      challenge operator identically.  This is the theorem, not an instance.
      (Requires sympy; skipped with a clear message if unavailable.)

  (3) High-precision evaluation of p_n/q_n against pi + e.
      (Requires mpmath; falls back to Fraction + float if unavailable.)

Usage:  python3 verify.py
"""

from fractions import Fraction
import math
import sys

# ---------------------------------------------------------------- coefficients

def coeffs(n):
    """The five coefficients c_0(n), ..., c_4(n) of the challenge recurrence,
    transcribed directly from the problem statement."""
    c0 = -n**3 + 2*n**2 + 7*n + 3
    c1 = (n + 2) * (2*n**4 + n**3 - 26*n**2 - 48*n - 19)
    c2 = (n + 2) * (n**6 + 9*n**5 + 8*n**4 - 87*n**3 - 249*n**2 - 234*n - 68)
    c3 = (n + 1)**2 * (n + 2) * (2*n**5 + 3*n**4 - 13*n**3 - 21*n**2 + 4)
    c4 = -n**3 * (n + 1)**2 * (n + 2) * (n**3 + n**2 - 8*n - 11)
    return c0, c1, c2, c3, c4


P_INIT = [1, 1, 20, 296]   # p_{-3}, p_{-2}, p_{-1}, p_0
Q_INIT = [1, 0, 4, 48]     # q_{-3}, q_{-2}, q_{-1}, q_0

# ------------------------------------------------------- the two order-2 systems

def lambert(N):
    """A_m, B_m for m = -1 .. N, keyed by m."""
    A = {-1: 1, 0: 1}
    B = {-1: 0, 0: 1}
    for m in range(1, N + 1):
        A[m] = (2*m + 1) * A[m-1] + m*m * A[m-2]
        B[m] = (2*m + 1) * B[m-1] + m*m * B[m-2]
    return A, B


def derangements(N):
    D = {0: 1, 1: 0}
    for m in range(2, N + 1):
        D[m] = (m - 1) * (D[m-1] + D[m-2])
    return D

# ------------------------------------------------------------------- check (1)

def forward_solve(init, N):
    """Solve the challenge recurrence over Q from initial data alone."""
    u = {-3: Fraction(init[0]), -2: Fraction(init[1]),
         -1: Fraction(init[2]),  0: Fraction(init[3])}
    for n in range(1, N + 1):
        c0, c1, c2, c3, c4 = coeffs(n)
        if c0 == 0:
            raise AssertionError("leading coefficient vanished at n=%d" % n)
        u[n] = -(c1*u[n-1] + c2*u[n-2] + c3*u[n-3] + c4*u[n-4]) / c0
    return u


def check_closed_forms(N=40):
    A, B = lambert(N + 4)
    D = derangements(N + 4)
    q = forward_solve(Q_INIT, N)
    p = forward_solve(P_INIT, N)

    bad = []
    for n in range(-3, N + 1):
        q_cf = A[n+2] * D[n+3]
        p_cf = 4 * B[n+2] * D[n+3] + A[n+2] * math.factorial(n+3)
        if q[n] != q_cf:
            bad.append(("q", n, q[n], q_cf))
        if p[n] != p_cf:
            bad.append(("p", n, p[n], p_cf))

    integral = all(v.denominator == 1 for v in q.values()) and \
               all(v.denominator == 1 for v in p.values())

    print("(1) forward solve vs. closed forms, n = -3 .. %d" % N)
    print("    all values integral            : %s" % integral)
    print("    closed forms match everywhere  : %s" % (not bad))
    if bad:
        for item in bad[:5]:
            print("      MISMATCH", item)
        return False
    print("    first new values: q_1 = %d, q_2 = %d, p_1 = %d, p_2 = %d"
          % (q[1], q[2], p[1], p[2]))
    return True

# ------------------------------------------------------------------- check (2)

def check_tensor_symbolic(NMAX=22):
    try:
        import sympy as sp
    except ImportError:
        print("(2) SKIPPED: sympy not installed "
              "(pip install sympy to run the symbolic check)")
        return None

    x0, x1, y0, y1 = sp.symbols('x0 x1 y0 y1')
    X = {-1: x0, 0: x1}
    for m in range(1, NMAX + 4):
        X[m] = sp.expand((2*m + 1) * X[m-1] + m*m * X[m-2])
    Y = {0: y0, 1: y1}
    for m in range(2, NMAX + 4):
        Y[m] = sp.expand((m - 1) * (Y[m-1] + Y[m-2]))

    def u(k):
        return sp.expand(X[k+2] * Y[k+3])

    ok = True
    for n in range(1, NMAX + 1):
        c0, c1, c2, c3, c4 = coeffs(n)
        r = sp.expand(c0*u(n) + c1*u(n-1) + c2*u(n-2) + c3*u(n-3) + c4*u(n-4))
        if r != 0:
            ok = False
            print("      NONZERO RESIDUAL at n=%d: %s" % (n, r))
            break

    print("(2) tensor theorem with four FREE initial values, n = 1 .. %d" % NMAX)
    print("    residual identically zero      : %s" % ok)
    return ok

# ------------------------------------------------------------------- check (3)

def check_limit(N=50):
    A, B = lambert(N + 4)
    D = derangements(N + 4)

    def q(n):
        return A[n+2] * D[n+3]

    def p(n):
        return 4 * B[n+2] * D[n+3] + A[n+2] * math.factorial(n+3)

    print("(3) convergence of p_n / q_n to pi + e")
    try:
        import mpmath as mp
        mp.mp.dps = 60
        target = mp.pi + mp.e
        print("    pi + e = %s" % mp.nstr(target, 45))
        for n in (10, 20, 30, 40, N):
            r = mp.mpf(p(n)) / mp.mpf(q(n))
            print("    n = %2d   p/q = %s   err = %s"
                  % (n, mp.nstr(r, 25), mp.nstr(r - target, 6)))
        err = abs(mp.mpf(p(N)) / mp.mpf(q(N)) - target)
        return err < mp.mpf(10) ** (-35)
    except ImportError:
        print("    mpmath not installed; falling back to double precision")
        target = math.pi + math.e
        for n in (5, 10, 15):
            r = p(n) / q(n)
            print("    n = %2d   p/q = %.15f   err = %.3e" % (n, r, r - target))
        return abs(p(15) / q(15) - target) < 1e-12

# ------------------------------------------------------------------------ main

def main():
    print(__doc__.strip().splitlines()[0])
    print("=" * 68)
    r1 = check_closed_forms()
    print()
    r2 = check_tensor_symbolic()
    print()
    r3 = check_limit()
    print()
    print("=" * 68)
    results = [("closed forms", r1), ("tensor theorem", r2), ("limit", r3)]
    for name, r in results:
        status = "PASS" if r else ("SKIP" if r is None else "FAIL")
        print("  %-16s %s" % (name, status))
    failed = [n for n, r in results if r is False]
    if failed:
        print("\nFAILED: %s" % ", ".join(failed))
        sys.exit(1)
    print("\nAll executed checks passed.")


if __name__ == "__main__":
    main()
