"""The two Mobius relations (R12), (R21), checked against the six known rows.

With I11, I20, I22 already proved in Lean and I10 the only row left, these two
relations reduce I12 and I21 to I10:

    I12 = I22 - I10 + I20 - (13/20) Z2^2
    I21 = I22 + 2 I10 - 2 I20 - (1/2) I11 + (37/40) Z2^2

Verified here as exact identities over Q in the basis
(Li4(1/2), L^4, L^2 Z2, L Z3, Z2^2), using the independently established row
values.  A relation that failed this check would be wrong, and it is much cheaper
to find that out here than in Lean.
"""
from fractions import Fraction as F

NAMES = ["Li4(1/2)", "L^4", "L^2*Z2", "L*Z3", "Z2^2"]
ROW = {
    "I10": (F(-2), F(-1, 12), F(-1),    F(-7, 4),  F(1, 10)),
    "I11": (F(0),  F(0),      F(0),     F(-7, 2),  F(3, 4)),
    "I12": (F(-6), F(-1, 4),  F(3),     F(-21, 4), F(9, 5)),
    "I20": (F(-2), F(-1, 12), F(1, 2),  F(-7, 4),  F(1, 4)),
    "I21": (F(-6), F(-1, 4),  F(-3, 2), F(-7, 2),  F(51, 20)),
    "I22": (F(-6), F(-1, 4),  F(3, 2),  F(-21, 4), F(23, 10)),
}
Z22 = (F(0), F(0), F(0), F(0), F(1))


def lin(*terms):
    out = [F(0)] * 5
    for c, v in terms:
        for i in range(5):
            out[i] += F(c) * v[i]
    return tuple(out)


checks = [
    ("R12", lin((1, ROW["I12"]), (-1, ROW["I22"]), (1, ROW["I10"]), (-1, ROW["I20"])),
     lin((F(-13, 20), Z22))),
    ("R21", lin((1, ROW["I21"]), (-1, ROW["I22"]), (-2, ROW["I10"]), (2, ROW["I20"]),
                (F(1, 2), ROW["I11"])),
     lin((F(37, 40), Z22))),
    ("I12 from R12", lin((1, ROW["I22"]), (-1, ROW["I10"]), (1, ROW["I20"]),
                         (F(-13, 20), Z22)), ROW["I12"]),
    ("I21 from R21", lin((1, ROW["I22"]), (2, ROW["I10"]), (-2, ROW["I20"]),
                         (F(-1, 2), ROW["I11"]), (F(37, 40), Z22)), ROW["I21"]),
]

ok = True
for name, lhs, rhs in checks:
    good = lhs == rhs
    ok &= good
    mark = "OK" if good else "MISMATCH"
    print(f"{name:<14} {mark}")
    if not good:
        for n, a, b in zip(NAMES, lhs, rhs):
            if a != b:
                print(f"    {n}: {a} vs {b}")
print("\nALL EXACT" if ok else "\nFAILED")


# ---------------------------------------------------------------------------
# The two relations as SINGLE integrals, and why they are cheap.
#
#   R12 = int_0^1 W0 (H1 - H2) (1/(2-t) + 1/t) dt                  = -(13/20) Z2^2
#   R21 = int_0^1 [ W0 H2 (1/(1-t) - 1/(2-t))
#                   - 2 W0 (H1 - H2)/t
#                   + (1/2) W0 H1/(1-t) ] dt                       =  (37/40) Z2^2
#
# Both answers are PURE Z2^2: no Li4(1/2), no log^4 2, no log2*zeta3.  The
# weight-four content cancels inside the Mobius structure, which is why these two
# integrals are much lighter than any single row -- their values need no
# polylogarithm machinery at all.
# ---------------------------------------------------------------------------

def _single_integral_check():
    from mpmath import mp, mpf, log, pi, polylog, quad, exp, nstr
    mp.dps = 30
    Z2 = pi ** 2 / 6
    W0 = lambda t: Z2 - 2 * polylog(2, t / 2) - log(t / 2) ** 2
    H1 = lambda t: -log(1 - t)
    H2 = lambda t: -log(1 - t / 2)

    def I(f):  # kill both endpoint singularities
        a = quad(lambda v: f(exp(-v)) * exp(-v), [log(2), 3, 12, 40])
        b = quad(lambda v: f(1 - exp(-v)) * exp(-v), [log(2), 3, 12, 40])
        return a + b

    r12 = I(lambda t: W0(t) * (H1(t) - H2(t)) * (1 / (2 - t) + 1 / t))
    r21 = I(lambda t: (W0(t) * H2(t) * (1 / (1 - t) - 1 / (2 - t))
                       - 2 * W0(t) * (H1(t) - H2(t)) / t
                       + mpf(1) / 2 * W0(t) * H1(t) / (1 - t)))
    print(f"R12 single integral  {nstr(r12, 20):>24}   target {nstr(-mpf(13) / 20 * Z2 ** 2, 20)}")
    print(f"R21 single integral  {nstr(r21, 20):>24}   target {nstr(mpf(37) / 40 * Z2 ** 2, 20)}")


if __name__ == "__main__":
    _single_integral_check()
