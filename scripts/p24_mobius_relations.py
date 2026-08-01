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
