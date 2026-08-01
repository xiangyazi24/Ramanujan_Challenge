"""Weight-4 half-argument Euler sums: do they close in the target basis?

This is the go/no-go test for the SERIES route to P2.4's six endpoint rows.
Route (A), through the antiderivative
    E(a) = -log^3 a log(1-a) - 3 log^2 a Li2(a) + 6 log a Li3(a) - 6 Li4(a),
needs the derivative of Li4, which neither the repo nor Mathlib v4.29 has.
Route (B) integrates termwise on [0,1/2]; the 2^{-n} factors turn the answers
into half-argument Euler sums.  Route (B) is only viable if every such sum lies
in span_Q{Li4(1/2), log^4 2, log^2 2 * zeta2, log2 * zeta3, zeta2^2}.

Verdict: all of them do.  Run at 60 dps; every relation found with coefficients
of denominator <= 40.  The first row agrees coefficient-by-coefficient with the
independent analytic derivation in /tmp/gpt/ccdex/Q6407.md eq (3.9).
"""
from fractions import Fraction as F
from mpmath import mp, mpf, pi, log, zeta, polylog, pslq, nstr, nsum, inf

mp.dps = 60
L, Z2, Z3 = log(2), pi**2 / 6, zeta(3)
A4 = polylog(4, mpf(1) / 2)
BASIS = [A4, L**4, L**2 * Z2, L * Z3, Z2**2]
NAMES = ["Li4(1/2)", "log^4 2", "log^2 2*Z2", "log2*Z3", "Z2^2"]


def H(n, p=1):
    return sum(mpf(1) / mpf(j) ** p for j in range(1, n + 1))


SUMS = {
    "sum H_n/(2^n n^3)":       lambda: nsum(lambda n: H(int(n)) / (2 ** int(n) * mpf(n) ** 3), [1, inf]),
    "sum H_n^2/(2^n n^2)":     lambda: nsum(lambda n: H(int(n)) ** 2 / (2 ** int(n) * mpf(n) ** 2), [1, inf]),
    "sum H_n^{(2)}/(2^n n^2)": lambda: nsum(lambda n: H(int(n), 2) / (2 ** int(n) * mpf(n) ** 2), [1, inf]),
    "sum H_n^{(3)}/(2^n n)":   lambda: nsum(lambda n: H(int(n), 3) / (2 ** int(n) * mpf(n)), [1, inf]),
    "sum 1/(2^n n^4)":         lambda: A4,
}

ok = True
for name, f in SUMS.items():
    v = f()
    rel = pslq([v] + BASIS, maxcoeff=10**8, maxsteps=8000)
    if not rel or rel[0] == 0:
        print(f"{name:<26} {nstr(v, 18):>20}  *** NOT IN BASIS -- route (B) is a trap ***")
        ok = False
        continue
    c0 = rel[0]
    coeffs = [F(int(-x), int(c0)).limit_denominator(10**6) for x in rel[1:]]
    body = " + ".join(f"({a})*{b}" for a, b in zip(coeffs, NAMES) if a != 0)
    print(f"{name:<26} {nstr(v, 18):>20}  =  {body}")
print("\nROUTE (B) VIABLE" if ok else "\nROUTE (B) BLOCKED")
