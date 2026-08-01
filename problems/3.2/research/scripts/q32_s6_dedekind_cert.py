#!/usr/bin/env python3
"""Dedekind cycle-type S6 certificates for the codegree polynomial f_x
(fourth verification leg of LEMMA-CODEGREE / Q6455 step 3.3; method
transplanted from cron's CRON_Mh_galois.py three-criteria approach).

Logic: for a good specialization (x, y0), Gal(f_x - y0 / Q) is a subgroup of
the generic-point monodromy group, which is contained in S_6. Certifying
S_6 at ONE specialization therefore proves the generic group equals S_6.
Criteria: irreducible mod q (transitive); (1,5) pattern (prime 5-cycle,
5 > 6/2, Jordan => primitive); single-even-part-2 pattern (power gives a
transposition); primitive + transposition => S_6.
"""
from sympy import Poly, factor_list, primerange, symbols

W = symbols("W")


def cert(xval, y0):
    X = xval + 5
    f = 34 * W**3 - 51 * W**4 + 27 * W**5 - X * W**6 - y0
    got = {"trans": False, "prim5": False, "transposition": False}
    for q in primerange(5, 400):
        if X % q == 0:
            continue
        fl = factor_list(Poly(f, W, modulus=q))
        degs = sorted(p.degree() for p, m in fl[1] for _ in range(m))
        if sum(degs) != 6:
            continue
        if degs == [6]:
            got["trans"] = True
        if degs == [1, 5]:
            got["prim5"] = True
        if [d for d in degs if d % 2 == 0] == [2]:
            got["transposition"] = True
        if all(got.values()):
            return True, q
    return all(got.values()), None


if __name__ == "__main__":
    import random

    random.seed(7)
    ok_all = True
    for xval in [7, -3, 1, 12, 100, -17]:
        y0 = random.randint(2, 50)
        ok, _ = cert(xval, y0)
        print(f"x={xval}, y0={y0}: S6 certificate {'OK' if ok else 'INCOMPLETE'}")
        ok_all &= ok
    print("ALL S6 CERTIFIED" if ok_all else "INCOMPLETE")
