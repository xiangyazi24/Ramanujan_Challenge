#!/usr/bin/env sage
"""Minimal Sage/Singular ideal-lift API probe."""

R.<x, y> = PolynomialRing(QQ, order="degrevlex")
I = R.ideal(x + y, y)
target = x
sI = singular(I)
s_target = singular(target)
lifted = singular.lift(sI, s_target)
print("TYPE", type(lifted))
print("TEXT", lifted)
try:
    print("SAGE", lifted.sage())
except Exception as error:
    print("SAGE_ERROR", repr(error))
