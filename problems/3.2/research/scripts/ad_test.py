#!/usr/bin/env python3
"""a_d = lim (r+d)^3 A_d(r) = -5*A_{d-1}(-d) + A_{d-2}(-d); closed form hunt.
Compare with d^3 * (Apery companion a-seq /6) and d^3*b_{d-1} patterns."""
from fractions import Fraction as F

def coefs(x):
    den = F((x+1)**3)
    return F((2*x+1)*(17*x*x+17*x+5),1)/den, F(-(x**3),1)/den

def A_at(k, base):
    if k == 0: return F(1)
    Aprev = F(1)
    a1, b1 = coefs(base); Acur = a1
    for j in range(1, k):
        aj, bj = coefs(base + j)
        Acur, Aprev = aj*Acur + bj*Aprev, Acur
    return Acur

# Apery b and companion atilde (atilde_0=0, atilde_1=6, same recurrence, rational)
b = [F(1), F(5)]; at = [F(0), F(6)]
for n in range(1, 25):
    P = F((2*n+1)*(17*n*n+17*n+5))
    b.append((P*b[n] - n**3*b[n-1]) / F((n+1)**3))
    at.append((P*at[n] - n**3*at[n-1]) / F((n+1)**3))

print(f"{'d':>3} {'a_d':>34} {'a_d/(d^3)':>26} {'vs b_{d-1}, atilde_{d-1}':>30}")
for d in range(2, 16):
    ad = F(-5)*A_at(d-1, -d) + A_at(d-2, -d)
    ratio = ad / F(d**3)
    print(f"{d:>3} {str(ad)[:34]:>34} {str(ratio)[:26]:>26}   b={b[d-1]}  at={at[d-1]}")
