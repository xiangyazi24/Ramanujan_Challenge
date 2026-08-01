#!/usr/bin/env python3
"""Spot-check polar Lucas law: p^3 a_{mp+r} ≡ a_m b_r (mod p) in Z_(p)."""
from fractions import Fraction as F
def seqs(N):
    a = [F(0), F(6)]; b = [F(1), F(5)]
    for n in range(1, N):
        P = F((2*n+1)*(17*n*n+17*n+5))
        a.append((P*a[n] - n**3*a[n-1]) / F((n+1)**3))
        b.append((P*b[n] - n**3*b[n-1]) / F((n+1)**3))
    return a, b
def vp_mod(x, p):
    # reduce rational x in Z_(p) to F_p
    num, den = x.numerator, x.denominator
    while num % p == 0 and den % p == 0: num //= p; den //= p
    assert den % p != 0, "not p-integral"
    return num * pow(den % p, p-2, p) % p
ok = bad = 0
for p in (7, 11, 13, 17):
    a, b = seqs(3*p + 5)
    for m in range(0, 3):
        for r in range(0, p):
            if m*p + r >= len(a): continue
            lhs = a[m*p+r] * p**3
            try:
                L = vp_mod(lhs, p)
            except AssertionError:
                bad += 1; print("non-integral", p, m, r); continue
            R = vp_mod(a[m], p) * (b[r].numerator % p) % p
            if L == R: ok += 1
            else: bad += 1; print("MISMATCH", p, m, r, L, R)
print(f"polar Lucas p^3 a_(mp+r) ≡ a_m b_r: {ok} ok, {bad} fail")
