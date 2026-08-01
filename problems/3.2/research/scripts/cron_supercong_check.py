#!/usr/bin/env python3
"""Verify cron-Fable's new supercongruence: b_{p+r} ≡ 5 b_r + 10 p D_r (mod p^2),
D_r = sum_k C(r,k)^2 C(r+k,k)^2 (H_{r+k} - H_{r-k}). Then extract W_p(r) =
(b_{p+r} - 5b_r - 10pD_r)/p^2 mod p and probe its structure (Apery-operator residual)."""
from fractions import Fraction as F
from math import comb

def harmonic(n):
    return sum(F(1, j) for j in range(1, n+1))

def apery_b(N):
    b = [1, 5]
    for n in range(1, N):
        num = (2*n+1)*(17*n*n+17*n+5)*b[n] - n**3*b[n-1]
        q, rm = divmod(num, (n+1)**3); assert rm == 0
        b.append(q)
    return b

def D(r):
    s = F(0)
    for k in range(r+1):
        s += comb(r,k)**2 * comb(r+k,k)**2 * (harmonic(r+k) - harmonic(r-k))
    return s

ok = bad = 0
Wtab = {}
for p in (7, 11, 13, 17, 19, 23, 29):
    b = apery_b(2*p + 3)
    W = []
    for r in range(0, (p+1)//2):
        Dr = D(r)
        assert Dr.denominator % p != 0
        lhs = b[p+r]
        rhs = 5*b[r] + 10*p*Dr
        diff = F(lhs) - rhs
        # check mod p^2: diff = p^2 * w with w p-integral
        num, den = diff.numerator, diff.denominator
        assert den % p != 0
        if num % (p*p) == 0 or (num * pow(den, -1, p*p)) % (p*p) == 0:
            ok += 1
            w = (num // 1)  # compute w = diff/p^2 mod p
            val = (num * pow(den, -1, p**3)) % p**3
            assert val % (p*p) == 0
            W.append(val // (p*p) % p)
        else:
            bad += 1; print("FAIL", p, r)
    Wtab[p] = W
print(f"supercongruence b_(p+r) = 5b_r + 10pD_r mod p^2: {ok} verified, {bad} failed")
# probe W_p structure: apply Apery operator L[w](r) = (r+1)^3 w_{r+1} - P(r) w_r + r^3 w_{r-1} mod p
for p in (13, 17, 19):
    W = Wtab[p]
    res = []
    for r in range(1, len(W)-1):
        v = (pow(r+1,3,p)*W[r+1] - (2*r+1)*(17*r*r+17*r+5)%p*W[r] + pow(r,3,p)*W[r-1]) % p
        res.append(v)
    print(f"p={p}: W_p = {W[:8]}...  Apery-operator residual = {res[:8]}...")
