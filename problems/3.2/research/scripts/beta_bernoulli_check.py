#!/usr/bin/env python3
"""Independent verification: beta_p = (b_p - 5)/p^3 = -(14/3) B_{p-3} (mod p)."""
from fractions import Fraction as F

def bernoulli_upto(N):
    B = [F(1)]
    for m in range(1, N+1):
        s = F(0)
        for k in range(m):
            from math import comb
            s += comb(m+1, k) * B[k]
        B.append(-s / (m+1))
    return B

def apery_b(N):
    b = [1, 5]
    for n in range(1, N):
        num = (2*n+1)*(17*n*n+17*n+5)*b[n] - n**3*b[n-1]
        q, r = divmod(num, (n+1)**3); assert r == 0
        b.append(q)
    return b

def sieve(n):
    s = bytearray([1])*(n+1); s[0]=s[1]=0
    for i in range(2,int(n**.5)+1):
        if s[i]: s[i*i::i]=bytearray(len(s[i*i::i]))
    return [i for i in range(n+1) if s[i]]

primes = [p for p in sieve(140) if p >= 7]
B = bernoulli_upto(max(primes))
b = apery_b(max(primes)+2)
ok = bad = 0
for p in primes:
    beta = ((b[p] - 5) // p**3) % p
    Bp3 = B[p-3]
    assert Bp3.denominator % p != 0
    rhs = (-14) * pow(3, p-2, p) % p * (Bp3.numerator % p) % p * pow(Bp3.denominator % p, p-2, p) % p
    if beta == rhs % p: ok += 1
    else:
        bad += 1
        print("MISMATCH p=", p, beta, rhs)
print(f"beta_p = -(14/3) B_(p-3) mod p: {ok} verified, {bad} failed (p in [7, 140])")
# Wolstenholme check: p=16843 would need big b — skip; check p=7 special: beta_7?
p = 7
print("p=7: beta =", ((b[7]-5)//343) % 7, " (cron: =0, coefficient degeneracy 14=2*7)")
