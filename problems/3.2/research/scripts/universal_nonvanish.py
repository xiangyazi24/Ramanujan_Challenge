#!/usr/bin/env python3
"""Universal nonvanishing sweep: flag ANY (p, h, k) with Psi_{h,k} == 0 identically mod p.
p <= 20000, 1<=h<k<=30. Screen: 40 random eval points; full-degree confirm on suspects."""
import random, sys
PMAX = 20000; KMAX = 30

def sieve(n):
    s = bytearray([1])*(n+1); s[0]=s[1]=0
    for i in range(2,int(n**.5)+1):
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(n+1) if s[i]]

def AB_eval(d, x, p):
    Ap, Bp = 1, 0; Ac = Bc = None
    for j in range(d):
        xx = (x + j) % p
        den = pow(xx+1, 3, p)
        if den == 0: return None, None
        di = pow(den, -1, p)
        a = (2*xx+1)*(17*xx*xx+17*xx+5) % p * di % p
        be = (-(xx**3)) % p * di % p
        if j == 0: Ac, Bc = a, be
        else:
            Ac, Ap = (a*Ac + be*Ap) % p, Ac
            Bc, Bp = (a*Bc + be*Bp) % p, Bc
    return Ac, Bc

rng = random.Random(7)
suspects = []
for p in sieve(PMAX):
    if p <= 2*KMAX+3: continue
    for h in range(1, KMAX):
        for k in range(h+1, KMAX+1):
            zero = True
            pts = tries = 0
            while pts < 30 and tries < 300:
                x = rng.randrange(k+1, p - k - 1); tries += 1
                A1, B1 = AB_eval(h, x, p)
                A2, B2 = AB_eval(k, x, p)
                if A1 is None or A2 is None: continue
                pts += 1
                if ((1-A1)*B2 - (1-A2)*B1) % p:
                    zero = False; break
            if zero and pts >= 30:
                suspects.append((p, h, k))
print("UNIVERSAL NONVANISHING suspects (Psi==0 mod p):", len(suspects), suspects[:20])
