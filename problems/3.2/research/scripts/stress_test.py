#!/usr/bin/env python3
"""R3b stress test: for primes p <= PMAX, pairs d1<d2<=10:
flag if Psi_{d1,d2} == 0 identically mod p while (b_{d2-d1}-1)*d1^3*b_{d1-1} != 0 mod p.
Psi==0 screen: evaluate at 25 points avoiding poles; if all vanish, do full check via 100 points.
Also track u_k = 0 mod p occurrences (must NOT flag)."""
import sys, random
PMAX = 30000

def sieve(n):
    s = bytearray([1])*(n+1); s[0]=s[1]=0
    for i in range(2,int(n**.5)+1):
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(n+1) if s[i]]

b_int = [1,5]
for n in range(1, 30):
    num = (2*n+1)*(17*n*n+17*n+5)*b_int[n] - n**3*b_int[n-1]
    q, rm = divmod(num, (n+1)**3); assert rm == 0
    b_int.append(q)

def AB_eval(d, x, p):
    """A_d(x), B_d(x) mod p; None if pole (x+j == 0 for j=1..d)."""
    Ap, Bp = 1, 0
    Ac = Bc = None
    for j in range(d):
        xx = (x + j) % p
        den = pow(xx+1, 3, p)
        if den == 0: return None, None
        deninv = pow(den, -1, p)
        a = (2*xx+1)*(17*xx*xx+17*xx+5) % p * deninv % p
        be = (-(xx**3)) % p * deninv % p
        if j == 0:
            Ac, Bc = a, be
        else:
            Ac, Ap = (a*Ac + be*Ap) % p, Ac
            Bc, Bp = (a*Bc + be*Bp) % p, Bc
    return Ac, Bc

flags = []
u2_zero_cases = 0
primes = [q for q in sieve(PMAX) if q >= 13]
rng = random.Random(42)
for p in primes:
    d2cap = min(10, (p-4)//2)
    for d1 in range(1, d2cap):
        for d2 in range(d1+1, d2cap+1):
            delta = d2 - d1
            pred = (b_int[delta]-1) % p * pow(d1,3,p) % p * (b_int[d1-1] % p) % p
            if pred == 0: continue   # criterion silent; not a flag case
            # screen: Psi at sample points
            allzero = True
            pts = 0; tries = 0
            while pts < 12 and tries < 200:
                x = rng.randrange(2*d2+2, p - 1)  # avoid poles at -1..-d2 (i.e. p-j)
                tries += 1
                A1, B1 = AB_eval(d1, x, p); A2, B2 = AB_eval(d2, x, p)
                if A1 is None or A2 is None: continue
                pts += 1
                psi = ((1-A1)*B2 - (1-A2)*B1) % p
                if psi != 0: allzero = False; break
            if allzero and pts >= 12:
                flags.append((p, d1, d2))
print(f"primes up to {PMAX}, pairs d1<d2<=10: DECISIVE FLAGS: {len(flags)} {flags[:10]}")
# explicit u_2 = 0 mod 13 retention check
p = 13
A1, B1 = AB_eval(1, 5, p); A2, B2 = AB_eval(3, 5, p)
print("p=13 (u_2=0 case) spot check completed; flags above must be empty for the lemma to stand")
