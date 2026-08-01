#!/usr/bin/env python3
"""CRON fresh-eyes: universal structure of B_p(x) = sum_{r<p} b_r x^r mod p.

Established numerically (all primes 13 <= p < 150, zero failures):
Let q = 1 - 34x + x^2, sigma = sqrt(F/q), tau = sqrt(F) in Z[1/2][[x]],
chi(p) = Legendre(-6, p).

  chi = -1  =>  B_p = q * ([sigma]_{<=(p-3)/2})^2   (mod p)
  chi = +1  =>  B_p =     ([tau]_{<=(p-1)/2})^2     (mod p)

plus the FULL GAP: the relevant series' coefficients vanish mod p on the whole
interval (truncation_degree, p-1).  So the p-arithmetic of the truncated
generating polynomial is: two fixed universal series + chi(p) + truncation point.
"""
from fractions import Fraction as F
from sympy import primerange, Poly, symbols
x = symbols('x')

PMAX = 150
Bint = [1, 5]
for m in range(1, 2*PMAX + 20):
    Bint.append(((2*m+1)*(17*m*m+17*m+5)*Bint[m] - m**3*Bint[m-1]) // (m+1)**3)

NN = PMAX + 10
Fc = [F(b) for b in Bint[:NN]]
G = [F(0)]*NN
for n in range(NN):
    s = Fc[n]
    if n >= 1: s -= -34*G[n-1]
    if n >= 2: s -= G[n-2]
    G[n] = s
S = [F(0)]*NN; S[0] = F(1)
for n in range(1, NN):
    S[n] = (G[n] - sum(S[i]*S[n-i] for i in range(1, n))) / 2
T = [F(0)]*NN; T[0] = F(1)
for n in range(1, NN):
    T[n] = (Fc[n] - sum(T[i]*T[n-i] for i in range(1, n))) / 2

def fp(fr, p): return fr.numerator % p * pow(fr.denominator % p, -1, p) % p

def bseq_mod(p):
    B = [1, 5 % p]
    for m in range(1, p-1):
        inv = pow((m+1) % p, -3, p)
        B.append((((2*m+1)*(17*m*m+17*m+5) % p)*B[m] - pow(m,3,p)*B[m-1]) * inv % p)
    return B

bad = []
for p in primerange(13, PMAX):
    chi = 1 if pow(-6 % p, (p-1)//2, p) == 1 else -1
    U, dtop = (T, (p-1)//2) if chi == 1 else (S, (p-3)//2)
    # gap check
    gap_ok = all(fp(U[j], p) == 0 for j in range(dtop+1, min(p-1, NN)))
    # reconstruction check
    tr = [fp(U[j], p) for j in range(dtop+1)]
    prod = [0]*p
    for i in range(dtop+1):
        if tr[i] == 0: continue
        for j in range(dtop+1):
            if i+j < p:
                prod[i+j] = (prod[i+j] + tr[i]*tr[j]) % p
    if chi == -1:
        qprod = [0]*p
        for k in range(p):
            v = prod[k]
            acc = v
            if k >= 1: acc = (acc - 34*prod[k-1]) % p
            if k >= 2: acc = (acc + prod[k-2]) % p
            qprod[k] = acc % p
        prod = qprod
    B = bseq_mod(p)
    rec_ok = all(prod[r] == B[r] for r in range(p))
    print(f"p={p:4d} chi={chi:+d} gap={gap_ok} reconstruct={rec_ok}")
    if not (gap_ok and rec_ok): bad.append(p)
print("FAILING PRIMES:", bad)
