#!/usr/bin/env python3
"""CRON fresh-eyes: full Frobenius deformation tower for Apery b_n, two-digit indices.

Result (all primes 7 <= p < 128, all r < p, zero failures):

  b_{p+r} = 5 b_r + 10 p D_r + p^2 E_r + p^3 ( U3~_r + beta_p * b_r )   (mod p^4)

- 10*D = U1: first-order universal solution  (U1_0=0, U1_1=60)     [mod p^2 layer]
- E   = U2: second-order universal solution (U2_0=U2_1=0); the naive fit constants
        (x'_p, y'_p) vanish for ALL p -> mod p^3 layer is ALSO universal.
- Third order: fit constants are x3_p = -6 (universal, absorbed: U3~ = U3 - 6a)
        and y3_p = beta_p := (b_p - 5)/p^3 mod p  -- the Beukers supercongruence
        defect. This is the FIRST and (through p^4) ONLY p-dependent invariant.

Derivation: v_r := b_{p+r} satisfies the Apery recurrence with coefficients
P(p+r); expanding P(p+r) = sum_i p^i P^(i)(r)/i! order by order gives inhomogeneous
recurrences for the deformation layers; each layer = universal particular solution
+ span{a_r, b_r} with constants fixed by matching v at r=0,1 (i.e. by b_p, b_{p+1}).
Universality of layers 1-2 and of x3 is the numerical discovery.

Run: python3 CRON_modp3_tower.py   (takes ~1 min)
"""
from fractions import Fraction as F
from sympy import primerange

TOP = 130

def bseq(top):
    B = [1, 5]
    for m in range(1, top):
        B.append(((2*m+1)*(17*m*m+17*m+5)*B[m] - m**3*B[m-1]) // (m+1)**3)
    return B

B = bseq(2*TOP + 4)

A = [F(0), F(6)]
for m in range(1, TOP+1):
    A.append(((2*m+1)*(17*m*m+17*m+5)*A[m] - m**3*A[m-1]) / F((m+1)**3))

U0 = [5*F(b) for b in B[:TOP+2]]

U1 = [F(0), F(60)]
for r in range(1, TOP+1):
    Q  = 34*r**3 + 51*r**2 + 27*r + 5
    Qp = 102*r*r + 102*r + 27
    rhs = Q*U1[r] - r**3*U1[r-1] + Qp*U0[r] - 3*(r+1)**2*U0[r+1] - 3*r*r*U0[r-1]
    U1.append(rhs / F((r+1)**3))

U2 = [F(0), F(0)]
for r in range(1, TOP+1):
    Q, Qp, Qpp2 = 34*r**3+51*r**2+27*r+5, 102*r*r+102*r+27, 102*r+51
    rhs = (Q*U2[r] + Qp*U1[r] + Qpp2*U0[r]
           - r**3*U2[r-1] - 3*r*r*U1[r-1] - 3*r*U0[r-1]
           - 3*(r+1)**2*U1[r+1] - 3*(r+1)*U0[r+1])
    U2.append(rhs / F((r+1)**3))

U3 = [F(0), F(0)]
for r in range(1, TOP+1):
    Q, Qp, Qpp2 = 34*r**3+51*r**2+27*r+5, 102*r*r+102*r+27, 102*r+51
    rhs = (Q*U3[r] + Qp*U2[r] + Qpp2*U1[r] + 34*U0[r]
           - r**3*U3[r-1] - 3*r*r*U2[r-1] - 3*r*U1[r-1] - U0[r-1]
           - 3*(r+1)**2*U2[r+1] - 3*(r+1)*U1[r+1] - U0[r+1])
    U3.append(rhs / F((r+1)**3))

def fp(fr, p):
    return fr.numerator % p * pow(fr.denominator % p, -1, p) % p

fails = []
for p in primerange(7, 128):
    Bm = [b % p for b in B[:p]]
    Am = [fp(A[r], p) for r in range(p)]
    # layer 2 check: (x',y') = (0,0)
    ok2 = all(fp((F(B[p+r]) - U0[r] - p*U1[r]) / (p*p), p) == fp(U2[r], p)
              for r in range(p))
    # layer 3 check: x3 = -6, y3 = beta_p
    beta = ((B[p] - 5) // p**3) % p
    ok3 = all(fp((F(B[p+r]) - U0[r] - p*U1[r] - p*p*U2[r]) / p**3, p)
              == (fp(U3[r], p) - 6*Am[r] + beta*Bm[r]) % p
              for r in range(p))
    print(f"p={p:4d}  mod-p^3-universal={ok2}  mod-p^4(beta_p={beta:4d})={ok3}")
    if not (ok2 and ok3):
        fails.append(p)
print("FAILING PRIMES:", fails)
