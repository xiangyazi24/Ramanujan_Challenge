#!/usr/bin/env python3
"""Verify: Psi_{h,k}(r*) == 0 mod p at r* = (p-1-h-k)/2 (mod p), exactly when h == k (mod 2).
Test all pairs h<k<=12 at p = 101, 199, 4001."""
def AB(d, rr, p):
    Ap_, Bp_ = 1, 0; Ac = Bc = None
    for j in range(d):
        x = (rr + j) % p
        den = pow((x+1) % p, 3, p)
        if den == 0: return None, None
        di = pow(den, p-2, p)
        a = (2*x+1)*(17*x*x+17*x+5) % p * di % p
        be = (-(x**3)) % p * di % p
        if j == 0: Ac, Bc = a, be
        else:
            Ac, Ap_ = (a*Ac + be*Ap_) % p, Ac
            Bc, Bp_ = (a*Bc + be*Bp_) % p, Bc
    return Ac, Bc
for p in (101, 199, 4001):
    good = bad = wrongparity_zero = 0
    inv2 = pow(2, p-2, p)
    for h in range(1, 12):
        for k in range(h+1, 13):
            rstar = (p-1-h-k) * inv2 % p
            A1, B1 = AB(h, rstar, p); A2, B2 = AB(k, rstar, p)
            if A1 is None or A2 is None: continue
            psi = ((1-A1)*B2 - (1-A2)*B1) % p
            if (h-k) % 2 == 0:
                if psi == 0: good += 1
                else: bad += 1
            else:
                if psi == 0: wrongparity_zero += 1
    print(f"p={p}: same-parity pairs with Psi(r*)=0: {good}, FAILURES: {bad}, odd-parity accidental zeros: {wrongparity_zero}")
