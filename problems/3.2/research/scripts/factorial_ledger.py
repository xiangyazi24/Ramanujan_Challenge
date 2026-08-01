#!/usr/bin/env python3
"""(1) Verify universal even-lag type-II: A_d(-(d+1)/2)=1, B_d(-(d+1)/2)=0 for even d.
(2) Factorial-moment ledger: R_2, R_3, C_2 (doubleton classes), R_ge3 over sigma classes;
prediction: R_3/R_2 ~ H/(3p) (doubletons dominate => recursion attacks negligible tail)."""
from collections import defaultdict

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

# (1) even-lag type-II verification
p = 4001
inv2 = pow(2, p-2, p)
ok = bad = 0
for d in range(2, 41, 2):
    rd = (-(d+1)) * inv2 % p
    A, B = AB(d, rd, p)
    if A == 1 and B == 0: ok += 1
    else: bad += 1; print("TYPE-II FAIL", d, A, B)
print(f"universal even-lag type-II centers: {ok} verified, {bad} failed (p={p}, d=2..40 even)")

# (2) factorial ledger
for p in (1009, 4001):
    H = 32
    inv2 = pow(2, p-2, p)
    R2 = R3 = C2 = Rge3 = 0
    for r in range(p):
        Ap_, Bp_ = 1, 0; Ac = Bc = None
        vals = []
        ok2 = True
        for d in range(1, H+1):
            x = (r + d - 1) % p
            if (x+1) % p == 0: ok2 = False; break
            den = pow((x+1) % p, 3, p); di = pow(den, p-2, p)
            a = (2*x+1)*(17*x*x+17*x+5) % p * di % p
            be = (-(x**3)) % p * di % p
            if d == 1: Ac, Bc = a, be
            else:
                Ac, Ap_ = (a*Ac + be*Ap_) % p, Ac
                Bc, Bp_ = (a*Bc + be*Bp_) % p, Bc
            if Bc: vals.append((d, (1-Ac)*pow(Bc, p-2, p) % p))
        if not ok2: continue
        cls = defaultdict(int)
        for d, s in vals: cls[s] += 1
        for s, m in cls.items():
            if m >= 2:
                R2 += m*(m-1)//2
                if m == 2: C2 += 1
                else: Rge3 += m*(m-1)//2
                if m >= 3: R3 += m*(m-1)*(m-2)//6
    print(f"p={p} H={H}: R2={R2} (R2/H^2={R2/H**2:.3f})  doubleton classes C2={C2} "
          f"(share {C2/max(1,R2):.3f})  R>=3 pair-mass {Rge3} ({Rge3/max(1,R2):.3f})  "
          f"R3={R3}  R3/R2={R3/max(1,R2):.4f}  pred H/(3p)={H/(3*p):.4f}")
