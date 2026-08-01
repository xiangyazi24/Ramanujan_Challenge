#!/usr/bin/env python3
"""Decisive rho measurement (R5's experiment):
For sample primes and dyadic H: R_H = total F_p-roots of Psi_{h,k} over pairs h<k<=H;
T_H = actual fiber triples (r, r+h, r+k same value, k<=H); wastage T_H/R_H;
report R_H/H^2 and R_H/H^3. rho=2 <=> R_H/H^2 bounded."""
import sys
from collections import defaultdict

def is_prime(n):
    if n < 2: return False
    d, s = n-1, 0
    while d % 2 == 0: d //= 2; s += 1
    for a in (2,3,5,7,11,13,17,19,23,29,31,37):
        if a >= n: continue
        x = pow(a, d, n)
        if x in (1, n-1): continue
        for _ in range(s-1):
            x = x*x % n
            if x == n-1: break
        else: return False
    return True

def run(p, HMAX):
    # transfer arrays A_d(r), B_d(r) for all r in [0,p), d <= HMAX, computed columnwise:
    # For each r, iterate d: A[d] etc. That's p*HMAX work with modular inverse per (r,d): precompute inverses.
    inv = [0]*p
    for i in range(1, p): inv[i] = pow(i, p-2, p)
    invcube = [0]*p
    for i in range(1, p): invcube[i] = inv[i]*inv[i]%p*inv[i]%p
    # coefficients a(x), beta(x)
    aa = [0]*p; bb = [0]*p
    for x in range(p):
        d3 = invcube[(x+1) % p] if (x+1) % p else None
        if d3 is None: aa[x] = bb[x] = None; continue
        aa[x] = (2*x+1)*(17*x*x+17*x+5) % p * d3 % p
        bb[x] = (-(x**3)) % p * d3 % p
    Rcount = defaultdict(int)   # H -> root incidences with k <= H
    for r in range(p):
        # iterate d, track (A_d, B_d) and detect sigma-collisions among type-I lags
        Ap_, Bp_ = 1, 0
        Ac = Bc = None
        sig = {}
        ok = True
        vals = []
        for d in range(1, HMAX+1):
            x = (r + d - 1) % p
            if aa[x] is None: ok = False; break
            if d == 1:
                Ac, Bc = aa[x], bb[x]
            else:
                Ac, Ap_ = (aa[x]*Ac + bb[x]*Ap_) % p, Ac
                Bc, Bp_ = (aa[x]*Bc + bb[x]*Bp_) % p, Bc
            if Bc % p:
                s = (1 - Ac) % p * pow(Bc, p-2, p) % p
                vals.append((d, s))
        if not ok: continue
        seen = defaultdict(list)
        for d, s in vals: seen[s].append(d)
        for s, ds in seen.items():
            if len(ds) >= 2:
                m = len(ds)
                # each pair (h,k) contributes a root at r; assign to H-bucket by k
                for i in range(m):
                    for j in range(i+1, m):
                        k = ds[j]
                        Rcount[k] += 1
    # actual triples
    b = [1,5]
    for n in range(1, p-1):
        b.append(((2*n+1)*(17*n*n+17*n+5)%p*b[n] - pow(n,3,p)*b[n-1]) * invcube[(n+1)%p] % p)
    fib = defaultdict(list)
    for r, v in enumerate(b): fib[v].append(r)
    Tcount = defaultdict(int)
    for v, rs in fib.items():
        L = len(rs)
        for i in range(L):
            for j in range(i+1, L):
                for l in range(j+1, L):
                    k = rs[l] - rs[i]
                    if k <= HMAX: Tcount[k] += 1
    Hs = [4, 8, 16, 32] if HMAX >= 32 else [4, 8, HMAX]
    out = []
    for H in Hs:
        RH = sum(c for k, c in Rcount.items() if k <= H)
        TH = sum(c for k, c in Tcount.items() if k <= H)
        out.append((H, RH, round(RH/H**2, 2), round(RH/H**3, 3), TH))
    return out

for p in (1009, 4001, 16001, 40009):
    if not is_prime(p): p += 2
    res = run(p, 32)
    print(f"p={p}: " + "  ".join(f"[H={H}: R={R} R/H^2={r2} R/H^3={r3} T={T}]" for H,R,r2,r3,T in res))
