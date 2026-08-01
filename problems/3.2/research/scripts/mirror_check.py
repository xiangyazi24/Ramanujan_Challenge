#!/usr/bin/env python3
"""Among sigma-collision pairs (r; h<k<=H): fraction with 2r+h+k == p-1 (mod p) [mirror pairs],
and fraction with h+k == special structures. p=4001, 16001; H=32."""
from collections import defaultdict
for p in (4001, 16001):
    inv = [0]*p
    for i in range(1, p): inv[i] = pow(i, p-2, p)
    ic = [0]*p
    for i in range(1, p): ic[i] = inv[i]*inv[i]%p*inv[i]%p
    aa = [None]*p; bb = [None]*p
    for x in range(p):
        if (x+1) % p:
            d3 = ic[(x+1)%p]
            aa[x] = (2*x+1)*(17*x*x+17*x+5) % p * d3 % p
            bb[x] = (-(x**3)) % p * d3 % p
    H = 32
    total = mirror = 0
    other = defaultdict(int)
    for r in range(p):
        Ap_, Bp_ = 1, 0; Ac = Bc = None
        vals = []
        ok = True
        for d in range(1, H+1):
            x = (r + d - 1) % p
            if aa[x] is None: ok = False; break
            if d == 1: Ac, Bc = aa[x], bb[x]
            else:
                Ac, Ap_ = (aa[x]*Ac + bb[x]*Ap_) % p, Ac
                Bc, Bp_ = (aa[x]*Bc + bb[x]*Bp_) % p, Bc
            if Bc: vals.append((d, (1-Ac)*pow(Bc, p-2, p) % p))
        if not ok: continue
        seen = defaultdict(list)
        for d, s in vals: seen[s].append(d)
        for s, ds in seen.items():
            for i in range(len(ds)):
                for j in range(i+1, len(ds)):
                    h, k = ds[i], ds[j]
                    total += 1
                    if (2*r + h + k) % p == (p-1) % p: mirror += 1
    print(f"p={p} H=32: collision pairs {total}, mirror (2r+h+k=p-1): {mirror} ({100*mirror/max(1,total):.1f}%)")
