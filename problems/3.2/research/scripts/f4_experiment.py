#!/usr/bin/env python3
"""R14's decisive experiment: actual cross-prime factorial moments F_2,F_3,F_4 from the Z_p census
+ reflection-preserving null model. Data: data_zp_pairs.bin (uint32 p, uint32 r pairs, LE)."""
import struct, random, os
from collections import defaultdict

# repo root = four levels up from problems/3.2/research/scripts/
REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..'))

pairs = []
with open(os.path.join(REPO, 'problems/3.2/data_zp_pairs.bin'), 'rb') as f:
    data = f.read()
n_rec = len(data)//8
for i in range(n_rec):
    p, r = struct.unpack_from('<II', data, i*8)
    pairs.append((p, r))
print(f"records: {n_rec}, p range: {pairs[0][0]}..{pairs[-1][0]}")

Zp = defaultdict(list)
for p, r in pairs: Zp[p].append(r)

def moments(N, zdict):
    H = defaultdict(int); L = 0
    for p, zs in zdict.items():
        if p <= N//2 or p > 2*N: continue
        for z in zs:
            if z < p:
                n = p + z
                if N < n <= 2*N: H[n] += 1; L += 1
    F = {1: L, 2: 0, 3: 0, 4: 0}
    mx = 0
    for n, h in H.items():
        mx = max(mx, h)
        F[2] += h*(h-1); F[3] += h*(h-1)*(h-2); F[4] += h*(h-1)*(h-2)*(h-3)
    return L, F, mx

def null_model(zdict, rng):
    out = {}
    for p, zs in zdict.items():
        # keep |Z_p| and reflection pairing: sample |Z|/2 random z (paired with p-1-z); odd size -> one self-paired center (p-1)/2
        m = len(zs); new = []
        k = m // 2
        for _ in range(k):
            z = rng.randrange(0, p)
            new += [z, p-1-z]
        if m % 2: new.append((p-1)//2)
        out[p] = new
    return out

import math
rng = random.Random(99)
for N in (100000, 300000, 800000, 3000000):
    L, F, mx = moments(N, Zp)
    if L == 0: continue
    R = {k: (N**(k-1)*F[k]/L**k if L else 0) for k in (2,3,4)}
    T4 = (math.log(N)**4)*F[4]/N**4
    # null comparison (5 samples for speed)
    nulls = []
    for s in range(5):
        Ln, Fn, _ = moments(N, null_model(Zp, rng))
        nulls.append((Fn[2], Fn[3], Fn[4]))
    nf2 = sum(x[0] for x in nulls)/5; nf3 = sum(x[1] for x in nulls)/5; nf4 = sum(x[2] for x in nulls)/5
    print(f"N={N}: L={L} maxH={mx} F2={F[2]} F3={F[3]} F4={F[4]} | R2={R[2]:.2f} R3={R[3]:.2f} R4={R[4]:.2f} "
          f"T4={T4:.2e} | null(F2,F3,F4)=({nf2:.1f},{nf3:.1f},{nf4:.1f})")
