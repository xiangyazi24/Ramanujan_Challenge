#!/usr/bin/env python3
"""|Z_p| growth statistics: calibrate the vertical zero-density target.

Z_p = {0<=r<p : b_r == 0 mod p}, Apery. Proved bound |Z_p| << p^{2/3}.
Question: empirical truth — bounded? log p? p^{1/2}?
Prints running stats every ~5s (long-run discipline).
"""
import time, math
from sympy import primerange

def zp_size(p):
    cnt = 0
    bprev, bcur = 1 % p, 5 % p
    if bprev == 0: cnt += 1
    if bcur == 0: cnt += 1
    for n in range(1, p - 1):
        num = ((34*n*n*n + 51*n*n + 27*n + 5) * bcur - n*n*n * bprev) % p
        bnext = num * pow(pow(n+1, 3, p), p - 2, p) % p
        if bnext == 0: cnt += 1
        bprev, bcur = bcur, bnext
    return cnt

P = 30000
t0 = time.time(); last = t0
sizes = []
maxrec = []
for p in primerange(7, P):
    s = zp_size(p)
    sizes.append((p, s))
    if not maxrec or s > maxrec[-1][1]:
        maxrec.append((p, s))
    now = time.time()
    if now - last > 5:
        n = len(sizes)
        mean = sum(x[1] for x in sizes)/n
        print(f"[progress] p={p} n={n} mean|Z|={mean:.3f} max={maxrec[-1]} {now-t0:.0f}s", flush=True)
        last = now

n = len(sizes)
vals = [s for _, s in sizes]
mean = sum(vals)/n
var = sum(v*v for v in vals)/n - mean*mean
zero_frac = sum(1 for v in vals if v == 0)/n
print(f"\nFINAL p<{P}: n={n} mean|Z_p|={mean:.4f} var={var:.4f} P(|Z|=0)={zero_frac:.4f} (Poisson(la) would give e^-la={math.exp(-mean):.4f})")
from collections import Counter
c = Counter(vals)
print("distribution:", dict(sorted(c.items())))
print("records (p, |Z_p|):", maxrec)
# growth comparison at records
for p, s in maxrec:
    print(f"  record p={p}: |Z|={s}  p^(1/2)={p**0.5:.1f}  p^(2/3)={p**(2/3):.1f}  log p={math.log(p):.1f}")
# second moment: sum |Z_p|^2 vs sum |Z_p| (pair-collision intensity)
s1 = sum(vals); s2 = sum(v*v for v in vals)
print(f"S1=sum|Z|={s1} S2=sum|Z|^2={s2} S2/S1={s2/s1:.3f} (bounded => variance harmless for GARQI)")
