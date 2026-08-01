#!/usr/bin/env python3
"""Fresh-eyes independent scan for P3.2 pointwise conjecture.

For primes 7 <= p < N compute Z_p = {r in [0,p): p | b_r} via the
division-free recurrence on c_r = (r!)^3 b_r:
    c_{r+1} = (2r+1)(17r^2+17r+5) c_r - r^6 c_{r-1}   (mod p)
(c_r = 0 <=> b_r = 0 mod p since r! is a unit for r < p).

Then H(n) = #{p in (n/2,n] : p | b_n} = #{p : n-p in Z_p}  (Gessel digits).

Outputs: max/record statistics of H(n) (pointwise-relevant), |Z_p| stats,
symmetry check, class-mod-8 breakdown, record structure dump.
"""
import sys, time, json
from math import log

N = int(sys.argv[1]) if len(sys.argv) > 1 else 30000

# ---- sieve ----
sieve = bytearray([1]) * N
sieve[0:2] = b'\x00\x00'
for i in range(2, int(N**0.5) + 1):
    if sieve[i]:
        sieve[i*i::i] = bytearray(len(sieve[i*i::i]))
primes = [i for i in range(7, N) if sieve[i]]

# ---- sanity: b_r integers for small r ----
B = [1, 5]
for n in range(1, 12):
    B.append(((2*n+1)*(17*n*n+17*n+5)*B[n] - n**3*B[n-1]) // (n+1)**3)
# known: 1,5,73,1445,33001,819005,...
assert B[2] == 73 and B[3] == 1445 and B[4] == 33001, B[:5]

# ---- main loop: Z_p ----
t0 = time.time()
ZP = []          # list of (p, tuple_of_zeros)
size_hist = {}
size_by_class = {1: [], 3: [], 5: [], 7: []}
asym_violations = 0
for p in primes:
    cp = 1        # c_0
    cc = 5 % p    # c_1
    zs = []
    for r in range(1, p):
        if cc == 0:
            zs.append(r)
        rr = r * r % p
        r3 = rr * r % p
        r6 = r3 * r3 % p
        coef = (2*r + 1) * ((17*rr + 17*r + 5) % p) % p
        cp, cc = cc, (coef * cc - r6 * cp) % p
    # sanity vs true b_r mod p for small r
    if p > 1000 and p == primes[-1]:
        pass
    ZP.append((p, tuple(zs)))
    k = len(zs)
    size_hist[k] = size_hist.get(k, 0) + 1
    size_by_class[p % 8].append(k)
    for z in zs:
        if (p - 1 - z) not in zs:
            asym_violations += 1
t1 = time.time()

# extra sanity: verify c-recurrence against integer b for one largish p
p = primes[-1]
for r in range(2, 12):
    fact3 = 1
    for i in range(1, r+1):
        fact3 = fact3 * pow(i, 3, p) % p
    # recompute c_r mod p independently
    cp, cc = 1, 5 % p
    for q in range(1, r):
        rr = q*q % p
        r6 = (rr*q % p)**2 % p
        coef = (2*q+1) * ((17*rr + 17*q + 5) % p) % p
        cp, cc = cc, (coef*cc - r6*cp) % p
    assert cc == fact3 * (B[r] % p) % p, (p, r)

# ---- H(n) ----
H = [0] * N
for p, zs in ZP:
    for z in zs:
        n = p + z
        if n < N:
            H[n] += 1

nmin = 1000
tail = H[nmin:]
maxH = max(tail)
hist = {}
for h in tail:
    hist[h] = hist.get(h, 0) + 1

records = sorted((h, n) for n, h in enumerate(H) if n >= nmin and h >= maxH - 1)[-15:]

print(f"N={N}  primes={len(primes)}  Zp-scan time={t1-t0:.1f}s")
print(f"reflection-symmetry violations in Z_p: {asym_violations} (expect 0)")
tot = sum(k*v for k, v in size_hist.items())
print(f"mean |Z_p| = {tot/len(primes):.4f}   (expect ~1.00)")
print("|Z_p| distribution:", dict(sorted(size_hist.items())),)
import statistics
print("mean |Z_p| by p mod 8:",
      {c: round(statistics.mean(v), 4) for c, v in size_by_class.items()})
print(f"\nH(n) stats for n in [{nmin},{N}):")
print("H histogram:", dict(sorted(hist.items())))
print(f"max H = {maxH}")
print("\nrecords (H, n):")
for h, n in records:
    bad = []
    for p, zs in ZP:
        if p <= n < 2*p and (n - p) in zs:
            bad.append(p)
    info = [(p, n - p, round((n-p)/p, 3), p % 8) for p in bad]
    print(f"  n={n} H={h}: (p, z=n-p, z/p, p%8) = {info}")

# element structure: distribution of z/(p-1) in 10 bins, distance to midpoint
bins = [0]*10
mid_hits = 0
tot_z = 0
for p, zs in ZP:
    for z in zs:
        tot_z += 1
        bins[min(9, int(10*z/(p-1)))] += 1
        if 2*z == p - 1:
            mid_hits += 1
print(f"\nz/(p-1) decile distribution (symmetric by reflection): {bins}")
print(f"z at exact midpoint (p-1)/2: {mid_hits}  (supersingular-type)")
