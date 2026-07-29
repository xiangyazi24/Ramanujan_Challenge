#!/usr/bin/env python3
"""Compute Z(p) for all primes up to 10^6 (or argv[1]).
Optimized: tight inner loop, array storage, progress reporting."""
import sys
import time
from math import isqrt
from array import array

def sieve_primes(N):
    is_prime = bytearray(N + 1)
    for i in range(2, N + 1):
        is_prime[i] = 1
    for i in range(2, isqrt(N) + 1):
        if is_prime[i]:
            for j in range(i*i, N + 1, i):
                is_prime[j] = 0
    return [p for p in range(2, N + 1) if is_prime[p]]

def count_zeros(p):
    """Count #{j in [0,p) : b_j ≡ 0 (mod p)} using recurrence mod p."""
    if p < 5:
        return 0
    b_prev = 1  # b[0]
    b_curr = 5 % p  # b[1]
    z = 0
    if b_prev == 0: z += 1
    if b_curr == 0: z += 1

    for n in range(1, p - 1):
        coeff = (34*n*n*n + 51*n*n + 27*n + 5) % p
        n3 = n * n * n % p
        np1_3 = (n + 1) ** 3 % p
        if np1_3 == 0:
            # p | (n+1), so n = p-1. This shouldn't happen in range [1, p-2].
            break
        inv = pow(np1_3, p - 2, p)
        b_next = (coeff * b_curr - n3 * b_prev) % p * inv % p
        b_prev = b_curr
        b_curr = b_next
        if b_curr == 0:
            z += 1
    return z

PMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 1000000
print(f"Computing Z(p) for primes 5 <= p <= {PMAX}...", flush=True)
t0 = time.time()

primes = sieve_primes(PMAX)
primes5 = [p for p in primes if p >= 5]
nprimes = len(primes5)
print(f"  {nprimes} primes to process", flush=True)

zp_data = []
zp_hist = {}
total_z = 0
max_z = 0
max_z_p = 0

report_interval = max(1, nprimes // 100)
for idx, p in enumerate(primes5):
    z = count_zeros(p)
    zp_data.append((p, z))
    zp_hist[z] = zp_hist.get(z, 0) + 1
    total_z += z
    if z > max_z:
        max_z = z
        max_z_p = p
        print(f"  NEW MAX Z(p) = {z} at p = {p}", flush=True)

    if (idx + 1) % report_interval == 0:
        elapsed = time.time() - t0
        pct = 100 * (idx + 1) / nprimes
        avg = total_z / (idx + 1)
        rate = (idx + 1) / elapsed
        eta = (nprimes - idx - 1) / rate
        print(f"  [{pct:5.1f}%] p={p}, avg Z={avg:.3f}, max Z={max_z} (p={max_z_p}), "
              f"elapsed {elapsed:.0f}s, ETA {eta:.0f}s", flush=True)

elapsed = time.time() - t0
count = len(zp_data)
print(f"\nDone in {elapsed:.1f}s ({count} primes)")

print(f"\n=== Z(p) histogram (p <= {PMAX}) ===")
for z in sorted(zp_hist.keys()):
    pct = 100 * zp_hist[z] / count
    print(f"  Z(p) = {z:3d}: {zp_hist[z]:6d} primes ({pct:6.2f}%)")

print(f"\n=== Summary ===")
print(f"  Total primes (>= 5): {count}")
print(f"  Mean Z(p): {total_z/count:.4f}")
print(f"  Max Z(p): {max_z} at p = {max_z_p}")
print(f"  Primes with Z(p)=0: {zp_hist.get(0,0)} ({100*zp_hist.get(0,0)/count:.1f}%)")

# Poisson comparison
import math
poi_half = [math.exp(-0.5) * 0.5**k / math.factorial(k) for k in range(10)]
print(f"\n=== Poisson(1/2) pair-count comparison ===")
print(f"  {'Z(p)':>5}  {'Observed':>10}  {'Poi(1/2)':>10}")
for k in range(0, min(max_z + 1, 20), 2):
    obs = zp_hist.get(k, 0) / count
    if k == 0:
        poi = poi_half[0]
    elif k == 1:
        poi = 0  # pairs
    else:
        poi = poi_half[k // 2]
    print(f"  {k:5d}  {obs:10.4f}  {poi:10.4f}")

# Average by range
print(f"\n=== Average Z(p) by range ===")
boundaries = [5, 100, 1000, 10000, 100000, 500000, 1000000, 10000000]
for i in range(len(boundaries) - 1):
    lo, hi = boundaries[i], boundaries[i + 1]
    if lo >= PMAX:
        break
    hi = min(hi, PMAX + 1)
    subset = [z for p, z in zp_data if lo <= p < hi]
    if subset:
        print(f"  p in [{lo:>7d}, {hi:>7d}): avg={sum(subset)/len(subset):.3f}, "
              f"max={max(subset)}, n={len(subset)}")

# Output notable primes
print(f"\n=== Primes with Z(p) >= 6 ===")
notable = [(p, z) for p, z in zp_data if z >= 6]
for p, z in notable[:50]:
    print(f"  p={p}: Z(p)={z}")
if len(notable) > 50:
    print(f"  ... ({len(notable)} total)")

# Power-law fit
from math import log
xs = [log(p) for p, z in zp_data if z > 0]
ys = [log(z) for p, z in zp_data if z > 0]
if len(xs) > 100:
    n_pts = len(xs)
    x_mean = sum(xs) / n_pts
    y_mean = sum(ys) / n_pts
    ss_xy = sum((xs[i] - x_mean) * (ys[i] - y_mean) for i in range(n_pts))
    ss_xx = sum((xs[i] - x_mean) ** 2 for i in range(n_pts))
    alpha = ss_xy / ss_xx
    C = math.exp(y_mean - alpha * x_mean)
    print(f"\n=== Power-law fit ===")
    print(f"  Z(p) ~ {C:.3f} * p^{alpha:.4f}")
    if abs(alpha) < 0.1:
        print(f"  => Z(p) = O(1) — supports Hypothesis Z")
