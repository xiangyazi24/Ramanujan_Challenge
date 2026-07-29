#!/usr/bin/env python3
"""Compute Z(p) = #{j < p : b_j ≡ 0 (mod p)} for primes p up to 10^4.

Uses the recurrence mod p (O(p) per prime) to avoid exact big-integer arithmetic.
For p | (n+1), falls back to exact computation for that index.
"""
import sys
import time
from math import isqrt

def sieve_primes(N):
    if N < 2:
        return []
    is_prime = [True] * (N + 1)
    is_prime[0] = is_prime[1] = False
    for i in range(2, isqrt(N) + 1):
        if is_prime[i]:
            for j in range(i*i, N + 1, i):
                is_prime[j] = False
    return [p for p in range(2, N + 1) if is_prime[p]]

def zp_mod_recurrence(p):
    """Compute b_0,...,b_{p-1} mod p using the recurrence.
    When p | (n+1), we need the exact value --- use a fallback."""
    b = [0] * p
    b[0] = 1
    if p == 2:
        return b
    b[1] = 5 % p

    # For n where p | (n+1), we need exact b_{n+1}.
    # Pre-compute these using the combinatorial formula mod p.
    # b_{p-1} = sum_{k=0}^{p-1} C(p-1,k)^2 C(p-1+k,k)^2
    # By Lucas, C(p-1,k) ≡ (-1)^k (mod p), so C(p-1,k)^2 ≡ 1.
    # C(p-1+k,k) = C(p-1+k, p-1). By Vandermonde/Lucas:
    # C(p-1+k, k) ≡ (-1)^k (mod p) [since C(p-1+k,k) = C(p-1,0)*C(k,k) for digit split].
    # Wait: p-1+k < 2p when k < p. Base-p digits of p-1+k:
    #   if k <= p-1: p-1+k = 1*p + (k-1) when k >= 1, or 0*p + (p-1) when k=0.
    #   Actually p-1+k for k=0 is p-1 (one digit), for k >= 1: p-1+k = p + (k-1),
    #   so digits are (k-1, 1). And k digits are (k, 0) for k < p.
    #   C(p-1+k, k) ≡ C(1,0)*C(k-1,k) for k >= 1... C(k-1,k) = 0 for k >= 1.
    # Hmm that can't be right since b_{p-1} is not zero in general.
    # Let me just compute b_{p-1} mod p exactly when needed.

    for n in range(1, p - 1):
        coeff = (34*n**3 + 51*n**2 + 27*n + 5) % p
        n3 = pow(n, 3, p)
        den = pow(n + 1, 3, p)
        if den == 0:
            # p | (n+1), i.e., n = p-1. But we only go up to n = p-2.
            # Actually n ranges from 1 to p-2, so n+1 from 2 to p-1, never = p.
            # So this shouldn't happen for the range we need.
            # But just in case:
            b[n+1] = compute_b_mod_p_direct(n+1, p)
            continue
        den_inv = pow(den, p - 2, p)
        num = (coeff * b[n] - n3 * b[n-1]) % p
        b[n+1] = (num * den_inv) % p

    return b

def compute_b_mod_p_direct(j, p):
    """Compute b_j mod p using the sum formula."""
    from math import comb
    s = 0
    for k in range(j + 1):
        s += pow(comb(j, k), 2, p) * pow(comb(j + k, k), 2, p)
        s %= p
    return s

PMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 10000
print(f"Computing Z(p) for primes up to {PMAX}...")
sys.stdout.flush()
t0 = time.time()

primes = sieve_primes(PMAX)
print(f"  {len(primes)} primes to process")
sys.stdout.flush()

zp_data = []
zp_hist = {}
total_z = 0
count = 0

for idx, p in enumerate(primes):
    if p < 5:
        continue

    b_mod = zp_mod_recurrence(p)
    zeros = [j for j in range(p) if b_mod[j] == 0]
    zp = len(zeros)
    zp_data.append((p, zp))
    zp_hist[zp] = zp_hist.get(zp, 0) + 1
    total_z += zp
    count += 1

    if (idx + 1) % 200 == 0:
        elapsed = time.time() - t0
        print(f"  p={p} ({idx+1}/{len(primes)}), elapsed {elapsed:.1f}s, "
              f"avg Z(p)={total_z/count:.3f}", flush=True)

elapsed = time.time() - t0
print(f"\nDone in {elapsed:.1f}s")

# Summary
print(f"\n=== Z(p) histogram (p ≤ {PMAX}) ===")
for z in sorted(zp_hist.keys()):
    pct = 100 * zp_hist[z] / count
    print(f"  Z(p) = {z}: {zp_hist[z]} primes ({pct:.1f}%)")

print(f"\n=== Z(p) statistics ===")
print(f"  Total primes (≥5): {count}")
print(f"  Mean Z(p): {total_z/count:.4f}")
print(f"  Max Z(p): {max(z for _, z in zp_data)}")
print(f"  Primes with Z(p)=0: {zp_hist.get(0,0)} ({100*zp_hist.get(0,0)/count:.1f}%)")

# Average by range
print(f"\n=== Average Z(p) by range ===")
ranges = [(5,100), (100,500), (500,1000), (1000,2000), (2000,5000), (5000,10000)]
for lo, hi in ranges:
    subset = [z for p, z in zp_data if lo <= p < hi]
    if subset:
        print(f"  p ∈ [{lo}, {hi}): avg={sum(subset)/len(subset):.3f}, "
              f"max={max(subset)}, n={len(subset)}")

# Power-law fit
print(f"\n=== Power-law fit ===")
from math import log
xs = [log(p) for p, z in zp_data if z > 0]
ys = [log(z) for p, z in zp_data if z > 0]
if len(xs) > 10:
    n_pts = len(xs)
    x_mean = sum(xs) / n_pts
    y_mean = sum(ys) / n_pts
    ss_xy = sum((xs[i] - x_mean) * (ys[i] - y_mean) for i in range(n_pts))
    ss_xx = sum((xs[i] - x_mean) ** 2 for i in range(n_pts))
    alpha = ss_xy / ss_xx
    C = 2.718281828**(y_mean - alpha * x_mean)
    print(f"  Z(p) ~ {C:.3f} * p^{alpha:.4f}")
    if abs(alpha) < 0.1:
        print(f"  => Z(p) = O(1) — supports Hypothesis Z")

# Poisson comparison
print(f"\n=== Poisson model comparison ===")
z0_frac = zp_hist.get(0, 0) / count
import math
print(f"  P(Z=0) observed: {z0_frac:.4f}")
print(f"  P(Z=0) Poisson(mean=0.5): {math.exp(-0.5):.4f}")
print(f"  (Each prime has ~(p+1)/2 independent values due to symmetry)")

# Verify symmetry b_j ≡ b_{p-1-j} (mod p)
print(f"\n=== Symmetry verification (sample) ===")
sym_ok = 0
sym_total = 0
for p in primes[:50]:
    if p < 5:
        continue
    b_mod = zp_mod_recurrence(p)
    for j in range(p // 2):
        sym_total += 1
        if b_mod[j] == b_mod[p - 1 - j]:
            sym_ok += 1
print(f"  b_j ≡ b_{{p-1-j}} (mod p): {sym_ok}/{sym_total} = {100*sym_ok/sym_total:.1f}%")

# Output data for selected large primes
print(f"\n=== Z(p) for largest primes ===")
for p, z in zp_data[-20:]:
    print(f"  p={p}: Z(p)={z}")
