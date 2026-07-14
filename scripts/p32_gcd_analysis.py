#!/usr/bin/env python3
"""P3.2: Deep analysis of gcd(d_n a_n, d_n b_n) for Apéry sequences.

Computes exact gcd, p-adic valuations, and density of 'bad' primes.
"""
from math import gcd, log, isqrt, comb
from fractions import Fraction
from collections import defaultdict
import time

def sieve_primes(N):
    is_prime = [True] * (N + 1)
    is_prime[0] = is_prime[1] = False
    for i in range(2, isqrt(N) + 1):
        if is_prime[i]:
            for j in range(i*i, N + 1, i):
                is_prime[j] = False
    return [p for p in range(2, N + 1) if is_prime[p]]

def lcm_cube(n):
    """Compute lcm(1,...,n)^3."""
    result = 1
    primes = sieve_primes(n)
    for p in primes:
        pk = p
        while pk <= n:
            pk *= p
        result *= (pk // p)
    return result ** 3

def vp(n, p):
    """p-adic valuation of integer n."""
    if n == 0:
        return float('inf')
    v = 0
    while n % p == 0:
        n //= p
        v += 1
    return v

def apery_sequences(N):
    """Compute a_n, b_n as exact Fractions for n=0,...,N using the recurrence."""
    a = [Fraction(0)] * (N + 1)
    b = [Fraction(0)] * (N + 1)
    a[0] = Fraction(0)
    a[1] = Fraction(6)
    b[0] = Fraction(1)
    b[1] = Fraction(5)
    for n in range(1, N):
        coeff_n = 34*n**3 + 51*n**2 + 27*n + 5
        a[n+1] = (Fraction(coeff_n) * a[n] - Fraction(n**3) * a[n-1]) / Fraction((n+1)**3)
        b[n+1] = (Fraction(coeff_n) * b[n] - Fraction(n**3) * b[n-1]) / Fraction((n+1)**3)
    return a, b

NMAX = 200
print(f"Computing Apéry sequences to n={NMAX}...")
t0 = time.time()
a, b = apery_sequences(NMAX)
t1 = time.time()
print(f"  Done in {t1-t0:.1f}s")

# Verify integrality of d_n*a_n, d_n*b_n for small n
print("\n=== Integrality check ===")
for n in [1, 2, 5, 10, 20]:
    dn = lcm_cube(n)
    da = dn * a[n]
    db = dn * b[n]
    assert da.denominator == 1, f"d_{n}*a_{n} not integer!"
    assert db.denominator == 1, f"d_{n}*b_{n} not integer!"
print("  d_n*a_n, d_n*b_n are integers for all tested n.")

# Compute gcd(d_n*a_n, d_n*b_n)
print("\n=== GCD analysis ===")
print(f"{'n':>4s} {'log(gcd)/n':>10s} {'gcd':>20s} {'prime factors':>40s}")
print("-" * 80)

gcd_data = []
for n in range(1, NMAX + 1):
    dn = lcm_cube(n)
    da = int(dn * a[n])
    db = int(dn * b[n])
    if da == 0:
        g = abs(db)
    else:
        g = gcd(abs(da), abs(db))

    ratio = log(g) / n if g > 1 else 0.0
    gcd_data.append((n, g, ratio))

    if n <= 20 or n % 10 == 0:
        # Factor the gcd for display
        g_copy = g
        factors = []
        primes = sieve_primes(min(g_copy, 10000))
        for p in primes:
            if p * p > g_copy:
                break
            e = 0
            while g_copy % p == 0:
                g_copy //= p
                e += 1
            if e > 0:
                factors.append(f"{p}^{e}" if e > 1 else str(p))
        if g_copy > 1:
            factors.append(str(g_copy))
        fstr = " * ".join(factors) if factors else "1"
        print(f"{n:4d} {ratio:10.4f} {g:20d} {fstr:>40s}")

# Analyze density of "bad" primes
print("\n=== Bad prime density ===")
print("For each n, count primes p <= n with p | gcd(d_n a_n, d_n b_n)")
print(f"{'n':>4s} {'#primes<=n':>10s} {'#bad_primes':>12s} {'ratio':>8s} {'bad primes':>30s}")
print("-" * 70)

for n in [10, 20, 30, 50, 70, 100, 150, 200]:
    if n > NMAX:
        break
    dn = lcm_cube(n)
    da = int(dn * a[n])
    db = int(dn * b[n])
    g = gcd(abs(da), abs(db)) if da != 0 else abs(db)

    primes_n = sieve_primes(n)
    bad_primes = [p for p in primes_n if g % p == 0]
    ratio = len(bad_primes) / len(primes_n) if primes_n else 0

    bad_str = str(bad_primes[:10])
    if len(bad_primes) > 10:
        bad_str += "..."
    print(f"{n:4d} {len(primes_n):10d} {len(bad_primes):12d} {ratio:8.4f} {bad_str:>30s}")

# p-adic analysis: for each prime p, track v_p(gcd) as function of n
print("\n=== p-adic valuations of gcd (for small primes) ===")
for p in [2, 3, 5, 7, 11, 13]:
    vals = []
    for n in range(1, min(NMAX + 1, 101)):
        dn = lcm_cube(n)
        da = int(dn * a[n])
        db = int(dn * b[n])
        if da == 0:
            v = vp(abs(db), p) if db != 0 else 0
        else:
            g = gcd(abs(da), abs(db))
            v = vp(g, p)
        vals.append(v)
    max_v = max(vals)
    nonzero_count = sum(1 for v in vals if v > 0)
    print(f"  p={p:2d}: max v_p = {max_v}, nonzero at {nonzero_count}/100 indices, max v_p/n = {max(vals[i]/(i+1) for i in range(len(vals))):.4f}")

# Track log(gcd)/n trend
print("\n=== Trend of log(gcd)/n ===")
# Moving average over windows of 20
window = 20
for start in range(0, NMAX - window + 1, window):
    avg = sum(gcd_data[start + i][2] for i in range(window)) / window
    print(f"  n={start+1:3d}-{start+window:3d}: avg log(gcd)/n = {avg:.4f}")

# Check: proportion of n where gcd = 1
count_trivial = sum(1 for n, g, _ in gcd_data if g == 1)
print(f"\n  gcd = 1 for {count_trivial}/{NMAX} values of n ({100*count_trivial/NMAX:.1f}%)")

# Large prime analysis: for primes p in (n/2, n], check if p | b_n
print("\n=== Large prime divisibility of b_n ===")
for n in [50, 100, 150, 200]:
    if n > NMAX:
        break
    bn = int(b[n].numerator)  # b_n is integer
    primes_half = [p for p in sieve_primes(n) if p > n // 2]
    divides = [p for p in primes_half if bn % p == 0]
    print(f"  n={n:3d}: {len(primes_half)} primes in (n/2,n], {len(divides)} divide b_n = {100*len(divides)/max(1,len(primes_half)):.1f}%")
    if divides:
        print(f"         divisors: {divides[:10]}")

print(f"\nTotal time: {time.time()-t0:.1f}s")
