#!/usr/bin/env python3
"""P3.2: Extended GCD analysis to n=500 with deeper p-adic structure.

Key questions:
1. Does log(gcd)/n -> 0? What's the rate?
2. How does the bad prime count grow with n?
3. For primes p > n/2, what fraction divides b_n? (Expected O(1))
4. Lucas-type analysis: does p | b_n depend on n mod p?
"""
from math import gcd, log, isqrt, comb
from fractions import Fraction
import time
import sys

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

def lcm_power(n):
    """Compute lcm(1,...,n) as integer."""
    result = 1
    primes = sieve_primes(n)
    for p in primes:
        pk = p
        while pk <= n:
            pk *= p
        result *= (pk // p)
    return result

def vp(n, p):
    if n == 0:
        return float('inf')
    v = 0
    n = abs(n)
    while n % p == 0:
        n //= p
        v += 1
    return v

NMAX = 500
print(f"Computing Apéry sequences to n={NMAX}...")
sys.stdout.flush()
t0 = time.time()

# Use integer arithmetic where possible for speed.
# b_n satisfies: (n+1)^3 b_{n+1} = (34n^3+51n^2+27n+5) b_n - n^3 b_{n-1}
# b_n is always integer, so we can use pure int arithmetic.

b = [0] * (NMAX + 1)
b[0] = 1
b[1] = 5
for n in range(1, NMAX):
    coeff = 34*n**3 + 51*n**2 + 27*n + 5
    num = coeff * b[n] - n**3 * b[n-1]
    den = (n+1)**3
    assert num % den == 0, f"b_{n+1} not integer! n={n}"
    b[n+1] = num // den

t1 = time.time()
print(f"  b_n computed in {t1-t0:.1f}s")
sys.stdout.flush()

# a_n needs Fraction arithmetic (it's rational)
a = [Fraction(0)] * (NMAX + 1)
a[0] = Fraction(0)
a[1] = Fraction(6)
for n in range(1, NMAX):
    coeff = 34*n**3 + 51*n**2 + 27*n + 5
    a[n+1] = (Fraction(coeff) * a[n] - Fraction(n**3) * a[n-1]) / Fraction((n+1)**3)
    if n % 100 == 0:
        print(f"  a_n: n={n}...", flush=True)

t2 = time.time()
print(f"  a_n computed in {t2-t1:.1f}s")
sys.stdout.flush()

# Compute d_n = lcm(1,...,n)^3 and gcd(d_n a_n, d_n b_n)
print(f"\nComputing gcds...")
sys.stdout.flush()

gcd_data = []
bad_prime_counts = []
all_primes = sieve_primes(NMAX)

for n in range(1, NMAX + 1):
    dn3 = lcm_power(n) ** 3
    da = int(dn3 * a[n])
    db = dn3 * b[n]  # already integer

    if da == 0:
        g = abs(db)
    else:
        g = gcd(abs(da), abs(db))

    ratio = log(g) / n if g > 1 else 0.0
    gcd_data.append((n, g, ratio))

    if n % 50 == 0:
        # Count bad primes
        primes_n = [p for p in all_primes if p <= n]
        bad = [p for p in primes_n if g % p == 0]
        bad_prime_counts.append((n, len(primes_n), len(bad)))
        print(f"  n={n}: log(gcd)/n = {ratio:.4f}, bad primes: {len(bad)}/{len(primes_n)} = {len(bad)/len(primes_n):.3f}", flush=True)

t3 = time.time()
print(f"  GCDs computed in {t3-t2:.1f}s")

# === Detailed output ===

# 1. Trend of log(gcd)/n
print("\n=== log(gcd)/n trend (moving average, window=20) ===")
window = 20
for start in range(0, NMAX - window + 1, window):
    avg = sum(gcd_data[start + i][2] for i in range(window)) / window
    print(f"  n={start+1:3d}-{start+window:3d}: {avg:.4f}")

# Larger windows for smoothing
print("\n=== log(gcd)/n trend (window=50) ===")
window = 50
for start in range(0, NMAX - window + 1, window):
    avg = sum(gcd_data[start + i][2] for i in range(window)) / window
    print(f"  n={start+1:3d}-{start+window:3d}: {avg:.4f}")

# 2. Bad prime density
print("\n=== Bad prime density ===")
print(f"{'n':>4s} {'#primes':>8s} {'#bad':>6s} {'ratio':>8s} {'bad primes (large>n/2)':>40s}")
for n in [10, 20, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500]:
    if n > NMAX:
        break
    g = gcd_data[n-1][1]
    primes_n = [p for p in all_primes if p <= n]
    bad = [p for p in primes_n if g % p == 0]
    large_bad = [p for p in bad if p > n // 2]
    print(f"{n:4d} {len(primes_n):8d} {len(bad):6d} {len(bad)/max(1,len(primes_n)):8.3f} {str(large_bad):>40s}")

# 3. Large prime analysis
print("\n=== Large primes dividing b_n ===")
print("For each n, count primes p in (n/2, n] that divide b_n")
for n in range(50, NMAX + 1, 50):
    primes_half = [p for p in all_primes if n // 2 < p <= n]
    divides = [p for p in primes_half if b[n] % p == 0]
    print(f"  n={n:3d}: {len(divides)}/{len(primes_half)} primes divide b_n ({100*len(divides)/max(1,len(primes_half)):.1f}%) = {divides[:8]}")

# 4. Fit log(gcd)/n to various decay models
print("\n=== Decay model fitting ===")
# Test: log(gcd)/n ~ C/n^alpha
# Take log: log(log(gcd)/n) ~ log(C) - alpha * log(n)
import statistics
xs = []
ys = []
for n, g, r in gcd_data:
    if n >= 20 and r > 0:
        xs.append(log(n))
        ys.append(log(r))

if len(xs) > 10:
    # Linear regression: y = a + b*x
    n_pts = len(xs)
    x_mean = sum(xs) / n_pts
    y_mean = sum(ys) / n_pts
    ss_xy = sum((xs[i] - x_mean) * (ys[i] - y_mean) for i in range(n_pts))
    ss_xx = sum((xs[i] - x_mean) ** 2 for i in range(n_pts))
    slope = ss_xy / ss_xx
    intercept = y_mean - slope * x_mean
    print(f"  log(log(gcd)/n) ~ {intercept:.3f} + ({slope:.3f}) * log(n)")
    print(f"  => log(gcd)/n ~ C * n^({slope:.3f}), C = e^{intercept:.3f} = {2.718281828**intercept:.4f}")
    print(f"  => log(gcd)/n decays as n^{slope:.3f}")
    if slope < -0.5:
        print(f"  This is consistent with log(gcd)/n -> 0 (and hence gcd = e^{{o(n)}})")
    else:
        print(f"  Decay is slow; need more data or different model")

# 5. p-adic structure: for each prime p, is v_p(gcd) periodic in n?
print("\n=== p-adic periodicity check ===")
for p in [2, 3, 5, 7, 11, 13, 17, 19, 23]:
    vals = []
    for n in range(1, min(NMAX + 1, 201)):
        g = gcd_data[n-1][1]
        vals.append(vp(g, p))

    # Check if v_p is eventually periodic
    max_v = max(vals)
    nonzero = sum(1 for v in vals if v > 0)
    # Check period p-1
    period_match = 0
    period_total = 0
    if len(vals) > 2 * p:
        for i in range(p, len(vals)):
            period_total += 1
            if vals[i] == vals[i - p]:
                period_match += 1
    pct = period_match / max(1, period_total) * 100
    print(f"  p={p:2d}: max v_p={max_v}, nonzero {nonzero}/{len(vals)}, period-{p} match: {pct:.0f}%")

# 6. The "denominator budget" analysis
# For each prime p in (sqrt(n), n], v_p(d_n) = 3.
# v_p(gcd) = v_p(d_n) + min(v_p(a_n), v_p(b_n)) where v_p(d_n) = 3 for these primes.
# "Bad" means v_p(a_n) > -3 for this prime, i.e., a_n doesn't use full denominator.
print("\n=== Denominator budget analysis (primes in (sqrt(n), n]) ===")
for n in [100, 200, 300, 400, 500]:
    if n > NMAX:
        break
    dn3 = lcm_power(n) ** 3
    da = int(dn3 * a[n])
    db = dn3 * b[n]
    g = gcd(abs(da), abs(db)) if da != 0 else abs(db)

    sqrt_n = isqrt(n)
    medium_primes = [p for p in all_primes if sqrt_n < p <= n]
    full_budget = 0  # primes where a_n uses full denominator
    partial = 0  # "bad" primes
    for p in medium_primes:
        vp_g = vp(g, p)
        if vp_g == 0:
            full_budget += 1
        else:
            partial += 1
    print(f"  n={n:3d}: {len(medium_primes)} primes in (sqrt(n),n], full budget: {full_budget}, partial ('bad'): {partial} ({100*partial/max(1,len(medium_primes)):.1f}%)")

print(f"\nTotal time: {time.time()-t0:.1f}s")
