#!/usr/bin/env python3
"""Verify the denominator connection lemma:
For p in (n/2, n] with p >= 5:
  sum_{k=p}^n c_k ≡ 4 b_{n-p} (mod p)
where c_k = C(n,k)^2 C(n+k,k)^2.

Also verify that v_p(G_n) = 0 iff p does not divide b_{n-p}.
"""
from math import comb, gcd, isqrt
from fractions import Fraction

def apery_b(N):
    b = [0] * (N + 1)
    b[0] = 1
    if N >= 1:
        b[1] = 5
    for n in range(1, N):
        coeff = 34*n**3 + 51*n**2 + 27*n + 5
        num = coeff * b[n] - n**3 * b[n-1]
        den = (n+1)**3
        assert num % den == 0
        b[n+1] = num // den
    return b

def apery_a(N):
    a = [Fraction(0)] * (N + 1)
    a[0] = Fraction(0)
    a[1] = Fraction(6)
    for n in range(1, N):
        coeff = 34*n**3 + 51*n**2 + 27*n + 5
        a[n+1] = (Fraction(coeff) * a[n] - Fraction(n**3) * a[n-1]) / Fraction((n+1)**3)
    return a

def lcm_cube(n):
    primes = sieve_primes(n)
    result = 1
    for p in primes:
        pk = p
        while pk <= n:
            pk *= p
        result *= (pk // p)
    return result ** 3

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

def vp(n, p):
    if n == 0:
        return float('inf')
    v = 0
    n = abs(n)
    while n % p == 0:
        n //= p
        v += 1
    return v

NMAX = 300
print(f"Computing sequences to n={NMAX}...")
b = apery_b(NMAX)
a = apery_a(NMAX)
primes = sieve_primes(NMAX)

# Test 1: sum_{k=p}^n c_k ≡ 4 b_{n-p} (mod p)
print("\n=== Denominator connection: sum c_k ≡ 4 b_r (mod p) ===")
matches = 0
total = 0
failures = []

for n in range(10, NMAX + 1, 10):
    for p in primes:
        if p < 5 or p <= n // 2 or p > n:
            continue
        r = n - p
        # Compute sum_{k=p}^n c_k mod p
        tail_sum = 0
        for k in range(p, n + 1):
            ck = pow(comb(n, k), 2, p) * pow(comb(n + k, k), 2, p) % p
            tail_sum = (tail_sum + ck) % p
        expected = (4 * b[r]) % p
        total += 1
        if tail_sum == expected:
            matches += 1
        else:
            failures.append((n, p, r, tail_sum, expected))

print(f"  {matches}/{total} matches ({100*matches/max(1,total):.1f}%)")
if failures:
    for n, p, r, got, exp in failures[:5]:
        print(f"  FAIL: n={n}, p={p}, r={r}: sum={got}, 4b_r={exp}")

# Test 2: v_p(G_n) = 0 iff p ∤ b_{n-p}
print("\n=== v_p(G_n) = 0 iff p ∤ b_{n-p} ===")
matches2 = 0
total2 = 0
for n in range(10, min(NMAX + 1, 201)):
    dn = lcm_cube(n)
    da = int(dn * a[n])
    db = dn * b[n]
    if da == 0:
        continue
    g = gcd(abs(da), abs(db))
    for p in primes:
        if p < 5 or p <= n // 2 or p > n:
            continue
        r = n - p
        vp_g = vp(g, p)
        br_div_p = (b[r] % p == 0)
        total2 += 1
        # Lemma: v_p(G_n) = 0 iff p ∤ b_r
        if (vp_g == 0) == (not br_div_p):
            matches2 += 1
        else:
            print(f"  MISMATCH: n={n}, p={p}, r={r}: v_p(G)={vp_g}, p|b_r={br_div_p}")

print(f"  {matches2}/{total2} matches ({100*matches2/max(1,total2):.1f}%)")

# Test 3: b_n ≡ 5 b_r (mod p) for p in (n/2, n]
print("\n=== b_n ≡ 5 b_r (mod p) for p in (n/2, n] ===")
m3 = 0
t3 = 0
for n in range(10, min(NMAX + 1, 501)):
    for p in primes:
        if p < 5 or p <= n // 2 or p > n:
            continue
        r = n - p
        t3 += 1
        if b[n] % p == (5 * b[r]) % p:
            m3 += 1
        else:
            print(f"  FAIL: n={n}, p={p}, r={r}: b_n%p={b[n]%p}, 5*b_r%p={(5*b[r])%p}")

print(f"  {m3}/{t3} matches ({100*m3/max(1,t3):.1f}%)")
print("\nDone.")
