#!/usr/bin/env python3
"""Verify Lucas-type congruences for Apéry numbers b_n = sum C(n,k)^2 C(n+k,k)^2.

Key question: does b_n ≡ prod b_{n_i} (mod p) where n = sum n_i p^i?
If yes, this is the cornerstone for proving gcd = e^{o(n)}.
"""
from math import comb, gcd
import sys

def apery_b(N):
    """Compute b_0,...,b_N via recurrence."""
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

def base_p_digits(n, p):
    """Return base-p digits of n, least significant first."""
    if n == 0:
        return [0]
    digits = []
    while n > 0:
        digits.append(n % p)
        n //= p
    return digits

NMAX = 500
print(f"Computing Apéry numbers b_0,...,b_{NMAX}...")
b = apery_b(NMAX)

# === Test 1: Multiplicative Lucas congruence ===
# b_n ≡ prod b_{n_i} (mod p) where n = sum n_i p^i?
print("\n=== Test: b_n ≡ prod b_{n_i} (mod p) ===")
for p in [5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]:
    matches = 0
    total = 0
    failures = []
    for n in range(p, NMAX + 1):
        digits = base_p_digits(n, p)
        if any(d >= len(b) for d in digits):
            continue
        prod_b = 1
        for d in digits:
            prod_b *= b[d]
        total += 1
        if (b[n] - prod_b) % p == 0:
            matches += 1
        else:
            if len(failures) < 3:
                failures.append((n, b[n] % p, prod_b % p, digits))

    pct = 100 * matches / total if total > 0 else 0
    print(f"  p={p:2d}: {matches}/{total} matches ({pct:.1f}%)")
    if failures:
        for n, bn_mod, prod_mod, digits in failures[:2]:
            print(f"    FAIL n={n}: b_n≡{bn_mod}, prod≡{prod_mod}, digits={digits}")

# === Test 2: Two-digit case b_{mp+j} ≡ b_m * b_j (mod p) ===
# This is the simplest version; maybe there's a sign correction.
print("\n=== Test: b_{mp+j} ≡ ± b_m * b_j (mod p) ===")
for p in [5, 7, 11, 13]:
    signs_plus = 0
    signs_minus = 0
    neither = 0
    total = 0
    for m in range(0, min(50, NMAX // p)):
        for j in range(0, p):
            n = m * p + j
            if n > NMAX:
                break
            total += 1
            val = b[n] % p
            prod_val = (b[m] * b[j]) % p
            neg_prod = (-b[m] * b[j]) % p
            if val == prod_val:
                signs_plus += 1
            elif val == neg_prod:
                signs_minus += 1
            else:
                neither += 1

    print(f"  p={p:2d}: +matches={signs_plus}, -matches={signs_minus}, neither={neither} (total={total})")

# === Test 3: Maybe b_{mp+j} ≡ (-1)^m * b_m * b_j (mod p)? ===
print("\n=== Test: b_{mp+j} ≡ (-1)^? * b_m * b_j (mod p) ===")
for p in [5, 7, 11, 13]:
    print(f"  p={p}:")
    for j in range(min(p, 5)):
        for m in range(min(10, NMAX // p)):
            n = m * p + j
            if n > NMAX:
                break
            val = b[n] % p
            prod_val = (b[m] * b[j]) % p
            ratio = None
            if prod_val != 0 and val != 0:
                # Find r such that val = r * prod_val (mod p)
                for r in range(p):
                    if (r * prod_val) % p == val:
                        ratio = r
                        break
            elif prod_val == 0 and val == 0:
                ratio = 0
            print(f"    m={m}, j={j}, n={n}: b_n%p={val}, b_m*b_j%p={prod_val}, ratio={ratio}")

# === Test 4: Beukers supercongruence b_{mp} ≡ b_m (mod p^3) ===
print("\n=== Test: b_{mp} ≡ b_m (mod p^3) [Beukers] ===")
for p in [5, 7, 11, 13]:
    for m in range(0, min(20, NMAX // p)):
        n = m * p
        if n > NMAX:
            break
        diff = b[n] - b[m]
        if diff == 0:
            vp = float('inf')
        else:
            vp = 0
            d = abs(diff)
            while d % p == 0:
                d //= p
                vp += 1
        status = "✓" if vp >= 3 else "✗"
        if m <= 5 or vp < 3:
            print(f"  p={p:2d}, m={m:2d}: v_p(b_{n} - b_{m}) = {vp} {status}")

# === Test 5: For n=200 bad primes 139, 181: why does p | b_n? ===
print("\n=== Bad prime analysis: why p | b_n ===")
for n, p in [(200, 139), (200, 181), (300, 191), (300, 227), (500, 271)]:
    if n > NMAX:
        continue
    digits = base_p_digits(n, p)
    bn_mod = b[n] % p
    digit_vals = [(d, b[d] % p) for d in digits]
    print(f"  n={n}, p={p}: digits_base_p={digits}, b_n%p={bn_mod}")
    for d, bv in digit_vals:
        print(f"    digit {d}: b_{d}%{p} = {bv}")
    prod = 1
    for d, bv in digit_vals:
        prod = (prod * bv) % p
    print(f"    product of b_{{digit}} mod p = {prod}")

print("\nDone.")
