#!/usr/bin/env python3
"""Count Z(p) = #{j in [0,p) : b_j ≡ 0 (mod p)} for Apéry numbers.

If Z(p) = O(1) or O(sqrt(p)), the proof of gcd = e^{o(n)} follows.
"""
from math import isqrt

def apery_b_modp(p, maxj=None):
    """Compute b_0,...,b_{maxj} mod p using the recurrence mod p."""
    if maxj is None:
        maxj = p - 1
    b = [0] * (maxj + 1)
    b[0] = 1
    if maxj >= 1:
        b[1] = 5 % p
    for n in range(1, maxj):
        coeff = (34*n**3 + 51*n**2 + 27*n + 5) % p
        # Need to invert (n+1)^3 mod p
        den = pow(n + 1, 3, p)
        if den == 0:
            # p | (n+1), need to handle carefully
            # Use the formula: (n+1)^3 b_{n+1} = coeff*b_n - n^3*b_{n-1}
            # If p | (n+1), then LHS ≡ 0 mod p, so coeff*b_n ≡ n^3*b_{n-1} mod p
            # But we can't determine b_{n+1} mod p from this alone!
            # Actually (n+1)^3 b_{n+1} is an integer, so b_{n+1} mod p
            # requires computing the full integer and reducing.
            # For simplicity, fall back to exact integer computation for this n.
            # Actually, we can resolve this: if p | (n+1), then b_{n+1} mod p
            # requires knowing the p-adic valuation of the numerator.
            # Let's just mark as unknown and skip.
            b[n+1] = -1  # sentinel
            continue
        if b[n] == -1 or (n >= 1 and b[n-1] == -1):
            b[n+1] = -1
            continue
        den_inv = pow(den, p - 2, p)
        num = (coeff * b[n] - pow(n, 3, p) * b[n-1]) % p
        b[n+1] = (num * den_inv) % p
    return b

def apery_b_exact(N):
    """Compute exact Apéry numbers b_0,...,b_N."""
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

# Compute exact b_n for reference
print("Computing exact Apéry numbers to verify...", flush=True)
b_exact = apery_b_exact(1000)

primes = sieve_primes(1000)
print(f"\n{'p':>5s} {'Z(p)':>6s} {'Z/p':>8s} {'Z/sqrt(p)':>10s} {'zeros (first few)':>30s}")
print("-" * 65)

zp_data = []
for p in primes:
    if p < 5:
        continue
    # Count zeros of b_j mod p for j = 0,...,p-1 using exact values
    zeros = []
    for j in range(p):
        if b_exact[j] % p == 0:
            zeros.append(j)
    zp = len(zeros)
    zp_data.append((p, zp))

    if p <= 100 or p in [101, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199,
                          211, 223, 227, 229, 233, 239, 241, 251, 271, 277, 281, 283, 293,
                          307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397,
                          401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499,
                          503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601,
                          607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691,
                          701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797,
                          809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887,
                          907, 911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991, 997]:
        zero_str = str(zeros[:8])
        if len(zeros) > 8:
            zero_str += "..."
        sqrt_p = p**0.5
        print(f"{p:5d} {zp:6d} {zp/p:8.4f} {zp/sqrt_p:10.4f} {zero_str:>30s}")

# Summary statistics
print("\n=== Z(p) growth rate analysis ===")
import statistics

# Check if Z(p) is O(1), O(log p), O(sqrt(p)), or O(p)
from math import log, sqrt

# Fit Z(p) ~ C * p^alpha
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
    print(f"Power law fit: Z(p) ~ {C:.3f} * p^{alpha:.3f}")
    if alpha < 0.1:
        print(f"  Z(p) is essentially O(1) — EXCELLENT for the proof!")
    elif alpha < 0.55:
        print(f"  Z(p) is O(p^{alpha:.2f}) ≈ O(sqrt(p)) — sufficient for the proof")
    elif alpha < 0.9:
        print(f"  Z(p) is O(p^{alpha:.2f}) — might be borderline")
    else:
        print(f"  Z(p) is O(p^{alpha:.2f}) ≈ O(p) — bad for the proof")

# Average Z(p) in ranges
print("\nAverage Z(p) by range:")
ranges = [(5, 50), (50, 100), (100, 200), (200, 300), (300, 500), (500, 700), (700, 1000)]
for lo, hi in ranges:
    subset = [z for p, z in zp_data if lo <= p < hi]
    if subset:
        print(f"  p ∈ [{lo}, {hi}): avg Z(p) = {sum(subset)/len(subset):.2f}, max = {max(subset)}, count = {len(subset)}")

# Z(p) = 0 cases
zero_cases = [p for p, z in zp_data if z == 0]
print(f"\nPrimes with Z(p) = 0 (b_j never ≡ 0 mod p for j < p): {len(zero_cases)}/{len(zp_data)}")
if zero_cases[:20]:
    print(f"  First few: {zero_cases[:20]}")
