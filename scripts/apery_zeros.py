#!/usr/bin/env python3
"""Problem 3.2: Compute Apéry number zeros mod primes.

For each prime p, compute ν_p = |{0 ≤ r < p : b_r ≡ 0 (mod p)}|.

Also verify:
1. Reflection: b_{p-1-r} ≡ b_r (mod p)
2. No adjacent zeros
3. ν_p ≤ (p-1)/2
"""

def apery_seq_mod(p, N=None):
    """Compute b_0, b_1, ..., b_{N-1} mod p.
    b_n = sum_{k=0}^n C(n,k)^2 C(n+k,k)^2.
    Uses the recurrence: (n+1)^3 b_{n+1} = (34n^3+51n^2+27n+5) b_n - n^3 b_{n-1}
    """
    if N is None:
        N = p
    b = [0] * N
    b[0] = 1
    if N > 1:
        b[1] = 5
    for n in range(1, N-1):
        P = (34*n**3 + 51*n**2 + 27*n + 5) % p
        rhs = (P * b[n] - pow(n, 3, p) * b[n-1]) % p
        inv = pow((n+1)**3 % p, -1, p) if (n+1) % p != 0 else None
        if inv is not None:
            b[n+1] = (rhs * inv) % p
        else:
            b[n+1] = 0  # p | (n+1)^3, need special handling
    return b

def count_zeros(b, p):
    """Count zeros of b mod p."""
    return sum(1 for x in b[:p] if x % p == 0)

def check_reflection(b, p):
    """Check b_{p-1-r} ≡ b_r (mod p) for all 0 ≤ r < p."""
    for r in range(p):
        if b[r] % p != b[p-1-r] % p:
            return False, r
    return True, -1

def check_no_adjacent_zeros(b, p):
    """Check that no two consecutive b_r are both 0 mod p."""
    for r in range(p-1):
        if b[r] % p == 0 and b[r+1] % p == 0:
            return False, r
    return True, -1

# ---- Main computation ----
print("=== Apéry number zeros mod primes ===")
print(f"{'p':>5} {'ν_p':>5} {'ν_p/p':>8} {'reflection':>10} {'no adj':>7} {'zeros':>40}")
print("-" * 85)

nu_data = []
for p in range(5, 200):
    # Check if p is prime
    if p < 4:
        continue
    is_prime = True
    for d in range(2, int(p**0.5)+1):
        if p % d == 0:
            is_prime = False
            break
    if not is_prime:
        continue

    b = apery_seq_mod(p)
    nu = count_zeros(b, p)
    refl_ok, refl_fail = check_reflection(b, p)
    adj_ok, adj_fail = check_no_adjacent_zeros(b, p)

    zeros = [r for r in range(p) if b[r] % p == 0]
    zeros_str = str(zeros[:10]) + ("..." if len(zeros) > 10 else "")

    print(f"{p:>5} {nu:>5} {nu/p:>8.4f} {'OK' if refl_ok else f'FAIL@{refl_fail}':>10} {'OK' if adj_ok else f'FAIL@{adj_fail}':>7} {zeros_str}")

    nu_data.append((p, nu))

# ---- Statistics ----
print("\n=== Statistics ===")
print(f"{'Range':>10} {'avg ν_p':>8} {'max ν_p':>8} {'avg ν_p/p':>10}")

ranges = [(5, 50), (50, 100), (100, 200)]
for lo, hi in ranges:
    subset = [(p, nu) for p, nu in nu_data if lo <= p < hi]
    if subset:
        avg_nu = sum(nu for _, nu in subset) / len(subset)
        max_nu = max(nu for _, nu in subset)
        avg_ratio = sum(nu/p for p, nu in subset) / len(subset)
        print(f"[{lo},{hi}){'>':>5} {avg_nu:>8.2f} {max_nu:>8} {avg_ratio:>10.4f}")

# ---- Larger primes for density trend ----
print("\n=== Larger primes ===")
for p in [251, 509, 1009, 2003, 5003]:
    b = apery_seq_mod(p)
    nu = count_zeros(b, p)
    print(f"p={p}: ν_p={nu}, ν_p/p={nu/p:.4f}, ν_p/√p={nu/p**0.5:.4f}")
