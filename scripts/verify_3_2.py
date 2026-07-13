#!/usr/bin/env python3
"""Problem 3.2: Computational verification of gcd(d_n a_n, d_n b_n) = e^{o(n)}.
Apery sequences for zeta(3) irrationality."""
from math import gcd, log, lcm
from functools import reduce

def compute_apery(N):
    """Compute Apery sequences a_n, b_n using exact rational arithmetic."""
    from fractions import Fraction as F
    # (n+1)^3 u_{n+1} - (34n^3+51n^2+27n+5) u_n + n^3 u_{n-1} = 0
    a = [F(0), F(6)]
    b = [F(1), F(5)]
    for n in range(1, N):
        coeff_n = 34*n**3 + 51*n**2 + 27*n + 5
        coeff_prev = n**3
        coeff_next = (n+1)**3
        a_next = F(coeff_n * a[-1] - coeff_prev * a[-2], coeff_next)
        b_next = F(coeff_n * b[-1] - coeff_prev * b[-2], coeff_next)
        a.append(a_next)
        b.append(b_next)
    return a, b

def compute_lcm_cubed(n):
    """Compute d_n = lcm(1,...,n)^3."""
    if n == 0:
        return 1
    L = 1
    for i in range(1, n+1):
        L = lcm(L, i)
    return L**3

def main():
    N = 200
    print(f"=== Problem 3.2: gcd(d_n*a_n, d_n*b_n) for n=1,...,{N} ===")

    a, b = compute_apery(N)

    # Verify first few values
    print(f"a: {a[:8]}")
    print(f"b: {b[:8]}")

    print(f"\n{'n':>4} {'log(gcd)/n':>12} {'gcd factors':>30}")
    print("-" * 50)

    max_ratio = 0
    for n in range(1, min(N+1, 101)):
        d_n = compute_lcm_cubed(n)
        da = d_n * a[n]
        db = d_n * b[n]

        # Check they're integers (they should be by Apery's theorem)
        assert da.denominator == 1, f"d_n*a_n not integer at n={n}"
        assert db.denominator == 1, f"d_n*b_n not integer at n={n}"
        da = abs(int(da))
        db = abs(int(db))
        g = gcd(da, db)

        ratio = log(g) / n if n > 0 and g > 0 else 0
        if ratio > max_ratio:
            max_ratio = ratio

        if n <= 20 or n % 10 == 0:
            # Factor g for small values
            if g < 10**15:
                print(f"{n:4d} {ratio:12.6f} {g}")
            else:
                print(f"{n:4d} {ratio:12.6f} (large)")

    print(f"\nMax log(gcd)/n over n=1..{min(N, 100)}: {max_ratio:.6f}")
    print("For e^{o(n)}: this ratio should → 0 as n → ∞")

    # Detailed prime analysis for small primes
    print(f"\n=== p-adic valuations of gcd ===")
    primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]
    print(f"{'n':>4}", end="")
    for p in primes:
        print(f" v_{p:>2}", end="")
    print()

    for n in range(1, 51):
        d_n = compute_lcm_cubed(n)
        da = abs(int(d_n * a[n]))
        db = abs(int(d_n * b[n]))
        g = gcd(da, db)

        print(f"{n:4d}", end="")
        for p in primes:
            v = 0
            temp = g
            while temp > 0 and temp % p == 0:
                v += 1
                temp //= p
            print(f" {v:4d}", end="")
        print()

if __name__ == "__main__":
    main()
