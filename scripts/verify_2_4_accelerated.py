#!/usr/bin/env python3
"""Problem 2.4: Accelerated verification using Euler-Maclaurin + Richardson.
The double sum converges as O(1/m), too slow for direct summation.
Strategy: compute the RECURRENCE for the outer sum's tail and accelerate."""
from mpmath import mp, mpf, binomial, polylog, log, zeta, matrix, lu_solve
from fractions import Fraction
import time

mp.dps = 80

def compute_exact_Am(m_max):
    """Compute A_m = sum_{k=0}^m (m choose k)^2 H_k^2 exactly as Fractions."""
    from math import comb
    H = [Fraction(0)]
    for k in range(1, m_max + 2):
        H.append(H[-1] + Fraction(1, k))

    results = []
    for m in range(m_max + 1):
        Am = Fraction(0)
        for k in range(m + 1):
            Am += Fraction(comb(m, k)**2) * H[k]**2
        results.append(Am)
    return results

def main():
    print("=== Problem 2.4: Accelerated verification ===")
    start = time.time()

    # RHS
    rhs = (20*polylog(4, mpf(1)/2) + mpf(5)/6 * log(2)**4 + 10*zeta(2)
           - mpf(65)/9 * zeta(2)**2 - log(2)**2 * (12 + 5*zeta(2))
           + mpf(1)/2 * zeta(3) + log(2) * (mpf(35)/2 * zeta(3) - 16))
    print(f"RHS = {rhs}")

    # Compute exact A_m for m=0..300
    M_MAX = 300
    print(f"\nComputing exact A_m for m=0..{M_MAX}...")
    Am_exact = compute_exact_Am(M_MAX)
    elapsed = time.time() - start
    print(f"Done in {elapsed:.1f}s")

    # Compute g_m = A_m / ((m+1)^2 * C(2m,m)) and partial sums
    partial_sums = []
    s = mpf(0)
    gm_values = []
    for m in range(M_MAX + 1):
        Am = mpf(Am_exact[m].numerator) / mpf(Am_exact[m].denominator)
        binom_2m = binomial(2*m, m)
        gm = Am / ((m+1)**2 * binom_2m)
        gm_values.append(gm)
        s += gm
        partial_sums.append(s)

    print(f"\nDirect partial sums:")
    for n in [50, 100, 150, 200, 250, 300]:
        diff = partial_sums[n] - rhs
        print(f"  S_{n:3d} - RHS = {mp.nstr(diff, 8)}")

    # Neville-Aitken / Richardson extrapolation using the tail structure
    # g_m ~ C/m as m → ∞ (from 1/((m+1)^2 * C(2m,m)) ~ sqrt(π m) * 4^{-m} * m^{-2})
    # Wait, C(2m,m) ~ 4^m / sqrt(π m), so 1/((m+1)^2 C(2m,m)) ~ sqrt(π m) * 4^{-m} / m^2
    # That's EXPONENTIALLY small! So g_m decays as 4^{-m} * poly(m) * A_m
    # And A_m grows roughly as C(2m,m) * (log m)^2 ~ 4^m * (log m)^2 / sqrt(m)
    # So g_m ~ (log m)^2 / m^{5/2}... that's polynomial decay
    # Actually let me check numerically
    print(f"\nDecay rate g_m:")
    for m in [50, 100, 200, 300]:
        if m <= M_MAX:
            print(f"  g_{m} = {mp.nstr(gm_values[m], 10)}, m^2*g_m = {mp.nstr(m**2*gm_values[m], 10)}")

    # Euler-Maclaurin acceleration: fit g_m ~ c1/m + c2/m^2 + c3/m^3 for large m
    # Then tail sum = c1*H_N + c2*zeta(2,N) + c3*zeta(3,N) + ...
    # Actually, better: use Wynn's epsilon algorithm on partial sums
    print(f"\nWynn epsilon acceleration:")
    # Take the last 30 partial sums
    sums = partial_sums[-30:]
    n = len(sums)
    e = [[mpf(0)] * (n + 1) for _ in range(n + 1)]
    for i in range(n):
        e[i][1] = sums[i]

    for j in range(2, n):
        for i in range(n - j):
            diff = e[i+1][j-1] - e[i][j-1]
            if abs(diff) < mpf(10)**(-mp.dps + 5):
                e[i][j] = e[i+1][j-1]
            else:
                e[i][j] = e[i][j-2] + 1/diff

    # Best estimates are in even columns
    for order in [2, 4, 6, 8, 10, 12]:
        if order < n:
            val = e[0][order]
            diff = val - rhs
            print(f"  order={order:2d}: diff = {mp.nstr(diff, 10)}")

    elapsed = time.time() - start
    print(f"\nTotal time: {elapsed:.1f}s")

if __name__ == "__main__":
    main()
