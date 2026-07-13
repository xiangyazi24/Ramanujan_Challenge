#!/usr/bin/env python3
"""High-precision verification of Problem 2.4.
Run on uisai2 for M=2000+ terms.
The double sum converges slowly — needs many terms."""
from mpmath import mp, mpf, binomial, polylog, log, zeta, inf
import sys
import time

# Use 100 digits (increase for more precision)
mp.dps = 100

def verify_2_4(M_max=2000):
    """Verify Problem 2.4 to high precision."""
    print(f"=== Problem 2.4: High-precision verification (M={M_max}) ===")
    print(f"Using {mp.dps} decimal digits")

    # RHS
    rhs = (20*polylog(4, mpf(1)/2) + mpf(5)/6 * log(2)**4 + 10*zeta(2)
           - mpf(65)/9 * zeta(2)**2 - log(2)**2 * (12 + 5*zeta(2))
           + mpf(1)/2 * zeta(3) + log(2) * (mpf(35)/2 * zeta(3) - 16))
    print(f"RHS = {rhs}")

    # Precompute harmonic numbers
    H = [mpf(0)] * (M_max + 1)
    for k in range(1, M_max + 1):
        H[k] = H[k-1] + mpf(1)/k

    # LHS: double sum
    lhs = mpf(0)
    start = time.time()
    for m in range(0, M_max + 1):
        inner = mpf(0)
        binom_2m_m = binomial(2*m, m)
        for k in range(0, m + 1):
            binom_mk = binomial(m, k)
            inner += binom_mk**2 * H[k]**2
        lhs += inner / ((m + 1)**2 * binom_2m_m)

        if m % 200 == 0 and m > 0:
            elapsed = time.time() - start
            diff = lhs - rhs
            print(f"  m={m}: diff = {diff:.6e}, elapsed = {elapsed:.1f}s")

    elapsed = time.time() - start
    print(f"\nLHS (M={M_max})  = {lhs}")
    print(f"RHS             = {rhs}")
    print(f"Difference      = {lhs - rhs}")
    print(f"Total time: {elapsed:.1f}s")

if __name__ == "__main__":
    M = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
    verify_2_4(M)
