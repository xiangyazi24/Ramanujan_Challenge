#!/usr/bin/env python3
"""
Search for the companion sequence p_n's structure.

Given: p₀ = -612218384750, p₁ = -9525021973931919/18100, p₂ = -29561828382772029/65380
And: p_n satisfies the same P2.7 recurrence as q_n.

Strategy:
1. Compute p_n to high precision using the recurrence
2. Compute the "reduced" companion r_n = p_n/q_n · denom - ζ(2)-related terms
3. Look for patterns in p_n that involve harmonic numbers, polylogarithms, etc.
4. Check: does p_n = Σ_k T₀(n,k) · [harmonic sum formula]?
"""
import mpmath
mpmath.mp.dps = 200

def binom(n, k):
    if k < 0 or k > n: return 0
    r = 1
    for i in range(k): r = r * (n - i) // (i + 1)
    return r

N = 30

# P2.7 recurrence
def A_p(n): return 1024*(2*n+5)**4*(2*n+7)**3*(2*n+9)**3*(946*n**2+6407*n+10860)
def B_p(n): return 128*(2*n+7)**3*(2*n+9)**3*(104060*n**6+1745370*n**5+12145238*n**4+44886481*n**3+92943995*n**2+102256019*n+46709052)
def C_p(n): return 16*(n+3)**4*(2*n+9)**3*(3784*n**5+57792*n**4+351019*n**3+1059230*n**2+1587211*n+944620)
def D_p(n): return (n+3)**4*(n+4)**6*(946*n**2+4515*n+5399)

# Compute q_n and p_n
q = [mpmath.mpf(0)] * N
q[0] = mpmath.mpf('-215040420000')
q[1] = mpmath.mpf('-167282265043404') / mpmath.mpf('905')
q[2] = mpmath.mpf('-964185327658080') / mpmath.mpf('6071')
for n in range(2, N-1):
    q[n+1] = (mpmath.mpf(B_p(n))/A_p(n)*q[n]
              - mpmath.mpf(C_p(n-1))/A_p(n-1)*q[n-1]
              + mpmath.mpf(D_p(n-2))/A_p(n-2)*q[n-2])

p = [mpmath.mpf(0)] * N
p[0] = mpmath.mpf('-612218384750')
p[1] = mpmath.mpf('-9525021973931919') / mpmath.mpf('18100')
p[2] = mpmath.mpf('-29561828382772029') / mpmath.mpf('65380')
for n in range(2, N-1):
    p[n+1] = (mpmath.mpf(B_p(n))/A_p(n)*p[n]
              - mpmath.mpf(C_p(n-1))/A_p(n-1)*p[n-1]
              + mpmath.mpf(D_p(n-2))/A_p(n-2)*p[n-2])

L = mpmath.zeta(2) + mpmath.zeta(3)
print(f"L = ζ(2)+ζ(3) = {float(L):.30f}")

# Verify convergence
print("\nConvergence p_n/q_n → ζ(2)+ζ(3):")
for n in range(0, min(20, N)):
    if abs(q[n]) > 1e-50:
        ratio = p[n] / q[n]
        digits = -mpmath.log10(abs(ratio - L)) if abs(ratio - L) > 0 else 999
        print(f"  n={n:2d}: p/q - L ≈ {float(ratio-L):+.6e}, digits={float(digits):.1f}")

# Error sequence
e = [p[n] - L * q[n] for n in range(N)]
print("\nError decay:")
for n in range(0, min(15, N)):
    if abs(e[n]) > 1e-200:
        print(f"  n={n:2d}: |e_n| = {float(abs(e[n])):.6e}")

# AESZ #209
a = [mpmath.mpf(0)] * N
for n in range(N):
    a[n] = sum(mpmath.mpf(binom(n,k)**2 * binom(n+k,n) * binom(n+2*k,n)) for k in range(n+1))

# Key ratio: q_n / a_n tells us the "gauge" between AESZ and P2.7
print("\n" + "="*70)
print("Structural analysis: q_n / a_n")
print("="*70)
for n in range(min(15, N)):
    if abs(a[n]) > 1e-50:
        r = q[n] / a[n]
        print(f"  n={n:2d}: q_n/a_n = {float(r):+.10e}")

# The ratio q_n/a_n should encode the gauge transformation
# Let's factor out 64^{-n} since q_n ~ (μ₀/64)^n while a_n ~ μ₀^n
print("\n  q_n / (a_n * 64^{-n}):")
for n in range(min(15, N)):
    if abs(a[n]) > 1e-50:
        r = q[n] * mpmath.mpf(64)**n / a[n]
        print(f"  n={n:2d}: {float(r):+.10e}")

# Now look at the companion p_n / a_n structure
print("\n" + "="*70)
print("p_n analysis")
print("="*70)

# Try: p_n = q_n · [something involving harmonic numbers of the first N terms of a]
# Look at d_n = L·q_n - p_n = -e_n  (which decays)
# and p_n/q_n ≈ L
# So p_n ≈ L·q_n + error

# What is the "denominator" structure?
# In Apéry, the denominators of p_n involve lcm(1,...,n)^s
# Let's look at the denominator of p_n in reduced form

from fractions import Fraction

# Using exact arithmetic for small n
p_exact = [None] * 5
p_exact[0] = Fraction(-612218384750, 1)
p_exact[1] = Fraction(-9525021973931919, 18100)
p_exact[2] = Fraction(-29561828382772029, 65380)

# For exact p_n computation with fractions:
def A_ex(n): return 1024*(2*n+5)**4*(2*n+7)**3*(2*n+9)**3*(946*n**2+6407*n+10860)
def B_ex(n): return 128*(2*n+7)**3*(2*n+9)**3*(104060*n**6+1745370*n**5+12145238*n**4+44886481*n**3+92943995*n**2+102256019*n+46709052)
def C_ex(n): return 16*(n+3)**4*(2*n+9)**3*(3784*n**5+57792*n**4+351019*n**3+1059230*n**2+1587211*n+944620)
def D_ex(n): return (n+3)**4*(n+4)**6*(946*n**2+4515*n+5399)

q_exact = [None] * 10
q_exact[0] = Fraction(-215040420000, 1)
q_exact[1] = Fraction(-167282265043404, 905)
q_exact[2] = Fraction(-964185327658080, 6071)

for n in range(2, 9):
    q_exact[n+1] = (Fraction(B_ex(n), A_ex(n)) * q_exact[n]
                    - Fraction(C_ex(n-1), A_ex(n-1)) * q_exact[n-1]
                    + Fraction(D_ex(n-2), A_ex(n-2)) * q_exact[n-2])

p_exact = [None] * 10
p_exact[0] = Fraction(-612218384750, 1)
p_exact[1] = Fraction(-9525021973931919, 18100)
p_exact[2] = Fraction(-29561828382772029, 65380)
for n in range(2, 9):
    p_exact[n+1] = (Fraction(B_ex(n), A_ex(n)) * p_exact[n]
                    - Fraction(C_ex(n-1), A_ex(n-1)) * p_exact[n-1]
                    + Fraction(D_ex(n-2), A_ex(n-2)) * p_exact[n-2])

print("Exact denominators:")
import math
for n in range(min(8, len(q_exact))):
    if q_exact[n] is not None:
        qd = q_exact[n].denominator
        pd = p_exact[n].denominator if p_exact[n] is not None else 0
        lcm_n = 1
        for k in range(1, n+1):
            lcm_n = lcm_n * k // math.gcd(lcm_n, k)
        print(f"  n={n}: q_denom={qd}, p_denom={pd}, lcm(1..n)={lcm_n}, lcm^2={lcm_n**2}, lcm^3={lcm_n**3}")

# Look for the denominators in terms of lcm
print("\nDenominator analysis:")
for n in range(min(8, len(q_exact))):
    if q_exact[n] is not None:
        qd = q_exact[n].denominator
        pd = p_exact[n].denominator if p_exact[n] is not None else 0
        # Factor denominators
        # Compute d_n = lcm(1,...,n+2)^3 or similar
        d1 = 1  # lcm(1,...,n+2)
        for k in range(1, n+3):
            d1 = d1 * k // math.gcd(d1, k)
        d2 = 1  # lcm(1,...,n+3)
        for k in range(1, n+4):
            d2 = d2 * k // math.gcd(d2, k)
        d3 = 1  # lcm(1,...,2n+3)
        for k in range(1, 2*n+4):
            d3 = d3 * k // math.gcd(d3, k)

        print(f"  n={n}: qd={qd}")
        print(f"    lcm(1..n+2)={d1}, d1 | qd? {qd % d1 == 0}")
        print(f"    qd/d1={qd//d1 if qd%d1==0 else 'N/A'}")
        if pd > 0:
            print(f"    pd={pd}")
            print(f"    lcm(1..n+2)={d1}, d1 | pd? {pd % d1 == 0}")

# The key check: does d_n³ · q_n ∈ Z for d_n = lcm(1,...,n)?
print("\nIntegrality check: d_n^s · q_n ∈ Z?")
for n in range(min(8, len(q_exact))):
    if q_exact[n] is not None:
        for s in range(1, 6):
            for offset in range(5):
                dd = 1
                for k in range(1, n+1+offset):
                    dd = dd * k // math.gcd(dd, k)
                val = q_exact[n] * dd**s
                if val.denominator == 1:
                    print(f"  n={n}: lcm(1..{n+offset})^{s} · q_n ∈ Z ✓ (= {val.numerator})")
                    break
            else:
                continue
            break
