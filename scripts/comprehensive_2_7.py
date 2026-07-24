#!/usr/bin/env python3
"""Problem 2.7: Comprehensive analysis.

1. Verify convergence to ζ(2)+ζ(3) with corrected initial conditions
2. Compute the remainder r_n = p_n - (ζ(2)+ζ(3))q_n
3. Analyze the remainder ratio to confirm Poincaré structure
4. Check if the recurrence is self-dual (palindromic)
5. Investigate the denominator structure of q_n (LCD, prime factorization)
"""
from mpmath import mp, mpf, zeta, log10, fabs, power, floor, sqrt

mp.dps = 300

def A(n):
    n = mpf(n)
    return 1024*(2*n+5)**4*(2*n+7)**3*(2*n+9)**3*(946*n**2+6407*n+10860)

def B(n):
    n = mpf(n)
    P6 = 104060*n**6+1745370*n**5+12145238*n**4+44886481*n**3+92943995*n**2+102256019*n+46709052
    return 128*(2*n+7)**3*(2*n+9)**3*P6

def C(n):
    n = mpf(n)
    P5 = 3784*n**5+57792*n**4+351019*n**3+1059230*n**2+1587211*n+944620
    return 16*(n+3)**4*(2*n+9)**3*P5

def D(n):
    n = mpf(n)
    return (n+3)**4*(n+4)**6*(946*n**2+4515*n+5399)

p = [mpf(-612218384750),
     mpf(-9525021973931919)/mpf(18100),
     mpf(-29561828382772029)/mpf(65380)]

q = [mpf(-215040420000),
     mpf(-167282265043404)/mpf(905),
     mpf(-964185327658080)/mpf(6071)]

target = zeta(2) + zeta(3)

N = 80
for n in range(2, N):
    p_next = B(n)/A(n)*p[n] - C(n-1)/A(n-1)*p[n-1] + D(n-2)/A(n-2)*p[n-2]
    q_next = B(n)/A(n)*q[n] - C(n-1)/A(n-1)*q[n-1] + D(n-2)/A(n-2)*q[n-2]
    p.append(p_next)
    q.append(q_next)

print("=== Problem 2.7: Comprehensive Analysis ===\n")
print(f"Target: ζ(2)+ζ(3) = {mp.nstr(target, 50)}")

# 1. Convergence verification
print("\n--- 1. Convergence verification ---")
for n in [5, 10, 20, 30, 50, 70]:
    ratio = p[n]/q[n]
    err = fabs(ratio - target)
    if err > 0:
        digits = -log10(err)
    else:
        digits = mpf('inf')
    print(f"  n={n:3d}: {mp.nstr(digits, 6)} digits match")

# 2. Remainder analysis
print("\n--- 2. Remainder ratio r_n/r_{n-1} ---")
remainders = [p[n] - target*q[n] for n in range(N)]
for n in [3, 5, 10, 15, 20, 30, 40, 50, 60, 70]:
    r = remainders[n]
    r_prev = remainders[n-1]
    if r_prev != 0:
        ratio_r = r / r_prev
        print(f"  n={n:3d}: r_n/r_{{n-1}} = {mp.nstr(ratio_r, 30)}")
    else:
        print(f"  n={n:3d}: r_{{n-1}} = 0 (precision exhausted)")

# The ratio should approach the COMPLEX Poincaré root.
# Since r_n is real but the subdominant roots are complex conjugate,
# the ratio oscillates. Let's check r_n/r_{n-2}
print("\n--- 3. Two-step ratio r_n/r_{n-2} ---")
for n in [4, 6, 10, 20, 30, 40, 50, 60]:
    r = remainders[n]
    r_prev2 = remainders[n-2]
    if r_prev2 != 0:
        ratio_r2 = r / r_prev2
        print(f"  n={n:3d}: r_n/r_{{n-2}} = {mp.nstr(ratio_r2, 20)}")
    else:
        print(f"  n={n:3d}: r_{{n-2}} = 0")

# |r±|² should be the limiting ratio of r_n/r_{n-2}
# |r±|² = 0.001054² ≈ 1.11e-6
print(f"\n  Expected |r±|² ≈ {mp.nstr(mpf('0.001054')**2, 8)}")

# 3. Self-duality check
# A recurrence is self-dual if a_k(n) = ±a_{m-k}(n+shift) for some shift.
# Check: is A_n D_n ∝ A_n D_n under some reflection?
print("\n--- 4. Coefficient symmetry check ---")
for n_test in [0, 1, 2, 5, 10]:
    n = mpf(n_test)
    an = A(n_test)
    dn = D(n_test)
    print(f"  n={n_test}: A_n/D_n = {mp.nstr(an/dn, 15)}")
    print(f"         A_n·D_n = {mp.nstr(an*dn, 15)}")

# 4. Denominator LCD structure
print("\n--- 5. Denominator structure ---")
from math import gcd
# Compute q_n as exact rationals using integer arithmetic
# q_n are given as p/q initially, let's track the LCD

# Actually, let's just check: are the q_n eventually integers when
# multiplied by some LCD?
# q_0 = -215040420000 (integer)
# q_1 = -167282265043404/905 (denominator 905 = 5·181)
# q_2 = -964185327658080/6071

# The LCD growth rate tells us about the "denominator obstruction"
# For Apéry-type proofs, the LCD usually grows like exp(c·n).
print("  q_0 = -215040420000 (integer)")
print("  q_1 = -167282265043404/905")
print("  q_2 = -964185327658080/6071")
print("  905 = 5 × 181")
print("  6071 = 7 × 867 + 2 = ...")
print(f"  6071 factors: ", end="")
n = 6071
for p_fac in range(2, 100):
    while n % p_fac == 0:
        print(f"{p_fac} ", end="")
        n //= p_fac
if n > 1:
    print(f"{n}")
else:
    print()

# 5. Growth rate of |q_n|
print("\n--- 6. Growth rate of |q_n| ---")
for n in [5, 10, 20, 30, 40, 50]:
    ratio = fabs(q[n]/q[n-1])
    print(f"  |q_{n}|/|q_{n-1}| = {mp.nstr(ratio, 15)}")

# Dominant Poincaré root should give asymptotic growth
# |q_n| ~ |r_0|^n · n^σ for some exponent σ
# log|q_n| ~ n·log|r_0| + σ·log(n) + const
# log|r_0| = log(0.8588) = -0.1523
# But q_n/q_{n-1} oscillates... The ratio approaches r_0 = 0.8588.
# Actually for GROWING q_n, the ratio q_{n+1}/q_n should approach
# the LARGEST root. But r_0 = 0.8588 < 1, so |q_n| → 0!
# No — the recurrence coefficients are NOT monic. The ratio
# B/A ~ 0.859 is the Poincaré dominant root. The actual growth
# of q_n includes the polynomial normalizing factors.

print(f"\n  Poincaré dominant root r_0 = 0.8588...")
print(f"  Since r_0 < 1, the NORMALIZED sequence decreases.")
print(f"  But the RAW sequence q_n includes polynomial growth factors.")
print(f"  The 'digits per step' = -log10(|r±|/r_0) = 2.91")
print(f"  Irrationality measure: μ ≤ 1 + log(r_0)/log|r±| = 1.022")

# 6. Check the Poincaré polynomial
print("\n--- 7. Poincaré polynomial verification ---")
# Leading coefficients of A, B, C, D as polynomials in n
# A_n = 1024·(2n+5)^4·(2n+7)^3·(2n+9)^3·(946n^2+...) → leading: 1024·16·8·8·946 n^12
# = 1024 · 2^4·n^4 · 2^3·n^3 · 2^3·n^3 · 946·n^2 = 1024·2^10·946·n^12
A_lead = 1024 * 2**10 * 946
# B_n → 128·2^3·n^3·2^3·n^3·104060·n^6 = 128·2^6·104060·n^12
B_lead = 128 * 2**6 * 104060
# C_n → 16·n^4·2^3·n^3·3784·n^5 = 16·2^3·3784·n^12
C_lead = 16 * 2**3 * 3784
# D_n → n^4·n^6·946·n^2 = 946·n^12
D_lead = 946

print(f"  A_∞ = {A_lead}")
print(f"  B_∞ = {B_lead}")
print(f"  C_∞ = {C_lead}")
print(f"  D_∞ = {D_lead}")
print(f"  Poincaré: {A_lead}c³ - {B_lead}c² + {C_lead}c - {D_lead} = 0")

from mpmath import polyroots
roots = polyroots([mpf(A_lead), -mpf(B_lead), mpf(C_lead), -mpf(D_lead)])
print(f"  Roots:")
for r in roots:
    print(f"    c = {mp.nstr(r, 20)}, |c| = {mp.nstr(abs(r), 15)}")

print(f"\n  Normalized: 4μ³ - 220μ² + 8μ - 1 = 0 (μ = 64c)")
for r in roots:
    mu = 64*r
    print(f"    μ = {mp.nstr(mu, 20)}")

# Verify: 4·55³ - 220·55² + 8·55 - 1
test = 4*55**3 - 220*55**2 + 8*55 - 1
print(f"\n  4·55³ - 220·55² + 8·55 - 1 = {test}")
# If the polynomial factors as 4(μ-55)(μ²+bμ+c)...
