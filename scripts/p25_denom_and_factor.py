#!/usr/bin/env python3
"""P2.5: Denominator analysis + operator factorization attempt.

1. Factor structure of Q̂_n denominators
2. Right-factor search for the order-3 recurrence
3. Check if Q̂_n has a decomposition matching the Poincaré splitting
"""
from fractions import Fraction
from math import gcd
import math

def M_entries(n):
    n = Fraction(n)
    m11 = (-2*n-5)*(n+3)**2*(136*n**4+1424*n**3+5548*n**2+9551*n+6141)
    m12 = 384*n**6+6384*n**5+44168*n**4+162698*n**3+336377*n**2+369933*n+169011
    m13 = -(480*n**4+4980*n**3+19210*n**2+32690*n+20730)
    m21 = (n+2)**2*(n+3)**2*(4*n+10)*(48*n**3+386*n**2+1017*n+879)
    m22 = (n+2)**2*(-272*n**5-3848*n**4-21732*n**3-61184*n**2-85761*n-47808)
    m23 = (n+2)**2*(320*n**3+2540*n**2+6610*n+5640)
    m31 = (-4*n-10)*(n+2)**2*(n+3)**2*(32*n**4+302*n**3+1037*n**2+1530*n+813)
    m32 = (n+2)**2*(192*n**6+2984*n**5+19116*n**4+64452*n**3+120256*n**2+117279*n+46476)
    m33 = (n+2)**2*(-16*n**5-408*n**4-2912*n**3-8884*n**2-12254*n-6240)
    return [[m11,m12,m13],[m21,m22,m23],[m31,m32,m33]]

def delta_H(n):
    n = Fraction(n)
    return Fraction(-2)*(n+2)**2*(n+3)**2*(2*n+5)*(2*n+7)**2

def MH_at(n):
    M = M_entries(n); d = delta_H(n)
    return [[M[i][j]/d for j in range(3)] for i in range(3)]

NMAX = 50
print("Computing CMF values...", flush=True)
q_row = [Fraction(33750), Fraction(-36000), Fraction(9000)]
p_row = [Fraction(30921), Fraction(-32972), Fraction(8240)]
cmf_q = []; cmf_p = []
# Also store all 3 columns for both Q and P rows
q_full = []; p_full = []
for N in range(NMAX):
    cmf_q.append(q_row[0]); cmf_p.append(p_row[0])
    q_full.append(list(q_row)); p_full.append(list(p_row))
    MH = MH_at(N)
    q_row = [sum(q_row[i]*MH[i][j] for i in range(3)) for j in range(3)]
    p_row = [sum(p_row[i]*MH[i][j] for i in range(3)) for j in range(3)]

# Denominator analysis
print("\n=== Denominator factorization of Q̂_n ===", flush=True)
def factorize_small(n):
    if n == 0: return {}
    if n < 0: n = -n
    factors = {}
    for p in range(2, min(n+1, 10000)):
        while n % p == 0:
            factors[p] = factors.get(p, 0) + 1
            n //= p
        if p*p > n: break
    if n > 1: factors[n] = 1
    return factors

for n in range(20):
    q = cmf_q[n]
    num, den = q.numerator, q.denominator
    if den != 1:
        df = factorize_small(den)
        print(f"  n={n}: den = {den} = {df}")
    else:
        print(f"  n={n}: integer Q̂={num}")

# Check: den(Q̂_n) divides lcm(1, 2, ..., 2n+1)^something?
print("\n=== LCM analysis ===", flush=True)
for n in range(1, 15):
    q = cmf_q[n]
    den = q.denominator
    lcm_val = 1
    for j in range(1, 2*n+2):
        lcm_val = lcm_val * j // gcd(lcm_val, j)
    # What power of lcm divides den?
    r = den
    power = 0
    while r > 1:
        g = gcd(r, lcm_val)
        if g == 1: break
        r //= g
        power += 1
    if r == 1:
        print(f"  n={n}: den | lcm(1..{2*n+1})^{power}")
    else:
        print(f"  n={n}: den has prime factor {factorize_small(r)} outside lcm(1..{2*n+1})")

# Now: operator factorization
# The order-3 recurrence operator L has Poincaré poly (ξ-1)(ξ²-34ξ+1)
# Try to find a right factor L₁ of order 1: L₁ = E - r(n)
# where r(n) is a rational function such that q̃_n satisfies L₁(q̃_n) = 0,
# i.e., q̃_{n+1} = r(n) · q̃_n

# If L = L₂ · L₁, then any solution of L₁ is a solution of L.
# The Poincaré root of L₁ must be 1 (from the factor (ξ-1)).

# A first-order recurrence with Poincaré root 1: q̃_{n+1}/q̃_n → 1.
# So r(n) → 1 as n → ∞.

# Approach: guess r(n) = P(n)/Q(n) where deg P = deg Q,
# leading coefficients equal (so ratio → 1), and
# check if the sequence b_n defined by b_{n+1} = r(n) · b_n, b_0 = 1
# satisfies the order-3 recurrence.

# Actually, the standard approach: if we HAVE the order-3 recurrence,
# we can factor it by finding a Poincaré-root-1 solution.
# One method: compute the ratio cmf_q[n+1] / cmf_q[n] and see if
# there's a companion sequence with ratio → 1.

# Let's compute the sequence ratios
print("\n=== Q̂_{n+1} / Q̂_n ===", flush=True)
for n in range(20):
    if cmf_q[n] != 0:
        r = cmf_q[n+1] / cmf_q[n]
        print(f"  n={n}: {float(r):.15f}")

# This ratio → (3+2√2)² ≈ 33.97... = 17+12√2
print(f"\n  (3+2√2)² = {float(17+12*2**0.5):.15f}")

# The D_n² sequence has ratio D_{n+1}²/D_n² → (3+2√2)² too.
# Let's compute the "refined ratio" Q̂_n / D_n² and see how it grows.

# Key idea: compute the LEFT quotient of L by (E - 1).
# If L = L₂ · (E - r(n)) where r(n) → 1, then we need to find r(n).
# Method: the recurrence Σ ell_j(n) q_{n+j} = 0 can be written as
# ell_0(n) q_n + ell_1(n) q_{n+1} + ell_2(n) q_{n+2} + ell_3(n) q_{n+3} = 0
#
# If there's a first-order right factor E - r(n), then there's a solution
# s_n of s_{n+1} = r(n) s_n. And r(n) = s_{n+1}/s_n where s_n satisfies L.
#
# To find r(n), we need to find a solution with ratio → 1.
# The general solution of L is c₁ s_n^{(1)} + c₂ s_n^{(2)} + c₃ s_n^{(3)}
# where s_n^{(i)} ~ λ_i^n with λ₁ = 1, λ₂ = 17+12√2, λ₃ = 17-12√2.
# The bounded solution s_n^{(1)} with λ₁ = 1 gives r(n) → 1.

# But we only have Q̂_n (which grows like λ₂^n). To find s_n^{(1)},
# we need to use the recurrence BACKWARDS with appropriate initial conditions.

# Better approach: solve the recurrence from LARGE n backwards.
# At large n, s_n^{(1)} is the SUBDOMINANT solution.
# Running the recurrence FORWARD amplifies λ₂, losing the λ₁ component.
# Running BACKWARD amplifies the INVERSE of the smallest root (λ₃^{-1}),
# which helps isolate s^{(1)}.

# Actually, for a 3-term situation, let me try to use the CMF structure.
# The 3-column CMF naturally gives THREE independent sequences.
# Column 0, 1, 2 of the Q-row are three solutions of the same recurrence.

print("\n=== Three columns of Q-row as solutions ===", flush=True)
for n in range(8):
    print(f"  n={n}: Q={[float(q_full[n][j]) for j in range(3)]}")

# Check which column has ratio → 1 (bounded)
print("\n=== Column ratios ===", flush=True)
for col in range(3):
    print(f"\n  Column {col}:")
    for n in range(15):
        if q_full[n][col] != 0:
            r = q_full[n+1][col] / q_full[n][col]
            print(f"    n={n}: ratio = {float(r):.10f}")

# Column 2 might be the bounded one (Poincaré root 1).
# Let's check which column grows slowest.
print("\n=== Column magnitudes ===", flush=True)
for n in [0, 5, 10, 15, 20]:
    mags = [abs(float(q_full[n][j])) for j in range(3)]
    print(f"  n={n}: |Q_j| = {[f'{m:.6e}' for m in mags]}")

# Compute the determinant-like quantity to check independence
print("\n=== Wronskian-like check ===", flush=True)
for n in range(5):
    W = [[q_full[n+i][j] for j in range(3)] for i in range(3)]
    det = (W[0][0]*(W[1][1]*W[2][2]-W[1][2]*W[2][1])
          -W[0][1]*(W[1][0]*W[2][2]-W[1][2]*W[2][0])
          +W[0][2]*(W[1][0]*W[2][1]-W[1][1]*W[2][0]))
    print(f"  n={n}: det = {float(det):.6e}")

print("\nDone.")
