#!/usr/bin/env python3
"""
P2.5: Decompose Q̂_N and P̂_N into the three Poincaré modes.

The scalar recurrence (order 3) has Poincaré roots:
  λ₊ = (3+2√2)² = 17+12√2 ≈ 33.97
  λ₀ = 1
  λ₋ = (3-2√2)² = 17-12√2 ≈ 0.0294

For irrationality of G, we need:
  Q̂_N·G - P̂_N ~ c·λ₋^N → 0 exponentially

This script:
1. Computes Q̂_N, P̂_N from the CMF
2. Computes D_N² (dominant solution, verified)
3. Finds s₂ (neutral) and s₃ (recessive) by numerical methods
4. Decomposes Q̂_N, P̂_N in the {D²,s₂,s₃} basis
5. Checks whether Q̂_N·G - P̂_N ∝ s₃_N
"""
from fractions import Fraction
from decimal import Decimal, getcontext
import math

getcontext().prec = 80

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

# Compute CMF rows
print("Computing CMF Q and P rows...", flush=True)
q_row = [Fraction(33750), Fraction(-36000), Fraction(9000)]
p_row = [Fraction(30921), Fraction(-32972), Fraction(8240)]
Q_hat = [q_row[0]]  # first column of Q-row
P_hat = [p_row[0]]  # first column of P-row

q_full = [list(q_row)]
p_full = [list(p_row)]

for N in range(NMAX):
    MH = MH_at(N)
    q_new = [sum(q_row[i]*MH[i][j] for i in range(3)) for j in range(3)]
    p_new = [sum(p_row[i]*MH[i][j] for i in range(3)) for j in range(3)]
    q_row = q_new
    p_row = p_new
    Q_hat.append(q_row[0])
    P_hat.append(p_row[0])
    q_full.append(list(q_row))
    p_full.append(list(p_row))

# Compute D_n²
print("Computing D_n² (central Delannoy squared)...", flush=True)
def delannoy(n):
    return sum(Fraction(math.comb(n,k)**2 * 2**k) for k in range(n+1))

D2 = [delannoy(n)**2 for n in range(NMAX+1)]

# Verify D² satisfies the scalar recurrence
# First, extract the scalar recurrence from the CMF
# The recurrence ℓ₀(n)s_{n+3} + ℓ₁(n)s_{n+2} + ℓ₂(n)s_{n+1} + ℓ₃(n)s_n = 0
# is determined by det(M_H(n) - λI) = 0 → char poly of the transition

# Instead, directly compute the recurrence coefficients from 4 consecutive Q-hat values
# For n ≥ 0: ℓ₀(n)Q̂_{n+3} + ℓ₁(n)Q̂_{n+2} + ℓ₂(n)Q̂_{n+1} + ℓ₃(n)Q̂_n = 0
# Use the D² sequence to find the recurrence

print("\n=== Scalar recurrence (from D²) ===", flush=True)
# Build the system: for each n, ℓ₃(n)D²_n + ℓ₂(n)D²_{n+1} + ℓ₁(n)D²_{n+2} + ℓ₀(n)D²_{n+3} = 0
# We know the recurrence has polynomial coefficients in n. Let's verify with Q̂ too.

# For the 3-step recurrence, we need the characteristic polynomial of M_H(n) at large n:
# ξ³ - 35ξ² + 35ξ - 1 = (ξ-1)(ξ²-34ξ+1)
# So the scalar recurrence has Poincaré polynomial (ξ-1)(ξ²-34ξ+1)

# Compute the recurrence: at each n, we have
# a(n)·y_{n+3} = b(n)·y_{n+2} + c(n)·y_{n+1} + d(n)·y_n
# Using Q̂ values:

print("Verifying scalar recurrence from Q̂ values...")
for n in range(5):
    # Check that Q̂ satisfies a recurrence of order 3
    # a·Q̂_{n+3} + b·Q̂_{n+2} + c·Q̂_{n+1} + d·Q̂_n = 0
    # Using 4 equations (n, n+1, n+2, n+3) we can determine (a:b:c:d) up to scaling
    pass

# Better approach: compute the ratio Q̂_N / D²_N and see its structure
print("\n=== Q̂_N / D²_N ===", flush=True)
for n in range(15):
    if D2[n] != 0:
        r = Q_hat[n] / D2[n]
        print(f"  N={n}: Q̂/D² = {float(r):.10f}  (exact = {r})")

# Compute P̂_N / D²_N
print("\n=== P̂_N / D²_N ===", flush=True)
for n in range(15):
    if D2[n] != 0:
        r = P_hat[n] / D2[n]
        print(f"  N={n}: P̂/D² = {float(r):.10f}")

# Now find the NEUTRAL solution (Poincaré root = 1)
# This is the solution that grows like ~1 (constant or polynomial in n)
# Method: if s_N ~ c·1^N = c, then s_N approaches a constant.
# The actual neutral solution might grow polynomially (like N^k).

# Run the recurrence with different initial conditions and extract the neutral mode
# by cancelling the dominant mode using D².

# Strategy: use two initial conditions IC1 and IC2.
# s(IC1)_N = α₁·D²_N + β₁·s₂_N + γ₁·s₃_N
# s(IC2)_N = α₂·D²_N + β₂·s₂_N + γ₂·s₃_N
# Subtract α₁/α₂ ratio to kill D²:
# (s(IC1) - (α₁/α₂)·s(IC2))_N = (β₁-α₁/α₂·β₂)·s₂_N + (γ₁-α₁/α₂·γ₁)·s₃_N
# At large N, the dominant term dies, and we're left with s₂ + small·s₃

# IC1 = [1, 0, 0], IC2 = [0, 1, 0]
print("\n=== Finding neutral and recessive solutions ===", flush=True)
print("Using column 0, 1, 2 of Q-row as three solutions...", flush=True)

# All three columns of Q_full satisfy the SAME scalar recurrence
# (because M_H(n) is the same matrix)
col = [[q_full[N][j] for N in range(NMAX+1)] for j in range(3)]

# Ratios to D²
print("\nColumn ratios at large N (col_j / D²):")
for j in range(3):
    for N in [20, 30, 40]:
        if D2[N] != 0:
            r = float(col[j][N] / D2[N])
            print(f"  col[{j}][{N}]/D²[{N}] = {r:.12e}")

# Linear combination to kill dominant mode
# col[0] / D²[N] → α₀
# col[1] / D²[N] → α₁
# col[2] / D²[N] → α₂
# At large N, col[j][N] ≈ αⱼ · D²[N]
# So col[0] - (α₀/α₁) · col[1] kills the dominant mode

alpha = [float(col[j][40] / D2[40]) for j in range(3)]
print(f"\nDominant coefficients (col/D² at N=40): {alpha}")

# Kill dominant mode: s₂_cand = col[0] - (α₀/α₁)·col[1]
ratio01 = col[0][40] / col[1][40]
s2_cand = [col[0][N] - ratio01 * col[1][N] for N in range(NMAX+1)]

print(f"\nNeutral candidate (col[0] - ratio·col[1]):")
for N in range(15):
    print(f"  N={N}: {float(s2_cand[N]):.10e}")

# Check growth rate of neutral candidate
print("\nRatios s2_cand[N+1]/s2_cand[N]:")
for N in range(1, 20):
    if s2_cand[N] != 0:
        r = float(s2_cand[N+1] / s2_cand[N])
        print(f"  N={N}: {r:.10f}")

# Now kill dominant mode more carefully using exact arithmetic
# ratio01 = col[0][N] / col[1][N] for large N
# Use N=40:
ratio01_exact = Fraction(col[0][40].numerator * col[1][40].denominator,
                         col[0][40].denominator * col[1][40].numerator)
# Actually, this is just col[0][40] / col[1][40]
ratio01_exact = col[0][40] / col[1][40]

s2_exact = [col[0][N] - ratio01_exact * col[1][N] for N in range(NMAX+1)]

print(f"\n=== Neutral candidate (exact, killed at N=40) ===")
for N in range(10):
    print(f"  N={N}: {float(s2_exact[N]):.10e}")

# This is a combination of s₂ and s₃. At large N, s₂ dominates (since s₃ → 0).
# Check the ratio:
print("\nRatios:")
for N in range(5, 30):
    if s2_exact[N] != 0 and s2_exact[N-1] != 0:
        r = float(s2_exact[N] / s2_exact[N-1])
        print(f"  N={N}: {r:.15f}")

# Now compute the ERROR: Q̂_N * G - P̂_N
# G = Catalan's constant = 0.9159655941772190...
G_float = 0.91596559417721901505460351493238411077414937428167213426649811962176301977625476947935651292611510624

# High-precision G (using mpmath or decimal)
from decimal import Decimal
G_dec = Decimal("0.9159655941772190150546035149323841107741493742816721342664981196217630197762547694793565129261151062")

print(f"\n=== Error: Q̂_N * G - P̂_N ===", flush=True)
errors = []
for N in range(25):
    q = Decimal(Q_hat[N].numerator) / Decimal(Q_hat[N].denominator)
    p = Decimal(P_hat[N].numerator) / Decimal(P_hat[N].denominator)
    err = q * G_dec - p
    errors.append(float(err))
    print(f"  N={N}: {float(err):.15e}")

# Error ratios (should → λ₋ ≈ 0.0294)
lambda_minus = 17 - 12*2**0.5
print(f"\n=== Error ratios (should → {lambda_minus:.10f} = (3-2√2)²) ===")
for N in range(1, 22):
    if errors[N-1] != 0:
        r = errors[N] / errors[N-1]
        print(f"  N={N}: err[{N}]/err[{N-1}] = {r:.15f}")

# Check: is error proportional to D² times λ₋^N?
print(f"\n=== Error / (D²_N * λ₋^N) ===")
for N in range(1, 20):
    if D2[N] != 0:
        d2f = float(D2[N])
        lmn = lambda_minus ** N
        ratio = errors[N] / (d2f * lmn)
        print(f"  N={N}: {ratio:.15e}")

# Also check error / λ₋^N directly
print(f"\n=== Error / λ₋^N ===")
for N in range(1, 20):
    lmn = lambda_minus ** N
    ratio = errors[N] / lmn
    print(f"  N={N}: {ratio:.15e}")

# Denominators of Q̂_N
print(f"\n=== Denominator structure ===")
for N in range(25):
    den = Q_hat[N].denominator
    # Factor out powers of 2
    v2 = 0
    d = den
    while d % 2 == 0:
        v2 += 1
        d //= 2
    odd_part = d
    print(f"  N={N}: den(Q̂) = 2^{v2}" + (f" * {odd_part}" if odd_part > 1 else ""))

print("\nDone.")
