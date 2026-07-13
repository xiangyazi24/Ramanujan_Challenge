#!/usr/bin/env python3
"""Verify ChatGPT's claimed decomposition of Problem 2.3.
Claim: p_n and q_n are built from Lambert continuants (A_m, B_m)
and derangement numbers D_m, factorial F_m = m!."""
from mpmath import mp, mpf, pi, e as euler_e, fac
mp.dps = 50

# Lambert continuants: A_{-1}=1, A_0=0, A_m = m*A_{m-1} + A_{m-2}
# (same for B: B_{-1}=0, B_0=1)
def lambert(N):
    A = [mpf(1), mpf(0)]  # A[-1], A[0]
    B = [mpf(0), mpf(1)]  # B[-1], B[0]
    for m in range(1, N+1):
        A.append(m * A[-1] + A[-2])
        B.append(m * B[-1] + B[-2])
    return A, B  # indexed: A[0]=A_{-1}, A[1]=A_0, A[m+1]=A_m

# Derangement numbers: D_0=1, D_1=0, D_m = (m-1)(D_{m-1} + D_{m-2})
def derangements(N):
    D = [mpf(1), mpf(0)]
    for m in range(2, N+1):
        D.append((m-1) * (D[-1] + D[-2]))
    return D  # D[m] = !m

A, B = lambert(20)
D = derangements(20)

print("=== Lambert continuants (A_m/B_m → π) ===")
for m in range(5, 15):
    ratio = A[m+1] / B[m+1]
    print(f"  A_{m}/B_{m} = {ratio} (diff from π: {ratio - pi})")

print("\n=== Derangement D_m/m! → 1/e ===")
for m in range(5, 15):
    ratio = D[m] / fac(m)
    print(f"  D_{m}/{m}! = {ratio} (diff from 1/e: {ratio - 1/euler_e})")

# Now try to match the problem's initial values
# The ChatGPT proof claims (formulas 5.1 and 5.2):
# p_n = A_{n+2} * (n+3)! + B_{n+2} * D_{n+3}
# q_n = B_{n+2} * (n+3)! + A_{n+2} * D_{n+3}
# Wait, actually I'm not sure of the exact formulas because LaTeX was stripped.
# Let me try both orderings and check initial values.

print("\n=== Testing decomposition formulas ===")
print("Target: p_{-3}=1, p_{-2}=1, p_{-1}=20, p_0=296")
print("Target: q_{-3}=1, q_{-2}=0, q_{-1}=4, q_0=48")

# A[m+1] = A_m, so A_{n+2} = A[n+3]
# B[m+1] = B_m, so B_{n+2} = B[n+3]
# D_{n+3} = D[n+3]
# (n+3)! = fac(n+3)

# Formula 1: p_n = A_{n+2} * (n+3)! + B_{n+2} * D_{n+3}
print("\nFormula 1: p_n = A_{n+2}*(n+3)! + B_{n+2}*D_{n+3}")
for n in [-3, -2, -1, 0]:
    # A_{n+2}: need A[-1], A[0], A[1], A[2] = A[0], A[1], A[2], A[3] in our array
    Am = A[n+3]  # A_{n+2}
    Bm = B[n+3]  # B_{n+2}
    Dm = D[n+3]  # D_{n+3}
    Fm = fac(n+3)  # (n+3)!
    p = Am * Fm + Bm * Dm
    print(f"  n={n}: A_{n+2}={Am}, B_{n+2}={Bm}, D_{n+3}={Dm}, (n+3)!={Fm} -> p={p}")

print("\nFormula 1: q_n = B_{n+2}*(n+3)! + A_{n+2}*D_{n+3}")
for n in [-3, -2, -1, 0]:
    Am = A[n+3]
    Bm = B[n+3]
    Dm = D[n+3]
    Fm = fac(n+3)
    q = Bm * Fm + Am * Dm
    print(f"  n={n}: q={q}")

# Try swapped formulas
print("\nFormula 2: p_n = A_{n+2}*(n+3)! - B_{n+2}*D_{n+3}")
for n in [-3, -2, -1, 0]:
    Am = A[n+3]
    Bm = B[n+3]
    Dm = D[n+3]
    Fm = fac(n+3)
    p = Am * Fm - Bm * Dm
    print(f"  n={n}: p={p}")

# Try with different index shifts
print("\n=== Trying different index shifts ===")
# Maybe: p_n = A_{n+3}*(n+3)! + B_{n+3}*D_{n+3}
print("\nShift+1: p_n = A_{n+3}*(n+3)! + B_{n+3}*D_{n+3}")
for n in [-3, -2, -1, 0]:
    if n+4 < len(A):
        Am = A[n+4]
        Bm = B[n+4]
        Dm = D[n+3]
        Fm = fac(n+3)
        p = Am * Fm + Bm * Dm
        print(f"  n={n}: p={p}")
