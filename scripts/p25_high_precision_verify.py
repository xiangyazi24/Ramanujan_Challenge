#!/usr/bin/env python3
"""
P2.5: High-precision verification that P_{N,j}/Q_{N,j} → G
using mpmath with 300-digit precision.

Also verifies the Birkhoff error structure:
|P/Q - G| ~ C · n^{-3} · ρ^N where ρ = 17-12√2 ≈ 0.0294
"""
from mpmath import mp, mpf, matrix, log10, fabs, sqrt, catalan

mp.dps = 350  # 350-digit precision

def M_entries(n):
    n = mpf(n)
    m11 = (-2*n-5)*(n+3)**2*(136*n**4+1424*n**3+5548*n**2+9551*n+6141)
    m12 = 384*n**6+6384*n**5+44168*n**4+162698*n**3+336377*n**2+369933*n+169011
    m13 = -(480*n**4+4980*n**3+19210*n**2+32690*n+20730)
    m21 = (n+2)**2*(n+3)**2*(4*n+10)*(48*n**3+386*n**2+1017*n+879)
    m22 = (n+2)**2*(-272*n**5-3848*n**4-21732*n**3-61184*n**2-85761*n-47808)
    m23 = (n+2)**2*(320*n**3+2540*n**2+6610*n+5640)
    m31 = (-4*n-10)*(n+2)**2*(n+3)**2*(32*n**4+302*n**3+1037*n**2+1530*n+813)
    m32 = (n+2)**2*(192*n**6+2984*n**5+19116*n**4+64452*n**3+120256*n**2+117279*n+46476)
    m33 = (n+2)**2*(-16*n**5-408*n**4-2912*n**3-8884*n**2-12254*n-6240)
    return matrix([[m11,m12,m13],[m21,m22,m23],[m31,m32,m33]])

G = catalan  # mpmath's built-in Catalan constant
rho = 17 - 12*sqrt(2)
c_plus = mpf(16)*(17 + 12*sqrt(2))

print(f"G = {mp.nstr(G, 50)}")
print(f"ρ = {mp.nstr(rho, 50)}")
print(f"|c₊| = {mp.nstr(c_plus, 20)}")
print(f"log₁₀(1/ρ) = {mp.nstr(log10(1/rho), 10)} digits/step")

# Compute CMF products
A = matrix([[30921, -32972, 8240],
            [33750, -36000, 9000]])

print("\nComputing CMF products...", flush=True)
NMAX = 120
product = matrix([[1,0,0],[0,1,0],[0,0,1]])  # identity
ratios = []
errors = []

for N in range(NMAX):
    M = M_entries(N)
    product = product * M
    AM = A * product

    for j in range(3):
        P_j = AM[0, j]
        Q_j = AM[1, j]
        if Q_j != 0:
            ratio = P_j / Q_j
            err = fabs(ratio - G)
            if j == 0:
                ratios.append(ratio)
                errors.append(err)

    if (N+1) % 20 == 0:
        digits = -log10(errors[-1]) if errors[-1] > 0 else mp.dps
        print(f"  N={N+1}: {mp.nstr(digits, 6)} digits match G", flush=True)

# Final verification: all three columns at N=NMAX
AM = A * product
print(f"\n=== Final verification at N={NMAX} ===")
for j in range(3):
    P_j = AM[0, j]
    Q_j = AM[1, j]
    ratio = P_j / Q_j
    err = fabs(ratio - G)
    digits = -log10(err) if err > 0 else mp.dps
    print(f"  Column {j+1}: P/Q - G has {mp.nstr(digits, 6)} correct digits")

# Verify Birkhoff error structure: |P/Q - G| ~ C · n^{-3} · ρ^N
print(f"\n=== Birkhoff error structure ===")
print(f"Testing: |P/Q - G| · |c₊/c₀|^N · N^3 → const")
print(f"  |c₊/c₀| = |c₊|/16 = {mp.nstr(c_plus/16, 15)}")

for N in range(10, min(NMAX, 80), 5):
    if errors[N] > 0:
        # |P/Q - G| should be ~ C * ρ^N * N^{-3}
        # So |P/Q - G| * ρ^{-N} * N^3 should → C
        adjusted = errors[N] * rho**(-N) * mpf(N+1)**3
        print(f"  N={N+1}: |err|·ρ^{{-N}}·N³ = {mp.nstr(adjusted, 20)}")

# Digits per step
print(f"\n=== Digits gained per step ===")
for N in range(20, min(NMAX, 80), 10):
    if errors[N] > 0 and errors[N-1] > 0:
        dps = log10(errors[N-1]/errors[N])
        print(f"  N={N}: {mp.nstr(dps, 8)} digits/step (theory: {mp.nstr(log10(1/rho), 8)})")

print("\nDone.")
