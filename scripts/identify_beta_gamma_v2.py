#!/usr/bin/env python3
"""
Problem 2.5: Identify β, γ using PSLQ with moderate N (before singularity).

β ≈ -5.770780163555853629...
γ ≈ 8.325475924022431191...
"""
from mpmath import mp, mpf, matrix, sqrt, pi, catalan, nstr, fabs, log, identify, pslq

mp.dps = 80

def M_mp(n):
    n = mpf(n)
    m11 = (-2*n-5)*(n+3)**2 * (136*n**4 + 1424*n**3 + 5548*n**2 + 9551*n + 6141)
    m12 = 384*n**6 + 6384*n**5 + 44168*n**4 + 162698*n**3 + 336377*n**2 + 369933*n + 169011
    m13 = -480*n**4 - 4980*n**3 - 19210*n**2 - 32690*n - 20730
    m21 = (n+2)**2*(n+3)**2*(4*n+10)*(48*n**3 + 386*n**2 + 1017*n + 879)
    m22 = (n+2)**2*(-272*n**5 - 3848*n**4 - 21732*n**3 - 61184*n**2 - 85761*n - 47808)
    m23 = (n+2)**2*(320*n**3 + 2540*n**2 + 6610*n + 5640)
    m31 = (-4*n-10)*(n+2)**2*(n+3)**2*(32*n**4 + 302*n**3 + 1037*n**2 + 1530*n + 813)
    m32 = (n+2)**2*(192*n**6 + 2984*n**5 + 19116*n**4 + 64452*n**3 + 120256*n**2 + 117279*n + 46476)
    m33 = (n+2)**2*(-16*n**5 - 408*n**4 - 2912*n**3 - 8884*n**2 - 12254*n - 6240)
    return matrix([[m11, m12, m13], [m21, m22, m23], [m31, m32, m33]])

# Compute q(N)
N_MAX = 35
I3 = matrix([[1,0,0],[0,1,0],[0,0,1]])
prod = I3.copy()
q = [mpf(1)]
for N in range(N_MAX + 3):
    prod = prod * M_mp(N)
    q.append(prod[0,0])

D = [mpf(1), mpf(3)]
E = [mpf(0), mpf(1)]
for n in range(1, N_MAX + 10):
    D.append((mpf(3*(2*n+1)) * D[n] - mpf(n) * D[n-1]) / mpf(n+1))
    E.append((mpf(3*(2*n+1)) * E[n] - mpf(n) * E[n-1]) / mpf(n+1))

# Use N=25 (good balance of precision and conditioning)
N = 25
Phi = matrix([
    [D[N]**2, D[N]*E[N], E[N]**2],
    [D[N+1]**2, D[N+1]*E[N+1], E[N+1]**2],
    [D[N+2]**2, D[N+2]*E[N+2], E[N+2]**2]
])
b_vec = matrix([q[N], q[N+1], q[N+2]])
r = Phi**(-1) * b_vec
beta = r[1] / r[0]
gamma = r[2] / r[0]

print(f"β = {nstr(beta, 60)}")
print(f"γ = {nstr(gamma, 60)}")

s2 = sqrt(2)

# PSLQ: find integer relations
print("\n=== PSLQ integer relation detection ===")
print(f"pslq([1, β, β²]): {pslq([1, beta, beta**2])}")
print(f"pslq([1, γ, γ²]): {pslq([1, gamma, gamma**2])}")
print(f"pslq([1, √2, β]): {pslq([1, s2, beta])}")
print(f"pslq([1, √2, γ]): {pslq([1, s2, gamma])}")

# Check β + γ and β·γ
bpg = beta + gamma
bxg = beta * gamma
print(f"\nβ + γ = {nstr(bpg, 50)}")
print(f"β · γ = {nstr(bxg, 50)}")
disc = bpg**2 - 4*bxg
print(f"Disc = (β+γ)² - 4βγ = {nstr(disc, 50)}")

print(f"\npslq([1, β+γ, (β+γ)²]): {pslq([1, bpg, bpg**2])}")
print(f"pslq([1, β·γ, (β·γ)²]): {pslq([1, bxg, bxg**2])}")
print(f"pslq([1, √2, β+γ]): {pslq([1, s2, bpg])}")
print(f"pslq([1, √2, β·γ]): {pslq([1, s2, bxg])}")
print(f"pslq([1, √2, disc]): {pslq([1, s2, disc])}")

# Check combined with G (Catalan's constant)
G = catalan
print(f"\npslq([1, G, β]): {pslq([1, G, beta])}")
print(f"pslq([1, G, γ]): {pslq([1, G, gamma])}")
print(f"pslq([1, π, β]): {pslq([1, pi, beta])}")
print(f"pslq([1, G/π, β]): {pslq([1, G/pi, beta])}")

# Try quadratic over Z[√2]
print("\n=== Quadratic over Q(√2) ===")
print(f"pslq([1, √2, β, β√2, β²]): {pslq([1, s2, beta, beta*s2, beta**2])}")
print(f"pslq([1, √2, γ, γ√2, γ²]): {pslq([1, s2, gamma, gamma*s2, gamma**2])}")

# Also try: maybe β = p/q + r/s * √2 for small integers
# Test: β + a + b√2 = 0?
# i.e., β = -a - b√2
print("\n=== Is β = a + b√2? ===")
result = pslq([beta, 1, s2])
print(f"pslq([β, 1, √2]): {result}")
if result:
    # result[0]*β + result[1] + result[2]*√2 = 0
    # β = -(result[1] + result[2]*√2) / result[0]
    print(f"  → β = ({-result[1]} + {-result[2]}√2) / {result[0]}")
    reconstructed = (-result[1] - result[2]*s2) / result[0]
    print(f"  → β = {nstr(reconstructed, 30)}")
    print(f"  → diff = {nstr(fabs(beta - reconstructed), 10)}")

result = pslq([gamma, 1, s2])
print(f"pslq([γ, 1, √2]): {result}")
if result:
    print(f"  → γ = ({-result[1]} + {-result[2]}√2) / {result[0]}")
    reconstructed = (-result[1] - result[2]*s2) / result[0]
    print(f"  → γ = {nstr(reconstructed, 30)}")
    print(f"  → diff = {nstr(fabs(gamma - reconstructed), 10)}")

# Also try: is β related to E/D ratio?
ED = E[30] / D[30]
print(f"\nE/D at N=30: {nstr(ED, 30)}")
print(f"pslq([β, 1, E/D]): {pslq([beta, 1, ED])}")
print(f"pslq([γ, 1, E/D]): {pslq([gamma, 1, ED])}")
print(f"pslq([β, 1, E/D, (E/D)²]): {pslq([beta, 1, ED, ED**2])}")

# Try log(3+2√2) or log(17+12√2)
from mpmath import ln
L1 = ln(3+2*s2)
L2 = ln(17+12*s2)
print(f"\npslq([β, 1, ln(3+2√2)]): {pslq([beta, 1, L1])}")
print(f"pslq([β, 1, ln(17+12√2)]): {pslq([beta, 1, L2])}")

# Try mpmath identify
print(f"\nmpmath identify(β): {identify(float(beta))}")
print(f"mpmath identify(γ): {identify(float(gamma))}")
print(f"mpmath identify(β+γ): {identify(float(bpg))}")

# Compute the gauge ratio r₀(N)/r₀(N-1) at various N
print("\n\n=== Gauge ratio analysis ===")
gauge_ratios_exact = []
for N in range(1, 25):
    Phi_N = matrix([
        [D[N]**2, D[N]*E[N], E[N]**2],
        [D[N+1]**2, D[N+1]*E[N+1], E[N+1]**2],
        [D[N+2]**2, D[N+2]*E[N+2], E[N+2]**2]
    ])
    Phi_prev = matrix([
        [D[N-1]**2, D[N-1]*E[N-1], E[N-1]**2],
        [D[N]**2, D[N]*E[N], E[N]**2],
        [D[N+1]**2, D[N+1]*E[N+1], E[N+1]**2]
    ])
    r_N = (Phi_N**(-1)) * matrix([q[N], q[N+1], q[N+2]])
    r_prev = (Phi_prev**(-1)) * matrix([q[N-1], q[N], q[N+1]])

    if r_prev[0] != 0:
        gr = r_N[0] / r_prev[0]
        gauge_ratios_exact.append((N, gr))
        # Normalize by the expected Pochhammer structure
        # g(N)/g(N-1) should involve Pochhammer products
        # From the recurrence, the dominant solution has g ~ ∏ Poincaré root
        # Actually, try dividing by (-16)^1 * N^7 * Pochhammer factors

        # Try normalizing by c₀(N-1)/c₃(N-1) (the full coefficient ratio)
        if N <= 10:
            print(f"N={N}: gauge_ratio = {nstr(gr, 25)}")

print("\nDone.")
