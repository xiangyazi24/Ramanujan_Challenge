#!/usr/bin/env python3
"""
Problem 2.5: Identify the connection constants β, γ to high precision.

q(N) = gauge(N) · (D_N² + β·D_N·E_N + γ·E_N²)

β = lim_{N→∞} r₁(N)/r₀(N), γ = lim_{N→∞} r₂(N)/r₀(N)

where r(N) = Φ(N)⁻¹ · [q(N), q(N+1), q(N+2)]^T

Use high-precision mpmath to identify β, γ as algebraic numbers.
"""
from mpmath import mp, mpf, matrix, sqrt, pi, catalan, nstr, fabs, log, identify
import sys

mp.dps = 100

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

# Compute q(N) in high precision
N_MAX = 60
I3 = matrix([[1,0,0],[0,1,0],[0,0,1]])
prod = I3.copy()
q = [mpf(1)]
for N in range(N_MAX + 3):
    prod = prod * M_mp(N)
    q.append(prod[0,0])

# Delannoy in high precision
D = [mpf(1), mpf(3)]
E = [mpf(0), mpf(1)]
for n in range(1, N_MAX + 10):
    D.append((mpf(3*(2*n+1)) * D[n] - mpf(n) * D[n-1]) / mpf(n+1))
    E.append((mpf(3*(2*n+1)) * E[n] - mpf(n) * E[n-1]) / mpf(n+1))

def v1(n): return D[n]**2
def v2(n): return D[n]*E[n]
def v3(n): return E[n]**2

print("=== High-precision β, γ identification ===\n")

# Compute r(N) via Casorati at large N for convergence
for N in [20, 30, 40, 50]:
    Phi = matrix([
        [v1(N), v2(N), v3(N)],
        [v1(N+1), v2(N+1), v3(N+1)],
        [v1(N+2), v2(N+2), v3(N+2)]
    ])
    b = matrix([q[N], q[N+1], q[N+2]])
    r = Phi**(-1) * b

    beta = r[1] / r[0]
    gamma = r[2] / r[0]
    print(f"N={N}: β = {nstr(beta, 50)}")
    print(f"       γ = {nstr(gamma, 50)}")
    print()

# Use N=50 for best precision
N = 50
Phi = matrix([
    [v1(N), v2(N), v3(N)],
    [v1(N+1), v2(N+1), v3(N+1)],
    [v1(N+2), v2(N+2), v3(N+2)]
])
b_vec = matrix([q[N], q[N+1], q[N+2]])
r = Phi**(-1) * b_vec
beta = r[1] / r[0]
gamma = r[2] / r[0]

print(f"β = {nstr(beta, 80)}")
print(f"γ = {nstr(gamma, 80)}")

# Try to identify
print("\n=== Attempting identification ===")
s2 = sqrt(2)

# Check various expressions
candidates_beta = {
    "-(3+2√2)": -(3+2*s2),
    "-17/3": mpf(-17)/3,
    "-23/4": mpf(-23)/4,
    "-(11+8√2)/3": -(11+8*s2)/3,
    "-2(1+√2)²": -2*(1+s2)**2,
    "-(5+4√2)/2": -(5+4*s2)/2,
    "-(2+√2)²": -(2+s2)**2,
    "-(1+√2)²": -(1+s2)**2,
    "-3(1+√2)": -3*(1+s2),
    "-(7+5√2)/2": -(7+5*s2)/2,
    "-(9+4√2)/2": -(9+4*s2)/2,
    "-(4+3√2)": -(4+3*s2),
    "3-3(1+√2)²": 3-3*(1+s2)**2,
    "1-2(1+√2)²": 1-2*(1+s2)**2,
    "3/(1-√2)-3": 3/(1-s2)-3,
    "-(53+38√2)/11": -(53+38*s2)/11,
    "-(8+3√2)": -(8+3*s2),
    "-(6+√2)": -(6+s2),
}

candidates_gamma = {
    "(3+2√2)²/4": (3+2*s2)**2/4,
    "25/3": mpf(25)/3,
    "(5+4√2)/2": (5+4*s2)/2,
    "(2+√2)²": (2+s2)**2,
    "2(1+√2)²": 2*(1+s2)**2,
    "(1+√2)²+2": (1+s2)**2+2,
    "3(1+√2)": 3*(1+s2),
    "(7+5√2)/2": (7+5*s2)/2,
    "(9+4√2)/2": (9+4*s2)/2,
    "(11+8√2)/3": (11+8*s2)/3,
    "(4+3√2)": 4+3*s2,
    "1+3(1+√2)": 1+3*(1+s2),
    "2+2(1+√2)²": 2+2*(1+s2)**2,
    "(8+3√2)": 8+3*s2,
    "(6+√2)": 6+s2,
    "5+2√2": 5+2*s2,
}

print("\nβ candidates (diff from target):")
for name, val in sorted(candidates_beta.items(), key=lambda x: fabs(beta - x[1])):
    diff = fabs(beta - val)
    if diff < 0.1:
        print(f"  {name} = {nstr(val, 20)}, diff = {nstr(diff, 10)}")

print("\nγ candidates (diff from target):")
for name, val in sorted(candidates_gamma.items(), key=lambda x: fabs(gamma - x[1])):
    diff = fabs(gamma - val)
    if diff < 0.1:
        print(f"  {name} = {nstr(val, 20)}, diff = {nstr(diff, 10)}")

# Try mpmath identify
print(f"\nmpmath identify(β): {identify(beta)}")
print(f"mpmath identify(γ): {identify(gamma)}")
print(f"mpmath identify(β+γ): {identify(beta + gamma)}")
print(f"mpmath identify(β*γ): {identify(beta * gamma)}")
print(f"mpmath identify(β/γ): {identify(beta / gamma)}")

# Check if β, γ are roots of a quadratic over Q
# t² - (β+γ)t + βγ = 0
bpg = beta + gamma
bxg = beta * gamma
print(f"\nβ + γ = {nstr(bpg, 50)}")
print(f"β · γ = {nstr(bxg, 50)}")
print(f"Discriminant = (β+γ)² - 4βγ = {nstr(bpg**2 - 4*bxg, 50)}")

# Maybe they're related to initial CMF data
# The initial CMF matrix A has rows (30921, -32972, 8240) and (33750, -36000, 9000)
# Maybe β, γ come from these?
print("\n=== Check CMF initial data connection ===")
A = matrix([[30921, -32972, 8240],
            [33750, -36000, 9000]])
print(f"A[0,:] / A[0,0] = {nstr(A[0,1]/A[0,0], 15)}, {nstr(A[0,2]/A[0,0], 15)}")
print(f"A[1,:] / A[1,0] = {nstr(A[1,1]/A[1,0], 15)}, {nstr(A[1,2]/A[1,0], 15)}")

# Also try β-related expressions with the convergent ratio
# E/D at infinity: E_N/D_N → C₊'/C₊
# Compute numerically
ED_ratio = E[50] / D[50]
print(f"\nE_N/D_N at N=50: {nstr(ED_ratio, 30)}")
print(f"β + γ·(E/D)² = {nstr(beta + gamma * ED_ratio**2, 20)}")
print(f"β·(E/D) + γ·(E/D)² = {nstr(beta*ED_ratio + gamma*ED_ratio**2, 20)}")

# Try integer relation: find a, b, c such that a + b·β + c·β² = 0
from mpmath import pslq
print(f"\npslq([1, β, β²]): {pslq([1, beta, beta**2])}")
print(f"pslq([1, γ, γ²]): {pslq([1, gamma, gamma**2])}")
print(f"pslq([1, β, γ]): {pslq([1, beta, gamma])}")
print(f"pslq([1, β, β², β³]): {pslq([1, beta, beta**2, beta**3])}")
print(f"pslq([1, √2, β]): {pslq([1, s2, beta])}")
print(f"pslq([1, √2, γ]): {pslq([1, s2, gamma])}")

print("\nDone.")
