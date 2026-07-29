#!/usr/bin/env python3
"""
P2.5: Quantify the convergence rate of ε_k = g(k)/f(k) - G.

From the previous check: g(k)/f(k) → G but ≠ C_k.
The convergence appears exponential: ε_k ~ C · r^k.
Compute ε_k to high precision and find r.

If |ε_k| ≤ C·r^k with r < 1, then since weights w_k(N) = f(k)F(N,k)/Q̂_N
concentrate near k ≈ αN with exponential tails:
  E_N/Q̂_N = Σ ε_k w_k → 0
proving L = G.

KEY: the ratio r must be identified. If r = ρ = 17-12√2 = 1/ξ_+,
this connects to the Poincaré structure of the k-recurrence.
"""
from mpmath import mp, mpf, catalan, sqrt, log, log10, pi, fabs
from fractions import Fraction as Fr
from math import comb

mp.dps = 100

def M_entries(n):
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
    return -2*(n+2)**2*(n+3)**2*(2*n+5)*(2*n+7)**2

def delannoy_summand(N, k):
    if k < 0 or k > N:
        return Fr(0)
    return Fr(2**k * comb(2*k, k) * comb(N, k) * comb(N+k, k))

KMAX = 40

# Compute CMF trajectories (exact fractions)
rows = {
    'e1': [Fr(1), Fr(0), Fr(0)],
    'e2': [Fr(0), Fr(1), Fr(0)],
    'e3': [Fr(0), Fr(0), Fr(1)],
}
history = {key: [[v[j] for j in range(3)]] for key, v in rows.items()}

for N in range(KMAX):
    M = M_entries(N)
    d = Fr(delta_H(N))
    MH = [[Fr(M[i][j]) / d for j in range(3)] for i in range(3)]
    for key in rows:
        r = rows[key]
        new_r = [sum(r[i]*MH[i][k] for i in range(3)) for k in range(3)]
        rows[key] = new_r
        history[key].append([new_r[j] for j in range(3)])

q = [Fr(33750), Fr(-36000), Fr(9000)]
p = [Fr(30921), Fr(-32972), Fr(8240)]

Q_vals = []
P_vals = []
for N in range(KMAX + 1):
    Q_N = sum(q[j] * history['e' + str(j+1)][N][0] for j in range(3))
    P_N = sum(p[j] * history['e' + str(j+1)][N][0] for j in range(3))
    Q_vals.append(Q_N)
    P_vals.append(P_N)

# Decompose in Delannoy basis
def decompose(vals, KMAX):
    coeffs = []
    for K in range(KMAX + 1):
        rhs = vals[K]
        for k in range(K):
            rhs -= coeffs[k] * delannoy_summand(K, k)
        bKK = delannoy_summand(K, K)
        coeffs.append(rhs / bKK)
    return coeffs

f_coeffs = decompose(Q_vals, KMAX)
g_coeffs = decompose(P_vals, KMAX)

G = catalan

print("ε_k = g(k)/f(k) - G  (high precision)")
print("=" * 80)

epsilons = []
ratios = []
for k in range(KMAX + 1):
    if f_coeffs[k] == 0:
        continue
    ratio_exact = g_coeffs[k] / f_coeffs[k]
    # Convert to mpf for ε_k computation
    ratio_mp = mpf(ratio_exact.numerator) / mpf(ratio_exact.denominator)
    eps_k = ratio_mp - G
    epsilons.append((k, eps_k))

    abs_eps = fabs(eps_k)
    if abs_eps > 0:
        log_abs = float(mp.log10(abs_eps))
    else:
        log_abs = -100

    sign = '+' if eps_k > 0 else '-'

    if k >= 1 and len(epsilons) >= 2:
        prev_eps = epsilons[-2][1]
        if prev_eps != 0:
            r = eps_k / prev_eps
            ratios.append((k, float(r)))
            print(f"  k={k:2d}: ε = {sign}{mp.nstr(abs_eps, 20):>30s}  "
                  f"log₁₀|ε| = {log_abs:8.3f}  ε_k/ε_{k-1} = {float(r):+.10f}")
        else:
            print(f"  k={k:2d}: ε = {sign}{mp.nstr(abs_eps, 20):>30s}  "
                  f"log₁₀|ε| = {log_abs:8.3f}")
    else:
        print(f"  k={k:2d}: ε = {sign}{mp.nstr(abs_eps, 20):>30s}  "
              f"log₁₀|ε| = {log_abs:8.3f}")

print()

# Analyze the ratio ε_{k+1}/ε_k
print("=== Ratio analysis ε_{k+1}/ε_k ===")
rho = 17 - 12*float(sqrt(mpf(2)))  # ≈ 0.02944
minus_rho = -rho
print(f"ρ = 17 - 12√2 = {rho:.10f}")
print(f"-1/16 = {-1/16:.10f}")
print(f"c₀/c₁ = -16/544 = {-16/544:.10f}")

# Fit exponential: ε_k = A · r^k
# From consecutive ratios:
if len(ratios) >= 5:
    late_ratios = [r for k, r in ratios if k >= 8]
    if late_ratios:
        avg_ratio = sum(late_ratios) / len(late_ratios)
        print(f"\nAverage ratio for k ≥ 8: {avg_ratio:.15f}")
        print(f"= 1/{1/avg_ratio:.10f}" if avg_ratio != 0 else "")

        # Is it -ρ?
        print(f"\nCompare with known quantities:")
        print(f"  ρ = {rho:.15f}")
        print(f"  -ρ = {-rho:.15f}")
        print(f"  ratio = {avg_ratio:.15f}")
        print(f"  ratio/(-ρ) = {avg_ratio/(-rho):.15f}")
        print(f"  -1/ξ₊ = {-1/(17+12*float(sqrt(mpf(2)))):.15f}")
        print(f"  ratio/(-1/ξ₊) = {avg_ratio/(-1/(17+12*float(sqrt(mpf(2))))):.15f}")

# Check: does ε_k satisfy a recurrence?
print("\n=== Recurrence check for ε_k ===")
print("If ε_k satisfies: a₃ε_{k+3} + a₂ε_{k+2} + a₁ε_{k+1} + a₀ε_k = 0")
print("then we can solve for the Poincaré roots.")

eps_vals = [e for _, e in epsilons]
if len(eps_vals) >= 10:
    # Try order 3 with constant coefficients
    # Solve a₃ε_{k+3} + a₂ε_{k+2} + a₁ε_{k+1} + a₀ε_k = 0
    # for a₀, a₁, a₂ (set a₃ = 1)
    # Use k = 5, 6, 7 (avoiding small k)
    from mpmath import matrix, lu_solve

    for start_k in [5, 10, 15, 20]:
        if start_k + 5 < len(eps_vals):
            A = matrix(3, 3)
            b = matrix(3, 1)
            for row in range(3):
                k = start_k + row
                A[row, 0] = eps_vals[k]
                A[row, 1] = eps_vals[k+1]
                A[row, 2] = eps_vals[k+2]
                b[row] = -eps_vals[k+3]
            try:
                x = lu_solve(A, b)
                a0, a1, a2 = float(x[0]), float(x[1]), float(x[2])

                # Verify
                residuals = []
                for k in range(start_k, min(start_k + 8, len(eps_vals) - 3)):
                    res = eps_vals[k+3] + a2*eps_vals[k+2] + a1*eps_vals[k+1] + a0*eps_vals[k]
                    residuals.append(float(fabs(res) / fabs(eps_vals[k+3])) if eps_vals[k+3] != 0 else 0)

                print(f"  Starting from k={start_k}: a₀={a0:.10f}, a₁={a1:.10f}, a₂={a2:.10f}")
                print(f"    Relative residuals: {[f'{r:.2e}' for r in residuals]}")

                # Characteristic polynomial: λ³ + a₂λ² + a₁λ + a₀ = 0
                import numpy as np
                roots = np.roots([1, a2, a1, a0])
                print(f"    Characteristic roots: {roots}")
            except:
                print(f"  k={start_k}: singular system")

# Also check: are the ratios ε_{k+1}/ε_k approaching a specific value?
print("\n=== Ratio convergence ===")
for k in range(max(0, len(epsilons)-15), len(epsilons)-1):
    idx, eps = epsilons[k]
    idx_next, eps_next = epsilons[k+1]
    if eps != 0:
        r = eps_next / eps
        print(f"  ε_{idx_next}/ε_{idx} = {mp.nstr(r, 30)}")

print("\nDone.")
