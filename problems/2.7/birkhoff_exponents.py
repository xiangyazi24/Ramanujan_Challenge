#!/usr/bin/env python3
"""
Compute Birkhoff exponents for both the P2.7 and Zudilin monic recurrences.

For a monic order-3 recurrence:
  u_{n+3} = α(n)u_{n+2} + β(n)u_{n+1} + γ(n)u_n

with α(n) = α₀ + α₁/n + O(1/n²), etc., the Birkhoff exponent for root ρ_j is:
  σ_j = (α₁ρ_j² + β₁ρ_j + γ₁) / (ρ_j P'(ρ_j))

where P(ν) = ν³ - α₀ν² - β₀ν - γ₀ and P'(ν) = 3ν² - 2α₀ν - β₀.

KEY CLAIM: σ_j^Z = -3/2 for all j (Zudilin), σ_j^P = 0 for all j (P2.7 scaled).
"""
from fractions import Fraction as F
from mpmath import mp, mpf, polyroots, fabs

mp.dps = 50

# === P2.7 recurrence coefficients ===
def A_c(n): return 1024*(2*n+5)**4*(2*n+7)**3*(2*n+9)**3*(946*n*n+6407*n+10860)
def B_c(n): return 128*(2*n+7)**3*(2*n+9)**3*(104060*n**6+1745370*n**5+12145238*n**4+44886481*n**3+92943995*n**2+102256019*n+46709052)
def C_c(n): return 16*(n+3)**4*(2*n+9)**3*(3784*n**5+57792*n**4+351019*n**3+1059230*n**2+1587211*n+944620)
def D_c(n): return (n+3)**4*(n+4)**6*(946*n*n+4515*n+5399)

# P2.7 monic coefficients (after 64^n scaling):
# α_P(n) = 64·B(n+2)/A(n+2)
# β_P(n) = -64²·C(n+1)/A(n+1)
# γ_P(n) = 64³·D(n)/A(n)
def alpha_P(n):
    return F(64) * F(B_c(n+2), A_c(n+2))

def beta_P(n):
    return F(-64**2) * F(C_c(n+1), A_c(n+1))

def gamma_P(n):
    return F(64**3) * F(D_c(n), A_c(n))

# === Zudilin recurrence coefficients ===
def QZ(n): return 946*n**2 - 731*n + 153
def MZ(n): return 104060*n**6+127710*n**5+12788*n**4-34525*n**3-8482*n**2+3298*n+1071
def NZ(n): return 3784*n**5-1032*n**4-1925*n**3+853*n**2+328*n-184
def RZ(n): return 946*n**2+1161*n+368

# Zudilin monic coefficients:
# α_Z(n) = P(n) / (Q(n)(2n+1)(n+1)³)
# β_Z(n) = -n·S(n) / (Q(n)(2n+1)(n+1)³)
# γ_Z(n) = R(n)·n·(n-1)³ / (2Q(n)(2n+1)(n+1)³)
def alpha_Z(n):
    d = 2 * QZ(n) * (2*n+1) * (n+1)**3
    return F(2 * MZ(n), d)

def beta_Z(n):
    d = 2 * QZ(n) * (2*n+1) * (n+1)**3
    return F(-2 * n * NZ(n), d)

def gamma_Z(n):
    d = 2 * QZ(n) * (2*n+1) * (n+1)**3
    return F(RZ(n) * n * (n-1)**3, d)

print("=== Verify leading coefficients (n → ∞) ===")
for n_val in [100, 1000, 10000]:
    aP = float(alpha_P(n_val))
    bP = float(beta_P(n_val))
    gP = float(gamma_P(n_val))
    aZ = float(alpha_Z(n_val))
    bZ = float(beta_Z(n_val))
    gZ = float(gamma_Z(n_val))
    print(f"n={n_val}:")
    print(f"  P2.7:   α={aP:.10f}, β={bP:.10f}, γ={gP:.10f}")
    print(f"  Zudilin: α={aZ:.10f}, β={bZ:.10f}, γ={gZ:.10f}")

print("\n=== Extract 1/n coefficients: coeff_1 = n·(coeff(n) - coeff_0) ===")
for n_val in [1000, 10000, 100000]:
    n = n_val
    # P2.7
    a1P = float(n * (alpha_P(n) - 55))
    b1P = float(n * (beta_P(n) - (-2)))
    g1P = float(n * (gamma_P(n) - F(1,4)))
    # Zudilin
    a1Z = float(n * (alpha_Z(n) - 55))
    b1Z = float(n * (beta_Z(n) - (-2)))
    g1Z = float(n * (gamma_Z(n) - F(1,4)))
    print(f"n={n_val}:")
    print(f"  P2.7:    α₁ ≈ {a1P:.6f}, β₁ ≈ {b1P:.6f}, γ₁ ≈ {g1P:.6f}")
    print(f"  Zudilin: α₁ ≈ {a1Z:.6f}, β₁ ≈ {b1Z:.6f}, γ₁ ≈ {g1Z:.6f}")

print("\n=== Exact 1/n coefficients (algebraic, large n) ===")
# For Zudilin, analytically verified:
# α₁^Z = -165/2, β₁^Z = 6, γ₁^Z = -9/8
a1Z_exact = F(-165, 2)
b1Z_exact = F(6)
g1Z_exact = F(-9, 8)
print(f"Zudilin (analytic): α₁ = {a1Z_exact} = {float(a1Z_exact):.6f}, β₁ = {b1Z_exact}, γ₁ = {g1Z_exact} = {float(g1Z_exact):.6f}")

# Extract P2.7 coefficients from n=10000
n = 10000
a1P_num = float(n * (alpha_P(n) - 55))
b1P_num = float(n * (beta_P(n) + 2))
g1P_num = float(n * (gamma_P(n) - F(1,4)))
print(f"P2.7 (numerical at n=10000): α₁ ≈ {a1P_num:.10f}, β₁ ≈ {b1P_num:.10f}, γ₁ ≈ {g1P_num:.10f}")

# Try to find exact P2.7 α₁ by extrapolation
# n·(α(n) - 55) should converge to α₁ as n → ∞
# Also compute n²·(α(n) - 55) to see if α₁ = 0 and the next term is α₂/n²
for n_val in [1000, 10000, 100000]:
    a1P = float(n_val * (alpha_P(n_val) - 55))
    a2P = float(n_val**2 * (alpha_P(n_val) - 55))
    print(f"  n={n_val}: n(α-55) = {a1P:.10f}, n²(α-55) = {a2P:.6f}")

print("\n=== Poincaré roots ===")
# 4ν³ - 220ν² + 8ν - 1 = 0, or ν³ - 55ν² + 2ν - 1/4 = 0
roots = polyroots([1, -55, 2, -F(1,4)])
for i, r in enumerate(roots):
    print(f"  ρ_{i} = {r}")

print("\n=== Birkhoff exponents ===")
# σ_j = (α₁ρ² + β₁ρ + γ₁) / (ρ P'(ρ))
# P'(ν) = 3ν² - 110ν + 2
for i, rho in enumerate(roots):
    Pp = 3*rho**2 - 110*rho + 2
    # Zudilin
    num_Z = float(a1Z_exact) * rho**2 + float(b1Z_exact) * rho + float(g1Z_exact)
    sigma_Z = num_Z / (rho * Pp)
    # P2.7 (using numerical values)
    num_P = a1P_num * rho**2 + b1P_num * rho + g1P_num
    sigma_P = num_P / (rho * Pp)
    print(f"  ρ_{i} = {rho}")
    print(f"    σ_Z = {sigma_Z}  (should be -1.5)")
    print(f"    σ_P = {sigma_P}  (should be 0)")

print("\n=== Algebraic verification: σ_Z = -3/2 ===")
print("Need: (-165/2)ρ² + 6ρ - 9/8 = (-3/2)·ρ·(3ρ²-110ρ+2)")
print("RHS = (-3/2)(3ρ³ - 110ρ² + 2ρ)")
print("    = (-9/2)ρ³ + 165ρ² - 3ρ")
print("Using ρ³ = 55ρ² - 2ρ + 1/4:")
print("    = (-9/2)(55ρ² - 2ρ + 1/4) + 165ρ² - 3ρ")
print("    = -495/2·ρ² + 9ρ - 9/8 + 165ρ² - 3ρ")
print("    = (-495/2 + 330/2)ρ² + 6ρ - 9/8")
print("    = (-165/2)ρ² + 6ρ - 9/8 = LHS  ✓")

print("\n=== Verify gauge growth rate numerically ===")
# Compute the gauge G̃(n) acting on the error vector, and check growth
from fractions import Fraction as F

N = 30
def zudilin_terms(init, N):
    u = list(init)
    for n in range(2, N):
        d = 2 * QZ(n) * (2*n+1) * (n+1)**3
        nxt = F(2*MZ(n)) * u[n] + F(-2*n*NZ(n)) * u[n-1] + F(RZ(n)*n*(n-1)**3) * u[n-2]
        u.append(nxt / F(d))
    return u

def p27_terms(init, N):
    u = list(init)
    for n in range(2, N):
        nxt = F(B_c(n), A_c(n)) * u[n] + F(-C_c(n-1), A_c(n-1)) * u[n-1] + F(D_c(n-2), A_c(n-2)) * u[n-2]
        u.append(nxt)
    return u

b  = zudilin_terms([F(1), F(7), F(163)], N)
bt = zudilin_terms([F(0), F(23,2), F(2145,8)], N)
btt = zudilin_terms([F(0), F(17,2), F(3135,16)], N)
q = p27_terms([F(-215040420000), F(-167282265043404, 905), F(-964185327658080, 6071)], N)
p = p27_terms([F(-612218384750), F(-9525021973931919, 18100), F(-29561828382772029, 65380)], N)

# m_n = b̃_n + b̃̃_n
m = [bt[n] + btt[n] for n in range(N)]
qhat = [F(64)**n * q[n] for n in range(N)]
phat = [F(64)**n * p[n] for n in range(N)]

# Use high-precision ζ(2)+ζ(3) for error computation
from mpmath import zeta
mp.dps = 80
L = zeta(2) + zeta(3)

print("n | |ε_n|=|m_n-L·b_n| | |ê_n|=|p̂_n-L·q̂_n| | ratio ê/ε | ε·n^{3/2} | ê/|ρ₁|^n")
rho1_abs = abs(roots[1])
for n in range(3, 25):
    eps_n = abs(mpf(m[n].numerator)/mpf(m[n].denominator) - L * mpf(b[n].numerator)/mpf(b[n].denominator))
    ehat_n = abs(mpf(phat[n].numerator)/mpf(phat[n].denominator) - L * mpf(qhat[n].numerator)/mpf(qhat[n].denominator))
    if eps_n > 0 and ehat_n > 0:
        ratio = ehat_n / eps_n
        eps_scaled = eps_n * mpf(n)**1.5
        ehat_over_rho1 = ehat_n / rho1_abs**n
        print(f"{n:2d} | {float(eps_n):12.4e} | {float(ehat_n):12.4e} | {float(ratio):12.4e} | {float(eps_scaled):12.4e} | {float(ehat_over_rho1):12.4e}")

print("\n=== Check: ratio ê_n/ε_n should grow as n^{3/2} ===")
print("If σ_Z=-3/2 and σ_P=0, gauge ~ n^{3/2}, so |ê|/|ε| ~ n^{3/2}")
for n in range(5, 22):
    eps_n = abs(mpf(m[n].numerator)/mpf(m[n].denominator) - L * mpf(b[n].numerator)/mpf(b[n].denominator))
    ehat_n = abs(mpf(phat[n].numerator)/mpf(phat[n].denominator) - L * mpf(qhat[n].numerator)/mpf(qhat[n].denominator))
    if eps_n > 0 and ehat_n > 0:
        ratio = ehat_n / eps_n
        ratio_over_n32 = ratio / mpf(n)**1.5
        print(f"  n={n:2d}: |ê|/|ε| = {float(ratio):12.4e}, / n^{3/2} = {float(ratio_over_n32):12.4e}")
