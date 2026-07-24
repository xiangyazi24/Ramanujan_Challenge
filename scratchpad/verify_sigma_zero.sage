"""
Verify σ=0 numerically by checking q_n / μ₀^n → constant.
If σ≠0, the ratio would grow/decay as n^{-σ}.

Also compute the adjoint bracket numerically.
"""
from sage.all import *

# Recurrence coefficients
def A_raw(n): return 1024*(2*n+5)**4*(2*n+7)**3*(2*n+9)**3*(946*n**2+6407*n+10860)
def B_raw(n): return 128*(2*n+7)**3*(2*n+9)**3*(104060*n**6+1745370*n**5+12145238*n**4+44886481*n**3+92943995*n**2+102256019*n+46709052)
def C_raw(n): return 16*(n+3)**4*(2*n+9)**3*(3784*n**5+57792*n**4+351019*n**3+1059230*n**2+1587211*n+944620)
def D_raw(n): return (n+3)**4*(n+4)**6*(946*n**2+4515*n+5399)

# Initial conditions
p = [QQ(-612218384750), QQ(-9525021973931919)/QQ(18100), QQ(-29561828382772029)/QQ(65380)]
q = [QQ(-215040420000), QQ(-167282265043404)/QQ(905), QQ(-964185327658080)/QQ(6071)]

N_MAX = 300
for _ in range(N_MAX - 2):
    n = len(p) - 1
    p_new = B_raw(n)/A_raw(n) * p[n] - C_raw(n-1)/A_raw(n-1) * p[n-1] + D_raw(n-2)/A_raw(n-2) * p[n-2]
    q_new = B_raw(n)/A_raw(n) * q[n] - C_raw(n-1)/A_raw(n-1) * q[n-1] + D_raw(n-2)/A_raw(n-2) * q[n-2]
    p.append(QQ(p_new))
    q.append(QQ(q_new))

print(f"Computed {len(p)} terms")

# Poincaré root
RR500 = RealField(500)
R_poly = PolynomialRing(QQ, 'x')
x = R_poly.gen()
P_cubic = 4*x**3 - 220*x**2 + 8*x - 1
CC500 = ComplexField(500)
roots = P_cubic.change_ring(CC500).roots(multiplicities=False)
mu0_full = max(roots, key=lambda z: z.real())
mu0 = RR500(mu0_full.real())
mu_poincare = mu0 / 64  # actual Poincaré root

print(f"\nμ₀ (root of 4x³-220x²+8x-1) = {float(mu0):.15f}")
print(f"μ_Poincaré = μ₀/64 = {float(mu_poincare):.15f}")
print(f"q_1/q_0 = {float(RR500(q[1])/RR500(q[0])):.15f}")
print(f"q_10/q_9 = {float(RR500(q[10])/RR500(q[9])):.15f}")
print(f"q_100/q_99 = {float(RR500(q[100])/RR500(q[99])):.15f}")

# Test: q_n / μ_Poincaré^n should converge to a constant if σ=0
print("\n=== Test σ=0: q_n / μ^n should converge ===")
print(f"{'n':>5} | {'q_n/μ^n':>30} | {'ratio(n)/ratio(n-1)':>20}")
prev_ratio = None
for nn in [10, 20, 30, 50, 80, 100, 150, 200, 250, 300]:
    if nn >= len(q):
        break
    ratio = RR500(q[nn]) / mu_poincare**nn
    if prev_ratio is not None and prev_ratio != 0:
        change = ratio / prev_ratio
        print(f"{nn:5d} | {float(ratio):30.15e} | {float(change):20.15f}")
    else:
        print(f"{nn:5d} | {float(ratio):30.15e} | ---")
    prev_ratio = ratio

# If σ≠0, q_n/μ^n ~ C·n^{-σ}, so log(ratio_n/ratio_m) / log(n/m) → -σ
print("\n=== Extract effective σ from consecutive ratios ===")
ratios = []
for nn in range(50, 300):
    if nn >= len(q):
        break
    ratios.append((nn, RR500(q[nn]) / mu_poincare**nn))

for i in range(0, len(ratios)-50, 50):
    n1, r1 = ratios[i]
    n2, r2 = ratios[i+50]
    if r1 != 0 and r2 != 0:
        effective_sigma = -RR500(log(abs(r2/r1))) / RR500(log(RR500(n2)/RR500(n1)))
        print(f"  n={n1}→{n2}: effective σ = {float(effective_sigma):.10f}")

# Compute L = ζ(2) + ζ(3)
from mpmath import mp, zeta as mpzeta
mp.dps = 200
L_mp = mpzeta(2) + mpzeta(3)
L = RR500(L_mp)

# Compute the error ratio more carefully
print(f"\n=== Error e_n = p_n - L*q_n asymptotic ===")
mu_sub = min(roots, key=lambda z: abs(z))  # subdominant root
mu_sub_poincare = mu_sub / 64
print(f"|μ±/64| = {float(abs(mu_sub_poincare)):.15f}")

for nn in [50, 100, 150, 200, 250]:
    if nn >= len(q):
        break
    e_n = RR500(p[nn]) - L * RR500(q[nn])
    rate = abs(e_n)**(QQ(1)/nn)
    # e_n / (μ_sub_poincare^n) should converge if σ=0 for subdominant too
    print(f"  n={nn}: |e_n|^(1/n) = {float(rate):.12f}, target |μ±/64| = {float(abs(mu_sub_poincare)):.12f}")

# Direct ratio p_n/q_n convergence
print(f"\n=== p_n/q_n convergence to L ===")
for nn in [10, 50, 100, 200, 300]:
    if nn >= len(q):
        break
    if q[nn] != 0:
        ratio = RR500(p[nn]) / RR500(q[nn])
        diff = ratio - L
        if diff != 0:
            digits = -float(log(abs(diff), 10))
            print(f"  n={nn}: p_n/q_n - L = {float(diff):.5e}, matching digits ≈ {digits:.0f}")

print("\nDone.")
