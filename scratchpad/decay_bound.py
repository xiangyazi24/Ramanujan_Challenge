"""
Compute e_n = p_n - (zeta(2)+zeta(3))·q_n for the P2.7 recurrence
and verify the decay rate.

If c_0(e) = 0, then |e_n|^{1/n} -> |mu_pm| ≈ 0.00105
If c_0(e) != 0, then |e_n|^{1/n} -> mu_0 ≈ 0.859

This computation uses the CORRECT monic convention from proof.tex.
"""
from mpmath import mp, mpf, polylog, zeta, log, fabs, power, pi

mp.dps = 600  # 600 digits of precision

# P2.7 coefficient polynomials (exact integer arithmetic)
def A(n):
    return (1024*(2*n+5)**4*(2*n+7)**3*(2*n+9)**3
            *(946*n**2 + 6407*n + 10860))

def B(n):
    return (128*(2*n+7)**3*(2*n+9)**3
            *(104060*n**6 + 1745370*n**5 + 12145238*n**4
              + 44886481*n**3 + 92943995*n**2
              + 102256019*n + 46709052))

def C(n):
    return (16*(n+3)**4*(2*n+9)**3
            *(3784*n**5 + 57792*n**4 + 351019*n**3
              + 1059230*n**2 + 1587211*n + 944620))

def D(n):
    return ((n+3)**4*(n+4)**6
            *(946*n**2 + 4515*n + 5399))

# Correct monic recurrence (proof.tex: each coeff / own A)
def alpha(n): return mpf(B(n)) / mpf(A(n))
def beta(n): return -mpf(C(n-1)) / mpf(A(n-1))
def gamma(n): return mpf(D(n-2)) / mpf(A(n-2))

# Initial values (exact rationals from proof.tex)
p0 = mpf(-612218384750)
p1 = mpf(-9525021973931919) / mpf(18100)
p2 = mpf(-29561828382772029) / mpf(65380)

q0 = mpf(-215040420000)
q1 = mpf(-167282265043404) / mpf(905)
q2 = mpf(-964185327658080) / mpf(6071)

L = zeta(2) + zeta(3)
print(f"L = zeta(2)+zeta(3) = {mp.nstr(L, 50)}")

# Compute e_n = p_n - L*q_n
e0 = p0 - L*q0
e1 = p1 - L*q1
e2 = p2 - L*q2

print(f"\ne[0] = {mp.nstr(e0, 30)}")
print(f"e[1] = {mp.nstr(e1, 30)}")
print(f"e[2] = {mp.nstr(e2, 30)}")

# Compute forward using the recurrence
N = 300
p = [mpf(0)]*(N+3)
q = [mpf(0)]*(N+3)

p[0], p[1], p[2] = p0, p1, p2
q[0], q[1], q[2] = q0, q1, q2

for n in range(2, N+2):
    p[n+1] = alpha(n)*p[n] + beta(n)*p[n-1] + gamma(n)*p[n-2]
    q[n+1] = alpha(n)*q[n] + beta(n)*q[n-1] + gamma(n)*q[n-2]

# Compute e_n and check decay rate
print("\n--- Decay analysis ---")
print(f"{'n':>5}  {'|e_n|':>20}  {'|e_n|^(1/n)':>20}  {'log|e_n|/n':>20}")
print("-" * 85)

decay_rates = []
for n in [5, 10, 20, 30, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300]:
    if n >= N+3:
        break
    en = p[n] - L*q[n]
    abs_en = fabs(en)
    if abs_en > 0:
        rate = power(abs_en, mpf(1)/n)
        log_rate = log(abs_en) / n
        print(f"{n:5d}  {mp.nstr(abs_en, 15):>20s}  {mp.nstr(rate, 15):>20s}  {mp.nstr(log_rate, 15):>20s}")
        if n >= 50:
            decay_rates.append((n, float(rate)))
    else:
        print(f"{n:5d}  {'ZERO':>20s}")

# Poincare roots for comparison
print("\n--- Poincare roots ---")
# Compute from the characteristic equation
# mu^3 = alpha_inf * mu^2 + beta_inf * mu + gamma_inf
# where alpha_inf = lim B(n)/A(n), etc.

# Leading coefficients
# A(n) ~ 1024 * 16 * 8 * 8 * 946 * n^12 = 1024 * 1024 * 946 * n^12
A_lead = 1024 * 16 * 8 * 8 * 946
# B(n) ~ 128 * 8 * 8 * 104060 * n^12
B_lead = 128 * 8 * 8 * 104060
# C(n) ~ 16 * 1 * 8 * 3784 * n^12
C_lead = 16 * 1 * 8 * 3784
# D(n) ~ 1 * 1 * 946 * n^12
D_lead = 946

alpha_inf = mpf(B_lead) / mpf(A_lead)
beta_inf = -mpf(C_lead) / mpf(A_lead)
gamma_inf = mpf(D_lead) / mpf(A_lead)

print(f"alpha_inf = {mp.nstr(alpha_inf, 30)}")
print(f"beta_inf  = {mp.nstr(beta_inf, 30)}")
print(f"gamma_inf = {mp.nstr(gamma_inf, 30)}")

# Solve mu^3 - alpha_inf*mu^2 - beta_inf*mu - gamma_inf = 0
from mpmath import polyroots
roots = polyroots([1, -alpha_inf, -beta_inf, -gamma_inf])
print("\nPoincare roots:")
for i, r in enumerate(roots):
    print(f"  mu_{i} = {mp.nstr(r, 20)}, |mu_{i}| = {mp.nstr(fabs(r), 20)}")

mu0 = max(fabs(r) for r in roots)
print(f"\nDominant |mu_0| = {mp.nstr(mu0, 20)}")
print(f"Subdominant |mu_pm| = {mp.nstr(min(fabs(r) for r in roots), 20)}")

# Verdict
if decay_rates:
    avg_rate = sum(r for _, r in decay_rates) / len(decay_rates)
    print(f"\nAverage decay rate (n>=50): {avg_rate:.10f}")
    print(f"Dominant root mu_0:         {float(mu0):.10f}")
    if avg_rate < float(mu0) * 0.1:
        print("\n*** DECAY RATE << mu_0 ***")
        print("*** This strongly implies c_0(e) = 0 ***")
    elif avg_rate < float(mu0):
        print("\n*** DECAY RATE < mu_0 (but close) ***")
    else:
        print("\n*** DECAY RATE >= mu_0: c_0(e) might be nonzero ***")

# Also compute ratio p_n/q_n to verify convergence to L
print("\n--- Convergence of p_n/q_n ---")
for n in [10, 50, 100, 200, 300]:
    if n >= N+3:
        break
    if q[n] != 0:
        ratio = p[n] / q[n]
        err = fabs(ratio - L)
        print(f"n={n:3d}: |p_n/q_n - L| = {mp.nstr(err, 15)}")
