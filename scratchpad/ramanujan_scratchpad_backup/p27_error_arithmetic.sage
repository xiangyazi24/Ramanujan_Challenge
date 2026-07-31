#!/usr/bin/env sage
"""
Study the arithmetic structure of e_n = p_n - (ζ(2)+ζ(3))·q_n.
Key questions:
1. What are the denominators of p_n and q_n?
2. What is the growth of |e_n| exactly?
3. Can we identify a "rational part" of e_n after clearing denominators?
4. Does e_n have a recognizable hypergeometric structure?
"""
from sage.all import *

# Recurrence coefficients (raw form, not monic)
def A_raw(n): return 1024*(2*n+5)**4*(2*n+7)**3*(2*n+9)**3*(946*n**2+6407*n+10860)
def B_raw(n): return 128*(2*n+7)**3*(2*n+9)**3*(104060*n**6+1745370*n**5+12145238*n**4+44886481*n**3+92943995*n**2+102256019*n+46709052)
def C_raw(n): return 16*(n+3)**4*(2*n+9)**3*(3784*n**5+57792*n**4+351019*n**3+1059230*n**2+1587211*n+944620)
def D_raw(n): return (n+3)**4*(n+4)**6*(946*n**2+4515*n+5399)

# Monic recurrence: u_{n+1} = (B_n/A_n) u_n - (C_{n-1}/A_{n-1}) u_{n-1} + (D_{n-2}/A_{n-2}) u_{n-2}
# for n >= 2.

# Initial conditions from proof.tex
p = [QQ(-612218384750), QQ(-9525021973931919)/QQ(18100), QQ(-29561828382772029)/QQ(65380)]
q = [QQ(-215040420000), QQ(-167282265043404)/QQ(905), QQ(-964185327658080)/QQ(6071)]

# Extend sequences using recurrence
N_MAX = 100
for _ in range(N_MAX - 2):
    n = len(p) - 1  # current last index
    p_new = B_raw(n)/A_raw(n) * p[n] - C_raw(n-1)/A_raw(n-1) * p[n-1] + D_raw(n-2)/A_raw(n-2) * p[n-2]
    q_new = B_raw(n)/A_raw(n) * q[n] - C_raw(n-1)/A_raw(n-1) * q[n-1] + D_raw(n-2)/A_raw(n-2) * q[n-2]
    p.append(QQ(p_new))
    q.append(QQ(q_new))

print(f"Computed {len(p)} terms of p_n and q_n")

# Verify q_n matches AESZ formula for first few terms
def inner_sum(n):
    return sum(binomial(n,k)**2 * binomial(n+k,n) * binomial(n+2*k,n) for k in range(n+1))

print("\n--- Verify q_n = C(2n,n) * inner_sum(n) ---")
for nn in range(10):
    q_formula = QQ(binomial(2*nn, nn)) * inner_sum(nn)
    if q[nn] == q_formula:
        pass  # OK
    else:
        print(f"  n={nn}: MISMATCH! q_rec = {q[nn]}, q_formula = {q_formula}")
        break
else:
    print("  First 10 terms match!")

# Study denominators
print("\n--- Denominator structure ---")
print(f"{'n':>4} | {'den(q_n)':>20} | {'den(p_n)':>20} | {'lcm':>20}")
for nn in [0, 1, 2, 3, 4, 5, 10, 20, 30, 50]:
    dq = q[nn].denominator()
    dp = p[nn].denominator()
    lc = lcm(dp, dq)
    print(f"{nn:4d} | {dq:20d} | {dp:20d} | {lc:20d}")

# Factor the denominators
print("\n--- Prime factorization of denominators ---")
for nn in [0, 1, 2, 3, 4, 5]:
    dq = q[nn].denominator()
    dp = p[nn].denominator()
    print(f"  q_{nn} den = {factor(dq)}")
    print(f"  p_{nn} den = {factor(dp)}")

# Compute e_n = p_n - L * q_n at high precision
PREC = 2000
RR = RealField(PREC)
L = RR.pi()**2/6 + RR(sum(QQ(1)/QQ(k**3) for k in range(1, 500)))
# Better: use mpmath for zeta
from mpmath import mp, zeta as mpzeta
mp.dps = 700
L_mp = mpzeta(2) + mpzeta(3)
L = RR(L_mp)

print(f"\nL = ζ(2)+ζ(3) = {str(L)[:50]}...")

# Compute e_n in high precision
print("\n--- Error sequence structure ---")
print(f"{'n':>4} | {'log10|e_n|':>15} | {'|e_n|^(1/n)':>15} | {'den(q_n)*e_n approx integer?':>30}")
for nn in [0, 1, 2, 3, 5, 10, 20, 30, 50, 80, 100]:
    e_n = RR(p[nn]) - L * RR(q[nn])
    if e_n != 0:
        log_e = float(log(abs(e_n), 10))
        if nn > 0:
            rate = float(abs(e_n)**(QQ(1)/nn))
        else:
            rate = float('inf')
        # Check if den(q_n) * e_n is close to integer
        dq = q[nn].denominator()
        dp = p[nn].denominator()
        lc = lcm(dp, dq)
        scaled = lc * e_n
        # How close to integer?
        dist = abs(scaled - round(float(scaled)))
        print(f"{nn:4d} | {log_e:15.3f} | {rate:15.10f} | dist={float(dist):.3e}")

# Compute v_n = g(n) * q_n (gauged sequence)
print("\n--- Gauged sequence v_n = g(n)*q_n ---")
g = lambda n: QQ(n)**2 + QQ(105)/22*QQ(n) + QQ(5399)/946
for nn in range(8):
    v_n = g(nn) * q[nn]
    print(f"  v_{nn} = {v_n}")
    print(f"    den(v_{nn}) = {factor(v_n.denominator())}")

# Study the dominant root ratio q_{n+1}/q_n
print("\n--- Dominant ratio convergence ---")
R_poly = PolynomialRing(QQ, 'x')
x = R_poly.gen()
P_cubic = 4*x**3 - 220*x**2 + 8*x - 1
CC = ComplexField(200)
roots = P_cubic.change_ring(CC).roots(multiplicities=False)
mu0 = max(roots, key=lambda z: z.real())
print(f"μ₀ = {float(mu0.real()):.15f}")
print(f"|μ±| = {float(abs(roots[1])):.15f}")
sigma0 = 24*(4*mu0 - 1)/(220*mu0**2 - 16*mu0 + 3)
print(f"σ₀ = {float(sigma0.real()):.15f}")

for nn in [10, 20, 50, 100]:
    if q[nn] != 0 and q[nn-1] != 0:
        ratio = RR(q[nn]) / RR(q[nn-1])
        err = abs(ratio - RR(mu0.real()))
        print(f"  n={nn}: q_n/q_{{n-1}} - μ₀ = {float(err):.3e}")

# Check: is e_n / (μ±^n * n^{-σ±}) converging?
mu_sub = roots[1]  # one of the subdominant complex roots
sigma_sub = 24*(4*mu_sub - 1)/(220*mu_sub**2 - 16*mu_sub + 3)
print(f"\nSubdominant: |μ±| = {float(abs(mu_sub)):.10f}, σ± = {CC(sigma_sub)}")

print("\n--- Subdominant coefficient extraction ---")
for nn in [20, 30, 50, 80, 100]:
    e_n = RR(p[nn]) - L * RR(q[nn])
    if e_n != 0:
        # e_n ~ c_+ μ_+^n n^{-σ_+} + c_- μ_-^n n^{-σ_-}
        # |e_n| ~ 2|c_+| |μ_+|^n n^{-Re(σ_+)} cos(...)
        rate_n = float(abs(e_n)**(QQ(1)/nn))
        expected = float(abs(mu_sub))
        print(f"  n={nn}: |e_n|^(1/n) = {rate_n:.10f}, expected |μ±| = {expected:.10f}")

print("\nDone.")
