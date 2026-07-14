#!/usr/bin/env python3
"""P2.7: Rigorous verification for the Poincaré-Birkhoff proof.
Compute e_n = p_n - (ζ(2)+ζ(3))*q_n to high precision,
verify it's recessive (two-step ratio → -|r_±|²)."""
from mpmath import mp, mpf, zeta, pi, nstr, log10, sqrt, polyroots, fabs

mp.dps = 500

# Target constant
L = zeta(2) + zeta(3)
print(f"L = ζ(2)+ζ(3) = {nstr(L, 80)}")

# Poincaré polynomial 4μ³ - 220μ² + 8μ - 1
roots = polyroots([mpf(4), mpf(-220), mpf(8), mpf(-1)])
roots_sorted = sorted(roots, key=lambda r: -abs(r))
mu0 = roots_sorted[0]
mu_plus = roots_sorted[1]
mu_minus = roots_sorted[2]
print(f"\nPoincaré polynomial 4μ³-220μ²+8μ-1 roots:")
print(f"  μ₀ = {nstr(mu0, 40)} (real)")
print(f"  μ₊ = {nstr(mu_plus, 40)}")
print(f"  μ₋ = {nstr(mu_minus, 40)}")
print(f"  |μ₊| = |μ₋| = {nstr(abs(mu_plus), 40)}")

r0 = mu0 / 64
rp = mu_plus / 64
rm = mu_minus / 64
print(f"\nScaled Poincaré roots (÷64):")
print(f"  r₀ = {nstr(r0, 40)}")
print(f"  |r₊| = |r₋| = {nstr(abs(rp), 40)}")
print(f"  Convergence rate |r₊/r₀| = {nstr(abs(rp)/abs(r0), 20)}")
print(f"  Two-step ratio -|r₊|² = {nstr(-abs(rp)**2, 20)}")

# Recurrence coefficients
def A(n):
    return 1024*(2*n+5)**4*(2*n+7)**3*(2*n+9)**3*(946*n**2+6407*n+10860)
def B(n):
    return 128*(2*n+7)**3*(2*n+9)**3*(104060*n**6+1745370*n**5+12145238*n**4+44886481*n**3+92943995*n**2+102256019*n+46709052)
def C(n):
    return 16*(n+3)**4*(2*n+9)**3*(3784*n**5+57792*n**4+351019*n**3+1059230*n**2+1587211*n+944620)
def D(n):
    return (n+3)**4*(n+4)**6*(946*n**2+4515*n+5399)

# Compute p_n, q_n to high precision
NMAX = 150
p = [mpf(0)] * (NMAX + 1)
q = [mpf(0)] * (NMAX + 1)

# Initial values (exact rationals)
p[0] = mpf(-612218384750)
p[1] = mpf(-9525021973931919) / mpf(18100)
p[2] = mpf(-29561828382772029) / mpf(65380)
q[0] = mpf(-215040420000)
q[1] = mpf(-167282265043404) / mpf(905)
q[2] = mpf(-964185327658080) / mpf(6071)

print(f"\nComputing p_n, q_n for n=0..{NMAX}...", flush=True)
for n in range(2, NMAX):
    An = mpf(A(n)); Bn = mpf(B(n))
    Cn1 = mpf(C(n-1)); An1 = mpf(A(n-1))
    Dn2 = mpf(D(n-2)); An2 = mpf(A(n-2))
    p[n+1] = (Bn/An)*p[n] - (Cn1/An1)*p[n-1] + (Dn2/An2)*p[n-2]
    q[n+1] = (Bn/An)*q[n] - (Cn1/An1)*q[n-1] + (Dn2/An2)*q[n-2]

# Compute error e_n = p_n - L*q_n
print(f"\n=== Error analysis ===")
e = [p[n] - L * q[n] for n in range(NMAX + 1)]

# Verify q_n ≠ 0
print("\nq_n ≠ 0 check:")
all_nonzero = True
for n in range(NMAX + 1):
    if q[n] == 0:
        print(f"  WARNING: q_{n} = 0!")
        all_nonzero = False
if all_nonzero:
    print(f"  q_n ≠ 0 for all n=0..{NMAX} ✓")

# Convergence of p_n/q_n
print("\nConvergence of p_n/q_n to L:")
print(f"{'n':>5}  {'|p_n/q_n - L|':>20}  {'digits':>8}  {'q_{n+1}/q_n':>20}")
for n in [0, 1, 2, 5, 10, 20, 30, 50, 70, 100, 130]:
    if n <= NMAX and q[n] != 0:
        ratio = p[n] / q[n]
        err = fabs(ratio - L)
        digits = -float(log10(err)) if err > 0 else 999
        qratio = q[n]/q[n-1] if n > 0 and q[n-1] != 0 else mpf(0)
        print(f"{n:5d}  {nstr(err, 8):>20s}  {digits:8.1f}  {nstr(qratio.real, 15):>20s}")

# Two-step error ratios
print(f"\nTwo-step error ratios e_n/e_{{n-2}}:")
print(f"  Expected: -|r₊|² = {nstr(-abs(rp)**2, 20)}")
print(f"{'n':>5}  {'e_n/e_{n-2}':>25}  {'|e_n|':>15}")
for n in range(4, min(NMAX+1, 102), 2):
    if e[n-2] != 0:
        ratio = e[n] / e[n-2]
        if n <= 20 or n % 10 == 0:
            print(f"{n:5d}  {nstr(ratio, 20):>25s}  {nstr(fabs(e[n]), 6):>15s}")

# Digits of agreement between ratio and -|r_pm|^2
print(f"\nDigits of agreement: e_n/e_{{n-2}} vs -|r₊|²")
expected = -abs(rp)**2
for n in [10, 20, 30, 50, 70, 100]:
    if n <= NMAX and e[n-2] != 0:
        ratio = e[n] / e[n-2]
        diff = fabs(ratio - expected)
        if diff > 0:
            digs = -float(log10(diff))
        else:
            digs = 999
        print(f"  n={n:3d}: {digs:.1f} digits")

# Verify: e_n with exact initial values
print(f"\n=== Exact Poincaré root verification ===")
# Product relation: μ₀·μ₊·μ₋ = 1/4 (from Vieta)
print(f"  μ₀·μ₊·μ₋ = {nstr(mu0*mu_plus*mu_minus, 20)} (should be 1/4)")
# Sum relation: μ₀+μ₊+μ₋ = 55
print(f"  μ₀+μ₊+μ₋ = {nstr(mu0+mu_plus+mu_minus, 20)} (should be 55)")
# Pairwise sum: μ₀μ₊+μ₀μ₋+μ₊μ₋ = 2
print(f"  μ₀μ₊+μ₀μ₋+μ₊μ₋ = {nstr(mu0*mu_plus+mu0*mu_minus+mu_plus*mu_minus, 20)} (should be 2)")

# Total matched digits
n_final = NMAX
ratio_final = p[n_final] / q[n_final]
err_final = fabs(ratio_final - L)
digits_final = -float(log10(err_final)) if err_final > 0 else 999
print(f"\nFinal: p_{n_final}/q_{n_final} matches ζ(2)+ζ(3) to {digits_final:.0f} digits")

# Explicit computation of |r_pm|^2 for the two-step ratio
print(f"\n|r₊|² = |μ₊|²/64² = {nstr(abs(mu_plus)**2/4096, 30)}")
print(f"1/(2μ₀^{{3/2}}·64) = {nstr(1/(2*mu0**mpf(1.5)*64), 30)}")
print(f"|r₊| from Vieta = {nstr(1/(2*abs(mu0)**(mpf(3)/2)*64), 30)}")

print(f"\nDone. ({mp.dps} digits precision)")
