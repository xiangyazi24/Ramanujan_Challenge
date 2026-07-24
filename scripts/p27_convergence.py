#!/usr/bin/env python3
"""P2.7: Compute many terms and verify convergence to ζ(2)+ζ(3).
Also compute the MOP step-line recurrence coefficients and check
if they match P2.7 after a pullback.
"""
import mpmath as mp
mp.mp.dps = 250

def A(n):
    return 1024*(2*n+5)**4*(2*n+7)**3*(2*n+9)**3*(946*n**2+6407*n+10860)

def B(n):
    return 128*(2*n+7)**3*(2*n+9)**3*(104060*n**6+1745370*n**5+12145238*n**4+44886481*n**3+92943995*n**2+102256019*n+46709052)

def C(n):
    return 16*(n+3)**4*(2*n+9)**3*(3784*n**5+57792*n**4+351019*n**3+1059230*n**2+1587211*n+944620)

def D(n):
    return (n+3)**4*(n+4)**6*(946*n**2+4515*n+5399)

# Initial values
p = [mp.mpf(-612218384750),
     mp.mpf(-9525021973931919)/mp.mpf(18100),
     mp.mpf(-29561828382772029)/mp.mpf(65380)]

q = [mp.mpf(-215040420000),
     mp.mpf(-167282265043404)/mp.mpf(905),
     mp.mpf(-964185327658080)/mp.mpf(6071)]

L = mp.zeta(2) + mp.zeta(3)
print(f"Target: ζ(2)+ζ(3) = {mp.nstr(L, 60)}")

# Recurrence: u_{n+1} = (B_n/A_n) u_n - (C_{n-1}/A_{n-1}) u_{n-1} + (D_{n-2}/A_{n-2}) u_{n-2}
# n >= 2
NMAX = 80
for n in range(2, NMAX):
    pnew = mp.mpf(B(n))/mp.mpf(A(n)) * p[-1] \
         - mp.mpf(C(n-1))/mp.mpf(A(n-1)) * p[-2] \
         + mp.mpf(D(n-2))/mp.mpf(A(n-2)) * p[-3]
    qnew = mp.mpf(B(n))/mp.mpf(A(n)) * q[-1] \
         - mp.mpf(C(n-1))/mp.mpf(A(n-1)) * q[-2] \
         + mp.mpf(D(n-2))/mp.mpf(A(n-2)) * q[-3]
    p.append(pnew)
    q.append(qnew)

print(f"\nComputed {len(p)} terms.")

# Convergence analysis
print(f"\n{'n':>4}  {'p_n/q_n - L':>25}  {'digits':>8}  {'ratio':>15}")
prev_err = None
for n in range(len(p)):
    ratio = p[n] / q[n]
    err = abs(ratio - L)
    digits = -mp.log10(err) if err > 0 else 999
    r_str = ""
    if prev_err and prev_err > 0 and err > 0:
        r = mp.log10(err) / mp.log10(prev_err)
        r_str = mp.nstr(r, 6)
    if n < 5 or n % 10 == 0 or n == len(p)-1:
        print(f"{n:4d}  {mp.nstr(err, 12):>25s}  {mp.nstr(digits, 6):>8s}  {r_str:>15s}")
    prev_err = err

# Check Poincaré roots
print(f"\n=== Poincaré analysis ===")
print(f"q_n/q_{'{n-1}'} for large n:")
for n in [30, 40, 50, 60, 70]:
    if n < len(q):
        r = q[n] / q[n-1]
        print(f"  n={n}: q_n/q_{{n-1}} = {mp.nstr(r, 20)}")

# Dominant root μ₀/64
import mpmath
roots = mpmath.polyroots([mp.mpf(-1), mp.mpf(8), mp.mpf(-220), mp.mpf(4)])
print(f"\nPoincaré roots of 4μ³-220μ²+8μ-1:")
for r in roots:
    print(f"  μ = {mp.nstr(r, 20)}, |μ| = {mp.nstr(abs(r), 20)}")

# Expected dominant ratio: μ₀
mu0 = max(roots, key=abs)
print(f"\nExpected q_n/q_{{n-1}} → μ₀ = {mp.nstr(mu0, 20)}")
print(f"Convergence rate: |μ_±/μ₀| = {mp.nstr(abs(roots[1])/abs(mu0), 20)}")

# Verify last ratio matches L
print(f"\n=== Final verification ===")
n_final = len(p) - 1
ratio_final = p[n_final] / q[n_final]
diff = ratio_final - L
digits_final = -mp.log10(abs(diff)) if diff != 0 else 999
print(f"p_{n_final}/q_{n_final} = {mp.nstr(ratio_final, 60)}")
print(f"ζ(2)+ζ(3)      = {mp.nstr(L, 60)}")
print(f"Digits of agreement: {mp.nstr(digits_final, 10)}")
