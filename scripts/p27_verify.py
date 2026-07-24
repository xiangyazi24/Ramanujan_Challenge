#!/usr/bin/env python3
"""P2.7: Verify the recurrence numerically with the large rational initial values.
Check convergence to ζ(2)+ζ(3)."""
from fractions import Fraction
from mpmath import mp, mpf, zeta, pi, nstr

mp.dps = 300

def A(n):
    return 1024*(2*n+5)**4*(2*n+7)**3*(2*n+9)**3*(946*n**2+6407*n+10860)

def B(n):
    return 128*(2*n+7)**3*(2*n+9)**3*(
        104060*n**6+1745370*n**5+12145238*n**4+44886481*n**3
        +92943995*n**2+102256019*n+46709052)

def C(n):
    return 16*(n+3)**4*(2*n+9)**3*(
        3784*n**5+57792*n**4+351019*n**3+1059230*n**2
        +1587211*n+944620)

def D(n):
    return (n+3)**4*(n+4)**6*(946*n**2+4515*n+5399)

# Large rational initial values from proof.tex
p = [None]*200
q = [None]*200
p[0] = Fraction(-612218384750)
p[1] = Fraction(-9525021973931919, 18100)
p[2] = Fraction(-29561828382772029, 65380)
q[0] = Fraction(-215040420000)
q[1] = Fraction(-167282265043404, 905)
q[2] = Fraction(-964185327658080, 6071)

print("Initial values:")
print(f"  p0 = {p[0]}")
print(f"  p1 = {float(p[1]):.6e}")
print(f"  p2 = {float(p[2]):.6e}")
print(f"  q0 = {q[0]}")
print(f"  q1 = {float(q[1]):.6e}")
print(f"  q2 = {float(q[2]):.6e}")

# Recurrence: u_{n+1} = (B_n/A_n) u_n - (C_{n-1}/A_{n-1}) u_{n-1} + (D_{n-2}/A_{n-2}) u_{n-2}
# for n >= 2
NMAX = 100
print(f"\nComputing p_n, q_n for n = 3..{NMAX}...")

for n in range(2, NMAX):
    An = A(n)
    Bn = B(n)
    Cn_1 = C(n-1)
    An_1 = A(n-1)
    Dn_2 = D(n-2)
    An_2 = A(n-2)

    p[n+1] = Fraction(Bn, An) * p[n] - Fraction(Cn_1, An_1) * p[n-1] + Fraction(Dn_2, An_2) * p[n-2]
    q[n+1] = Fraction(Bn, An) * q[n] - Fraction(Cn_1, An_1) * q[n-1] + Fraction(Dn_2, An_2) * q[n-2]

target = mpf(zeta(2)) + mpf(zeta(3))
print(f"\nTarget = ζ(2)+ζ(3) = {nstr(target, 50)}")

print(f"\nConvergence check:")
for n in [5, 10, 20, 30, 40, 50, 70]:
    if q[n] is not None and q[n] != 0:
        ratio = mpf(p[n].numerator) / mpf(p[n].denominator) / (mpf(q[n].numerator) / mpf(q[n].denominator))
        err = abs(ratio - target)
        if err > 0:
            digits = -float(mp.log10(err))
        else:
            digits = mp.dps
        print(f"  n={n:3d}: p_n/q_n - target ≈ {nstr(err, 6)}, digits ≈ {digits:.1f}")

# Poincaré analysis: check ratio q_{n+1}/q_n
print(f"\nPoincaré ratio q_{{n+1}}/q_n:")
for n in [10, 20, 30, 50, 70]:
    if q[n] is not None and q[n-1] is not None and q[n-1] != 0:
        r = mpf(q[n].numerator) / mpf(q[n].denominator) / (mpf(q[n-1].numerator) / mpf(q[n-1].denominator))
        print(f"  n={n:3d}: ratio = {nstr(r, 15)}")

# Check the Poincaré polynomial 4μ³-220μ²+8μ-1
# Dominant root μ₀ ≈ 54.964, so λ₀ = μ₀/64 ≈ 0.8588
print(f"\nPoincaré polynomial 4μ³-220μ²+8μ-1:")
from mpmath import polyroots
roots = polyroots([4, -220, 8, -1])
for i, r in enumerate(roots):
    print(f"  μ_{i} = {nstr(r, 20)}, λ_{i} = μ/64 = {nstr(r/64, 20)}")

# Error decay: check |p_n - target*q_n| / |p_{n-2} - target*q_{n-2}|
print(f"\nError decay (two-step ratio):")
errs = {}
for n in range(3, min(NMAX+1, 80)):
    if p[n] is not None and q[n] is not None:
        pn = mpf(p[n].numerator) / mpf(p[n].denominator)
        qn = mpf(q[n].numerator) / mpf(q[n].denominator)
        errs[n] = pn - target * qn

for n in [10, 20, 30, 40, 50, 60, 70]:
    if n in errs and n-2 in errs and errs[n-2] != 0:
        ratio = errs[n] / errs[n-2]
        print(f"  n={n:3d}: e_n/e_{{n-2}} = {nstr(ratio, 10)}")

print("\nDone.")
