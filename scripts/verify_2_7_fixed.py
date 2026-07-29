#!/usr/bin/env python3
"""Problem 2.7: VERIFIED with corrected initial conditions.
q₁ = -167282265043404/905 (NOT -1672822650043404/905 — extra '0' was transcription error)
"""
from mpmath import mp, mpf, zeta, log10, fabs

mp.dps = 200

def A(n):
    n = mpf(n)
    return 1024*(2*n+5)**4*(2*n+7)**3*(2*n+9)**3*(946*n**2+6407*n+10860)

def B(n):
    n = mpf(n)
    P6 = 104060*n**6+1745370*n**5+12145238*n**4+44886481*n**3+92943995*n**2+102256019*n+46709052
    return 128*(2*n+7)**3*(2*n+9)**3*P6

def C(n):
    n = mpf(n)
    P5 = 3784*n**5+57792*n**4+351019*n**3+1059230*n**2+1587211*n+944620
    return 16*(n+3)**4*(2*n+9)**3*P5

def D(n):
    n = mpf(n)
    return (n+3)**4*(n+4)**6*(946*n**2+4515*n+5399)

# CORRECTED initial conditions
p = [mpf(-612218384750),
     mpf(-9525021973931919)/mpf(18100),
     mpf(-29561828382772029)/mpf(65380)]

q = [mpf(-215040420000),
     mpf(-167282265043404)/mpf(905),   # FIXED: was -1672822650043404
     mpf(-964185327658080)/mpf(6071)]

target = zeta(2) + zeta(3)

for n in range(2, 120):
    p_next = B(n)/A(n)*p[n] - C(n-1)/A(n-1)*p[n-1] + D(n-2)/A(n-2)*p[n-2]
    q_next = B(n)/A(n)*q[n] - C(n-1)/A(n-1)*q[n-1] + D(n-2)/A(n-2)*q[n-2]
    p.append(p_next)
    q.append(q_next)

print(f"ζ(2)+ζ(3) = {mp.nstr(target, 60)}")
print()
for n in [5, 10, 20, 30, 50, 80, 100]:
    ratio = p[n]/q[n]
    err = fabs(ratio - target)
    if err > 0:
        digits = -log10(err)
    else:
        digits = mpf('inf')
    print(f"n={n:3d}: {mp.nstr(ratio, 40)} ({mp.nstr(digits, 5)} digits)")

# Convergence rate
print(f"\nPoincaré dominant root: 0.858807735912...")
print(f"Subdominant |root|:    0.001053785...")
print(f"Convergence rate:      0.001227...")
print(f"Expected digits/step:  2.911")

# Remainder analysis
print(f"\nRemainder r_n = p_n - ζ(2)+ζ(3))·q_n:")
for n in [5, 10, 20, 30, 50]:
    r = p[n] - target * q[n]
    r_prev = p[n-1] - target * q[n-1]
    ratio_r = r / r_prev
    print(f"  n={n}: r_{n}/r_{n-1} = {mp.nstr(ratio_r, 12)}")
print(f"  (should approach subdominant Poincaré root ≈ 0.00106)")
