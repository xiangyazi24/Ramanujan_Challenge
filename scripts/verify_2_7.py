#!/usr/bin/env python3
"""Problem 2.7: Verify convergence to ζ(2)+ζ(3) with arbitrary precision."""
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

p = [mpf(-612218384750), mpf(-9525021973931919)/mpf(18100), mpf(-29561828382772029)/mpf(65380)]
q = [mpf(-215040420000), mpf(-1672822650043404)/mpf(905), mpf(-964185327658080)/mpf(6071)]

target = zeta(2) + zeta(3)
print(f"ζ(2)+ζ(3) = {mp.nstr(target, 50)}")

for n in range(2, 120):
    p_next = B(n)*p[n]/A(n) - C(n-1)*p[n-1]/A(n-1) + D(n-2)*p[n-2]/A(n-2)
    q_next = B(n)*q[n]/A(n) - C(n-1)*q[n-1]/A(n-1) + D(n-2)*q[n-2]/A(n-2)
    p.append(p_next)
    q.append(q_next)

for n in [5, 10, 20, 30, 40, 50, 60, 80, 100]:
    ratio = p[n]/q[n]
    err = fabs(ratio - target)
    if err > 0:
        digits = -log10(err)
        print(f"  n={n}: {mp.nstr(ratio, 30)}, {mp.nstr(digits, 5)} digits")
    else:
        print(f"  n={n}: exact match!")
