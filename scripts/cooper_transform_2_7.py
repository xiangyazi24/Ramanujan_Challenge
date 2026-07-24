#!/usr/bin/env python3
"""
Problem 2.7: Cooper level-11 → Problem 2.7 via binomial transform.

Cooper's recurrence (level 11):
  (n+1)³T_{n+1} = 2(2n+1)(5n²+5n+2)T_n - 8n(7n²+1)T_{n-1} + 22n(2n-1)(n-1)T_{n-2}
  T_0=1, T_1=4, T_2=28, T_3=268

The shifted-binomial/even-section transform (from Q4793):
  (KT)_n = (1/256^n) Σ_{k=0}^{2n} C(2n,k)(-2)^{2n-k} T_k

Problem 2.7 recurrence:
  A_n u_{n+1} = B_n u_n - C_{n-1} u_{n-1} + D_{n-2} u_{n-2}

where:
  A_n = 1024(2n+5)⁴(2n+7)³(2n+9)³(946n²+6407n+10860)
  B_n = 128(2n+7)³(2n+9)³(104060n⁶+1745370n⁵+12145238n⁴+44886481n³+92943995n²+102256019n+46709052)
  C_n = 16(n+3)⁴(2n+9)³(3784n⁵+57792n⁴+351019n³+1059230n²+1587211n+944620)
  D_n = (n+3)⁴(n+4)⁶(946n²+4515n+5399)

Initial conditions:
  p_0 = -612218384750
  p_1 = -9525021973931919/18100
  p_2 = -29561828382772029/65380

  q_0 = -215040420000
  q_1 = -1672822650043404/905
  q_2 = -964185327658080/6071

Prove: lim p_n/q_n = ζ(2)+ζ(3)
"""
from fractions import Fraction as F
from mpmath import mp, mpf, zeta, pi
import math

mp.dps = 80

# Cooper sequence
def cooper_seq(N_max):
    T = [F(1), F(4)]
    for n in range(1, N_max):
        t3 = 2*(2*n+1)*(5*n**2+5*n+2)*T[n] - 8*n*(7*n**2+1)*T[n-1]
        if n >= 2:
            t3 += 22*n*(2*n-1)*(n-1)*T[n-2]
        T.append(t3 / F((n+1)**3))
    return T

# Compute Cooper sequence
N_T = 200
T = cooper_seq(N_T)
print("Cooper T_0..T_5:", [T[n] for n in range(6)])

# Binomial transform
def binomial_transform(T, n):
    """(KT)_n = (1/256^n) Σ_{k=0}^{2n} C(2n,k)(-2)^{2n-k} T_k"""
    result = F(0)
    for k in range(2*n + 1):
        binom = F(math.comb(2*n, k))
        result += binom * F(-2)**(2*n - k) * T[k]
    return result / F(256)**n

print("\nBinomial transform (KT)_n:")
KT = []
for n in range(30):
    if 2*n < N_T:
        kt_n = binomial_transform(T, n)
        KT.append(kt_n)
        if n < 10:
            print(f"  (KT)_{n} = {kt_n} = {float(kt_n):.10f}")

# Problem 2.7 recurrence coefficients
def A(n):
    return 1024*(2*n+5)**4*(2*n+7)**3*(2*n+9)**3*(946*n**2+6407*n+10860)

def B(n):
    P6 = 104060*n**6+1745370*n**5+12145238*n**4+44886481*n**3+92943995*n**2+102256019*n+46709052
    return 128*(2*n+7)**3*(2*n+9)**3*P6

def C(n):
    P5 = 3784*n**5+57792*n**4+351019*n**3+1059230*n**2+1587211*n+944620
    return 16*(n+3)**4*(2*n+9)**3*P5

def D(n):
    return (n+3)**4*(n+4)**6*(946*n**2+4515*n+5399)

# Problem 2.7 initial conditions
p = [F(-612218384750), F(-9525021973931919, 18100), F(-29561828382772029, 65380)]
q = [F(-215040420000), F(-1672822650043404, 905), F(-964185327658080, 6071)]

# Generate more p,q using the recurrence
# u_{n+1} = (B_n/A_n) u_n - (C_{n-1}/A_{n-1}) u_{n-1} + (D_{n-2}/A_{n-2}) u_{n-2}
for n in range(2, 80):
    p_next = F(B(n)) * p[n] / F(A(n)) - F(C(n-1)) * p[n-1] / F(A(n-1)) + F(D(n-2)) * p[n-2] / F(A(n-2))
    q_next = F(B(n)) * q[n] / F(A(n)) - F(C(n-1)) * q[n-1] / F(A(n-1)) + F(D(n-2)) * q[n-2] / F(A(n-2))
    p.append(p_next)
    q.append(q_next)

# Verify convergence
target = float(zeta(2) + zeta(3))
print(f"\nζ(2)+ζ(3) = {target:.15f}")
print(f"ζ(2) = {float(zeta(2)):.15f}")
print(f"ζ(3) = {float(zeta(3)):.15f}")

print("\nConvergence p_n/q_n → ζ(2)+ζ(3):")
for n in [5, 10, 20, 30, 40, 50]:
    ratio = float(p[n]) / float(q[n])
    err = abs(ratio - target)
    if err > 0:
        digits = -math.log10(err)
    else:
        digits = float('inf')
    print(f"  n={n}: p/q = {ratio:.15f}, digits = {digits:.1f}")

# Check if KT matches q_n (or some rescaling)
print("\n=== Matching KT with Problem 2.7 denominators ===")
if len(KT) > 0 and KT[0] != 0:
    for n in range(min(8, len(KT))):
        if KT[n] != 0 and q[n] != 0:
            ratio_kq = q[n] / KT[n]
            print(f"  n={n}: q_n/KT_n = {float(ratio_kq):.6e}")

# Also check if KT_n satisfies the Problem 2.7 recurrence directly
print("\n=== Checking if KT satisfies Problem 2.7 recurrence ===")
for n in range(2, min(15, len(KT)-1)):
    # A_n * KT_{n+1} - B_n * KT_n + C_{n-1} * KT_{n-1} - D_{n-2} * KT_{n-2}
    residual = A(n) * KT[n+1] - B(n) * KT[n] + C(n-1) * KT[n-1] - D(n-2) * KT[n-2]
    print(f"  n={n}: residual = {float(residual):.6e}")

# Poincaré analysis
print("\n=== Poincaré roots ===")
# All four coefficients have degree 12
# Leading terms: A ~ 1024·2⁴·2³·2³·946 · n^{12} = 1024·16·8·8·946·n^{12}
A_lead = 1024 * 16 * 8 * 8 * 946
B_lead = 128 * 8 * 8 * 104060
C_lead = 16 * 1 * 8 * 3784  # (n+3)⁴ ~ n⁴, (2n+9)³ ~ 8n³, leading of P5 = 3784n⁵ → deg = 4+3+5=12
D_lead = 1 * 1 * 946  # (n+3)⁴·(n+4)⁶·(946n²) → deg = 4+6+2=12
print(f"  A leading: {A_lead}")
print(f"  B leading: {B_lead}")
print(f"  C leading: {C_lead}")
print(f"  D leading: {D_lead}")

# Poincaré polynomial: A_lead · c³ - B_lead · c² + C_lead · c - D_lead = 0
import numpy as np
coeffs = [A_lead, -B_lead, C_lead, -D_lead]
roots = np.roots(coeffs)
print(f"  Poincaré roots: {roots}")
print(f"  |roots|: {[abs(r) for r in roots]}")

# Convergence rate
dominant = max(abs(r) for r in roots)
subdominant = sorted(abs(r) for r in roots)[1]
rate = subdominant / dominant
print(f"  Convergence rate: {rate:.6f}")
print(f"  Digits/step: {-math.log10(rate):.3f}")
