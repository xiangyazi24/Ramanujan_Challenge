#!/usr/bin/env python3
"""P2.5: Compare the neutral Legendre integral N_n from Q4869
with the actual CMF error ê_n = G·Q̂_n - P̂_n.

If they are proportional, the proportionality constant closes the proof.
"""
import mpmath as mp
from fractions import Fraction

mp.mp.dps = 80
sqrt2 = mp.sqrt(2)
Lam = 3 + 2*sqrt2
rho = Lam**(-2)
G = mp.catalan

# --- Neutral integral N_n (Q4869 formula) ---
def legendre_P_vals(N):
    P = [mp.mpf(1)]
    if N == 0: return P
    P.append(mp.mpf(3))
    for n in range(1, N):
        P.append((3*(2*n+1)*P[n] - n*P[n-1])/(n+1))
    return P

def legendre_Q_vals(N):
    Q = [mp.log(2)/2]
    if N == 0: return Q
    Q.append(3*Q[0] - 1)
    for n in range(1, N):
        Q.append((3*(2*n+1)*Q[n] - n*Q[n-1])/(n+1))
    return Q

NMAX = 25
Pvals = legendre_P_vals(NMAX + 2)
Qvals = legendre_Q_vals(NMAX + 2)

def neutral_n(n):
    return (Pvals[n] - Pvals[n+1]/Lam) * (Qvals[n] - Lam*Qvals[n+1])

# --- CMF error ê_n = G·Q̂_n - P̂_n ---
def M_entries(n):
    m11 = (-2*n-5)*(n+3)**2*(136*n**4+1424*n**3+5548*n**2+9551*n+6141)
    m12 = 384*n**6+6384*n**5+44168*n**4+162698*n**3+336377*n**2+369933*n+169011
    m13 = -480*n**4-4980*n**3-19210*n**2-32690*n-20730
    m21 = (n+2)**2*(n+3)**2*(4*n+10)*(48*n**3+386*n**2+1017*n+879)
    m22 = (n+2)**2*(-272*n**5-3848*n**4-21732*n**3-61184*n**2-85761*n-47808)
    m23 = (n+2)**2*(320*n**3+2540*n**2+6610*n+5640)
    m31 = (-4*n-10)*(n+2)**2*(n+3)**2*(32*n**4+302*n**3+1037*n**2+1530*n+813)
    m32 = (n+2)**2*(192*n**6+2984*n**5+19116*n**4+64452*n**3+120256*n**2+117279*n+46476)
    m33 = (n+2)**2*(-16*n**5-408*n**4-2912*n**3-8884*n**2-12254*n-6240)
    return [[mp.mpf(m11),mp.mpf(m12),mp.mpf(m13)],
            [mp.mpf(m21),mp.mpf(m22),mp.mpf(m23)],
            [mp.mpf(m31),mp.mpf(m32),mp.mpf(m33)]]

def H(nn):
    if nn == 0: return mp.mpf(1)
    val = mp.mpf(1)
    for k in range(nn):
        val *= mp.mpf(-16) * mp.mpf(k+2)**2 * mp.mpf(k+3)**2 * mp.mpf(2*k+5)/2 * (mp.mpf(2*k+7)/2)**2
    return val

# Compute Q̂_n, P̂_n
row_q = [mp.mpf(33750), mp.mpf(-36000), mp.mpf(9000)]
row_p = [mp.mpf(30921), mp.mpf(-32972), mp.mpf(8240)]

qhat = [row_q[0]]
phat = [row_p[0]]

for n in range(NMAX):
    M = M_entries(n)
    new_q = [sum(row_q[k]*M[k][j] for k in range(3)) for j in range(3)]
    new_p = [sum(row_p[k]*M[k][j] for k in range(3)) for j in range(3)]
    row_q = new_q
    row_p = new_p
    h = H(n+1)
    qhat.append(row_q[0] / h)
    phat.append(row_p[0] / h)

print("=== Comparing N_n (neutral Legendre) with ê_n (CMF error) ===")
print(f"{'n':>3}  {'N_n':>20}  {'ê_n':>20}  {'ê_n/N_n':>25}  {'n^3·ê_n':>20}")
for n in range(NMAX + 1):
    Nn = neutral_n(n)
    en = G * qhat[n] - phat[n]

    ratio_str = ""
    if abs(Nn) > 1e-50:
        ratio = en / Nn
        ratio_str = mp.nstr(ratio, 15)

    n3en = (n**3 * en) if n > 0 else en
    print(f"{n:3d}  {mp.nstr(Nn, 12):>20s}  {mp.nstr(en, 12):>20s}  {ratio_str:>25s}  {mp.nstr(n3en, 12):>20s}")

# Check if the ratio converges
print("\n=== Ratio ê_n / N_n for larger n ===")
for n in [5, 10, 15, 20, 25]:
    if n <= NMAX:
        Nn = neutral_n(n)
        en = G * qhat[n] - phat[n]
        if abs(Nn) > 1e-50:
            ratio = en / Nn
            print(f"  n={n}: ê_n/N_n = {mp.nstr(ratio, 30)}")

# Also check: n^3 * ê_n limit
print("\n=== n^3 · ê_n (should approach C_0) ===")
for n in range(1, NMAX + 1):
    en = G * qhat[n] - phat[n]
    val = n**3 * en
    print(f"  n={n:2d}: n^3·ê_n = {mp.nstr(val, 20)}")

# Predicted limit: 1/(16*sqrt2) for N_n
print(f"\nPredicted n^3 limit of N_n: {mp.nstr(1/(16*sqrt2), 20)}")
