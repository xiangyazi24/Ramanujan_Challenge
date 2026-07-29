#!/usr/bin/env python3
"""P2.5: Compute exact Q̂_n and D_n², check their ratio for polynomial/rational structure.
Also test if the neutral mode is a rational multiple of P_n(3)Q_n(3)."""
from fractions import Fraction
from mpmath import mp, mpf, log, nstr

mp.dps = 200

# Correct globally-normalized coefficients
c_coeffs = [
    [-170972650800, -826494925500, -1792449886332, -2317972607944, -2000297648936,
     -1219354055500, -541255279788, -177419351856, -43002662976, -7620091136,
     -960400960, -81589760, -4190208, -98304],
    [8781630505200, 38850314624124, 78557994908508, 96136040496551, 79442239242197,
     46814452218572, 20241514501104, 6502490145168, 1552168938336, 271943188864,
     33995217088, 2871763456, 146952192, 3440640],
    [-21132458248680, -87529225645944, -165451256319618, -189073879129764, -145809619841418,
     -80164318460172, -32338316008004, -9694892892592, -2160716677664, -353683596544,
     -41340724928, -3268370944, -156684288, -3440640],
    [587448626688, 2442715444224, 4635428285664, 5317694979920, 4116150568664,
     2270943978716, 919036676572, 276298241680, 61721801728, 10120470656,
     1184128064, 93632000, 4485120, 98304],
]

def eval_c(i, n):
    val = Fraction(0)
    nk = Fraction(1)
    for coeff in c_coeffs[i]:
        val += coeff * nk
        nk *= n
    return val

# CMF matrix entries
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
    return [[m11,m12,m13],[m21,m22,m23],[m31,m32,m33]]

def H(n):
    """Pochhammer normalization factor H_n = (-16)^n (2)_n^2 (3)_n^2 (5/2)_n (7/2)_n^2"""
    if n == 0:
        return Fraction(1)
    val = Fraction(1)
    for k in range(n):
        f = Fraction(-16) * Fraction(k+2)**2 * Fraction(k+3)**2 * Fraction(2*k+5, 2) * Fraction(2*k+7, 2)**2
        val *= f
    return val

# Compute Q̂_n = Q_{N,0} / H_n via CMF
# Q_{N,0} is the first component of (33750, -36000, 9000) M(0) M(1) ... M(N-1) divided by H_N
NMAX = 30
q_row = [Fraction(33750), Fraction(-36000), Fraction(9000)]

qhat = [None] * (NMAX + 1)
qhat[0] = Fraction(33750)

# Accumulate Q_{n,0} = row * M(0) * M(1) * ... * M(n-1), component 0, divided by H_n
row = [Fraction(33750), Fraction(-36000), Fraction(9000)]
for n in range(NMAX):
    M = M_entries(n)
    new_row = [Fraction(0)] * 3
    for j in range(3):
        for k in range(3):
            new_row[j] += Fraction(row[k]) * Fraction(M[k][j])
    row = new_row
    h = H(n + 1)
    qhat[n + 1] = Fraction(row[0], h) if h != 0 else None

# Delannoy numbers D_n = P_n(3)
D = [Fraction(0)] * (NMAX + 1)
D[0] = Fraction(1)
D[1] = Fraction(3)
for n in range(1, NMAX):
    D[n+1] = (Fraction(6*n+3) * D[n] - Fraction(n) * D[n-1]) / Fraction(n+1)

# Q_n(3) - Legendre second kind
Qleg = [mpf(0)] * (NMAX + 1)
Qleg[0] = log(2) / 2
Qleg[1] = 3 * Qleg[0] - 1
for n in range(1, NMAX):
    Qleg[n+1] = (mpf(6*n+3) * Qleg[n] - mpf(n) * Qleg[n-1]) / mpf(n+1)

print("Q̂_n / D_n² ratios:")
print("="*80)
for n in range(NMAX + 1):
    if qhat[n] is not None and D[n] != 0:
        dn2 = D[n]**2
        ratio = qhat[n] / dn2
        print(f"  n={n:3d}: Q̂_n/D_n² = {ratio}")

# Check if Q̂_n / D_n² is a polynomial in n
print("\n\nDifferences of Q̂_n / D_n²:")
ratios = []
for n in range(NMAX + 1):
    if qhat[n] is not None and D[n] != 0:
        ratios.append(qhat[n] / D[n]**2)
    else:
        ratios.append(None)

# First differences
diffs1 = []
for n in range(len(ratios)-1):
    if ratios[n] is not None and ratios[n+1] is not None:
        diffs1.append(ratios[n+1] - ratios[n])
    else:
        diffs1.append(None)

print("First differences:")
for i, d in enumerate(diffs1[:15]):
    if d is not None:
        print(f"  Δ[{i}] = {d}")

# Second differences
diffs2 = []
for n in range(len(diffs1)-1):
    if diffs1[n] is not None and diffs1[n+1] is not None:
        diffs2.append(diffs1[n+1] - diffs1[n])
    else:
        diffs2.append(None)

print("Second differences:")
for i, d in enumerate(diffs2[:15]):
    if d is not None:
        print(f"  Δ²[{i}] = {d}")

# Check Q̂_n against P_n(3)*Q_n(3) (mixed mode)
print("\n\nQ̂_n vs D_n² and P_n Q_n:")
from mpmath import catalan
G = catalan
for n in range(min(15, NMAX+1)):
    if qhat[n] is not None:
        q_mp = mpf(qhat[n].numerator) / mpf(qhat[n].denominator)
        d_mp = mpf(int(D[n]))
        pq = d_mp * Qleg[n]
        print(f"  n={n:3d}: Q̂_n={nstr(q_mp,12)}, D_n²={int(D[n]**2)}, PQ={nstr(pq,12)}, Q̂_n/D_n²={nstr(q_mp/d_mp**2,12)}")

# Verify the recurrence is satisfied
print("\n\nRecurrence verification:")
for n in range(min(NMAX - 2, 20)):
    if all(qhat[n+j] is not None for j in range(4)):
        res = sum(eval_c(j, Fraction(n)) * qhat[n+j] for j in range(4))
        print(f"  n={n:3d}: residual = {res}")

print("\nDone.")
