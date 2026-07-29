#!/usr/bin/env python3
"""P2.5: Find the Ore intertwiner between the CMF and Delannoy-square modules.

If Q̂_n/D_n² satisfies a simple recurrence, this gives the module map.
Also try: does Q̂_n satisfy a SPECIFIC higher-order recurrence related to
the integrated-K module?
"""
from fractions import Fraction
from mpmath import mp, mpf, matrix, nstr, rf, catalan

mp.dps = 200

def M_exact_int(n):
    """Return M(n) as a list of lists of Python ints."""
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

def H_frac(n):
    """H_n as exact Fraction."""
    result = Fraction(1)
    for k in range(n):
        result *= Fraction(-16)
        result *= Fraction(k+2)**2
        result *= Fraction(k+3)**2
        result *= Fraction(2*k+5, 2)
        result *= Fraction(2*k+7, 2)**2
    return result

def legendre_P_int(n, x=3):
    """Legendre P_n(3) using exact integer arithmetic."""
    if n == 0: return 1
    if n == 1: return x
    p0, p1 = 1, x
    for k in range(1, n):
        p2 = ((2*k+1)*x*p1 - k*p0) // (k+1)
        p0, p1 = p1, p2
    return p1

# Compute Q_{N,0} as exact integers
print("Computing Q_{N,0} for N = 0..70...")
q = [33750, -36000, 9000]
p = [30921, -32972, 8240]
Q0_seq = [33750]
P0_seq = [30921]

for N in range(70):
    M = M_exact_int(N)
    new_q = [0, 0, 0]
    new_p = [0, 0, 0]
    for j in range(3):
        for k in range(3):
            new_q[j] += q[k] * M[k][j]
            new_p[j] += p[k] * M[k][j]
    q = new_q
    p = new_p
    Q0_seq.append(q[0])
    P0_seq.append(p[0])

print(f"  Done. {len(Q0_seq)} terms.")

# Compute Q̂_n = Q_{N,0}/H_n and D_n² = P_n(3)²
print("\nComputing Q̂_n/D_n² ratios...")
ratios = []
Q_hat = []
D2_vals = []

for n in range(50):
    Hn = H_frac(n)
    Dn = legendre_P_int(n, 3)
    D2 = Dn * Dn

    if Hn != 0 and D2 != 0:
        q_hat = Fraction(Q0_seq[n]) / Hn
        ratio = q_hat / D2
        Q_hat.append(q_hat)
        D2_vals.append(Fraction(D2))
        ratios.append(ratio)
    else:
        Q_hat.append(Fraction(Q0_seq[n]))
        D2_vals.append(Fraction(D2))
        ratios.append(None)

# Check: is the ratio Q̂_n/D_n² a polynomial in n?
print("\nQ̂_n/D_n² values:")
for n in range(min(15, len(ratios))):
    if ratios[n] is not None:
        print(f"  n={n:2d}: ratio = {float(ratios[n]):.6f}")

# Compute differences of the ratio
print("\nFirst differences:")
diffs = [ratios[i+1] - ratios[i] for i in range(min(20, len(ratios)-1)) if ratios[i] is not None and ratios[i+1] is not None]
for i, d in enumerate(diffs[:10]):
    print(f"  Δ({i}) = {float(d):.6f}")

print("\nSecond differences:")
ddiffs = [diffs[i+1] - diffs[i] for i in range(len(diffs)-1)]
for i, d in enumerate(ddiffs[:10]):
    print(f"  Δ²({i}) = {float(d):.6f}")

print("\nThird differences:")
dddiffs = [ddiffs[i+1] - ddiffs[i] for i in range(len(ddiffs)-1)]
for i, d in enumerate(dddiffs[:10]):
    print(f"  Δ³({i}) = {float(d):.10f}")

print("\nFourth differences:")
d4 = [dddiffs[i+1] - dddiffs[i] for i in range(len(dddiffs)-1)]
for i, d in enumerate(d4[:10]):
    print(f"  Δ⁴({i}) = {float(d):.15f}")

# Check if the differences stabilize (indicating a polynomial)
print("\n5th through 8th differences (last 3 values):")
dk = d4
for order in range(5, 9):
    dk = [dk[i+1] - dk[i] for i in range(len(dk)-1)]
    if len(dk) >= 3:
        print(f"  Δ^{order}: {[float(x) for x in dk[-3:]]}")

# Now try: does the ratio a_n satisfy a first-order RECURRENCE?
# i.e., a_{n+1}/a_n = r(n) for some rational function r(n)?
print("\n" + "="*60)
print("Ratio a_{n+1}/a_n:")
print("="*60)
for n in range(min(15, len(ratios)-1)):
    if ratios[n] is not None and ratios[n+1] is not None and ratios[n] != 0:
        r = ratios[n+1] / ratios[n]
        print(f"  n={n:2d}: a_{n+1}/a_n = {float(r):.15f}")

# Try: a_n - (αn + β) for best linear fit
# a_n ≈ 147094·(n-1) + ... seems wrong. Let me check the actual values
print("\n" + "="*60)
print("Checking if a_n = α·n + β + γ/n + ...")
print("="*60)
# Use n=10, 20, 30 to fit α, β, γ
if len(ratios) > 30:
    a10 = float(ratios[10])
    a20 = float(ratios[20])
    a30 = float(ratios[30])
    print(f"  a_10 = {a10:.6f}")
    print(f"  a_20 = {a20:.6f}")
    print(f"  a_30 = {a30:.6f}")

    # Linear fit: a_n = α·n + β
    alpha = (a30 - a10) / 20
    beta = a10 - 10*alpha
    print(f"  Linear fit: α = {alpha:.6f}, β = {beta:.6f}")
    for n in [5, 15, 25, 35, 40]:
        if n < len(ratios) and ratios[n] is not None:
            pred = alpha*n + beta
            actual = float(ratios[n])
            print(f"    n={n:2d}: predicted={pred:.2f}, actual={actual:.2f}, diff={actual-pred:.6f}")

# Try: a_n = (P(n) + Q(n)·G) / R(n) for polynomials P, Q, R
# The key insight: G should appear because the CMF converges to G
print("\n" + "="*60)
print("Testing a_n - linear = O(1/n) correction")
print("="*60)
# Compute b_n = (a_n - α·n - β) · n for large n
if len(ratios) > 30:
    for n in [10, 15, 20, 25, 30, 35, 40, 45]:
        if n < len(ratios) and ratios[n] is not None:
            correction = (float(ratios[n]) - alpha*n - beta) * n
            print(f"  n={n:2d}: (a_n - αn - β)·n = {correction:.6f}")

# CRITICAL: check if the EXACT ratio a_n = Q̂_n / D_n² is a rational number!
print("\n" + "="*60)
print("Exact rationality of Q̂_n / D_n²:")
print("="*60)
for n in range(min(10, len(ratios))):
    if ratios[n] is not None:
        # ratios[n] is already a Fraction
        r = ratios[n]
        print(f"  n={n}: Q̂/D² = {r.numerator}/{r.denominator}")
        # Check denominator: is it a power of 2?
        d = r.denominator
        while d % 2 == 0:
            d //= 2
        print(f"          denom/2^k = {d}")
