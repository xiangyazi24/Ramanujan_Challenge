#!/usr/bin/env python3
"""P2.5: Find the neutral invariant line of the normalized CMF.

For each n, compute the balanced matrix B̄(n) and find its eigenvalue closest to 1.
The corresponding eigenvector v(n) should be a rational function of n.
If v(n) is found, the block-extension certificate exists.
"""
from mpmath import mp, mpf, matrix, sqrt, eigsy
from fractions import Fraction

mp.dps = 60

sqrt2 = sqrt(2)

def M_exact(n):
    m11 = (-2*n-5)*(n+3)**2*(136*n**4+1424*n**3+5548*n**2+9551*n+6141)
    m12 = 384*n**6+6384*n**5+44168*n**4+162698*n**3+336377*n**2+369933*n+169011
    m13 = -480*n**4-4980*n**3-19210*n**2-32690*n-20730
    m21 = (n+2)**2*(n+3)**2*(4*n+10)*(48*n**3+386*n**2+1017*n+879)
    m22 = (n+2)**2*(-272*n**5-3848*n**4-21732*n**3-61184*n**2-85761*n-47808)
    m23 = (n+2)**2*(320*n**3+2540*n**2+6610*n+5640)
    m31 = (-4*n-10)*(n+2)**2*(n+3)**2*(32*n**4+302*n**3+1037*n**2+1530*n+813)
    m32 = (n+2)**2*(192*n**6+2984*n**5+19116*n**4+64452*n**3+120256*n**2+117279*n+46476)
    m33 = (n+2)**2*(-16*n**5-408*n**4-2912*n**3-8884*n**2-12254*n-6240)
    return matrix([[m11,m12,m13],[m21,m22,m23],[m31,m32,m33]])

def delta(n):
    return mpf(-2) * (n+2)**2 * (n+3)**2 * (2*n+5) * (2*n+7)**2

def balanced_matrix(n):
    """B̄(n) = D_n^{-1} M(n) D_{n+1} / δ(n)"""
    M = M_exact(n)
    d = delta(n)
    # D_n = diag(1, n+1, (n+1)²), D_{n+1} = diag(1, n+2, (n+2)²)
    n1 = mpf(n+1)
    n2 = mpf(n+2)
    B = matrix(3, 3)
    for i in range(3):
        for j in range(3):
            scale_inv = [1, 1/n1, 1/n1**2][i]
            scale_next = [1, n2, n2**2][j]
            B[i,j] = M[i,j] * scale_inv * scale_next / d
    return B

# Compute neutral eigenvalue and eigenvector for various n
print("="*80)
print("Neutral eigenvalue and eigenvector of B̄(n)")
print("="*80)
print(f"{'n':>3} {'λ_neutral':>20} {'v1':>20} {'v2':>20} {'v3':>20}")
print("-"*85)

neutral_data = []

for n in range(30):
    B = balanced_matrix(n)

    # Eigenvalues of B (general, not symmetric)
    # Use characteristic polynomial
    # det(B - λI) = 0
    # Compute eigenvalues numerically
    tr = B[0,0] + B[1,1] + B[2,2]
    cofsum = (B[0,0]*B[1,1] - B[0,1]*B[1,0]
            + B[0,0]*B[2,2] - B[0,2]*B[2,0]
            + B[1,1]*B[2,2] - B[1,2]*B[2,1])
    det = (B[0,0]*(B[1,1]*B[2,2]-B[1,2]*B[2,1])
         - B[0,1]*(B[1,0]*B[2,2]-B[1,2]*B[2,0])
         + B[0,2]*(B[1,0]*B[2,1]-B[1,1]*B[2,0]))

    # Eigenvalues are roots of λ³ - tr·λ² + cofsum·λ - det = 0
    # Find by numerical root-finding
    from mpmath import polyroots
    roots = polyroots([1, -tr, cofsum, -det])

    # Find the root closest to 1
    neutral_idx = min(range(3), key=lambda i: abs(roots[i] - 1))
    lam_n = roots[neutral_idx]

    # Find eigenvector: (B - λI)v = 0
    C = matrix(3, 3)
    for i in range(3):
        for j in range(3):
            C[i,j] = B[i,j] - (lam_n if i==j else 0)

    # Null vector: use cofactor method
    # v = cross product of any two rows of C
    row0 = [C[0,0], C[0,1], C[0,2]]
    row1 = [C[1,0], C[1,1], C[1,2]]
    v1 = row0[1]*row1[2] - row0[2]*row1[1]
    v2 = row0[2]*row1[0] - row0[0]*row1[2]
    v3 = row0[0]*row1[1] - row0[1]*row1[0]

    # Normalize by v3 (or largest component)
    if abs(v3) > 1e-50:
        v1n, v2n, v3n = v1/v3, v2/v3, mpf(1)
    elif abs(v1) > 1e-50:
        v1n, v2n, v3n = mpf(1), v2/v1, v3/v1
    else:
        v1n, v2n, v3n = v1, v2, v3

    neutral_data.append((n, lam_n, v1n, v2n, v3n))
    if n <= 10 or n % 5 == 0:
        print(f"{n:3d} {mp.nstr(lam_n, 12):>20s} {mp.nstr(v1n, 12):>20s} "
              f"{mp.nstr(v2n, 12):>20s} {mp.nstr(v3n, 12):>20s}")

# Check if the eigenvector components are rational in n
# by testing v1(n) * known denominators
print()
print("="*80)
print("Rational function test for v1(n)")
print("="*80)

# v1(n) should be p(n)/q(n) for some polynomials p, q.
# Try: v1(n) = (an² + bn + c) / (dn² + en + f) or similar

# Collect values for Lagrange interpolation
from mpmath import mpf
v1_vals = [(d[0], d[2]) for d in neutral_data]

# Check: is v1(n) a polynomial? Compute differences.
diffs = [v1_vals[i+1][1] - v1_vals[i][1] for i in range(len(v1_vals)-1)]
ddiffs = [diffs[i+1] - diffs[i] for i in range(len(diffs)-1)]
dddiffs = [ddiffs[i+1] - ddiffs[i] for i in range(len(ddiffs)-1)]

print("First differences of v1(n):")
for i in range(min(10, len(diffs))):
    print(f"  Δv1({i}) = {mp.nstr(diffs[i], 15)}")

print("\nSecond differences of v1(n):")
for i in range(min(10, len(ddiffs))):
    print(f"  Δ²v1({i}) = {mp.nstr(ddiffs[i], 15)}")

print("\nThird differences of v1(n):")
for i in range(min(10, len(dddiffs))):
    print(f"  Δ³v1({i}) = {mp.nstr(dddiffs[i], 15)}")

# Also check: v1 * (n+1)^k for small k
print()
print("="*80)
print("Check if v1·(n+1)^k or v1·(2n+3)·(2n+5)·... is polynomial")
print("="*80)
for k in range(5):
    vals = [(d[0], d[2] * (d[0]+1)**k) for d in neutral_data[:15]]
    diffs_k = [vals[i+1][1] - vals[i][1] for i in range(len(vals)-1)]
    ddiffs_k = [diffs_k[i+1] - diffs_k[i] for i in range(len(diffs_k)-1)]
    dddiffs_k = [ddiffs_k[i+1] - ddiffs_k[i] for i in range(len(ddiffs_k)-1)]
    d4 = [dddiffs_k[i+1] - dddiffs_k[i] for i in range(len(dddiffs_k)-1)]
    d5 = [d4[i+1] - d4[i] for i in range(len(d4)-1)]

    # Check if k-th order differences stabilize (=polynomial of degree k)
    last_diffs = d5[-3:] if len(d5) >= 3 else d5
    max_ratio = max(abs(x) for x in last_diffs) if last_diffs else mpf(1)
    print(f"  k={k}: max |Δ⁵(v1·(n+1)^k)| ≈ {mp.nstr(max_ratio, 5)}")

# Check r₀(n) = λ_neutral(n): the neutral multiplier
print()
print("="*80)
print("Neutral multiplier r₀(n) analysis")
print("="*80)
for d in neutral_data[:15]:
    n, lam, v1, v2, v3 = d
    # r₀(n) should approach 1 as n → ∞
    # More precisely, r₀(n) → 1 with correction ~ -3/(n+c)
    # (since formal index is -3, giving n^{-3} decay)
    r0_corr = (lam - 1) * (n+1) if n > 0 else 0
    print(f"  n={n:2d}: r₀ = {mp.nstr(lam, 20)}, (r₀-1)·(n+1) = {mp.nstr(r0_corr, 15)}")
