#!/usr/bin/env python3
"""P2.5: Extract the neutral mode of the CMF error ê_n = G·Q̂_n - P̂_n.
Check if f₀(n+1)/f₀(n) is a rational function of n (hypergeometric solution)."""
from fractions import Fraction
from mpmath import mp, mpf, catalan, nstr, log10

mp.dps = 300

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

def H_mpf(n):
    if n == 0:
        return mpf(1)
    val = mpf(1)
    for k in range(n):
        f = mpf(-16) * mpf(k+2)**2 * mpf(k+3)**2 * (mpf(2*k+5)/2) * (mpf(2*k+7)/2)**2
        val *= f
    return val

NMAX = 120
G = catalan

# Compute both p and q rows via CMF
p_row = [mpf(30921), mpf(-32972), mpf(8240)]
q_row = [mpf(33750), mpf(-36000), mpf(9000)]

phat = [None] * (NMAX + 1)
qhat = [None] * (NMAX + 1)
phat[0] = mpf(30921)
qhat[0] = mpf(33750)

pr = [mpf(30921), mpf(-32972), mpf(8240)]
qr = [mpf(33750), mpf(-36000), mpf(9000)]

for n in range(NMAX):
    M = M_entries(n)
    new_pr = [mpf(0)] * 3
    new_qr = [mpf(0)] * 3
    for j in range(3):
        for k in range(3):
            new_pr[j] += pr[k] * mpf(M[k][j])
            new_qr[j] += qr[k] * mpf(M[k][j])
    pr = new_pr
    qr = new_qr
    h = H_mpf(n + 1)
    phat[n + 1] = pr[0] / h
    qhat[n + 1] = qr[0] / h

# Compute error ê_n = G Q̂_n - P̂_n
ehat = [None] * (NMAX + 1)
for n in range(NMAX + 1):
    if phat[n] is not None and qhat[n] is not None:
        ehat[n] = G * qhat[n] - phat[n]

print("Neutral mode analysis of ê_n = G·Q̂_n - P̂_n:")
print("="*80)

# Check ê_n behavior: should be ~ C₀ n⁻³ for large n
print("\n1. ê_n values and n³·ê_n (should converge to constant):")
for n in [1, 2, 5, 10, 20, 30, 50, 70, 100]:
    if ehat[n] is not None and ehat[n] != 0:
        n3e = mpf(n)**3 * ehat[n]
        print(f"  n={n:4d}: ê_n = {nstr(ehat[n], 15)}, n³·ê_n = {nstr(n3e, 15)}")

# Check ratio ê_(n+1)/ê_n — should approach 1 as (1-3/n + O(1/n²))
print("\n2. Successive ratios r(n) = ê_(n+1)/ê_n:")
for n in [3, 5, 10, 20, 30, 50, 70, 100]:
    if ehat[n] is not None and ehat[n+1] is not None and ehat[n] != 0:
        r = ehat[n+1] / ehat[n]
        # Subtract 1 and multiply by n to see the -3/n coefficient
        deviation = (r - 1) * mpf(n)
        print(f"  n={n:4d}: r(n) = {nstr(r, 20)}, (r-1)·n = {nstr(deviation, 15)}")

# Check if (r(n)-1)·n converges (it should approach -3 for formal index -3)
print("\n3. Testing if r(n) = 1 - 3/n + a/n² + ..., compute n²·(r(n) - 1 + 3/n):")
for n in [5, 10, 20, 30, 50, 70, 100]:
    if ehat[n] is not None and ehat[n+1] is not None and ehat[n] != 0:
        r = ehat[n+1] / ehat[n]
        correction = mpf(n)**2 * (r - 1 + mpf(3)/mpf(n))
        print(f"  n={n:4d}: n²·(r-1+3/n) = {nstr(correction, 15)}")

# Check if r(n) is EXACTLY rational: try r(n) = (n-2)³/(n+1)³ or similar
print("\n4. Testing specific rational forms:")
for n in [10, 20, 30, 50, 70, 100]:
    if ehat[n] is not None and ehat[n+1] is not None and ehat[n] != 0:
        r_actual = ehat[n+1] / ehat[n]

        # Test r(n) = n³/(n+1)³ = (n/(n+1))³ → formal index -3
        r_test1 = (mpf(n)/(mpf(n)+1))**3
        err1 = abs(r_actual - r_test1)

        # Test r(n) = (n-1)²(n-2)/((n+2)(n+1)²) — preserving degree, formal index -3
        r_test2 = mpf(n-1)**2 * mpf(n-2) / (mpf(n+2) * mpf(n+1)**2)
        err2 = abs(r_actual - r_test2)

        # Test r(n) = (n-1)(n-2)(n-3)/((n+1)(n+2)(n+3)) — (n-a)/(n+a) pattern
        r_test3 = mpf(n-1)*mpf(n-2)*mpf(n-3)/(mpf(n+1)*mpf(n+2)*mpf(n+3))
        err3 = abs(r_actual - r_test3)

        print(f"  n={n:4d}: |r - n³/(n+1)³| = {nstr(err1, 6)}, "
              f"|r - (n-1)²(n-2)/((n+2)(n+1)²)| = {nstr(err2, 6)}, "
              f"|r - (n-1)(n-2)(n-3)/((n+1)(n+2)(n+3))| = {nstr(err3, 6)}")

# More sophisticated: compute the expansion r(n) = 1 + a₁/n + a₂/n² + a₃/n³ + ...
print("\n5. Asymptotic expansion of r(n):")
# Use large n values to extract coefficients
ns = [60, 70, 80, 90, 100, 110]
from mpmath import matrix, lu_solve

A_mat = []
b_vec = []
for n in ns:
    if ehat[n] is not None and ehat[n+1] is not None and ehat[n] != 0:
        r = ehat[n+1] / ehat[n]
        row = [mpf(1)/mpf(n)**k for k in range(1, len(ns)+1)]
        A_mat.append(row)
        b_vec.append(r - 1)

if len(A_mat) >= 4:
    A = matrix(A_mat)
    b = matrix(b_vec)
    try:
        coeffs = lu_solve(A, b)
        print("  r(n) = 1 + Σ aₖ/nᵏ:")
        for k in range(min(len(coeffs), 8)):
            print(f"    a_{k+1} = {nstr(coeffs[k], 20)}")
    except:
        print("  LU solve failed")

# Also extract the TWO-step ratio (since subdominant is complex pair)
print("\n6. Two-step ratio ê_(n+2)/ê_n (more stable for complex subdominant):")
for n in [10, 20, 30, 50, 70, 100]:
    if ehat[n] is not None and ehat[n+2] is not None and ehat[n] != 0:
        r2 = ehat[n+2] / ehat[n]
        print(f"  n={n:4d}: ê_(n+2)/ê_n = {nstr(r2, 20)}")

# Check: is n³ ê_n itself a known sequence?
print("\n7. n³·ê_n for small n:")
for n in range(1, 30):
    if ehat[n] is not None:
        val = mpf(n)**3 * ehat[n]
        print(f"  n={n:3d}: n³·ê_n = {nstr(val, 30)}")

print("\nDone.")
