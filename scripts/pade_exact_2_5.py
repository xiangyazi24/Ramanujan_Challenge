#!/usr/bin/env python3
"""Exact Padé for r_int: 31 points, 31 unknowns. Check if coefficients are in Q or Q(sqrt2).

Use Lagrange-style interpolation: pick 31 evenly-spaced points in [80, 140],
solve the system r(n)*D(n) = N(n) exactly, check coefficient structure.
"""
from mpmath import mp, mpf, nstr, matrix as mp_matrix, lu_solve, sqrt
import time
mp.dps = 2000

def M_mat(n):
    n = mpf(n)
    m11 = (-2*n-5)*(n+3)**2 * (136*n**4 + 1424*n**3 + 5548*n**2 + 9551*n + 6141)
    m12 = 384*n**6 + 6384*n**5 + 44168*n**4 + 162698*n**3 + 336377*n**2 + 369933*n + 169011
    m13 = -480*n**4 - 4980*n**3 - 19210*n**2 - 32690*n - 20730
    m21 = (n+2)**2*(n+3)**2*(4*n+10)*(48*n**3 + 386*n**2 + 1017*n + 879)
    m22 = (n+2)**2*(-272*n**5 - 3848*n**4 - 21732*n**3 - 61184*n**2 - 85761*n - 47808)
    m23 = (n+2)**2*(320*n**3 + 2540*n**2 + 6610*n + 5640)
    m31 = (-4*n-10)*(n+2)**2*(n+3)**2*(32*n**4 + 302*n**3 + 1037*n**2 + 1530*n + 813)
    m32 = (n+2)**2*(192*n**6 + 2984*n**5 + 19116*n**4 + 64452*n**3 + 120256*n**2 + 117279*n + 46476)
    m33 = (n+2)**2*(-16*n**5 - 408*n**4 - 2912*n**3 - 8884*n**2 - 12254*n - 6240)
    return [[m11,m12,m13],[m21,m22,m23],[m31,m32,m33]]

N_max = 170
v1 = [mpf(1), mpf(0), mpf(0)]
v2 = [mpf(0), mpf(1), mpf(0)]
u1_vals = [v1[0]]; u2_vals = [v2[0]]
for N in range(N_max + 5):
    M = M_mat(N)
    v1_new = [sum(M[i][j]*v1[j] for j in range(3)) for i in range(3)]
    v2_new = [sum(M[i][j]*v2[j] for j in range(3)) for i in range(3)]
    v1 = v1_new; v2 = v2_new
    u1_vals.append(v1[0]); u2_vals.append(v2[0])

R = u1_vals[160] / u2_vals[160]
w_vals = [u1_vals[n] - R * u2_vals[n] for n in range(len(u1_vals))]
r_int = {}
for n in range(N_max + 2):
    if abs(w_vals[n]) > mpf('1e-1500') and abs(w_vals[n+1]) > mpf('1e-1500'):
        r_int[n] = w_vals[n+1] / w_vals[n]

# Use the FIRST few (small n) values where recessive contamination is the issue,
# and ALSO large n values. The best precision is at intermediate n.
# Let's check: use n=80..110 (31 points)
dP, dQ = 19, 12
n_unk = dP + dQ  # 31

eval_pts = list(range(80, 111))  # 31 points
assert len(eval_pts) == n_unk

# System: r(n) * D(n) - N(n) = 0
# D(n) = n^12 + d11*n^11 + ... + d0 (monic, 12 unknowns)
# N(n) = c19*n^19 + c18*n^18 + ... + c0 (20 unknowns)
# Total: 32 unknowns, but we have 31 equations.
# Fix D monic degree 12. Then:
# r(n) * (n^12 + sum d_j n^j) = sum c_k n^k
# Rewrite: sum c_k n^k - r(n) * sum d_j n^j = r(n) * n^12
# 32 unknowns (c0..c19, d0..d11), 31 equations.
# Need to fix one more. Fix c19 so that leading term is -16.
# Then: c19 = -16, and we have 31 unknowns and 31 equations.

mat = mp_matrix(31, 31)
rhs_vec = mp_matrix(31, 1)
for idx, nv_int in enumerate(eval_pts):
    nv = mpf(nv_int)
    rv = r_int[nv_int]
    col = 0
    # c0..c18 (19 unknowns for N)
    for k in range(19):
        mat[idx, col] = nv**k
        col += 1
    # d0..d11 (12 unknowns for D)
    for j in range(12):
        mat[idx, col] = -rv * nv**j
        col += 1
    # RHS: r(n)*n^12 - c19*n^19 = r(n)*n^12 + 16*n^19
    rhs_vec[idx, 0] = rv * nv**12 + mpf(16) * nv**19

print("Solving 31x31 exact Padé system...")
sol = lu_solve(mat, rhs_vec)

print("\n=== Numerator N(n) = -16*n^19 + c18*n^18 + ... + c0 ===")
c_coeffs = [sol[k, 0] for k in range(19)] + [mpf(-16)]  # c0..c18, c19=-16
print("\n=== Denominator D(n) = n^12 + d11*n^11 + ... + d0 ===")
d_coeffs = [sol[19 + j, 0] for j in range(12)] + [mpf(1)]  # d0..d11, d12=1

# Try to recognize as a + b*sqrt(2)
sqrt2 = mp.sqrt(2)

print("\nRecognizing N coefficients as a + b*sqrt(2):")
for k in range(20):
    val = c_coeffs[k]
    # Try: val = a + b*sqrt2 with a, b rational
    # Method: if val = a + b*sqrt2, then
    # val' = a - b*sqrt2 (conjugate)
    # a = (val + val')/2, b = (val - val')/(2*sqrt2)
    # But we don't know val'. Instead, use LLL/PSLQ.
    # Simple approach: try val = p/q + (r/s)*sqrt2
    # i.e., s*q*val = s*p + r*q*sqrt2
    # Use integer relation: find integers m1, m2, m3 such that
    # m1 + m2*val + m3*sqrt2 ≈ 0, i.e. val ≈ -m1/m2 - (m3/m2)*sqrt2

    # Use PSLQ on [1, val, sqrt2]
    from mpmath import pslq
    rel = pslq([1, val, sqrt2], maxcoeff=10**30, tol=mpf('1e-500'))
    if rel is not None:
        m1, m2, m3 = rel
        # m1 + m2*val + m3*sqrt2 = 0
        # val = (-m1 - m3*sqrt2) / m2
        a_part = mpf(-m1) / m2
        b_part = mpf(-m3) / m2
        # Check
        resid = abs(val - a_part - b_part * sqrt2)
        if resid < mpf('1e-100'):
            from fractions import Fraction
            fa = Fraction(-m1, m2).limit_denominator(10**15)
            fb = Fraction(-m3, m2).limit_denominator(10**15)
            print(f"  c[{k:>2}] = {fa} + ({fb})*√2   [resid={float(mp.log10(resid)):.0f}]")
        else:
            print(f"  c[{k:>2}] = {nstr(val, 40)} (PSLQ relation poor: resid {float(resid):.2e})")
    else:
        # Try just rational
        rel2 = pslq([1, val], maxcoeff=10**30, tol=mpf('1e-500'))
        if rel2 is not None:
            m1, m2 = rel2
            rat = Fraction(-m1, m2).limit_denominator(10**15)
            print(f"  c[{k:>2}] = {rat}  (rational)")
        else:
            print(f"  c[{k:>2}] = {nstr(val, 40)} (no relation found)")

print("\nRecognizing D coefficients as a + b*sqrt(2):")
for j in range(13):
    val = d_coeffs[j]
    from mpmath import pslq
    rel = pslq([1, val, sqrt2], maxcoeff=10**30, tol=mpf('1e-500'))
    if rel is not None:
        m1, m2, m3 = rel
        a_part = mpf(-m1) / m2
        b_part = mpf(-m3) / m2
        resid = abs(val - a_part - b_part * sqrt2)
        if resid < mpf('1e-100'):
            from fractions import Fraction
            fa = Fraction(-m1, m2).limit_denominator(10**15)
            fb = Fraction(-m3, m2).limit_denominator(10**15)
            print(f"  d[{j:>2}] = {fa} + ({fb})*√2   [resid={float(mp.log10(resid)):.0f}]")
        else:
            print(f"  d[{j:>2}] = {nstr(val, 40)} (PSLQ poor)")
    else:
        rel2 = pslq([1, val], maxcoeff=10**30, tol=mpf('1e-500'))
        if rel2 is not None:
            from fractions import Fraction
            m1, m2 = rel2
            rat = Fraction(-m1, m2).limit_denominator(10**15)
            print(f"  d[{j:>2}] = {rat}  (rational)")
        else:
            print(f"  d[{j:>2}] = {nstr(val, 40)} (no relation found)")

# Verify at points outside the fitting range
print("\nVerification at out-of-sample points:")
for n_test in [0, 1, 5, 10, 20, 50, 70, 120, 140, 150]:
    if n_test not in r_int:
        continue
    nv = mpf(n_test)
    N_val = sum(c_coeffs[k] * nv**k for k in range(20))
    D_val = sum(d_coeffs[j] * nv**j for j in range(13))
    r_rat = N_val / D_val
    rel = abs((r_rat - r_int[n_test]) / r_int[n_test])
    logr = float(mp.log10(rel)) if rel > 0 else -9999
    print(f"  n={n_test:>3}: log10(rel_diff) = {logr:.1f}")
