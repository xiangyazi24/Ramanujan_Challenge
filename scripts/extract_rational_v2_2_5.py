#!/usr/bin/env python3
"""Extract exact rational coefficients of r_int(n) = -16 * P19(n)/Q12(n).

The zero at n=160 is an artifact (R computed from n=160).
Use n=80..155 where intermediate mode dominates.
"""
from mpmath import mp, mpf, nstr, matrix as mp_matrix, lu_solve
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
u1_vals = [v1[0]]
u2_vals = [v2[0]]
for N in range(N_max + 5):
    M = M_mat(N)
    v1_new = [sum(M[i][j]*v1[j] for j in range(3)) for i in range(3)]
    v2_new = [sum(M[i][j]*v2[j] for j in range(3)) for i in range(3)]
    v1 = v1_new; v2 = v2_new
    u1_vals.append(v1[0]); u2_vals.append(v2[0])

R = u1_vals[160] / u2_vals[160]
w = [u1_vals[n] - R * u2_vals[n] for n in range(len(u1_vals))]
r_int = {}
for n in range(N_max + 2):
    if abs(w[n]) > mpf('1e-1500') and abs(w[n+1]) > mpf('1e-1500'):
        r_int[n] = w[n+1] / w[n]

# Solve for r_int(n) = -16 * P(n)/Q(n) where P has degree 19, Q has degree 12
# Rewrite: r_int(n) * Q(n) + 16 * P(n) = 0
# With monic P, Q: P = n^19 + a18*n^18 + ... + a0, Q = n^12 + b11*n^11 + ... + b0
# At each n: r(n) * (n^12 + sum b_j n^j) + 16 * (n^19 + sum a_k n^k) = 0
# => sum_k 16*a_k*n^k + sum_j r(n)*b_j*n^j = -16*n^19 - r(n)*n^12

dP, dQ = 19, 12
n_unk = dP + dQ

eval_pts = [n for n in range(80, 156) if n in r_int]
n_eq = len(eval_pts)
print(f"Using {n_eq} eval points from n={min(eval_pts)} to {max(eval_pts)} for {n_unk} unknowns")

mat = mp_matrix(n_eq, n_unk)
rhs = mp_matrix(n_eq, 1)
for idx, nv_int in enumerate(eval_pts):
    nv = mpf(nv_int)
    rv = r_int[nv_int]
    for k in range(dP):
        mat[idx, k] = mpf(16) * nv**k
    for j in range(dQ):
        mat[idx, dP + j] = rv * nv**j
    rhs[idx, 0] = mpf(-16) * nv**dP - rv * nv**dQ

ATA = mat.T * mat
ATb = mat.T * rhs
sol = lu_solve(ATA, ATb)

print("\n=== Coefficients ===")
print("\nNumerator P(n) = n^19 + a[18]*n^18 + ... + a[0]:")
a_coeffs = [sol[k, 0] for k in range(dP)]
for k in range(dP-1, -1, -1):
    val = a_coeffs[k]
    found = False
    for denom in range(1, 100001):
        numer = mp.nint(val * denom)
        if abs(val - mpf(numer)/denom) < mpf('1e-100'):
            # Simplify
            from math import gcd
            g = gcd(abs(int(numer)), denom)
            num_s, den_s = int(numer)//g, denom//g
            if den_s == 1:
                print(f"  a[{k:>2}] = {num_s}")
            else:
                print(f"  a[{k:>2}] = {num_s}/{den_s}")
            found = True
            break
    if not found:
        print(f"  a[{k:>2}] = {nstr(val, 60)} (no rational ≤ 100000)")

print("\nDenominator Q(n) = n^12 + b[11]*n^11 + ... + b[0]:")
b_coeffs = [sol[dP + j, 0] for j in range(dQ)]
for j in range(dQ-1, -1, -1):
    val = b_coeffs[j]
    found = False
    for denom in range(1, 100001):
        numer = mp.nint(val * denom)
        if abs(val - mpf(numer)/denom) < mpf('1e-100'):
            from math import gcd
            g = gcd(abs(int(numer)), denom)
            num_s, den_s = int(numer)//g, denom//g
            if den_s == 1:
                print(f"  b[{j:>2}] = {num_s}")
            else:
                print(f"  b[{j:>2}] = {num_s}/{den_s}")
            found = True
            break
    if not found:
        print(f"  b[{j:>2}] = {nstr(val, 60)} (no rational ≤ 100000)")

# Verify
print("\nVerification: r_rational(n) vs r_int(n)")
for n_test in [0, 1, 5, 10, 20, 50, 80, 100, 130, 150]:
    if n_test not in r_int:
        continue
    nv = mpf(n_test)
    P_val = nv**dP + sum(a_coeffs[k] * nv**k for k in range(dP))
    Q_val = nv**dQ + sum(b_coeffs[j] * nv**j for j in range(dQ))
    r_rat = mpf(-16) * P_val / Q_val
    rel_diff = abs((r_rat - r_int[n_test]) / r_int[n_test])
    log_rd = float(mp.log10(rel_diff)) if rel_diff > 0 else -9999
    print(f"  n={n_test:>3}: log10(rel_diff) = {log_rd:.1f}")
