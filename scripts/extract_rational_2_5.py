#!/usr/bin/env python3
"""Extract exact rational coefficients of r_int(n) = -16 * P19(n)/Q12(n).

From the reconstruction: dP=19, dQ=12 gives residual 10^{-1956} (machine precision).
r_int(n) = -16 * (n^19 + a18*n^18 + ... + a0) / (n^12 + b11*n^11 + ... + b0)

Strategy: compute r_int at 35 integer points with 2000-digit precision,
then solve the 31-unknown linear system exactly.
"""
from mpmath import mp, mpf, nstr, matrix as mp_matrix, lu_solve
import time, json
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

# Compute forward solutions
N_max = 170
v1 = [mpf(1), mpf(0), mpf(0)]
v2 = [mpf(0), mpf(1), mpf(0)]
u1_vals = [v1[0]]
u2_vals = [v2[0]]

print("Computing forward solutions...")
t0 = time.time()
for N in range(N_max + 5):
    M = M_mat(N)
    def mat_vec(M, v):
        return [sum(M[i][j]*v[j] for j in range(3)) for i in range(3)]
    v1 = mat_vec(M, v1)
    v2 = mat_vec(M, v2)
    u1_vals.append(v1[0])
    u2_vals.append(v2[0])
print(f"Done in {time.time()-t0:.1f}s")

R = u1_vals[160] / u2_vals[160]
w = [u1_vals[n] - R * u2_vals[n] for n in range(len(u1_vals))]
r_int = {}
for n in range(N_max + 2):
    if abs(w[n]) > mpf('1e-1500'):
        r_int[n] = w[n+1] / w[n]

# Now solve the exact system for dP=19, dQ=12
# r_int(n) * Q(n) = -16 * P(n)
# r_int(n) * (n^12 + b11*n^11 + ... + b0) = -16 * (n^19 + a18*n^18 + ... + a0)
# Unknowns: a0..a18 (19), b0..b11 (12) = 31 total
# At each evaluation point n_val:
#   -16*a0 - 16*a1*n - ... - 16*a18*n^18 + r*b0 + r*b1*n + ... + r*b11*n^11 = -r*n^12 - (-16)*n^19
#   i.e., sum_{k=0}^{18} (-16*n^k)*a_k + sum_{j=0}^{11} (r*n^j)*b_j = 16*n^19 - r*n^12

dP, dQ = 19, 12
n_unk = dP + dQ  # 31

# Use evaluation points n=130..165 (36 points, overdetermined)
eval_pts = list(range(130, 166))
n_eq = len(eval_pts)
print(f"\nSolving {n_unk} unknowns from {n_eq} equations...")

mat = mp_matrix(n_eq, n_unk)
rhs = mp_matrix(n_eq, 1)
for idx, nv_int in enumerate(eval_pts):
    nv = mpf(nv_int)
    rv = r_int[nv_int]
    for k in range(dP):
        mat[idx, k] = mpf(-16) * nv**k
    for j in range(dQ):
        mat[idx, dP + j] = rv * nv**j
    rhs[idx, 0] = mpf(16) * nv**dP - rv * nv**dQ

# Solve via normal equations
ATA = mat.T * mat
ATb = mat.T * rhs
sol = lu_solve(ATA, ATb)

# Print coefficients and try to recognize as rationals
print("\nNumerator P(n) = n^19 + a18*n^18 + ... + a0:")
a_coeffs = [sol[k, 0] for k in range(dP)]
for k in range(dP-1, -1, -1):
    # Try to recognize as rational p/q with small denominator
    val = a_coeffs[k]
    found_rational = False
    for denom in range(1, 10001):
        numer = mp.nint(val * denom)
        if abs(val - numer/denom) < mpf('1e-500'):
            print(f"  a[{k}] = {int(numer)}/{denom}")
            found_rational = True
            break
    if not found_rational:
        print(f"  a[{k}] = {nstr(val, 80)} (no small rational found)")

print("\nDenominator Q(n) = n^12 + b11*n^11 + ... + b0:")
b_coeffs = [sol[dP + j, 0] for j in range(dQ)]
for j in range(dQ-1, -1, -1):
    val = b_coeffs[j]
    found_rational = False
    for denom in range(1, 10001):
        numer = mp.nint(val * denom)
        if abs(val - numer/denom) < mpf('1e-500'):
            print(f"  b[{j}] = {int(numer)}/{denom}")
            found_rational = True
            break
    if not found_rational:
        print(f"  b[{j}] = {nstr(val, 80)} (no small rational found)")

# Verify: compute r_int(n) from the rational function at some points and compare
print("\nVerification: r_rational(n) vs r_int(n)")
for n_test in [0, 5, 10, 20, 50, 100, 130, 150]:
    if n_test not in r_int:
        continue
    nv = mpf(n_test)
    P_val = nv**dP + sum(a_coeffs[k] * nv**k for k in range(dP))
    Q_val = nv**dQ + sum(b_coeffs[j] * nv**j for j in range(dQ))
    r_rat = mpf(-16) * P_val / Q_val
    diff = abs(r_rat - r_int[n_test])
    log_diff = float(mp.log10(diff / abs(r_int[n_test]))) if diff > 0 else -9999
    print(f"  n={n_test:>3}: log10(rel_diff) = {log_diff:.1f}")

