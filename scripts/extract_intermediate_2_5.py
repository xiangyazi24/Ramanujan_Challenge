#!/usr/bin/env python3
"""Extract intermediate Poincaré mode ratio r_int(n) via gauge-and-peel.

Strategy:
1. Backward iteration on order-3 → get r_rec(n) [recessive mode]
2. Gauge by r_rec: dj(n) = cj(n) * prod r_rec(n+i)
3. Since (S-1) is right factor of gauged recurrence, divide it out → order-2
4. Backward iteration on order-2 → get σ(n) [smaller of two remaining modes]
   σ(n) should be the intermediate/recessive ratio ≈ 34.1
5. r_int(n) = σ(n) * r_rec(n)
"""
from fractions import Fraction
from mpmath import mp, mpf, nstr, polyroots

mp.dps = 200

# Rebuild exact recurrence
def M_mat_qq(n):
    n = Fraction(n)
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

A_row = [Fraction(33750), Fraction(-36000), Fraction(9000)]
T = [[Fraction(1 if i==j else 0) for j in range(3)] for i in range(3)]
q_vals = []
for N in range(120):
    q = sum(A_row[k] * T[k][0] for k in range(3))
    q_vals.append(q)
    T_new = [[Fraction(0)]*3 for _ in range(3)]
    M = M_mat_qq(N)
    for i in range(3):
        for j in range(3):
            for k in range(3):
                T_new[i][j] += T[i][k] * M[k][j]
    T = T_new

order = 3; degs = [28, 21, 14, 7]
n_unknowns = sum(d+1 for d in degs)
rows = []
for N in range(n_unknowns + 10):
    row = []
    for j in range(order + 1):
        for k in range(degs[j] + 1):
            row.append(Fraction(N)**k * q_vals[N + j])
    rows.append(row)
mat_g = [list(row) for row in rows]
n_rows = len(mat_g); n_cols = len(mat_g[0])
pivot_cols = []
for col in range(n_cols):
    pivot_row = None
    for row in range(len(pivot_cols), n_rows):
        if mat_g[row][col] != 0:
            pivot_row = row; break
    if pivot_row is None: continue
    mat_g[pivot_row], mat_g[len(pivot_cols)] = mat_g[len(pivot_cols)], mat_g[pivot_row]
    pr = len(pivot_cols); pivot_cols.append(col)
    piv = mat_g[pr][col]
    for row in range(n_rows):
        if row == pr: continue
        if mat_g[row][col] != 0:
            f = mat_g[row][col] / piv
            for c2 in range(n_cols):
                mat_g[row][c2] -= f * mat_g[pr][c2]
free_cols = [c for c in range(n_cols) if c not in pivot_cols]
null_vec = [Fraction(0)] * n_cols
null_vec[free_cols[0]] = Fraction(1)
for i in range(len(pivot_cols) - 1, -1, -1):
    pc = pivot_cols[i]
    val = -sum(mat_g[i][c] * null_vec[c] for c in range(n_cols) if c != pc) / mat_g[i][pc]
    null_vec[pc] = val
idx = 0; polys_qq = []
for j in range(order + 1):
    coeffs = null_vec[idx:idx + degs[j] + 1]
    idx += degs[j] + 1
    polys_qq.append(coeffs)

polys_mp = [[mpf(c.numerator)/mpf(c.denominator) for c in p] for p in polys_qq]

def eval_poly(coeffs, x):
    val = mpf(0)
    for c in reversed(coeffs):
        val = val * x + c
    return val

def cj(j, n):
    return eval_poly(polys_mp[j], mpf(n))

# Poincaré roots
print("Poincaré roots:")
c3_root = mpf(-272) + mpf(192)*mp.sqrt(2)
c1_root = mpf(-272) - mpf(192)*mp.sqrt(2)
c2_root = mpf(-16)
print(f"  c_dom = {nstr(c1_root, 20)}")
print(f"  c_int = {nstr(c2_root, 20)}")
print(f"  c_rec = {nstr(c3_root, 20)}")
print(f"  c_int/c_rec = {nstr(c2_root/c3_root, 20)}")

# STEP 1: Backward iteration for r_rec
N_max = 600
r_rec = {}
for N in range(N_max, N_max + 3):
    n = mpf(N)
    r_rec[N] = c3_root * n**7 * (1 + mpf(33)/2/n)  # approximate asymptotics

for n_val in range(N_max - 1, -1, -1):
    n = mpf(n_val)
    c0n, c1n, c2n, c3n = cj(0,n), cj(1,n), cj(2,n), cj(3,n)
    rn1, rn2 = r_rec[n_val+1], r_rec[n_val+2]
    denom = c1n + rn1*(c2n + c3n*rn2)
    r_rec[n_val] = -c0n / denom

print("\nStep 1: r_rec computed.")
print(f"  r_rec(0) = {nstr(r_rec[0], 20)}")
print(f"  r_rec(5) = {nstr(r_rec[5], 20)}")

# Verify: check functional equation
print("\nVerifying r_rec (functional equation residual):")
for n_val in [0, 5, 10, 50, 100]:
    n = mpf(n_val)
    F = cj(3,n)*r_rec[n_val]*r_rec[n_val+1]*r_rec[n_val+2] + cj(2,n)*r_rec[n_val]*r_rec[n_val+1] + cj(1,n)*r_rec[n_val] + cj(0,n)
    print(f"  n={n_val}: |F|/|c0| = {nstr(abs(F)/abs(cj(0,n)), 5)}")

# STEP 2: Gauge by r_rec → coefficients d_j(n)
# d0(n) = c0(n)
# d1(n) = c1(n) * r_rec(n)
# d2(n) = c2(n) * r_rec(n) * r_rec(n+1)
# d3(n) = c3(n) * r_rec(n) * r_rec(n+1) * r_rec(n+2)

def dj_coeff(j, n_val):
    n = mpf(n_val)
    cv = cj(j, n)
    prod = mpf(1)
    for i in range(j):
        prod *= r_rec[n_val + i]
    return cv * prod

# Verify (S-1) is right factor: d0+d1+d2+d3 = 0
print("\nStep 2: Verifying gauged recurrence has (S-1) right factor:")
for n_val in [0, 1, 5, 10, 50]:
    s = sum(dj_coeff(j, n_val) for j in range(4))
    d0 = dj_coeff(0, n_val)
    print(f"  n={n_val}: |d0+d1+d2+d3|/|d0| = {nstr(abs(s)/abs(d0), 5)}")

# STEP 3: Order-2 quotient L2' after dividing out (S-1)
# L2' = d3(n) S^2 + (d2(n)+d3(n)) S + (-d0(n))
# (equivalently a0 = d1+d2+d3 = -d0)

# STEP 4: Backward iteration on L2'
# a2(n)*w(n+2) + a1(n)*w(n+1) + a0(n)*w(n) = 0
# with a2 = d3, a1 = d2+d3, a0 = -d0
# Ratio σ(n) = w(n+1)/w(n):
# a2*σ*σ(n+1) + a1*σ + a0 = 0
# σ(n) = -a0(n) / [a1(n) + a2(n)*σ(n+1)]
# σ(n) = d0(n) / [(d2(n)+d3(n)) + d3(n)*σ(n+1)]

# Poincaré root of L2': the ratio r_int/r_rec ~ c_int/c_rec ≈ 34.1
target_sigma = c2_root / c3_root
print(f"\nStep 4: Expected σ ≈ r_int/r_rec = {nstr(target_sigma, 20)}")
print(f"  = -16/(-272+192√2) = 16/(272-192√2) × (272+192√2)/(272+192√2)")

# For large n, σ(n) → target_sigma (constant)
sigma = {}
for N in range(N_max, N_max + 2):
    sigma[N] = target_sigma  # start with constant asymptotic

for n_val in range(N_max - 1, -1, -1):
    d0n = dj_coeff(0, n_val)
    d2n = dj_coeff(2, n_val)
    d3n = dj_coeff(3, n_val)
    sigma[n_val] = d0n / ((d2n + d3n) + d3n * sigma[n_val + 1])

print("\nσ(n) for small n:")
for n_val in range(21):
    print(f"  σ({n_val}) = {nstr(sigma[n_val], 25)}")

# Convergence check
sigma2 = {}
for N in range(N_max, N_max + 2):
    sigma2[N] = target_sigma * (1 + mpf(1)/N)

for n_val in range(N_max - 1, -1, -1):
    d0n = dj_coeff(0, n_val)
    d2n = dj_coeff(2, n_val)
    d3n = dj_coeff(3, n_val)
    sigma2[n_val] = d0n / ((d2n + d3n) + d3n * sigma2[n_val + 1])

print("\nConvergence check (σ vs σ2):")
for n_val in [0, 1, 2, 5, 10, 50, 100]:
    diff = abs(sigma[n_val] - sigma2[n_val]) / abs(sigma[n_val])
    print(f"  n={n_val}: rel diff = {nstr(diff, 5)}")

# STEP 5: r_int(n) = σ(n) * r_rec(n)
print("\nr_int(n) = σ(n) * r_rec(n) for small n:")
for n_val in range(21):
    r_int = sigma[n_val] * r_rec[n_val]
    print(f"  r_int({n_val}) = {nstr(r_int, 25)}")

# Verify r_int satisfies the functional equation
print("\nVerifying r_int (functional equation residual):")
for n_val in [0, 1, 2, 5, 10, 20, 50]:
    r_int_n = sigma[n_val] * r_rec[n_val]
    r_int_n1 = sigma[n_val+1] * r_rec[n_val+1]
    r_int_n2 = sigma[n_val+2] * r_rec[n_val+2]
    n = mpf(n_val)
    F = cj(3,n)*r_int_n*r_int_n1*r_int_n2 + cj(2,n)*r_int_n*r_int_n1 + cj(1,n)*r_int_n + cj(0,n)
    print(f"  n={n_val}: |F|/|c0| = {nstr(abs(F)/abs(cj(0,n)), 5)}")

# Check r_int(n)/(-16*n^7) for asymptotic behavior
print("\nr_int(n)/(-16*n^7) for moderate n:")
for n_val in [5, 10, 20, 50, 100, 200]:
    r_int = sigma[n_val] * r_rec[n_val]
    n = mpf(n_val)
    ratio = r_int / (-16 * n**7)
    print(f"  n={n_val}: {nstr(ratio, 20)}")
