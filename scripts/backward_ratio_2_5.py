#!/usr/bin/env python3
"""Backward iteration to extract r(n) = h(n+1)/h(n) for the intermediate Poincaré mode.

r(n) = -c0(n) / [c1(n) + r(n+1)*(c2(n) + c3(n)*r(n+2))]

Start at large N with r(N) ~ -16*N^7*(1 + 33/(2N) + ...) and iterate backward.
Then reconstruct the rational function from integer evaluations.
"""
from fractions import Fraction
from mpmath import mp, mpf, nstr

mp.dps = 300  # high precision

# Rebuild exact recurrence coefficients
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

# Compute exact sequence and recurrence
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

# Gaussian elimination for exact recurrence coefficients
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

# Convert to mpf for evaluation
polys_mp = [[mpf(c.numerator)/mpf(c.denominator) for c in p] for p in polys_qq]

def eval_poly(coeffs, x):
    val = mpf(0)
    for c in reversed(coeffs):
        val = val * x + c
    return val

def c_j(j, n):
    return eval_poly(polys_mp[j], mpf(n))

print("Recurrence coefficients computed.")
print(f"Leading coefficients: c0={nstr(polys_mp[0][-1],5)}, c1={nstr(polys_mp[1][-1],5)}, c2={nstr(polys_mp[2][-1],5)}, c3={nstr(polys_mp[3][-1],5)}")

# Compute asymptotic expansion r(n) = -16*n^7 * (1 + s1/n + s2/n^2 + ...)
# We know s1 = 33/2 from ChatGPT Q4752
# Compute more terms by substituting into functional equation

# Method: substitute r(n) = -16*n^7*sum(s_k/n^k) into
# c3(n)*r(n)*r(n+1)*r(n+2) + c2(n)*r(n)*r(n+1) + c1(n)*r(n) + c0(n) = 0
# and solve order by order

# Instead, let's just start with first-order and iterate
N_start = 500
r_vals = {}
for N in range(N_start, N_start + 3):
    n = mpf(N)
    r_vals[N] = -16 * n**7 * (1 + mpf(33)/mpf(2)/n)

# Backward iteration
print(f"\nBackward iteration from N={N_start} to N=0...")
for n_val in range(N_start - 1, -1, -1):
    n = mpf(n_val)
    c0n = c_j(0, n)
    c1n = c_j(1, n)
    c2n = c_j(2, n)
    c3n = c_j(3, n)
    rn1 = r_vals[n_val + 1]
    rn2 = r_vals[n_val + 2]
    denom = c1n + rn1 * (c2n + c3n * rn2)
    if abs(denom) < mpf('1e-200'):
        print(f"  WARNING: denom near zero at n={n_val}")
        break
    r_vals[n_val] = -c0n / denom

print("Done.")

# Check convergence: run with a DIFFERENT starting point and compare
print("\nConvergence check: perturbing initial conditions...")
r_vals2 = {}
for N in range(N_start, N_start + 3):
    n = mpf(N)
    # Use different asymptotic (add s2 term perturbation)
    r_vals2[N] = -16 * n**7 * (1 + mpf(33)/mpf(2)/n + mpf(100)/n**2)

for n_val in range(N_start - 1, -1, -1):
    n = mpf(n_val)
    c0n = c_j(0, n)
    c1n = c_j(1, n)
    c2n = c_j(2, n)
    c3n = c_j(3, n)
    rn1 = r_vals2[n_val + 1]
    rn2 = r_vals2[n_val + 2]
    denom = c1n + rn1 * (c2n + c3n * rn2)
    if abs(denom) < mpf('1e-200'):
        break
    r_vals2[n_val] = -c0n / denom

# Compare
print("\nComparison (should agree at small n if converged):")
for n_val in [0, 1, 2, 3, 4, 5, 10, 20, 50, 100, 200]:
    if n_val in r_vals and n_val in r_vals2:
        diff = abs(r_vals[n_val] - r_vals2[n_val])
        reldiff = diff / abs(r_vals[n_val]) if abs(r_vals[n_val]) > 0 else diff
        print(f"  n={n_val}: r = {nstr(r_vals[n_val], 20)}, rel diff = {nstr(reldiff, 5)}")

# Print r(n) for small n
print("\nr(n) for n=0..20:")
for n_val in range(21):
    print(f"  r({n_val}) = {nstr(r_vals[n_val], 30)}")

# Check if r(n) / (-16 * n^7) converges for larger n
print("\nr(n) / (-16*n^7) for moderate n:")
for n_val in [5, 10, 20, 30, 50, 100]:
    n = mpf(n_val)
    ratio = r_vals[n_val] / (-16 * n**7)
    print(f"  n={n_val}: {nstr(ratio, 20)}")
