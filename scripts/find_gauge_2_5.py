#!/usr/bin/env python3
"""Problem 2.5: Find the gauge r(n) = -16*P(n) where P(n) = (n+a1)...(n+a7).

Uses the exact QQ recurrence coefficients and Newton's method.

Functional equation: c3(n)*r(n)*r(n+1)*r(n+2) + c2(n)*r(n)*r(n+1) + c1(n)*r(n) + c0(n) = 0
"""
from fractions import Fraction
from mpmath import mp, mpf, nstr, matrix as mp_matrix, lu_solve, eye, sqrt, log
mp.dps = 100

# First, recompute the exact QQ recurrence coefficients (same as modular_recurrence_2_5.py)
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

def mat_mul(A, B):
    C = [[Fraction(0)]*3 for _ in range(3)]
    for i in range(3):
        for j in range(3):
            for k in range(3):
                C[i][j] += A[i][k] * B[k][j]
    return C

A_row = [Fraction(33750), Fraction(-36000), Fraction(9000)]

print("Computing exact QQ terms...", flush=True)
N_terms = 120
T = [[Fraction(1 if i==j else 0) for j in range(3)] for i in range(3)]
q_vals = []
for N in range(N_terms):
    q = sum(A_row[k] * T[k][0] for k in range(3))
    q_vals.append(q)
    T_new = [[Fraction(0)]*3 for _ in range(3)]
    M = M_mat_qq(N)
    for i in range(3):
        for j in range(3):
            for k in range(3):
                T_new[i][j] += T[i][k] * M[k][j]
    T = T_new
    if N % 40 == 0:
        print(f"  N={N} done", flush=True)

# Build exact recurrence
print("Building recurrence...", flush=True)
order = 3; degs = [28, 21, 14, 7]
n_unknowns = sum(d+1 for d in degs)
rows = []
for N in range(n_unknowns + 10):
    row = []
    for j in range(order + 1):
        qNj = q_vals[N + j]
        for k in range(degs[j] + 1):
            row.append(Fraction(N)**k * qNj)
    rows.append(row)

# Gaussian elimination
mat = [list(row) for row in rows]
n_rows = len(mat); n_cols = len(mat[0])
pivot_cols = []
for col in range(n_cols):
    pivot_row = None
    for row in range(len(pivot_cols), n_rows):
        if mat[row][col] != 0:
            pivot_row = row; break
    if pivot_row is None: continue
    mat[pivot_row], mat[len(pivot_cols)] = mat[len(pivot_cols)], mat[pivot_row]
    pivot_row = len(pivot_cols); pivot_cols.append(col)
    piv = mat[pivot_row][col]
    for row in range(n_rows):
        if row == pivot_row: continue
        if mat[row][col] != 0:
            factor = mat[row][col] / piv
            for c in range(n_cols):
                mat[row][c] -= factor * mat[pivot_row][c]

free_cols = [c for c in range(n_cols) if c not in pivot_cols]
null_vec = [Fraction(0)] * n_cols
null_vec[free_cols[0]] = Fraction(1)
for i in range(len(pivot_cols) - 1, -1, -1):
    pc = pivot_cols[i]
    val = -sum(mat[i][c] * null_vec[c] for c in range(n_cols) if c != pc) / mat[i][pc]
    null_vec[pc] = val

idx = 0; polys_qq = []
for j in range(order + 1):
    coeffs = null_vec[idx:idx + degs[j] + 1]
    idx += degs[j] + 1
    polys_qq.append(coeffs)

print(f"Recurrence found, LCs: {[p[-1] for p in polys_qq]}", flush=True)

# Convert to mpf polynomials for Newton's method
def eval_poly(coeffs, x):
    """Evaluate polynomial at x (horner)."""
    val = mpf(0)
    for c in reversed(coeffs):
        val = val * x + mpf(c)
    return val

# Store polynomials as mpf coefficient lists
polys = [[mpf(c.numerator) / mpf(c.denominator) for c in p] for p in polys_qq]

def c_j(j, n):
    return eval_poly(polys[j], n)

# The functional equation:
# F(n; σ) = c3(n)*r(n)*r(n+1)*r(n+2) + c2(n)*r(n)*r(n+1) + c1(n)*r(n) + c0(n) = 0
# where r(n) = -16 * P(n) and P(n) = n^7 + σ1*n^6 + ... + σ7

def P_val(n, sigma):
    """P(n) = n^7 + σ1*n^6 + ... + σ7."""
    val = mpf(1)
    for i in range(7):
        val = val * n + sigma[i]
    return val

def r_val(n, sigma):
    return mpf(-16) * P_val(n, sigma)

def F_val(n, sigma):
    """Evaluate the functional equation at n with parameters sigma."""
    r0 = r_val(n, sigma)
    r1 = r_val(n+1, sigma)
    r2 = r_val(n+2, sigma)
    return c_j(3, n)*r0*r1*r2 + c_j(2, n)*r0*r1 + c_j(1, n)*r0 + c_j(0, n)

# Initial guess: from the structure of c0(N), likely parameters are
# half-integers like 1, 3/2, 2, 5/2, 3, 7/2, 4
# P(n) = (n+1)(n+3/2)(n+2)(n+5/2)(n+3)(n+7/2)(n+4) has σ from elementary symmetric functions
# Let me compute these σ values
from itertools import combinations
from functools import reduce
import operator

guess_roots = [mpf(1), mpf('1.5'), mpf(2), mpf('2.5'), mpf(3), mpf('3.5'), mpf(4)]
# Elementary symmetric functions σ_k = (-1)^k * e_k(roots)
# P(n) = (n+a1)...(n+a7) = n^7 + (Σai)n^6 + (Σ ai*aj)n^5 + ... + Π(ai)
# So σ1 = Σai, σ2 = Σ_{i<j} ai*aj, ..., σ7 = Πai

def elem_sym(roots):
    """Compute elementary symmetric functions."""
    n = len(roots)
    sigma = []
    for k in range(1, n+1):
        s = mpf(0)
        for combo in combinations(roots, k):
            s += reduce(operator.mul, combo, mpf(1))
        sigma.append(s)
    return sigma

sigma0 = elem_sym(guess_roots)
print(f"\nInitial guess (roots 1,3/2,...,4): σ = {[nstr(s,10) for s in sigma0]}")

# Test F at several n values
print("\nF values at initial guess:")
for n in range(7):
    f = F_val(n, sigma0)
    print(f"  F({n}) = {nstr(f, 15)}")

# Newton's method: solve F(n_i; σ) = 0 for n_i = 0,...,6
# Jacobian: ∂F(n_i)/∂σ_j computed by finite differences

def newton_step(sigma):
    """One Newton step. Returns new sigma and max residual."""
    residuals = [F_val(n, sigma) for n in range(7)]

    # Jacobian by central differences
    eps = mpf('1e-40')
    J = mp_matrix(7, 7)
    for j in range(7):
        sigma_plus = list(sigma)
        sigma_minus = list(sigma)
        sigma_plus[j] += eps
        sigma_minus[j] -= eps
        for i in range(7):
            fp = F_val(i, sigma_plus)
            fm = F_val(i, sigma_minus)
            J[i, j] = (fp - fm) / (2 * eps)

    # Solve J · delta = -residuals
    b = mp_matrix(7, 1)
    for i in range(7):
        b[i, 0] = -residuals[i]

    delta = lu_solve(J, b)

    new_sigma = [sigma[j] + delta[j, 0] for j in range(7)]
    max_res = max(abs(r) for r in residuals)
    return new_sigma, max_res

print("\nRunning Newton's method...")
sigma = list(sigma0)
for iteration in range(30):
    sigma, max_res = newton_step(sigma)
    print(f"  Iter {iteration}: max |F| = {nstr(max_res, 5)}")
    if max_res < mpf('1e-80'):
        break

# Extract the roots of P(n) from sigma
print(f"\nFinal σ = {[nstr(s, 30) for s in sigma]}")

# Find the roots of P(n) = n^7 + σ1*n^6 + ... + σ7
# Using numpy for root-finding
print("\nRoots of P(n):")
import numpy as np
coeffs_np = [1.0] + [float(s) for s in sigma]
roots = np.roots(coeffs_np)
for i, r in enumerate(sorted(roots, key=lambda x: x.real)):
    print(f"  a_{i+1} = {r.real:.15f}" + (f" + {r.imag:.15f}i" if abs(r.imag) > 1e-10 else ""))

# Check: are the roots close to simple fractions?
print("\nClosest simple fractions:")
for r in sorted(roots, key=lambda x: x.real):
    if abs(r.imag) < 1e-10:
        x = r.real
        # Check half-integers
        best = None; best_err = 1e10
        for num in range(-10, 30):
            for den in [1, 2, 3, 4, 5, 6]:
                frac = num/den
                err = abs(x - frac)
                if err < best_err:
                    best = f"{num}/{den}" if den > 1 else str(num)
                    best_err = err
        print(f"  {x:.15f} ≈ {best} (error {best_err:.3e})")

# Verify: check F at many more n values
print("\nVerification at n = 0,...,30:")
max_res = mpf(0)
for n in range(31):
    f = F_val(n, sigma)
    if abs(f) > max_res:
        max_res = abs(f)
print(f"  Max |F| over n=0..30: {nstr(max_res, 5)}")
