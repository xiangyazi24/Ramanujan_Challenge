#!/usr/bin/env python3
"""Problem 2.5: Compute the order-3 recurrence using modular arithmetic.

Same method as the Sage computation in recurrence_proof.txt:
1. Compute 200 terms of Q_N mod p (large prime)
2. Set up linear system for recurrence coefficients
3. Solve and use rational reconstruction

Then: check (S-1) factorization and find gauge.
"""
import sys
from fractions import Fraction

# Use a moderate prime for rational reconstruction
# The coefficients have degree up to 28, so ~29 rational coefficients per polynomial
# With Fraction, we can use multiple primes and CRT, or just exact QQ computation.

# Actually, let's just compute in QQ directly using Fraction (exact arithmetic).
# For 100 terms this should be feasible.

def M_mat_qq(n):
    """Return 3x3 matrix M(n) as list-of-lists of Fraction."""
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
    """3x3 matrix multiply."""
    C = [[Fraction(0)]*3 for _ in range(3)]
    for i in range(3):
        for j in range(3):
            for k in range(3):
                C[i][j] += A[i][k] * B[k][j]
    return C

def mat_vec_mul(A, v):
    """3x3 matrix × 3-vector."""
    return [sum(A[i][k]*v[k] for k in range(3)) for i in range(3)]

# Compute Q_N = (A · M(0) · M(1) · ... · M(N-1))_{1,0}
# A = [[30921, -32972, 8240], [33750, -36000, 9000]]
A_row = [Fraction(33750), Fraction(-36000), Fraction(9000)]  # Second row of A

print("Computing exact QQ terms of Q_N...", flush=True)
N_terms = 120  # enough for order 3, degree 28 recurrence
T = [[Fraction(1),Fraction(0),Fraction(0)],
     [Fraction(0),Fraction(1),Fraction(0)],
     [Fraction(0),Fraction(0),Fraction(1)]]
q_vals = []

for N in range(N_terms):
    # q_N = A_row · T · e_0
    q = sum(A_row[k] * T[k][0] for k in range(3))
    q_vals.append(q)
    T = mat_mul(T, M_mat_qq(N))
    if N % 20 == 0:
        num_digits = len(str(abs(q.numerator)))
        print(f"  N={N}: {num_digits} digits", flush=True)

print(f"Computed {len(q_vals)} terms", flush=True)

# Now find the recurrence: order 3, degree pattern (28, 21, 14, 7)
# c₀(N)·Q(N) + c₁(N)·Q(N+1) + c₂(N)·Q(N+2) + c₃(N)·Q(N+3) = 0
# c₃ has degree 7 → 8 unknowns
# c₂ has degree 14 → 15 unknowns
# c₁ has degree 21 → 22 unknowns
# c₀ has degree 28 → 29 unknowns
# Total: 74 unknowns, minus 1 normalization = 73 free parameters
# Need at least 73 equations (values of N)

order = 3
degs = [28, 21, 14, 7]  # degrees of c₀, c₁, c₂, c₃
n_unknowns = sum(d+1 for d in degs)  # 29+22+15+8 = 74
print(f"\nSetting up linear system: {n_unknowns} unknowns", flush=True)

# Build the linear system
# For each N, one equation: Σ_{j=0}^{3} c_j(N) · Q(N+j) = 0
# c_j(N) = Σ_{k=0}^{degs[j]} a_{j,k} · N^k
# So the equation becomes: Σ_j Σ_k a_{j,k} · N^k · Q(N+j) = 0

n_eqs = N_terms - order  # can use N = 0, ..., N_terms - order - 1
print(f"Available equations: {n_eqs}", flush=True)

# Build matrix A_sys and vector b (Ax = 0, normalize by setting one coeff = 1)
rows = []
for N in range(min(n_eqs, n_unknowns + 10)):
    row = []
    for j in range(order + 1):
        qNj = q_vals[N + j]
        for k in range(degs[j] + 1):
            row.append(Fraction(N)**k * qNj)
    rows.append(row)

print(f"Matrix size: {len(rows)} × {len(rows[0])}", flush=True)

# Gaussian elimination to find the null space
# Since we expect rank = n_unknowns - 1, the null space is 1-dimensional
import copy

mat = [list(row) for row in rows]
n_rows = len(mat)
n_cols = len(mat[0])

print("Running Gaussian elimination...", flush=True)
pivot_cols = []
for col in range(n_cols):
    # Find pivot
    pivot_row = None
    for row in range(len(pivot_cols), n_rows):
        if mat[row][col] != 0:
            pivot_row = row
            break
    if pivot_row is None:
        continue
    
    # Swap
    mat[pivot_row], mat[len(pivot_cols)] = mat[len(pivot_cols)], mat[pivot_row]
    pivot_row = len(pivot_cols)
    pivot_cols.append(col)
    
    # Eliminate
    piv = mat[pivot_row][col]
    for row in range(n_rows):
        if row == pivot_row:
            continue
        if mat[row][col] != 0:
            factor = mat[row][col] / piv
            for c in range(n_cols):
                mat[row][c] -= factor * mat[pivot_row][c]
    
    if len(pivot_cols) % 10 == 0:
        print(f"  {len(pivot_cols)} pivots found...", flush=True)

rank = len(pivot_cols)
null_dim = n_cols - rank
print(f"Rank = {rank}, null space dim = {null_dim}", flush=True)

if null_dim != 1:
    print(f"WARNING: Expected null dim = 1, got {null_dim}")
    print("Trying with fewer unknowns or checking data...")
    sys.exit(1)

# Extract the null vector
# Set the free variable (non-pivot column) to 1, back-substitute
free_cols = [c for c in range(n_cols) if c not in pivot_cols]
null_vec = [Fraction(0)] * n_cols
null_vec[free_cols[0]] = Fraction(1)

# Back-substitute
for i in range(rank - 1, -1, -1):
    pc = pivot_cols[i]
    val = -sum(mat[i][c] * null_vec[c] for c in range(n_cols) if c != pc) / mat[i][pc]
    null_vec[pc] = val

# Extract polynomial coefficients
print("\nExtracting recurrence coefficients...", flush=True)
idx = 0
polys = []
for j in range(order + 1):
    coeffs = null_vec[idx:idx + degs[j] + 1]
    idx += degs[j] + 1
    polys.append(coeffs)
    print(f"  c_{j}: degree {degs[j]}, LC = {coeffs[-1]}")

# Verify on a few terms
print("\nVerifying recurrence...", flush=True)
for N in [50, 60, 70, 80, 90]:
    if N + order < len(q_vals):
        val = Fraction(0)
        for j in range(order + 1):
            cj = sum(polys[j][k] * Fraction(N)**k for k in range(len(polys[j])))
            val += cj * q_vals[N + j]
        print(f"  N={N}: residual = {val}")

# Poincaré polynomial
print("\nPoincaré polynomial:")
lcs = [p[-1] for p in polys]
print(f"  Leading coefficients: c₀={lcs[0]}, c₁={lcs[1]}, c₂={lcs[2]}, c₃={lcs[3]}")
# Normalize
for i in range(4):
    print(f"  LC(c_{i})/LC(c_3) = {lcs[i]/lcs[3]}")

# Check (S-1) condition: sum of all coefficients = 0 for each N
# c₀(N) + c₁(N) + c₂(N) + c₃(N) = 0 for all N
# This means the sum polynomial = 0
print("\n=== Checking (S-1) factor condition ===")
sum_coeffs = [Fraction(0)] * (max(degs) + 1)
for j in range(order + 1):
    for k in range(len(polys[j])):
        sum_coeffs[k] += polys[j][k]

nonzero = [(k, c) for k, c in enumerate(sum_coeffs) if c != 0]
if not nonzero:
    print("  c₀+c₁+c₂+c₃ ≡ 0  ✓  (S-1) IS a right factor of the RAW recurrence!")
    has_factor = True
else:
    print(f"  c₀+c₁+c₂+c₃ ≠ 0 ({len(nonzero)} nonzero coefficients)")
    print(f"  First few: {nonzero[:5]}")
    has_factor = False

