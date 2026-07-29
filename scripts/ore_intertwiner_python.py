#!/usr/bin/env python3
"""
Problem 2.5: Compute Ore intertwiner between L_25 and Sym^2(Delannoy)
using pure Python (no Sage/ore_algebra needed).

Strategy:
1. Compute q_N from CMF matrix using exact integer arithmetic
2. Extract L_25 recurrence coefficients via linear algebra
3. Implement basic Ore algebra operations over Q(n)
4. Compute GCRD(L_25, L_D2)
"""
from fractions import Fraction as F
from functools import reduce

# ---- CMF matrix M(n) with exact integer entries ----
def M_int(n):
    """Return 3x3 matrix M(n) as list of lists of Python ints."""
    m11 = (-2*n-5)*(n+3)**2 * (136*n**4 + 1424*n**3 + 5548*n**2 + 9551*n + 6141)
    m12 = 384*n**6 + 6384*n**5 + 44168*n**4 + 162698*n**3 + 336377*n**2 + 369933*n + 169011
    m13 = -480*n**4 - 4980*n**3 - 19210*n**2 - 32690*n - 20730

    m21 = (n+2)**2*(n+3)**2*(4*n+10)*(48*n**3 + 386*n**2 + 1017*n + 879)
    m22 = (n+2)**2*(-272*n**5 - 3848*n**4 - 21732*n**3 - 61184*n**2 - 85761*n - 47808)
    m23 = (n+2)**2*(320*n**3 + 2540*n**2 + 6610*n + 5640)

    m31 = (-4*n-10)*(n+2)**2*(n+3)**2*(32*n**4 + 302*n**3 + 1037*n**2 + 1530*n + 813)
    m32 = (n+2)**2*(192*n**6 + 2984*n**5 + 19116*n**4 + 64452*n**3 + 120256*n**2 + 117279*n + 46476)
    m33 = (n+2)**2*(-16*n**5 - 408*n**4 - 2912*n**3 - 8884*n**2 - 12254*n - 6240)

    return [[m11, m12, m13],
            [m21, m22, m23],
            [m31, m32, m33]]

def mat_mul(A, B):
    """Multiply two 3x3 integer matrices."""
    n = len(A)
    return [[sum(A[i][k]*B[k][j] for k in range(n)) for j in range(n)] for i in range(n)]

def mat_vec(A, v):
    """Multiply 3x3 matrix by 3-vector."""
    return [sum(A[i][k]*v[k] for k in range(3)) for i in range(3)]

# ---- Compute scalar sequences q_N ----
print("Computing CMF scalar sequences with exact integer arithmetic...")

# Start with identity, accumulate M(0)*M(1)*...*M(N-1)
N_MAX = 120  # Need ~80 for degree-28 recurrence extraction
I3 = [[1,0,0],[0,1,0],[0,0,1]]

# Store all 3 columns of the product matrix for 3 independent solutions
# q[k][N] = [M(0)...M(N-1)]_{k,0} (first column)
prod = [row[:] for row in I3]  # deep copy
q = [[], [], []]  # q[row][N]

for N in range(N_MAX + 5):
    for k in range(3):
        q[k].append(prod[k][0])
    if N < N_MAX + 4:
        MN = M_int(N)
        prod = mat_mul(prod, MN)

print(f"Computed q[0..{len(q[0])-1}] (3 independent solutions)")
print(f"q[0][0..4] = {q[0][:5]}")
print(f"q[1][0..4] = {q[1][:5]}")
print(f"q[2][0..4] = {q[2][:5]}")
print(f"|q[0][50]| has {len(str(abs(q[0][50])))} digits")

# ---- Extract recurrence coefficients ----
# Order 3, degree pattern (28, 21, 14, 7)
# c₀(N) has degree 28, c₁(N) degree 21, c₂(N) degree 14, c₃(N) degree 7
# Total coefficients: 29 + 22 + 15 + 8 = 74

print("\n=== Extracting L_25 recurrence ===")

deg = [28, 21, 14, 7]
n_coeffs = [d+1 for d in deg]
total_unknowns = sum(n_coeffs)  # 74
print(f"Degree pattern: {deg}")
print(f"Coefficients per poly: {n_coeffs}")
print(f"Total unknowns: {total_unknowns}")

# Use q[0] (first row of product) as the representative sequence
# Build the linear system: for each N, sum_{i=0}^{3} c_i(N) * q[N+i] = 0
# where c_i(N) = sum_{m=0}^{deg[i]} a_{i,m} * N^m

# We need at least total_unknowns equations
# Each equation at N: sum_i sum_m a_{i,m} * N^m * q[0][N+i] = 0
# But this is homogeneous, so fix one coefficient (say the leading coeff of c_3)

# First try with Fraction for exact arithmetic
use_row = 0  # use q[0]
n_eqns = total_unknowns + 10  # overdetermined

# Build matrix A_sys[eq, unknown] = N^m * q[use_row][N+i]
A_sys = []
for eq_idx in range(n_eqns):
    N = eq_idx
    row = []
    for i in range(4):  # c_0, c_1, c_2, c_3
        for m in range(deg[i]+1):
            row.append(F(N**m * q[use_row][N+i]))
    A_sys.append(row)

print(f"System: {len(A_sys)} equations × {total_unknowns} unknowns (Fraction)")

# Solve via Gaussian elimination to find the null space
# Since system is homogeneous, find one nonzero solution
# Use the approach: fix a_{3, deg[3]} = 1 (leading coeff of c_3)

# Move the last column (a_{3,7}) to RHS
A_reduced = []
b_reduced = []
for row in A_sys:
    A_reduced.append(row[:-1])
    b_reduced.append(-row[-1])

n_vars = total_unknowns - 1  # 73

print("Solving with Gaussian elimination over Q...")

# Gaussian elimination
mat = [A_reduced[i][:] + [b_reduced[i]] for i in range(len(A_reduced))]
n_rows = len(mat)
n_cols = n_vars

pivot_row = 0
for col in range(n_cols):
    # Find pivot
    found = -1
    for r in range(pivot_row, n_rows):
        if mat[r][col] != 0:
            found = r
            break
    if found == -1:
        continue
    mat[found], mat[pivot_row] = mat[pivot_row], mat[found]
    # Eliminate
    pivot_val = mat[pivot_row][col]
    for r in range(n_rows):
        if r != pivot_row and mat[r][col] != 0:
            factor = mat[r][col] / pivot_val
            for c in range(n_cols + 1):
                mat[r][c] -= factor * mat[pivot_row][c]
    pivot_row += 1

# Back-substitute
solution = [F(0)] * n_vars
for r in range(min(pivot_row, n_rows)):
    # Find the pivot column
    pivot_col = -1
    for c in range(n_cols):
        if mat[r][c] != 0:
            pivot_col = c
            break
    if pivot_col >= 0:
        solution[pivot_col] = mat[r][n_cols] / mat[r][pivot_col]

# Add the fixed variable
solution.append(F(1))  # a_{3,7} = 1

# Extract polynomial coefficients
print("\nExtracted recurrence coefficients:")
idx = 0
polys = []
for i in range(4):
    coeffs = solution[idx:idx+deg[i]+1]
    polys.append(coeffs)
    idx += deg[i]+1
    # Print leading and constant terms
    print(f"  c_{i}(N): degree {deg[i]}, "
          f"leading = {coeffs[-1]}, constant = {coeffs[0]}")

# ---- Verify recurrence ----
print("\n=== Verifying recurrence ===")

def eval_poly(coeffs, N):
    """Evaluate polynomial with Fraction coefficients at integer N."""
    return sum(c * F(N)**m for m, c in enumerate(coeffs))

max_residual = F(0)
for N in range(min(80, len(q[use_row])-3)):
    val = sum(eval_poly(polys[i], N) * q[use_row][N+i] for i in range(4))
    if val != 0:
        print(f"  NONZERO residual at N={N}: {val}")
        max_residual = max(max_residual, abs(val))
        if N > 5:
            break

if max_residual == 0:
    print("  All residuals are exactly ZERO. Recurrence verified!")
else:
    print(f"  Max residual: {max_residual}")

# ---- Normalize: make c_3 monic (leading coeff = 1) ----
# Already done since we fixed a_{3,7} = 1

# ---- Poincaré polynomial check ----
# Leading coefficients of c_i give the Poincaré polynomial
# c_i(N) ~ alpha_i * N^{deg[i]} as N → infinity
# Poincaré eq: alpha_0 + alpha_1*c + alpha_2*c^2 + alpha_3*c^3 = 0
print("\n=== Poincaré polynomial ===")
alphas = [float(polys[i][-1]) for i in range(4)]
print(f"Leading coefficients: {alphas}")
print(f"Poincaré polynomial: {alphas[3]:.6g}·c³ + {alphas[2]:.6g}·c² + "
      f"{alphas[1]:.6g}·c + {alphas[0]:.6g}")

# Normalize to monic
a3 = polys[3][-1]
poincare = [polys[i][-1] / a3 for i in range(4)]
print(f"Monic: c³ + {float(poincare[2]):.4f}·c² + "
      f"{float(poincare[1]):.4f}·c + {float(poincare[0]):.4f}")

# ---- Now compute GCRD with L_D2 ----
# L_D2 (Sym^2 Delannoy):
# (n+3)^2(2n+3) U_{n+3} - (2n+5)(35n^2+140n+131) U_{n+2}
# + (2n+3)(35n^2+140n+131) U_{n+1} - (2n+5)(n+1)^2 U_n = 0

print("\n=== Ore intertwiner via numerical identification ===")

# Compute Delannoy numbers D_n
N_DEL = N_MAX + 10
D = [F(0)] * N_DEL
D[0] = F(1)
D[1] = F(3)
for n in range(1, N_DEL-1):
    D[n+1] = (F(3*(2*n+1)) * D[n] - F(n) * D[n-1]) / F(n+1)

# Second solution E_n (linearly independent)
E = [F(0)] * N_DEL
E[0] = F(0)
E[1] = F(1)
for n in range(1, N_DEL-1):
    E[n+1] = (F(3*(2*n+1)) * E[n] - F(n) * E[n-1]) / F(n+1)

print(f"D[0..5] = {[D[i] for i in range(6)]}")
print(f"E[0..5] = {[E[i] for i in range(6)]}")

# Sym^2 solutions: u1 = D^2, u2 = D*E, u3 = E^2
# Verify they satisfy L_D2
def L_D2_eval(u, n):
    """Evaluate L_D2[u](n)."""
    return ((n+3)**2*(2*n+3)*u[n+3]
            - (2*n+5)*(35*n**2+140*n+131)*u[n+2]
            + (2*n+3)*(35*n**2+140*n+131)*u[n+1]
            - (2*n+5)*(n+1)**2*u[n])

u1 = [D[n]**2 for n in range(N_DEL)]
u2 = [D[n]*E[n] for n in range(N_DEL)]
u3 = [E[n]**2 for n in range(N_DEL)]

print("\nVerifying Sym^2 solutions satisfy L_D2:")
for name, u in [("D^2", u1), ("D*E", u2), ("E^2", u3)]:
    ok = all(L_D2_eval(u, n) == 0 for n in range(min(20, N_DEL-3)))
    print(f"  {name}: {'OK' if ok else 'FAIL'}")

# ---- Find the intertwiner T ----
# T = t_0(n) + t_1(n)*S + t_2(n)*S^2 maps solutions of L_D2 to solutions of L_25
#
# For each N:
#   t_0(N)*u_j(N) + t_1(N)*u_j(N+1) + t_2(N)*u_j(N+2) = sum_k a_{jk} * q[k][N]
#
# for j = 0,1,2 (three Sym^2 solutions) and k = 0,1,2 (three L_25 solutions)
#
# If t_i(N) are rational functions, try polynomial ansatz first.
# At each N, we have a 3x3 system (for each j):
#
# [ u1(N)    u1(N+1)  u1(N+2) ] [t_0(N)]     [q0(N)  q1(N)  q2(N)] [a_00]
# [ u2(N)    u2(N+1)  u2(N+2) ] [t_1(N)]  =  [q0(N)  q1(N)  q2(N)] [a_10]
# [ u3(N)    u3(N+1)  u3(N+2) ] [t_2(N)]     [q0(N)  q1(N)  q2(N)] [a_20]
#
# Wait, this isn't right. Let me reformulate.
#
# For each j and each N:
# t_0(N)*u_j(N) + t_1(N)*u_j(N+1) + t_2(N)*u_j(N+2) = sum_k a_{jk} * q_k(N)
#
# This is 1 equation per (j, N), with unknowns t_0(N), t_1(N), t_2(N), a_{jk}.
# The a_{jk} are 9 constants, and t_i(N) = poly(N) with degree d gives 3*(d+1) unknowns.
#
# Let's solve it differently. At fixed N, given the 3x3 systems for j=0,1,2:
#
# U_matrix(N) * t_vec(N) = Q_matrix(N) * a_vec_j  (one column of A at a time)
#
# Actually, let me think of it as: define w_j(N) = T(u_j)(N) = t_0*u_j(N) + t_1*u_j(N+1) + t_2*u_j(N+2)
# This w_j must be a solution of L_25. The general solution of L_25 is spanned by q_0, q_1, q_2.
# So w_j = sum_k a_{jk} q_k.
#
# At each N: w_j(N) = sum_k a_{jk} q_k(N)
# And w_j(N) = t_0(N) u_j(N) + t_1(N) u_j(N+1) + t_2(N) u_j(N+2)
#
# Strategy: first fix A by initial conditions, then solve for t_i.
#
# Actually simpler: at each N, the 3 equations (j=0,1,2) can be written as:
# U(N) * t(N) = Q(N) * a_col
# where U(N) is the 3x3 matrix with rows [u_j(N), u_j(N+1), u_j(N+2)]
# and Q(N) is the 3x3 matrix with columns [q_0(N), q_1(N), q_2(N)]
# and a_col = [a_{j0}, a_{j1}, a_{j2}]^T ... no this still mixes things.
#
# Let me write it as a matrix equation:
# Define U(N) = 3x3 matrix with U[j,i] = u_{j+1}(N+i) for i=0,1,2
# Define Q(N) = 3x3 matrix with Q[j,k] = q_k(N) for k=0,1,2 (same row repeated)
#
# No wait. Let me just be concrete.
#
# At each N, we want:
# for j = 0: t_0*u1(N) + t_1*u1(N+1) + t_2*u1(N+2) = a_{00}*q0(N) + a_{01}*q1(N) + a_{02}*q2(N)
# for j = 1: t_0*u2(N) + t_1*u2(N+1) + t_2*u2(N+2) = a_{10}*q0(N) + a_{11}*q1(N) + a_{12}*q2(N)
# for j = 2: t_0*u3(N) + t_1*u3(N+1) + t_2*u3(N+2) = a_{20}*q0(N) + a_{21}*q1(N) + a_{22}*q2(N)
#
# In matrix form:
# [u1(N)  u1(N+1)  u1(N+2)]   [t_0]     [q0(N)  q1(N)  q2(N)]   [a_00  a_01  a_02]^T
# [u2(N)  u2(N+1)  u2(N+2)] * [t_1]  =  [q0(N)  q1(N)  q2(N)] * [a_10  a_11  a_12]^T  <-- NO
# [u3(N)  u3(N+1)  u3(N+2)]   [t_2]     [q0(N)  q1(N)  q2(N)]   [a_20  a_21  a_22]^T
#
# Hmm, the right-hand side should be:
# row j = sum_k a_{jk} q_k(N)
# So RHS vector = A * [q0(N), q1(N), q2(N)]^T where A is the 3x3 connection matrix.
#
# So: U(N) * t(N) = A * q(N)
# where U(N)[j,i] = u_{j+1}(N+i), q(N) = [q0(N), q1(N), q2(N)]^T, t(N) = [t_0(N), t_1(N), t_2(N)]^T
#
# This gives: t(N) = U(N)^{-1} * A * q(N)
#
# For t(N) to be a polynomial/rational function of N, the expression
# U(N)^{-1} * A * q(N) must simplify appropriately.
#
# Approach: compute t(N) = U(N)^{-1} * A * q(N) for many N values,
# and identify the rational function form.
# But we don't know A!
#
# Alternative: from two different N values (N=N1, N=N2):
# U(N1) * t(N1) = A * q(N1)
# U(N2) * t(N2) = A * q(N2)
#
# If t(N) = t_const (degree 0), then t(N1) = t(N2) = t, and:
# U(N1) * t = A * q(N1)
# U(N2) * t = A * q(N2)
# So A = U(N1) * t * q(N1)^{-1} ... but q(N) is a vector, not invertible.
#
# Let me use a different approach. Stack equations from multiple N values.

# For polynomial degree d, t_i(N) = sum_{m=0}^d c_{im} N^m
# Unknowns: c_{im} (3*(d+1) values) and a_{jk} (9 values)
# At each N, 3 equations.
# Total equations from N=0..K-1: 3*K
# Need 3*K >= 3*(d+1) + 9, i.e., K >= d+4

# The system is homogeneous (scale t and A together), so fix one unknown.

# Try d = 0 first (constant t_i, should fail), then increase
for d in range(15):
    n_t_coeffs = 3 * (d + 1)  # c_{00},...,c_{0d}, c_{10},...,c_{2d}
    n_a_coeffs = 9             # a_{00},...,a_{22}
    n_total = n_t_coeffs + n_a_coeffs
    K = d + 6  # number of N values to use (overdetermined)
    n_eqns = 3 * K

    if n_eqns < n_total:
        K = (n_total + 2) // 3 + 1
        n_eqns = 3 * K

    # Build the system
    # For each N in 0..K-1 and j in 0..2:
    # sum_m sum_i c_{im} * N^m * u_{j+1}(N+i) - sum_k a_{jk} * q_k(N) = 0

    sys_mat = []
    for N in range(1, K+1):  # start at N=1 to avoid potential singularity at N=0
        for j in range(3):
            row = [F(0)] * n_total
            u_seqs = [u1, u2, u3]
            uj = u_seqs[j]

            # t coefficients: c_{im} contributes N^m * u_j(N+i)
            for i in range(3):
                for m in range(d+1):
                    col = i * (d+1) + m
                    row[col] = F(N**m) * uj[N+i]

            # a coefficients: a_{jk} contributes -q_k(N)
            for k in range(3):
                col = n_t_coeffs + j * 3 + k
                row[col] = -F(q[k][N])

            sys_mat.append(row)

    # Fix last a coefficient (a_{22}) = 1 to break homogeneity
    # Move last column to RHS
    A_mat = [row[:-1] for row in sys_mat]
    b_vec = [-row[-1] for row in sys_mat]
    n_unknowns = n_total - 1

    # Solve via least squares (Gaussian elimination)
    # For exact Fraction arithmetic, this is exact
    aug = [A_mat[i][:] + [b_vec[i]] for i in range(len(A_mat))]
    n_r = len(aug)
    n_c = n_unknowns

    piv_row = 0
    pivots = []
    for col in range(n_c):
        found = -1
        for r in range(piv_row, n_r):
            if aug[r][col] != 0:
                found = r
                break
        if found == -1:
            continue
        aug[found], aug[piv_row] = aug[piv_row], aug[found]
        pivot_val = aug[piv_row][col]
        for r in range(n_r):
            if r != piv_row and aug[r][col] != 0:
                factor = aug[r][col] / pivot_val
                for c2 in range(n_c + 1):
                    aug[r][c2] -= factor * aug[piv_row][c2]
        pivots.append((piv_row, col))
        piv_row += 1

    # Check consistency: are remaining rows all zero?
    consistent = True
    for r in range(piv_row, n_r):
        if aug[r][n_c] != 0:
            consistent = False
            break

    if not consistent:
        if d <= 3:
            print(f"  d={d}: INCONSISTENT (no solution)")
        continue

    # Check rank
    rank = len(pivots)
    if rank < n_c:
        if d <= 3:
            print(f"  d={d}: UNDERDETERMINED (rank {rank} < {n_c})")
        continue

    # Extract solution
    sol = [F(0)] * n_c
    for pr, pc in reversed(pivots):
        sol[pc] = aug[pr][n_c] / aug[pr][pc]
    sol.append(F(1))  # the fixed a_{22}

    # Verify on additional data points
    verified = True
    for N in range(K+1, K+20):
        if N+2 >= len(u1) or N >= len(q[0]):
            break
        for j in range(3):
            uj = [u1, u2, u3][j]
            lhs = sum(sum(sol[i*(d+1)+m] * F(N**m) for m in range(d+1)) * uj[N+i] for i in range(3))
            rhs = sum(sol[n_t_coeffs + j*3 + k] * F(q[k][N]) for k in range(3))
            if lhs != rhs:
                verified = False
                break
        if not verified:
            break

    if verified:
        print(f"\n*** FOUND INTERTWINER at degree d={d} ***")

        # Extract t_i polynomials
        for i in range(3):
            coeffs = sol[i*(d+1):(i+1)*(d+1)]
            print(f"  t_{i}(N) coefficients: {[str(c) for c in coeffs]}")

        # Extract connection matrix A
        print("  Connection matrix A:")
        for j in range(3):
            row_a = sol[n_t_coeffs + j*3:n_t_coeffs + (j+1)*3]
            print(f"    row {j}: {[str(c) for c in row_a]}")

        break
    else:
        if d <= 5:
            print(f"  d={d}: solution found but does NOT verify on additional points")

print("\nDone.")
