#!/usr/bin/env python3
"""
Verify that q♯_k (k=0,1,2) all satisfy the same twisted recurrence.
Also try null-space approach for finding the Ore intertwiner.
"""
from fractions import Fraction as F
from functools import reduce
from math import gcd

def pochhammer_int(a_num, a_den, n):
    result = F(1)
    for k in range(n):
        result *= F(a_num + k * a_den, a_den)
    return result

def H(n):
    if n == 0:
        return F(1)
    neg16_n = F(-16)**n
    poch_2 = pochhammer_int(2, 1, n)
    poch_3 = pochhammer_int(3, 1, n)
    poch_5_2 = pochhammer_int(5, 2, n)
    poch_7_2 = pochhammer_int(7, 2, n)
    return neg16_n * poch_2**2 * poch_3**2 * poch_5_2 * poch_7_2**2

def M_int(n):
    m11 = (-2*n-5)*(n+3)**2 * (136*n**4 + 1424*n**3 + 5548*n**2 + 9551*n + 6141)
    m12 = 384*n**6 + 6384*n**5 + 44168*n**4 + 162698*n**3 + 336377*n**2 + 369933*n + 169011
    m13 = -480*n**4 - 4980*n**3 - 19210*n**2 - 32690*n - 20730
    m21 = (n+2)**2*(n+3)**2*(4*n+10)*(48*n**3 + 386*n**2 + 1017*n + 879)
    m22 = (n+2)**2*(-272*n**5 - 3848*n**4 - 21732*n**3 - 61184*n**2 - 85761*n - 47808)
    m23 = (n+2)**2*(320*n**3 + 2540*n**2 + 6610*n + 5640)
    m31 = (-4*n-10)*(n+2)**2*(n+3)**2*(32*n**4 + 302*n**3 + 1037*n**2 + 1530*n + 813)
    m32 = (n+2)**2*(192*n**6 + 2984*n**5 + 19116*n**4 + 64452*n**3 + 120256*n**2 + 117279*n + 46476)
    m33 = (n+2)**2*(-16*n**5 - 408*n**4 - 2912*n**3 - 8884*n**2 - 12254*n - 6240)
    return [[m11, m12, m13], [m21, m22, m23], [m31, m32, m33]]

def mat_mul(A, B):
    n = len(A)
    return [[sum(A[i][k]*B[k][j] for k in range(n)) for j in range(n)] for i in range(n)]

N_MAX = 60
I3 = [[1,0,0],[0,1,0],[0,0,1]]
prod_mat = [row[:] for row in I3]
q_raw = [[], [], []]
for N in range(N_MAX + 5):
    for k in range(3):
        q_raw[k].append(prod_mat[k][0])
    if N < N_MAX + 4:
        prod_mat = mat_mul(prod_mat, M_int(N))

H_vals = [H(n) for n in range(N_MAX + 5)]
q_tw = [[], [], []]
for k in range(3):
    for n in range(N_MAX + 5):
        q_tw[k].append(F(q_raw[k][n]) / H_vals[n])

# Twisted recurrence coefficients (extracted in previous script)
# ℓ₀ through ℓ₃, degree 13 each
# I'll re-extract them here
deg = [13, 13, 13, 13]
n_coeffs_per = [d+1 for d in deg]
total = sum(n_coeffs_per)

A_sys = []
for eq_idx in range(min(total + 10, N_MAX - 3)):
    N = eq_idx
    row = []
    for i in range(4):
        for m in range(deg[i]+1):
            row.append(F(N)**m * q_tw[0][N+i])
    A_sys.append(row)

A_mat = [row[:-1] for row in A_sys]
b_vec = [-row[-1] for row in A_sys]
n_vars = total - 1

aug = [A_mat[i][:] + [b_vec[i]] for i in range(len(A_mat))]
n_r = len(aug); n_c = n_vars

piv_row = 0
for col in range(n_c):
    found = -1
    for r in range(piv_row, n_r):
        if aug[r][col] != 0:
            found = r; break
    if found == -1: continue
    aug[found], aug[piv_row] = aug[piv_row], aug[found]
    pivot_val = aug[piv_row][col]
    for r in range(n_r):
        if r != piv_row and aug[r][col] != 0:
            fac = aug[r][col] / pivot_val
            for c2 in range(n_c + 1):
                aug[r][c2] -= fac * aug[piv_row][c2]
    piv_row += 1

solution = [F(0)] * n_c
for r in range(piv_row):
    pc = -1
    for c in range(n_c):
        if aug[r][c] != 0: pc = c; break
    if pc >= 0:
        solution[pc] = aug[r][n_c] / aug[r][pc]
solution.append(F(1))

poly_coeffs = []
idx = 0
for i in range(4):
    poly_coeffs.append(solution[idx:idx+deg[i]+1])
    idx += deg[i]+1

def eval_poly_F(coeffs, N):
    return sum(c * F(N)**m for m, c in enumerate(coeffs))

# Verify q♯_0 satisfies recurrence
print("=== Verify q♯_0 satisfies twisted recurrence ===")
for N in range(10):
    val = sum(eval_poly_F(poly_coeffs[i], N) * q_tw[0][N+i] for i in range(4))
    print(f"  N={N}: residual = {val}")

# Verify q♯_1 satisfies the SAME recurrence
print("\n=== Verify q♯_1 satisfies twisted recurrence ===")
for N in range(10):
    val = sum(eval_poly_F(poly_coeffs[i], N) * q_tw[1][N+i] for i in range(4))
    print(f"  N={N}: residual = {'0' if val == 0 else str(float(val))}")

# Verify q♯_2
print("\n=== Verify q♯_2 satisfies twisted recurrence ===")
for N in range(10):
    val = sum(eval_poly_F(poly_coeffs[i], N) * q_tw[2][N+i] for i in range(4))
    print(f"  N={N}: residual = {'0' if val == 0 else str(float(val))}")

# Also check: are q♯_0, q♯_1, q♯_2 linearly independent?
print("\n=== Linear independence of q♯_0, q♯_1, q♯_2 ===")
# Wronskian at n=1
W = [[q_tw[j][n] for n in range(1, 4)] for j in range(3)]
det_W = (W[0][0]*(W[1][1]*W[2][2] - W[1][2]*W[2][1])
       - W[0][1]*(W[1][0]*W[2][2] - W[1][2]*W[2][0])
       + W[0][2]*(W[1][0]*W[2][1] - W[1][1]*W[2][0]))
print(f"  Wronskian at n=1,2,3: {float(det_W):.6e}")
print(f"  (nonzero = independent)")

# Now try the null-space approach: find the Ore intertwiner as
# the null space of a BIGGER system where we don't fix any variable
print("\n=== Null-space approach for Ore intertwiner ===")
D = [F(1), F(3)]
E = [F(0), F(1)]
for n in range(1, N_MAX + 8):
    D.append((F(3*(2*n+1)) * D[n] - F(n) * D[n-1]) / F(n+1))
    E.append((F(3*(2*n+1)) * E[n] - F(n) * E[n-1]) / F(n+1))
u_seqs = [[D[n]**2 for n in range(N_MAX+8)],
           [D[n]*E[n] for n in range(N_MAX+8)],
           [E[n]**2 for n in range(N_MAX+8)]]

# Try ORDER-1 Ore intertwiner first (t₂ = 0)
# T = t₀(N) + t₁(N)·S
# T(u_j)(N) = t₀(N)·u_j(N) + t₁(N)·u_j(N+1) = Σ_k a_{jk}·q♯_k(N)
print("\n--- Testing ORDER 1 (T = t₀ + t₁S) ---")
for d in range(8):
    n_t_coeffs = 2 * (d + 1)  # only t₀ and t₁
    n_a_coeffs = 9
    n_total = n_t_coeffs + n_a_coeffs
    K = max(d + 6, (n_total + 2) // 3 + 1)
    n_eqns = 3 * K

    sys_mat = []
    for N in range(1, K+1):
        for j in range(3):
            row = [F(0)] * n_total
            uj = u_seqs[j]
            # t₀ coefficients
            for m in range(d+1):
                row[m] = F(N**m) * uj[N]
            # t₁ coefficients
            for m in range(d+1):
                row[(d+1) + m] = F(N**m) * uj[N+1]
            # a coefficients
            for k in range(3):
                row[n_t_coeffs + j*3 + k] = -q_tw[k][N]
            sys_mat.append(row)

    # Find null space by computing rank
    aug = [row[:] for row in sys_mat]
    n_r = len(aug); n_c = n_total

    piv_row = 0; pivots = []
    for col in range(n_c):
        found = -1
        for r in range(piv_row, n_r):
            if aug[r][col] != 0: found = r; break
        if found == -1: continue
        aug[found], aug[piv_row] = aug[piv_row], aug[found]
        pivot_val = aug[piv_row][col]
        for r in range(n_r):
            if r != piv_row and aug[r][col] != 0:
                fac = aug[r][col] / pivot_val
                for c2 in range(n_c):
                    aug[r][c2] -= fac * aug[piv_row][c2]
        pivots.append((piv_row, col)); piv_row += 1

    rank = len(pivots)
    null_dim = n_c - rank
    if null_dim > 0:
        print(f"  d={d}: rank={rank}/{n_c}, null_dim={null_dim}")

        # Find null vector: for each free variable, set it to 1, others to 0
        pivot_cols = [pc for _, pc in pivots]
        free_cols = [c for c in range(n_c) if c not in pivot_cols]

        for fc in free_cols[:1]:  # just first free variable
            # Set free variable fc to 1, solve for pivot variables
            sol = [F(0)] * n_c
            sol[fc] = F(1)
            for pr, pc in reversed(pivots):
                sol[pc] = -sum(aug[pr][c] * sol[c] for c in range(n_c) if c != pc) / aug[pr][pc]

            # Verify
            verified = True
            for N in range(K+1, K+20):
                if N+1 >= len(u_seqs[0]) or N >= len(q_tw[0]): break
                for j in range(3):
                    uj = u_seqs[j]
                    lhs = sum(sol[m] * F(N**m) * uj[N] for m in range(d+1))
                    lhs += sum(sol[(d+1)+m] * F(N**m) * uj[N+1] for m in range(d+1))
                    rhs = sum(sol[n_t_coeffs + j*3 + k] * q_tw[k][N] for k in range(3))
                    if lhs != rhs:
                        verified = False; break
                if not verified: break

            if verified:
                print(f"    *** ORDER-1 INTERTWINER FOUND at d={d}! ***")
                print(f"    t₀ coeffs: {[str(sol[m]) for m in range(d+1)]}")
                print(f"    t₁ coeffs: {[str(sol[(d+1)+m]) for m in range(d+1)]}")
                print(f"    Connection matrix:")
                for jj in range(3):
                    print(f"      row {jj}: {[str(sol[n_t_coeffs + jj*3 + kk]) for kk in range(3)]}")
            else:
                print(f"    null vector NOT verified on extra points")
    else:
        if d <= 4:
            print(f"  d={d}: rank={rank}/{n_c}, null_dim=0 (no solution)")

# Try ORDER-2 with null space
print("\n--- Testing ORDER 2 (T = t₀ + t₁S + t₂S²), null-space ---")
for d in range(8):
    n_t_coeffs = 3 * (d + 1)
    n_a_coeffs = 9
    n_total = n_t_coeffs + n_a_coeffs
    K = max(d + 6, (n_total + 2) // 3 + 1)
    n_eqns = 3 * K

    sys_mat = []
    for N in range(1, K+1):
        for j in range(3):
            row = [F(0)] * n_total
            uj = u_seqs[j]
            for i in range(3):
                for m in range(d+1):
                    row[i*(d+1) + m] = F(N**m) * uj[N+i]
            for k in range(3):
                row[n_t_coeffs + j*3 + k] = -q_tw[k][N]
            sys_mat.append(row)

    aug = [row[:] for row in sys_mat]
    n_r = len(aug); n_c = n_total
    piv_row = 0; pivots = []
    for col in range(n_c):
        found = -1
        for r in range(piv_row, n_r):
            if aug[r][col] != 0: found = r; break
        if found == -1: continue
        aug[found], aug[piv_row] = aug[piv_row], aug[found]
        pivot_val = aug[piv_row][col]
        for r in range(n_r):
            if r != piv_row and aug[r][col] != 0:
                fac = aug[r][col] / pivot_val
                for c2 in range(n_c):
                    aug[r][c2] -= fac * aug[piv_row][c2]
        pivots.append((piv_row, col)); piv_row += 1

    rank = len(pivots)
    null_dim = n_c - rank
    if null_dim > 0:
        print(f"  d={d}: rank={rank}/{n_c}, null_dim={null_dim}")

        pivot_cols = [pc for _, pc in pivots]
        free_cols = [c for c in range(n_c) if c not in pivot_cols]

        for fc in free_cols[:1]:
            sol = [F(0)] * n_c
            sol[fc] = F(1)
            for pr, pc in reversed(pivots):
                sol[pc] = -sum(aug[pr][c] * sol[c] for c in range(n_c) if c != pc) / aug[pr][pc]

            verified = True
            for N in range(K+1, K+20):
                if N+2 >= len(u_seqs[0]) or N >= len(q_tw[0]): break
                for j in range(3):
                    uj = u_seqs[j]
                    lhs = sum(sum(sol[i*(d+1)+m] * F(N**m) for m in range(d+1)) * uj[N+i] for i in range(3))
                    rhs = sum(sol[n_t_coeffs + j*3 + k] * q_tw[k][N] for k in range(3))
                    if lhs != rhs:
                        verified = False; break
                if not verified: break

            if verified:
                print(f"    *** ORDER-2 INTERTWINER FOUND at d={d}! ***")
                for i in range(3):
                    print(f"    t_{i} coeffs: {[str(sol[i*(d+1)+m]) for m in range(d+1)]}")
            else:
                print(f"    null vector NOT verified")
    else:
        if d <= 4:
            print(f"  d={d}: rank={rank}/{n_c}, null_dim=0")

print("\nDone.")
