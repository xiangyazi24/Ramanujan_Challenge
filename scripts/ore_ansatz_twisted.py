#!/usr/bin/env python3
"""
Problem 2.5: Find the Ore intertwiner between L₂₅♯ (fully twisted) and L_Sym²(Delannoy).

Uses polynomial ansatz for t_i(N) combined with constant connection matrix A.
The twist removes the (N!)^7 factor, making polynomial/rational solutions possible.

From Q4792: deg_∞ = 1. So:
  - polynomial ansatz: t_i(N) = linear polynomial (degree 1)
  - with P₆'(n) denominator: numerator degree 7
  - with det-normal-form denominator (deg 9): numerator degree 10

We try polynomial first (degrees 0 through 7), then rational with P₆' denominator.
"""
from fractions import Fraction as F
from functools import reduce
from math import gcd

# Pochhammer
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

# Compute Delannoy sequences
N_MAX = 80
D = [F(1), F(3)]
E = [F(0), F(1)]
for n in range(1, N_MAX + 8):
    D.append((F(3*(2*n+1)) * D[n] - F(n) * D[n-1]) / F(n+1))
    E.append((F(3*(2*n+1)) * E[n] - F(n) * E[n-1]) / F(n+1))

u1 = [D[n]**2 for n in range(N_MAX+8)]   # DD
u2 = [D[n]*E[n] for n in range(N_MAX+8)] # DE
u3 = [E[n]**2 for n in range(N_MAX+8)]   # EE

# CMF matrix
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

# Compute 3 independent twisted CMF solutions
print("Computing twisted CMF solutions...")
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

print(f"q♯_0[0..4] = {[float(q_tw[0][i]) for i in range(5)]}")

# ---- Polynomial ansatz for Ore intertwiner ----
# For polynomial degree d: t_i(N) = Σ_{m=0}^d c_{i,m} · N^m
# Connection: T(u_j)(N) = Σ_k a_{jk} · q♯_k(N)
# i.e. Σ_i Σ_m c_{im} N^m u_j(N+i) = Σ_k a_{jk} q♯_k(N)
# Unknowns: c_{im} (3(d+1)) + a_{jk} (9) = 3d + 12
# Each N gives 3 equations, so need K ≥ d + 4 values of N
# System is homogeneous up to one fixable variable.

u_seqs = [u1, u2, u3]

for d in range(12):
    n_t_coeffs = 3 * (d + 1)
    n_a_coeffs = 9
    n_total = n_t_coeffs + n_a_coeffs

    K = max(d + 6, (n_total + 2) // 3 + 1)
    n_eqns = 3 * K

    # Build system
    sys_mat = []
    for N in range(1, K+1):
        for j in range(3):
            row = [F(0)] * n_total
            uj = u_seqs[j]
            for i in range(3):
                for m in range(d+1):
                    col = i * (d+1) + m
                    row[col] = F(N**m) * uj[N+i]
            for k in range(3):
                col = n_t_coeffs + j * 3 + k
                row[col] = -q_tw[k][N]
            sys_mat.append(row)

    # Fix last variable = 1
    A_mat = [row[:-1] for row in sys_mat]
    b_vec = [-row[-1] for row in sys_mat]
    n_unknowns = n_total - 1

    # Gaussian elimination
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

    # Check consistency
    consistent = True
    for r in range(piv_row, n_r):
        if aug[r][n_c] != 0:
            consistent = False
            break

    if not consistent:
        if d <= 5:
            print(f"  d={d}: INCONSISTENT")
        continue

    rank = len(pivots)
    if rank < n_c:
        if d <= 5:
            print(f"  d={d}: UNDERDETERMINED (rank {rank}/{n_c})")
        continue

    # Extract solution
    sol = [F(0)] * n_c
    for pr, pc in reversed(pivots):
        sol[pc] = aug[pr][n_c] / aug[pr][pc]
    sol.append(F(1))

    # Verify on additional data points
    verified = True
    for N in range(K+1, K+30):
        if N+2 >= len(u1) or N >= len(q_tw[0]):
            break
        for j in range(3):
            uj = u_seqs[j]
            lhs = sum(sum(sol[i*(d+1)+m] * F(N**m) for m in range(d+1)) * uj[N+i] for i in range(3))
            rhs = sum(sol[n_t_coeffs + j*3 + k] * q_tw[k][N] for k in range(3))
            if lhs != rhs:
                verified = False
                break
        if not verified:
            break

    if verified:
        print(f"\n*** FOUND INTERTWINER at polynomial degree d={d} ***")
        for i in range(3):
            coeffs = sol[i*(d+1):(i+1)*(d+1)]
            print(f"  t_{i}(N) coefficients: {[str(c) for c in coeffs]}")
        print("  Connection matrix A:")
        for j in range(3):
            row_a = sol[n_t_coeffs + j*3:n_t_coeffs + (j+1)*3]
            print(f"    row {j}: {[str(c) for c in row_a]}")
        break
    else:
        if d <= 5:
            print(f"  d={d}: NOT VERIFIED on extra points")
        elif d == 11:
            print(f"  d={d}: still not verified")

else:
    print("\nNo polynomial Ore intertwiner found for d ≤ 11.")
    print("Trying RATIONAL intertwiner with P₆'(N) denominator...")

    # P₆'(N) = 3072N⁶ + 55680N⁵ + 414064N⁴ + 1615610N³ + 3483853N² + 3929280N + 1806156
    def P6prime(N):
        return (3072*N**6 + 55680*N**5 + 414064*N**4 + 1615610*N**3
                + 3483853*N**2 + 3929280*N + 1806156)

    # t_i(N) = p_i(N) / P₆'(N) where p_i has degree d_num
    # T(u_j)(N) = Σ_i [p_i(N)/P₆'(N)] u_j(N+i) = Σ_k a_{jk} q♯_k(N)
    # Multiply through by P₆'(N):
    # Σ_i p_i(N) u_j(N+i) = P₆'(N) · Σ_k a_{jk} q♯_k(N)

    for d_num in range(7, 12):
        n_p_coeffs = 3 * (d_num + 1)
        n_a_coeffs = 9
        n_total = n_p_coeffs + n_a_coeffs

        K = max(d_num + 6, (n_total + 2) // 3 + 1)
        n_eqns = 3 * K

        sys_mat = []
        for N in range(1, K+1):
            P6p_val = F(P6prime(N))
            for j in range(3):
                row = [F(0)] * n_total
                uj = u_seqs[j]
                for i in range(3):
                    for m in range(d_num+1):
                        col = i * (d_num+1) + m
                        row[col] = F(N**m) * uj[N+i]
                for k in range(3):
                    col = n_p_coeffs + j * 3 + k
                    row[col] = -P6p_val * q_tw[k][N]
                sys_mat.append(row)

        A_mat = [row[:-1] for row in sys_mat]
        b_vec = [-row[-1] for row in sys_mat]
        n_unknowns = n_total - 1

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

        consistent = True
        for r in range(piv_row, n_r):
            if aug[r][n_c] != 0:
                consistent = False
                break

        if not consistent:
            print(f"  d_num={d_num}: INCONSISTENT (P₆' denominator)")
            continue

        rank = len(pivots)
        if rank < n_c:
            print(f"  d_num={d_num}: UNDERDETERMINED (rank {rank}/{n_c})")
            continue

        sol = [F(0)] * n_c
        for pr, pc in reversed(pivots):
            sol[pc] = aug[pr][n_c] / aug[pr][pc]
        sol.append(F(1))

        verified = True
        for N in range(K+1, K+20):
            if N+2 >= len(u1) or N >= len(q_tw[0]):
                break
            P6p_val = F(P6prime(N))
            for j in range(3):
                uj = u_seqs[j]
                lhs = sum(sum(sol[i*(d_num+1)+m] * F(N**m) for m in range(d_num+1)) * uj[N+i] for i in range(3))
                rhs = P6p_val * sum(sol[n_p_coeffs + j*3 + k] * q_tw[k][N] for k in range(3))
                if lhs != rhs:
                    verified = False
                    break
            if not verified:
                break

        if verified:
            print(f"\n*** FOUND RATIONAL INTERTWINER: numerator deg {d_num} / P₆'(N) ***")
            for i in range(3):
                coeffs = sol[i*(d_num+1):(i+1)*(d_num+1)]
                print(f"  p_{i}(N) coefficients: {[str(c) for c in coeffs[:5]]}...")
            print("  Connection matrix A:")
            for j in range(3):
                row_a = sol[n_p_coeffs + j*3:n_p_coeffs + (j+1)*3]
                print(f"    row {j}: {[float(c) for c in row_a]}")

            # Factor the numerators with sympy
            from sympy import Symbol, factor as sfact, Rational as SR
            N_sym = Symbol('N')
            for i in range(3):
                coeffs = sol[i*(d_num+1):(i+1)*(d_num+1)]
                denoms = [c.denominator for c in coeffs if c != 0]
                if not denoms:
                    continue
                L = reduce(lambda a,b: a*b//gcd(a,b), denoms)
                p_int = [int(c * L) for c in coeffs]
                g_val = reduce(gcd, [abs(x) for x in p_int if x != 0]) if any(x != 0 for x in p_int) else 1
                p_int = [x // g_val for x in p_int]
                p_sym = sum(c * N_sym**m for m, c in enumerate(p_int))
                print(f"  Factored p_{i}(N): {sfact(p_sym)}")
            break
        else:
            print(f"  d_num={d_num}: NOT VERIFIED (P₆' denominator)")
    else:
        print("  No rational intertwiner with P₆' denominator found")

print("\nDone.")
