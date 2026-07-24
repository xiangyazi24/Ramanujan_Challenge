#!/usr/bin/env python3
"""
Problem 2.5: Find Ore intertwiner between GAUGED L_25 and Sym^2(Delannoy).

Key insight: L_25 has degree pattern (28,21,14,7). After extracting the
factorial gauge g_N = (-16)^N * (N!)^7, the gauged operator L'_25 has
constant degree 28 and Poincaré polynomial t^3 - 35t^2 + 35t - 1, matching
Sym^2(Delannoy) EXACTLY. Now the Ore intertwiner has Q(n) coefficients.
"""
from fractions import Fraction as F
from math import factorial

# ---- CMF matrix M(n) with exact integer entries ----
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

# ---- Compute CMF scalar sequences ----
print("Computing CMF scalar sequences...")
N_MAX = 90
I3 = [[1,0,0],[0,1,0],[0,0,1]]
prod = [row[:] for row in I3]
q_raw = [[], [], []]

for N in range(N_MAX + 5):
    for k in range(3):
        q_raw[k].append(prod[k][0])
    if N < N_MAX + 4:
        prod = mat_mul(prod, M_int(N))

print(f"Computed q_raw[0..{len(q_raw[0])-1}]")

# ---- Apply factorial gauge: t_N = Q_N / ((-16)^N * (N!)^7) ----
print("\nApplying factorial gauge g_N = (-16)^N * (N!)^7...")

t = [[], [], []]
for k in range(3):
    for N in range(len(q_raw[k])):
        gauge = ((-16)**N) * (factorial(N))**7
        t[k].append(F(q_raw[k][N], gauge))

print(f"Gauged t[0..4]:")
for k in range(3):
    print(f"  t[{k}]: {[str(t[k][N]) for N in range(5)]}")

# Check: gauged Poincaré should be (t-1)(t^2-34t+1)
# t_N ~ c_j * N^sigma, ratio t_{N+1}/t_N → lambda_j ∈ {1, 17±12√2}
print("\nRatio t[0][N+1]/t[0][N] for large N (should → one of {1, 33.97, 0.029}):")
for N in [20, 30, 40, 50]:
    if N+1 < len(t[0]) and t[0][N] != 0:
        r = float(t[0][N+1] / t[0][N])
        print(f"  N={N}: ratio = {r:.6f}")

# ---- Delannoy sequences ----
N_DEL = N_MAX + 10
D = [F(0)] * N_DEL
D[0] = F(1); D[1] = F(3)
for n in range(1, N_DEL-1):
    D[n+1] = (F(3*(2*n+1)) * D[n] - F(n) * D[n-1]) / F(n+1)

E = [F(0)] * N_DEL
E[0] = F(0); E[1] = F(1)
for n in range(1, N_DEL-1):
    E[n+1] = (F(3*(2*n+1)) * E[n] - F(n) * E[n-1]) / F(n+1)

u1 = [D[n]**2 for n in range(N_DEL)]
u2 = [D[n]*E[n] for n in range(N_DEL)]
u3 = [E[n]**2 for n in range(N_DEL)]

# ---- Find intertwiner T' (degree d polynomial coefficients) ----
# t_0(N)*u_j(N) + t_1(N)*u_j(N+1) + t_2(N)*u_j(N+2) = sum_k a_{jk} * t[k][N]

print("\n=== Searching for Ore intertwiner T' (gauged) ===")

for d in range(25):
    n_t_coeffs = 3 * (d + 1)
    n_a_coeffs = 9
    n_total = n_t_coeffs + n_a_coeffs

    K = max(d + 6, (n_total + 2) // 3 + 2)
    n_eqns = 3 * K

    if K + 3 >= len(u1) or K >= len(t[0]):
        print(f"  d={d}: not enough data (need K={K})")
        break

    sys_mat = []
    for N in range(1, K+1):
        for j in range(3):
            row = [F(0)] * n_total
            uj = [u1, u2, u3][j]
            for i in range(3):
                for m in range(d+1):
                    col = i * (d+1) + m
                    row[col] = F(N**m) * uj[N+i]
            for k in range(3):
                col = n_t_coeffs + j * 3 + k
                row[col] = -t[k][N]
            sys_mat.append(row)

    # Fix a_{22} = 1
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
                factor = F(aug[r][col], aug[piv_row][col]) if isinstance(aug[piv_row][col], int) else aug[r][col] / pivot_val
                for c2 in range(n_c + 1):
                    aug[r][c2] -= factor * aug[piv_row][c2]
        pivots.append((piv_row, col))
        piv_row += 1

    consistent = all(aug[r][n_c] == 0 for r in range(piv_row, n_r))
    rank = len(pivots)

    if not consistent:
        if d <= 5 or d % 5 == 0:
            print(f"  d={d}: INCONSISTENT")
        continue

    if rank < n_c:
        if d <= 5 or d % 5 == 0:
            print(f"  d={d}: UNDERDETERMINED (rank {rank} < {n_c})")
        continue

    # Extract solution
    sol = [F(0)] * n_c
    for pr, pc in reversed(pivots):
        sol[pc] = aug[pr][n_c] / aug[pr][pc]
    sol.append(F(1))

    # Verify on additional points
    verified = True
    n_verify = 0
    for N in range(K+1, min(K+30, len(u1)-2, len(t[0]))):
        for j in range(3):
            uj = [u1, u2, u3][j]
            lhs = sum(sum(sol[i*(d+1)+m] * F(N**m) for m in range(d+1)) * uj[N+i] for i in range(3))
            rhs = sum(sol[n_t_coeffs + j*3 + k] * t[k][N] for k in range(3))
            if lhs != rhs:
                verified = False
                break
        if not verified:
            break
        n_verify += 1

    if verified and n_verify >= 5:
        print(f"\n*** INTERTWINER FOUND at degree d={d} ***")
        print(f"    (verified on {n_verify} additional points)")

        for i in range(3):
            coeffs = sol[i*(d+1):(i+1)*(d+1)]
            # Try to identify common denominator
            denoms = [c.denominator for c in coeffs if c != 0]
            from math import gcd
            from functools import reduce
            if denoms:
                lcm_d = reduce(lambda a, b: a * b // gcd(a, b), denoms)
                int_coeffs = [int(c * lcm_d) for c in coeffs]
                g = reduce(gcd, [abs(x) for x in int_coeffs if x != 0], 0)
                if g > 0:
                    int_coeffs = [x // g for x in int_coeffs]
                print(f"  t_{i}(N) = (1/{lcm_d // g}) * ({' + '.join(f'{c}*N^{m}' if m > 0 else str(c) for m, c in enumerate(int_coeffs) if c != 0)})")
            else:
                print(f"  t_{i}(N) = 0")

        print("\n  Connection matrix A (maps Sym² basis to gauged L_25 basis):")
        for j in range(3):
            row_a = sol[n_t_coeffs + j*3:n_t_coeffs + (j+1)*3]
            print(f"    [{', '.join(str(c) for c in row_a)}]")

        break
    else:
        if d <= 5 or d % 5 == 0:
            msg = f"  d={d}: "
            if not verified:
                msg += "fails verification"
            elif n_verify < 5:
                msg += f"only {n_verify} verification points"
            print(msg)

print("\nDone.")
