#!/usr/bin/env python3
"""Problem 2.5: Search for the TWISTED rational coboundary matrix.

From Q4811: The untwisted coboundary M(n)·U(n+1) = U(n)·C(n) does NOT exist
(determinant obstruction: det(M(n)) ~ n^21, impossible for rational U).

The correct equation with rank-one twist:
M(n) · U(n+1) = s(n) · U(n) · C(n)
where s(n) = (-16)(n+1)^7.

Strategy: Compute M(n), s(n), C(n) numerically at many integer points n,
then solve for U(n) as polynomial matrix of degree d.
"""
from fractions import Fraction as F

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
    return [[F(m11), F(m12), F(m13)],
            [F(m21), F(m22), F(m23)],
            [F(m31), F(m32), F(m33)]]

def C_mat(n):
    """Companion matrix for Sym²(Delannoy) recurrence:
    (2n+3)(n+3)² u_{n+3} = (2n+5)(35n²+140n+131) u_{n+2}
                          - (2n+3)(35n²+140n+131) u_{n+1}
                          + (2n+5)(n+1)² u_n
    """
    a3 = F((2*n+3)*(n+3)**2)
    a2 = F((2*n+5)*(35*n**2 + 140*n + 131))
    a1 = F((2*n+3)*(35*n**2 + 140*n + 131))
    a0 = F((2*n+5)*(n+1)**2)
    return [[F(0), F(1), F(0)],
            [F(0), F(0), F(1)],
            [a0/a3, -a1/a3, a2/a3]]

def s(n):
    """Scalar twist: s(n) = -16 * (n+1)^7."""
    return F(-16 * (n+1)**7)

def mat_mul(A, B):
    return [[sum(A[i][k]*B[k][j] for k in range(3)) for j in range(3)] for i in range(3)]

def mat_scale(c, A):
    return [[c * A[i][j] for j in range(3)] for i in range(3)]

# The equation: M(n) · U(n+1) = s(n) · U(n) · C(n)
# So: U(n) = s(n)^{-1} · M(n) · U(n+1) · C(n)^{-1}
# Or: U(n+1) = M(n)^{-1} · s(n) · U(n) · C(n)

# If U(n) has polynomial entries of degree d, then U(n) has 9*(d+1) unknowns.
# For each n, the equation gives 9 polynomial equations in n.

# Strategy: evaluate at many integer n values and solve the linear system.
# U(n) = sum_{m=0}^d u_m * n^m where u_m are 3x3 constant matrices.
# Total unknowns: 9*(d+1).
# Each integer n gives 9 equations.
# Need at least (d+1) integer points.

print("=== Twisted coboundary search ===")
print("Equation: M(n) · U(n+1) = s(n) · U(n) · C(n)")
print(f"s(n) = -16 * (n+1)^7\n")

# First, check the determinant ratio to verify it telescopes.
print("--- Determinant ratio check ---")
for n in range(8):
    Mn = M_int(n)
    Cn = C_mat(n)
    sn = s(n)

    det_M = (Mn[0][0]*(Mn[1][1]*Mn[2][2]-Mn[1][2]*Mn[2][1])
            -Mn[0][1]*(Mn[1][0]*Mn[2][2]-Mn[1][2]*Mn[2][0])
            +Mn[0][2]*(Mn[1][0]*Mn[2][1]-Mn[1][1]*Mn[2][0]))

    det_C = Cn[2][0]  # det of companion = c0

    # det(M(n)) / (s(n)^3 * det(C(n))) should be det(U(n+1))/det(U(n))
    ratio = det_M / (sn**3 * det_C)
    print(f"  n={n}: det(M)/(s^3·det(C)) = {ratio}")

# Now search for U(n) = sum u_m n^m with polynomial entries of degree d.
print("\n--- Searching for polynomial U(n) ---")

for d in range(15):
    n_unknowns = 9 * (d + 1)
    # We need enough equations. Each n gives 9 equations.
    n_points_needed = (n_unknowns + 8) // 9 + 5  # some slack

    if n_points_needed > 50:
        print(f"  d={d}: need {n_points_needed} points, skipping")
        continue

    # Build system: for each n, and each entry (i,j):
    # sum_m [M(n) · u_m · (n+1)^m]_{ij} = s(n) · sum_m [u_m · n^m · C(n)]_{ij}
    # where u_m is a 3x3 matrix.

    # Rewrite: for entry (i,j) of the matrix equation:
    # sum_{k=0}^2 M_{ik}(n) * U_{kj}(n+1) = s(n) * sum_{k=0}^2 U_{ik}(n) * C_{kj}(n)
    # where U_{ab}(n) = sum_{m=0}^d u_{ab,m} * n^m

    # LHS: sum_k M_{ik}(n) * sum_m u_{km,j} * (n+1)^m
    # RHS: s(n) * sum_k (sum_m u_{ik,m} * n^m) * C_{kj}(n)

    # Each entry (i,j) at each n gives one equation.
    # Unknowns: u_{a,b,m} for a,b in {0,1,2}, m in {0,...,d}.
    # Index: a*3*(d+1) + b*(d+1) + m

    rows = []
    for n in range(n_points_needed):
        Mn = M_int(n)
        Cn = C_mat(n)
        sn = s(n)

        for i in range(3):
            for j in range(3):
                row = [F(0)] * n_unknowns

                # LHS: sum_k M_{ik}(n) * U_{kj}(n+1) = sum_k M_ik * sum_m u_{kj,m} (n+1)^m
                for k in range(3):
                    for m in range(d+1):
                        col = k * 3 * (d+1) + j * (d+1) + m
                        row[col] += Mn[i][k] * F((n+1)**m)

                # RHS: s(n) * sum_k U_{ik}(n) * C_{kj}(n) = s(n) * sum_k sum_m u_{ik,m} n^m * C_{kj}(n)
                for k in range(3):
                    for m in range(d+1):
                        col = i * 3 * (d+1) + k * (d+1) + m
                        row[col] -= sn * F(n**m) * Cn[k][j]

                rows.append(row)

    # Gaussian elimination
    mat = [row[:] for row in rows]
    n_r = len(mat)
    n_c = n_unknowns

    pivots = []
    piv_row = 0
    for col in range(n_c):
        found = -1
        for rr in range(piv_row, n_r):
            if mat[rr][col] != 0:
                found = rr
                break
        if found == -1:
            continue
        mat[found], mat[piv_row] = mat[piv_row], mat[found]
        pv = mat[piv_row][col]
        for rr in range(n_r):
            if rr != piv_row and mat[rr][col] != 0:
                f = mat[rr][col] / pv
                for cc in range(n_c):
                    mat[rr][cc] -= f * mat[piv_row][cc]
        pivots.append((piv_row, col))
        piv_row += 1

    rank = len(pivots)
    null_dim = n_c - rank

    # Check if the remaining rows are all zero
    consistent = all(all(mat[rr][cc] == 0 for cc in range(n_c)) for rr in range(piv_row, n_r))

    if not consistent:
        if d <= 5 or d % 5 == 0:
            print(f"  d={d}: INCONSISTENT (rank={rank})")
        continue

    if null_dim == 0:
        print(f"  d={d}: unique trivial solution only")
        continue

    print(f"  d={d}: rank={rank}, null_dim={null_dim} — NONTRIVIAL SOLUTION EXISTS")

    # Extract a null vector
    free_cols = [c for c in range(n_c) if c not in [p[1] for p in pivots]]
    sol = [F(0)] * n_c
    sol[free_cols[0]] = F(1)
    for pr, pc in reversed(pivots):
        val = F(0)
        for cc in range(pc + 1, n_c):
            val += mat[pr][cc] * sol[cc]
        sol[pc] = -val / mat[pr][pc]

    # Display the solution
    print(f"\n  Solution U(n) entries:")
    for a in range(3):
        for b in range(3):
            coeffs = [sol[a*3*(d+1) + b*(d+1) + m] for m in range(d+1)]
            nonzero = [m for m in range(d+1) if coeffs[m] != 0]
            if nonzero:
                print(f"    U[{a}][{b}](n) = ", end="")
                terms = []
                for m in nonzero:
                    c = coeffs[m]
                    if m == 0:
                        terms.append(str(c))
                    elif m == 1:
                        terms.append(f"({c})*n")
                    else:
                        terms.append(f"({c})*n^{m}")
                print(" + ".join(terms))
            else:
                print(f"    U[{a}][{b}](n) = 0")

    # Verify on additional points
    n_verify = 0
    for n in range(n_points_needed, n_points_needed + 20):
        Mn = M_int(n)
        Cn = C_mat(n)
        sn = s(n)

        # Evaluate U(n) and U(n+1)
        Un = [[sum(sol[a*3*(d+1)+b*(d+1)+m] * F(n**m) for m in range(d+1))
               for b in range(3)] for a in range(3)]
        Un1 = [[sum(sol[a*3*(d+1)+b*(d+1)+m] * F((n+1)**m) for m in range(d+1))
                for b in range(3)] for a in range(3)]

        LHS = mat_mul([[Mn[i][j] for j in range(3)] for i in range(3)], Un1)
        RHS = mat_mul(Un, [[Cn[i][j] for j in range(3)] for i in range(3)])
        RHS = mat_scale(sn, RHS)

        match = all(LHS[i][j] == RHS[i][j] for i in range(3) for j in range(3))
        if not match:
            print(f"\n  VERIFICATION FAILED at n={n}")
            break
        n_verify += 1

    print(f"\n  Verified on {n_verify} additional points.")

    # Check det(U(n)) at a few points
    print(f"\n  det(U(n)) at n=0,...,5:")
    for n in range(6):
        Un = [[sum(sol[a*3*(d+1)+b*(d+1)+m] * F(n**m) for m in range(d+1))
               for b in range(3)] for a in range(3)]
        det_U = (Un[0][0]*(Un[1][1]*Un[2][2]-Un[1][2]*Un[2][1])
                -Un[0][1]*(Un[1][0]*Un[2][2]-Un[1][2]*Un[2][0])
                +Un[0][2]*(Un[1][0]*Un[2][1]-Un[1][1]*Un[2][0]))
        print(f"    det(U({n})) = {det_U}")

    break

print("\nDone.")
