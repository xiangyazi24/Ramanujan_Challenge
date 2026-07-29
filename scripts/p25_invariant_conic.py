#!/usr/bin/env python3
"""P2.5: Search for the invariant conic J(n).

From Q4880: if the CMF is gauge-equivalent to Sym²(Delannoy), then
there exists a nonsingular symmetric J(n) satisfying:

  M̄(n) · J(n+1) · M̄(n)^T = α_n² · J(n)

where M̄(n) = D(n)^{-1} M_H(n) D(n+1), D(n) = diag(1, n+1, (n+1)²),
α_n = (n+1)/(n+2), and J(n) has 6 independent entries.

This is a LINEAR system. We parametrize J(n) with polynomial entries
of degree d and solve.
"""
from fractions import Fraction

def M_entries(n):
    n = Fraction(n)
    m11 = (-2*n-5)*(n+3)**2*(136*n**4+1424*n**3+5548*n**2+9551*n+6141)
    m12 = 384*n**6+6384*n**5+44168*n**4+162698*n**3+336377*n**2+369933*n+169011
    m13 = -(480*n**4+4980*n**3+19210*n**2+32690*n+20730)
    m21 = (n+2)**2*(n+3)**2*(4*n+10)*(48*n**3+386*n**2+1017*n+879)
    m22 = (n+2)**2*(-272*n**5-3848*n**4-21732*n**3-61184*n**2-85761*n-47808)
    m23 = (n+2)**2*(320*n**3+2540*n**2+6610*n+5640)
    m31 = (-4*n-10)*(n+2)**2*(n+3)**2*(32*n**4+302*n**3+1037*n**2+1530*n+813)
    m32 = (n+2)**2*(192*n**6+2984*n**5+19116*n**4+64452*n**3+120256*n**2+117279*n+46476)
    m33 = (n+2)**2*(-16*n**5-408*n**4-2912*n**3-8884*n**2-12254*n-6240)
    return [[m11,m12,m13],[m21,m22,m23],[m31,m32,m33]]

def delta_H(n):
    n = Fraction(n)
    return Fraction(-2)*(n+2)**2*(n+3)**2*(2*n+5)*(2*n+7)**2

def M_H(n):
    """Normalized matrix M_H(n) = M(n)/delta_H(n)"""
    M = M_entries(n)
    d = delta_H(n)
    return [[M[i][j] / d for j in range(3)] for i in range(3)]

def balanced_M(n):
    """M̄(n) = D(n)^{-1} · M_H(n) · D(n+1), where D(n) = diag(1, n+1, (n+1)²)"""
    MH = M_H(n)
    nn = Fraction(n)
    # D(n)^{-1} = diag(1, 1/(n+1), 1/(n+1)²)
    # D(n+1) = diag(1, n+2, (n+2)²)
    d_inv = [Fraction(1), Fraction(1, int(nn+1)), Fraction(1, int((nn+1)**2))]
    d_next = [Fraction(1), nn+2, (nn+2)**2]
    result = [[Fraction(0)]*3 for _ in range(3)]
    for i in range(3):
        for j in range(3):
            result[i][j] = d_inv[i] * MH[i][j] * d_next[j]
    return result

def mat_mul_3(A, B):
    C = [[Fraction(0)]*3 for _ in range(3)]
    for i in range(3):
        for j in range(3):
            for k in range(3):
                C[i][j] += A[i][k] * B[k][j]
    return C

def mat_transpose(A):
    return [[A[j][i] for j in range(3)] for i in range(3)]

# The 6 independent entries of a symmetric 3×3 matrix:
# J = [[j00, j01, j02],
#       [j01, j11, j12],
#       [j02, j12, j22]]
# We index them as: 0→(0,0), 1→(0,1), 2→(0,2), 3→(1,1), 4→(1,2), 5→(2,2)
SYM_INDICES = [(0,0), (0,1), (0,2), (1,1), (1,2), (2,2)]

def sym_to_matrix(entries):
    """Convert 6-vector to symmetric 3×3 matrix."""
    J = [[Fraction(0)]*3 for _ in range(3)]
    for idx, (i, j) in enumerate(SYM_INDICES):
        J[i][j] = entries[idx]
        J[j][i] = entries[idx]
    return J

def matrix_to_sym(J):
    """Extract 6-vector from symmetric 3×3 matrix."""
    return [J[i][j] for (i, j) in SYM_INDICES]

print("=== P2.5 Invariant Conic Search ===", flush=True)
print("Testing M̄(n)·J(n+1)·M̄(n)^T = α_n²·J(n)", flush=True)

# First, verify limiting matrices
print("\n--- Limiting behavior ---", flush=True)
Mbar_large = balanced_M(1000)
print(f"M̄(1000) ≈")
for row in Mbar_large:
    print(f"  [{', '.join(f'{float(x):.6f}' for x in row)}]")

# Expected A_∞ = [[17,-24,0],[-12,17,0],[8,-12,1]]
print(f"\nExpected A_∞ = [[17,-24,0],[-12,17,0],[8,-12,1]]")

# Try polynomial J(n) of degree d
# J_{ij}(n) = Σ_{k=0}^d c_{ij,k} · n^k
# Unknowns: 6*(d+1) coefficients
# For each test n, M̄(n)·J(n+1)·M̄(n)^T - α_n²·J(n) = 0 gives 6 equations
# (from the 6 independent entries of the resulting symmetric matrix)

for deg in range(8):
    n_unknowns = 6 * (deg + 1)
    n_train = n_unknowns + 8

    # Build the system
    rows = []
    rhs = []

    for nn in range(n_train):
        Mbar = balanced_M(nn)
        MbarT = mat_transpose(Mbar)
        alpha_sq = Fraction(nn+1, nn+2)**2

        for sym_idx, (si, sj) in enumerate(SYM_INDICES):
            row = [Fraction(0)] * n_unknowns

            for unk_idx in range(6):
                ui, uj = SYM_INDICES[unk_idx]

                for power in range(deg + 1):
                    col = unk_idx * (deg + 1) + power

                    # Contribution from J(n+1) at position (ui, uj):
                    # [M̄ · e_{ui,uj} · M̄^T]_{si,sj}
                    # where e_{ui,uj} is the symmetric matrix with 1 at (ui,uj) and (uj,ui)

                    # M̄ · e · M̄^T: entry (si,sj) =
                    # Σ_a Σ_b M̄[si][a] · e[a][b] · M̄^T[b][sj]
                    # = Σ_a Σ_b M̄[si][a] · e[a][b] · M̄[sj][b]

                    # e[a][b] = δ_{a,ui}δ_{b,uj} + δ_{a,uj}δ_{b,ui} (if ui≠uj)
                    #         = δ_{a,ui}δ_{b,ui} (if ui==uj)

                    contrib = Fraction(0)
                    if ui == uj:
                        contrib = Mbar[si][ui] * Mbar[sj][ui]
                    else:
                        contrib = (Mbar[si][ui] * Mbar[sj][uj] +
                                   Mbar[si][uj] * Mbar[sj][ui])

                    # This is the coefficient of J_{ui,uj}(n+1) in the (si,sj) equation
                    # J_{ui,uj}(n+1) has monomial (n+1)^power at this coefficient slot
                    # Minus α_n² · J_{ui,uj}(n) which has monomial n^power

                    n_frac = Fraction(nn)
                    jp_coeff = contrib * Fraction(nn + 1)**power  # from J(n+1)
                    jn_coeff = -alpha_sq * Fraction(nn)**power     # from -α²·J(n)

                    if si == sj:
                        # Diagonal: only one equation
                        row[col] += jp_coeff + jn_coeff
                    else:
                        # Off-diagonal: the (si,sj) entry equals (sj,si) entry
                        row[col] += jp_coeff + jn_coeff

            rows.append(row)
            rhs.append(Fraction(0))  # homogeneous system

    # Solve via Gaussian elimination (find null space dimension)
    m_eq = len(rows)
    A = [list(row) for row in rows]

    pivot_cols = []
    row_idx = 0
    for col in range(n_unknowns):
        found = -1
        for r in range(row_idx, m_eq):
            if A[r][col] != 0:
                found = r
                break
        if found == -1:
            continue
        A[row_idx], A[found] = A[found], A[row_idx]
        piv = A[row_idx][col]
        for j in range(n_unknowns):
            A[row_idx][j] /= piv
        for r in range(m_eq):
            if r == row_idx:
                continue
            if A[r][col] == 0:
                continue
            factor = A[r][col]
            for j in range(n_unknowns):
                A[r][j] -= factor * A[row_idx][j]
        pivot_cols.append(col)
        row_idx += 1

    rank = len(pivot_cols)
    nullity = n_unknowns - rank

    print(f"\ndeg={deg}: unknowns={n_unknowns}, equations={m_eq}, rank={rank}, nullity={nullity}", flush=True)

    if nullity > 0:
        # Find a basis vector for the null space
        free_cols = [c for c in range(n_unknowns) if c not in pivot_cols]

        # Set first free variable to 1, rest to 0
        x = [Fraction(0)] * n_unknowns
        x[free_cols[0]] = Fraction(1)

        # Back-substitute
        for pi in range(rank - 1, -1, -1):
            pc = pivot_cols[pi]
            val = Fraction(0)
            for j in range(n_unknowns):
                if j != pc:
                    val += A[pi][j] * x[j]
            x[pc] = -val

        # Display the solution
        print(f"  Solution (first null vector):")
        for unk_idx in range(6):
            ui, uj = SYM_INDICES[unk_idx]
            coeffs = x[unk_idx*(deg+1):(unk_idx+1)*(deg+1)]
            if any(c != 0 for c in coeffs):
                poly_str = " + ".join(f"({c})*n^{k}" for k, c in enumerate(coeffs) if c != 0)
                print(f"  J[{ui},{uj}](n) = {poly_str}")

        # Verify on holdout
        print(f"\n  Verifying on holdout n={n_train}..{n_train+4}:", flush=True)
        all_ok = True
        for nn in range(n_train, n_train + 5):
            # Reconstruct J(n) and J(n+1) from solution
            Jn_entries = [sum(x[idx*(deg+1)+k] * Fraction(nn)**k for k in range(deg+1)) for idx in range(6)]
            Jn1_entries = [sum(x[idx*(deg+1)+k] * Fraction(nn+1)**k for k in range(deg+1)) for idx in range(6)]

            Jn = sym_to_matrix(Jn_entries)
            Jn1 = sym_to_matrix(Jn1_entries)

            Mbar = balanced_M(nn)
            MbarT = mat_transpose(Mbar)
            alpha_sq = Fraction(nn+1, nn+2)**2

            # LHS = M̄ · J(n+1) · M̄^T
            LHS = mat_mul_3(mat_mul_3(Mbar, Jn1), MbarT)
            # RHS = α² · J(n)
            RHS = [[alpha_sq * Jn[i][j] for j in range(3)] for i in range(3)]

            max_err = max(abs(LHS[i][j] - RHS[i][j]) for i in range(3) for j in range(3))
            if max_err != 0:
                all_ok = False
                print(f"  n={nn}: FAIL, max |LHS-RHS| = {float(max_err):.3e}")
            else:
                print(f"  n={nn}: EXACT ✓")

        if all_ok:
            print(f"\n  *** INVARIANT CONIC FOUND at degree {deg} ***")

            # Check if J(n) is nonsingular at n=0
            Jn_entries = [sum(x[idx*(deg+1)+k] * Fraction(0)**k for k in range(deg+1)) for idx in range(6)]
            J0 = sym_to_matrix(Jn_entries)
            det_J0 = (J0[0][0]*(J0[1][1]*J0[2][2]-J0[1][2]*J0[2][1])
                      -J0[0][1]*(J0[1][0]*J0[2][2]-J0[1][2]*J0[2][0])
                      +J0[0][2]*(J0[1][0]*J0[2][1]-J0[1][1]*J0[2][0]))
            print(f"  det J(0) = {det_J0}")
            if det_J0 != 0:
                print(f"  J is nonsingular → CMF IS gauge-equivalent to Sym²(Delannoy)!")
            else:
                print(f"  J(0) is singular → degenerate conic, need further analysis")

            # Also compute J at n=∞ (leading coefficients)
            print(f"\n  J(n) leading coefficients (degree {deg}):")
            for unk_idx in range(6):
                ui, uj = SYM_INDICES[unk_idx]
                lc = x[unk_idx*(deg+1) + deg]
                print(f"    J[{ui},{uj}] leading: {lc}")

            break

print("\nDone.")
