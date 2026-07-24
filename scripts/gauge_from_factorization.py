#!/usr/bin/env python3
"""
Problem 2.5: Extract gauge R(N) = r₀(N+1)/r₀(N) via operator factorization.

Strategy:
1. Find the Sym²(Delannoy) scalar recurrence by guessing from values
2. Apply L_Sym² to q(N) to get g(N) = L_Sym²[q](N)
3. Check if g(N+1)/g(N) is a rational function of N
4. The gauge ratio R(N) is then determined by this factorization

If L_CMF = L₁ ∘ L_Sym² (left factorization), then:
- L_Sym²[q] = g satisfies the first-order L₁
- g(N+1)/g(N) = rational function of N
"""
from fractions import Fraction as F
import sys

# Delannoy numbers (exact)
D = [F(1), F(3)]
E = [F(0), F(1)]
NMAX = 50
for n in range(1, NMAX):
    D.append((F(3*(2*n+1)) * D[n] - F(n) * D[n-1]) / F(n+1))
    E.append((F(3*(2*n+1)) * E[n] - F(n) * E[n-1]) / F(n+1))

# Sym² basis: D², DE, E²
DD = [D[n]**2 for n in range(NMAX+1)]
DE = [D[n]*E[n] for n in range(NMAX+1)]
EE = [E[n]**2 for n in range(NMAX+1)]

# Step 1: Find the order-3 scalar recurrence for D(n)²
# a₀(n)·D(n)² + a₁(n)·D(n+1)² + a₂(n)·D(n+2)² + a₃(n)·D(n+3)² = 0
# where aᵢ(n) = ∑ c_{i,j} n^j, degree d polynomials

print("=== Finding Sym²(Delannoy) scalar recurrence ===")
print("Looking for order-3 recurrence for D(n)²...")

def find_recurrence(seq, order, max_deg):
    """Find recurrence of given order with polynomial coefficients of degree max_deg."""
    num_coeffs = (order + 1) * (max_deg + 1)
    num_eqs = len(seq) - order
    if num_eqs < num_coeffs:
        return None

    # Build the system: for each n, ∑_{i=0}^{order} ∑_{j=0}^{max_deg} c_{i,j} · n^j · seq[n+i] = 0
    rows = []
    for n in range(min(num_eqs, num_coeffs + 5)):  # slight overdetermination
        row = []
        for i in range(order + 1):
            for j in range(max_deg + 1):
                row.append(F(n)**j * seq[n + i])
        rows.append(row)

    # Gaussian elimination over Q
    m = len(rows)
    ncols = num_coeffs
    mat = [list(row) for row in rows]

    pivots = []
    col = 0
    for row_idx in range(m):
        if col >= ncols:
            break
        # Find pivot
        pivot_row = None
        for r in range(row_idx, m):
            if mat[r][col] != 0:
                pivot_row = r
                break
        if pivot_row is None:
            col += 1
            continue
        mat[row_idx], mat[pivot_row] = mat[pivot_row], mat[row_idx]
        pivot_val = mat[row_idx][col]
        for c in range(ncols):
            mat[row_idx][c] /= pivot_val
        for r in range(m):
            if r != row_idx and mat[r][col] != 0:
                factor = mat[r][col]
                for c in range(ncols):
                    mat[r][c] -= factor * mat[row_idx][c]
        pivots.append(col)
        col += 1

    # Find null space
    free_cols = [c for c in range(ncols) if c not in pivots]
    if not free_cols:
        return None

    # Take first free variable = 1, rest = 0
    solution = [F(0)] * ncols
    fc = free_cols[0]
    solution[fc] = F(1)
    for row_idx, pc in enumerate(pivots):
        solution[pc] = -mat[row_idx][fc]

    # Extract polynomials
    polys = []
    idx = 0
    for i in range(order + 1):
        coeffs = solution[idx:idx + max_deg + 1]
        polys.append(coeffs)
        idx += max_deg + 1

    return polys

def eval_poly(coeffs, n):
    return sum(c * F(n)**j for j, c in enumerate(coeffs))

# Try increasing degrees
for deg in range(1, 8):
    result = find_recurrence(DD, 3, deg)
    if result is not None:
        # Verify
        ok = True
        for n in range(30):
            val = sum(eval_poly(result[i], n) * DD[n+i] for i in range(4))
            if val != 0:
                ok = False
                break
        if ok:
            print(f"  Found recurrence with degree-{deg} coefficients!")
            # Also check DE and EE satisfy same recurrence
            ok_DE = all(sum(eval_poly(result[i], n) * DE[n+i] for i in range(4)) == 0 for n in range(30))
            ok_EE = all(sum(eval_poly(result[i], n) * EE[n+i] for i in range(4)) == 0 for n in range(30))
            print(f"  DE satisfies it: {ok_DE}")
            print(f"  EE satisfies it: {ok_EE}")

            # Print the polynomials
            for i in range(4):
                # Clear denominators
                denoms = [c.denominator for c in result[i] if c != 0]
                if denoms:
                    from math import lcm
                    from functools import reduce
                    L = reduce(lcm, denoms)
                    cleared = [int(c * L) for c in result[i]]
                    from math import gcd as mgcd
                    g = reduce(mgcd, [abs(x) for x in cleared if x != 0])
                    cleared = [x // g for x in cleared]
                else:
                    cleared = [0] * len(result[i])

                terms = []
                for j, c in enumerate(cleared):
                    if c == 0:
                        continue
                    if j == 0:
                        terms.append(str(c))
                    elif j == 1:
                        terms.append(f"{c}n")
                    else:
                        terms.append(f"{c}n^{j}")
                print(f"  a_{i}(n) = {' + '.join(terms)}")

            sym2_rec = result
            break
    else:
        print(f"  degree {deg}: no solution")
else:
    print("  Failed to find recurrence up to degree 7!")
    sys.exit(1)

# Step 2: CMF scalar q(N) via matrix product
print("\n=== Computing CMF scalar q(N) ===")

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

N_COMPUTE = 30
prod = [[1,0,0],[0,1,0],[0,0,1]]
q = [1]
for N in range(N_COMPUTE + 10):
    prod = mat_mul(prod, M_int(N))
    q.append(prod[0][0])
q_F = [F(x) for x in q]

# Step 3: Compute g(N) = L_Sym²[q](N) = ∑ a_k(N) · q(N+k)
print("\n=== Computing g(N) = L_Sym²[q](N) ===")
g_vals = []
for N in range(N_COMPUTE):
    g = sum(eval_poly(sym2_rec[k], N) * q_F[N+k] for k in range(4))
    g_vals.append(g)
    if N < 5:
        print(f"  g({N}) = {g}")

# Step 4: Check g(N+1)/g(N)
print("\n=== g(N+1)/g(N) ratios ===")
g_ratios = []
for N in range(min(len(g_vals) - 1, 20)):
    if g_vals[N] == 0:
        print(f"  N={N}: g(N)=0, skip")
        continue
    ratio = g_vals[N+1] / g_vals[N]
    g_ratios.append((N, ratio))
    if N < 15:
        nd = len(str(abs(ratio.numerator)))
        dd = len(str(abs(ratio.denominator)))
        print(f"  N={N}: g({N+1})/g({N}) = ({nd}d)/({dd}d) ≈ {float(ratio):.10e}")

# Step 5: Try to identify g(N+1)/g(N) as a rational function of N
# Use Thiele interpolation (continued fraction interpolation for rational functions)
print("\n=== Identifying g-ratio as rational function ===")

# Method: if R(N) = p(N)/q(N), then R(N)·q(N) - p(N) = 0
# Try degree (d_p, d_q) with d_p + d_q + 2 = #data points
# Start with small degrees and increase

def try_rational_fit(data_points, dp, dq):
    """Try to fit R(N) = p(N)/q(N) with deg p = dp, deg q = dq.
    data_points = [(n, R(n)), ...]
    Returns (p_coeffs, q_coeffs) or None.
    """
    # R(n) * q(n) = p(n) for each data point
    # ∑ p_j n^j - R(n) ∑ q_j n^j = 0
    # Fix q_0 = 1 (normalization)
    # Then: ∑_{j=0}^{dp} p_j n^j - R(n) - R(n) ∑_{j=1}^{dq} q_j n^j = 0
    # i.e., ∑_{j=0}^{dp} p_j n^j - ∑_{j=1}^{dq} q_j R(n) n^j = R(n)

    num_unknowns = (dp + 1) + dq  # p_0..p_dp, q_1..q_dq
    if len(data_points) < num_unknowns:
        return None

    rows = []
    rhs = []
    for n_val, r_val in data_points[:num_unknowns + 2]:
        row = []
        # p coefficients
        for j in range(dp + 1):
            row.append(F(n_val)**j)
        # q coefficients (q_1, ..., q_dq)
        for j in range(1, dq + 1):
            row.append(-r_val * F(n_val)**j)
        rows.append(row)
        rhs.append(r_val)

    # Solve with Gaussian elimination
    m = len(rows)
    ncols = num_unknowns
    mat = [list(rows[i]) + [rhs[i]] for i in range(m)]

    pivots = []
    col = 0
    for row_idx in range(min(m, ncols)):
        if col >= ncols:
            break
        pivot_row = None
        for r in range(row_idx, m):
            if mat[r][col] != 0:
                pivot_row = r
                break
        if pivot_row is None:
            col += 1
            continue
        mat[row_idx], mat[pivot_row] = mat[pivot_row], mat[row_idx]
        pivot_val = mat[row_idx][col]
        for c in range(ncols + 1):
            mat[row_idx][c] /= pivot_val
        for r in range(m):
            if r != row_idx and mat[r][col] != 0:
                factor = mat[r][col]
                for c in range(ncols + 1):
                    mat[r][c] -= factor * mat[row_idx][c]
        pivots.append(col)
        col += 1

    if len(pivots) < ncols:
        return None  # underdetermined

    solution = [mat[i][ncols] for i in range(ncols)]
    p_coeffs = solution[:dp+1]
    q_coeffs = [F(1)] + solution[dp+1:]

    # Verify on remaining data points
    for n_val, r_val in data_points[num_unknowns:]:
        p_val = sum(c * F(n_val)**j for j, c in enumerate(p_coeffs))
        q_val = sum(c * F(n_val)**j for j, c in enumerate(q_coeffs))
        if q_val == 0:
            continue
        predicted = p_val / q_val
        if predicted != r_val:
            return None

    return p_coeffs, q_coeffs

for dp in range(1, 25):
    for dq in range(max(0, dp-10), dp+1):
        result = try_rational_fit(g_ratios, dp, dq)
        if result is not None:
            p_coeffs, q_coeffs = result
            # Clear denominators for nice printing
            all_coeffs = p_coeffs + q_coeffs
            denoms = [c.denominator for c in all_coeffs if c != 0]
            if denoms:
                from math import lcm
                from functools import reduce
                L = reduce(lcm, denoms)
            else:
                L = 1

            p_int = [int(c * L) for c in p_coeffs]
            q_int = [int(c * L) for c in q_coeffs]

            # Simplify by GCD
            from math import gcd as mgcd
            all_nonzero = [abs(x) for x in p_int + q_int if x != 0]
            if all_nonzero:
                g_val = reduce(mgcd, all_nonzero)
                p_int = [x // g_val for x in p_int]
                q_int = [x // g_val for x in q_int]

            def format_poly(coeffs, var='N'):
                terms = []
                for j, c in enumerate(coeffs):
                    if c == 0:
                        continue
                    if j == 0:
                        terms.append(str(c))
                    elif j == 1:
                        if c == 1:
                            terms.append(var)
                        elif c == -1:
                            terms.append(f'-{var}')
                        else:
                            terms.append(f'{c}{var}')
                    else:
                        if c == 1:
                            terms.append(f'{var}^{j}')
                        elif c == -1:
                            terms.append(f'-{var}^{j}')
                        else:
                            terms.append(f'{c}{var}^{j}')
                return ' + '.join(terms).replace('+ -', '- ') if terms else '0'

            print(f"\n  FOUND: deg(p)={dp}, deg(q)={dq}")
            print(f"  g(N+1)/g(N) = [{format_poly(p_int)}] / [{format_poly(q_int)}]")

            # Factor using sympy
            from sympy import Symbol, Poly, factor as sym_factor
            N = Symbol('N')
            p_sym = sum(c * N**j for j, c in enumerate(p_int))
            q_sym = sum(c * N**j for j, c in enumerate(q_int))
            print(f"\n  Factored numerator: {sym_factor(p_sym)}")
            print(f"  Factored denominator: {sym_factor(q_sym)}")

            # Now use this to verify against gauge ratios
            print("\n=== Verifying against Casorati gauge ratios ===")

            # Compute gauge ratios from Casorati
            r0_vals = []
            for NN in range(20):
                Phi = [[DD[NN], DE[NN], EE[NN]],
                       [DD[NN+1], DE[NN+1], EE[NN+1]],
                       [DD[NN+2], DE[NN+2], EE[NN+2]]]
                b = [q_F[NN], q_F[NN+1], q_F[NN+2]]
                # Solve 3x3 exactly
                det = (Phi[0][0]*(Phi[1][1]*Phi[2][2]-Phi[1][2]*Phi[2][1])
                     - Phi[0][1]*(Phi[1][0]*Phi[2][2]-Phi[1][2]*Phi[2][0])
                     + Phi[0][2]*(Phi[1][0]*Phi[2][1]-Phi[1][1]*Phi[2][0]))
                if det == 0:
                    r0_vals.append(None)
                    continue
                det0 = (b[0]*(Phi[1][1]*Phi[2][2]-Phi[1][2]*Phi[2][1])
                      - Phi[0][1]*(b[1]*Phi[2][2]-Phi[1][2]*b[2])
                      + Phi[0][2]*(b[1]*Phi[2][1]-Phi[1][1]*b[2]))
                r0_vals.append(det0 / det)

            for NN in range(1, 15):
                if r0_vals[NN] is None or r0_vals[NN-1] is None or r0_vals[NN-1] == 0:
                    continue
                casorati_ratio = r0_vals[NN] / r0_vals[NN-1]

                # The gauge ratio R(N) where r₀(N+1)/r₀(N) = R(N)
                # ... this is NOT directly the g-ratio.
                # The relationship is more complex.
                if NN <= 8:
                    print(f"  N={NN}: Casorati ratio ≈ {float(casorati_ratio):.10e}")
                    print(f"         g-ratio at N={NN-1} ≈ {float(g_ratios[NN-1][1]) if NN-1 < len(g_ratios) else 'N/A':.10e}" if NN-1 < len(g_ratios) else "")

            break
    else:
        continue
    break
else:
    print("\n  No rational function found for g-ratio up to degree 24!")
    # Print more data
    print("\n  First 10 g-ratios:")
    for N, ratio in g_ratios[:10]:
        print(f"    N={N}: {float(ratio):.15e}")

print("\nDone.")
