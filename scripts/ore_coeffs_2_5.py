#!/usr/bin/env python3
"""
Problem 2.5: Compute the Ore intertwiner coefficients t₀(n), t₁(n), t₂(n).

From Φ_C(n) = T(n) · Φ_S(n), the first row of T(n) gives [t₀(n), t₁(n), t₂(n)].

Φ_C(n) = Casorati matrix of 3 independent CMF solutions.
Φ_S(n) = Casorati matrix of Sym² basis {D², DE, E²}.

The Ore intertwiner is T = t₀(n) + t₁(n)·S + t₂(n)·S² (rational in n).
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce

# Delannoy
D = [F(1), F(3)]
E = [F(0), F(1)]
for n in range(1, 40):
    D.append((F(3*(2*n+1)) * D[n] - F(n) * D[n-1]) / F(n+1))
    E.append((F(3*(2*n+1)) * E[n] - F(n) * E[n-1]) / F(n+1))

DD = [D[n]**2 for n in range(40)]
DE = [D[n]*E[n] for n in range(40)]
EE = [E[n]**2 for n in range(40)]

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

def mat_mul_int(A, B):
    n = len(A)
    return [[sum(A[i][k]*B[k][j] for k in range(n)) for j in range(n)] for i in range(n)]

# Compute P(n) = M(0)·M(1)·...·M(n-1) for n=0,...,N
N_MAX = 18
P_mats = [[[1,0,0],[0,1,0],[0,0,1]]]  # P(0) = I
for n in range(N_MAX + 5):
    P_mats.append(mat_mul_int(P_mats[-1], M_int(n)))

# Three CMF solutions (first, second, third columns of P(n)):
# q_j(n) = P(n)_{0,j}
def q(n, j):
    return F(P_mats[n][0][j])

# Build Φ_C(n) and Φ_S(n), compute T(n) = Φ_C · Φ_S⁻¹
def inv3(M):
    """Exact 3x3 inverse over Q."""
    a,b,c = M[0]
    d,e,f = M[1]
    g,h,i = M[2]
    det = a*(e*i-f*h) - b*(d*i-f*g) + c*(d*h-e*g)
    if det == 0:
        return None
    id = F(1) / det
    return [[(e*i-f*h)*id, (c*h-b*i)*id, (b*f-c*e)*id],
            [(f*g-d*i)*id, (a*i-c*g)*id, (c*d-a*f)*id],
            [(d*h-e*g)*id, (b*g-a*h)*id, (a*e-b*d)*id]]

def mat_mul_F(A, B):
    n = len(A)
    m = len(B[0])
    k = len(B)
    return [[sum(A[i][l]*B[l][j] for l in range(k)) for j in range(m)] for i in range(n)]

print("=== Ore intertwiner coefficients t₀(n), t₁(n), t₂(n) ===\n")

t_vals = []  # list of (t0, t1, t2) at each n
for n in range(N_MAX):
    # Φ_S(n): Casorati of Sym² basis
    Phi_S = [[DD[n], DE[n], EE[n]],
             [DD[n+1], DE[n+1], EE[n+1]],
             [DD[n+2], DE[n+2], EE[n+2]]]

    # Φ_C(n): Casorati of 3 CMF solutions
    Phi_C = [[q(n,0), q(n,1), q(n,2)],
             [q(n+1,0), q(n+1,1), q(n+1,2)],
             [q(n+2,0), q(n+2,1), q(n+2,2)]]

    Phi_S_inv = inv3(Phi_S)
    if Phi_S_inv is None:
        t_vals.append(None)
        continue

    T = mat_mul_F(Phi_C, Phi_S_inv)
    t0, t1, t2 = T[0]
    t_vals.append((t0, t1, t2))

    if n < 10:
        print(f"n={n}:")
        print(f"  t₀ = {float(t0):.10e}  ({len(str(abs(t0.numerator)))}d/{len(str(t0.denominator))}d)")
        print(f"  t₁ = {float(t1):.10e}  ({len(str(abs(t1.numerator)))}d/{len(str(t1.denominator))}d)")
        print(f"  t₂ = {float(t2):.10e}  ({len(str(abs(t2.numerator)))}d/{len(str(t2.denominator))}d)")

# Check: is t₀(n)/t₂(n) a rational function of n?
print("\n=== Ratios ===")
for n in range(min(N_MAX, 12)):
    if t_vals[n] is None:
        continue
    t0, t1, t2 = t_vals[n]
    if t2 != 0:
        r01 = t0 / t1 if t1 != 0 else None
        r02 = t0 / t2
        r12 = t1 / t2 if t1 != 0 else None
        print(f"n={n}: t₀/t₂ = {float(r02):.10f}, t₁/t₂ = {float(r12) if r12 else 'N/A':.10f}")

# Check: t₂(n+1)/t₂(n) — should be a rational function of n
print("\n=== t₂(n+1)/t₂(n) ratio ===")
for n in range(min(N_MAX-1, 15)):
    if t_vals[n] is None or t_vals[n+1] is None:
        continue
    t2_n = t_vals[n][2]
    t2_np1 = t_vals[n+1][2]
    if t2_n != 0:
        ratio = t2_np1 / t2_n
        print(f"  n={n}: t₂({n+1})/t₂({n}) = {float(ratio):.10e}  ({len(str(abs(ratio.numerator)))}d/{len(str(ratio.denominator))}d)")

# Also check t₀ ratio
print("\n=== t₀(n+1)/t₀(n) ratio ===")
for n in range(min(N_MAX-1, 15)):
    if t_vals[n] is None or t_vals[n+1] is None:
        continue
    t0_n = t_vals[n][0]
    t0_np1 = t_vals[n+1][0]
    if t0_n != 0:
        ratio = t0_np1 / t0_n
        print(f"  n={n}: t₀({n+1})/t₀({n}) = {float(ratio):.10e}")

# Try to identify t₂(n+1)/t₂(n) as a rational function
# by fitting p(n)/q(n) from the values
print("\n\n=== Rational function fit for t₂ ratio ===")

t2_ratios = []
for n in range(min(N_MAX-1, 15)):
    if t_vals[n] is None or t_vals[n+1] is None:
        continue
    t2_n = t_vals[n][2]
    t2_np1 = t_vals[n+1][2]
    if t2_n != 0:
        t2_ratios.append((n, t2_np1 / t2_n))

def try_rational_fit(data, dp, dq):
    """Fit R(n) = p(n)/q(n) to data points [(n, R(n)), ...]."""
    num_unknowns = (dp + 1) + dq  # p_0..p_dp, q_1..q_dq (q_0=1)
    if len(data) < num_unknowns:
        return None

    rows = []
    rhs = []
    for n_val, r_val in data[:num_unknowns + 2]:
        row = []
        for j in range(dp + 1):
            row.append(F(n_val)**j)
        for j in range(1, dq + 1):
            row.append(-r_val * F(n_val)**j)
        rows.append(row)
        rhs.append(r_val)

    # Solve
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
        pv = mat[row_idx][col]
        for c in range(ncols + 1):
            mat[row_idx][c] /= pv
        for r in range(m):
            if r != row_idx and mat[r][col] != 0:
                fac = mat[r][col]
                for c in range(ncols + 1):
                    mat[r][c] -= fac * mat[row_idx][c]
        pivots.append(col)
        col += 1

    if len(pivots) < ncols:
        return None

    solution = [mat[i][ncols] for i in range(ncols)]
    p_coeffs = solution[:dp+1]
    q_coeffs = [F(1)] + solution[dp+1:]

    # Verify
    for n_val, r_val in data[num_unknowns:]:
        p_val = sum(c * F(n_val)**j for j, c in enumerate(p_coeffs))
        q_val = sum(c * F(n_val)**j for j, c in enumerate(q_coeffs))
        if q_val == 0:
            continue
        if p_val / q_val != r_val:
            return None
    return p_coeffs, q_coeffs

for dp in range(1, 15):
    for dq in range(max(0, dp-8), dp+1):
        result = try_rational_fit(t2_ratios, dp, dq)
        if result is not None:
            p_coeffs, q_coeffs = result
            print(f"  FOUND t₂ ratio: deg(p)={dp}, deg(q)={dq}")
            from sympy import Symbol, factor as sym_factor
            N = Symbol('N')

            # Clear denominators
            all_c = p_coeffs + q_coeffs
            denoms = [c.denominator for c in all_c if c != 0]
            L = reduce(lambda a,b: a*b//gcd(a,b), denoms) if denoms else 1
            p_int = [int(c * L) for c in p_coeffs]
            q_int = [int(c * L) for c in q_coeffs]
            g_val = reduce(gcd, [abs(x) for x in p_int + q_int if x != 0])
            p_int = [x // g_val for x in p_int]
            q_int = [x // g_val for x in q_int]

            p_sym = sum(c * N**j for j, c in enumerate(p_int))
            q_sym = sum(c * N**j for j, c in enumerate(q_int))
            print(f"  Factored num: {sym_factor(p_sym)}")
            print(f"  Factored den: {sym_factor(q_sym)}")
            break
    else:
        continue
    break
else:
    print("  No fit found up to degree 14")

# Similarly for t₀ and t₁ ratios
print("\n=== Rational function fit for t₀ ratio ===")
t0_ratios = []
for n in range(min(N_MAX-1, 15)):
    if t_vals[n] is None or t_vals[n+1] is None:
        continue
    t0_n = t_vals[n][0]
    t0_np1 = t_vals[n+1][0]
    if t0_n != 0:
        t0_ratios.append((n, t0_np1 / t0_n))

for dp in range(1, 15):
    for dq in range(max(0, dp-8), dp+1):
        result = try_rational_fit(t0_ratios, dp, dq)
        if result is not None:
            p_coeffs, q_coeffs = result
            print(f"  FOUND t₀ ratio: deg(p)={dp}, deg(q)={dq}")
            from sympy import Symbol, factor as sym_factor
            N = Symbol('N')
            all_c = p_coeffs + q_coeffs
            denoms = [c.denominator for c in all_c if c != 0]
            L = reduce(lambda a,b: a*b//gcd(a,b), denoms) if denoms else 1
            p_int = [int(c * L) for c in p_coeffs]
            q_int = [int(c * L) for c in q_coeffs]
            g_val = reduce(gcd, [abs(x) for x in p_int + q_int if x != 0])
            p_int = [x // g_val for x in p_int]
            q_int = [x // g_val for x in q_int]
            p_sym = sum(c * N**j for j, c in enumerate(p_int))
            q_sym = sum(c * N**j for j, c in enumerate(q_int))
            print(f"  Factored num: {sym_factor(p_sym)}")
            print(f"  Factored den: {sym_factor(q_sym)}")
            break
    else:
        continue
    break
else:
    print("  No fit found up to degree 14")

print("\nDone.")
