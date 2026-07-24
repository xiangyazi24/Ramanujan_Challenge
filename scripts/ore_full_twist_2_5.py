#!/usr/bin/env python3
"""
Problem 2.5: Compute Ore intertwiner with FULL hypergeometric twist.

H_n = (-16)^n · (2)_n² · (3)_n² · (5/2)_n · (7/2)_n²
q♯(n) = q(n) / H_n  satisfies L₂₅♯ with degree pattern (13,13,13,13).

The intertwiner T maps Sym² solutions to L₂₅♯ solutions:
T(u_j)(n) = t₀(n)·u_j(n) + t₁(n)·u_j(n+1) + t₂(n)·u_j(n+2)

Since both operators now have the same Poincaré roots and are regular
at infinity, t_i should be RATIONAL FUNCTIONS of n.
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce

NMAX = 40

# Delannoy
D = [F(1), F(3)]
E = [F(0), F(1)]
for n in range(1, NMAX + 8):
    D.append((F(3*(2*n+1)) * D[n] - F(n) * D[n-1]) / F(n+1))
    E.append((F(3*(2*n+1)) * E[n] - F(n) * E[n-1]) / F(n+1))

DD = [D[n]**2 for n in range(NMAX+8)]
DE = [D[n]*E[n] for n in range(NMAX+8)]
EE = [E[n]**2 for n in range(NMAX+8)]

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

# Pochhammer
def pochhammer(a_num, a_den, n):
    result = F(1)
    for k in range(n):
        result *= F(a_num + k * a_den, a_den)
    return result

def H(n):
    if n == 0:
        return F(1)
    neg16_n = F(-16)**n
    poch_2 = pochhammer(2, 1, n)
    poch_3 = pochhammer(3, 1, n)
    poch_5_2 = pochhammer(5, 2, n)
    poch_7_2 = pochhammer(7, 2, n)
    return neg16_n * poch_2**2 * poch_3**2 * poch_5_2 * poch_7_2**2

# Compute CMF products
P_mats = [[[1,0,0],[0,1,0],[0,0,1]]]
for n in range(NMAX + 8):
    P_mats.append(mat_mul_int(P_mats[-1], M_int(n)))

def q(n, j):
    return F(P_mats[n][0][j])

# Twist factors
H_vals = [H(n) for n in range(NMAX + 8)]

def q_tw(n, j):
    return q(n, j) / H_vals[n]

def inv3(M):
    a,b,c = M[0]; d,e,f = M[1]; g,h,i_ = M[2]
    det = a*(e*i_-f*h) - b*(d*i_-f*g) + c*(d*h-e*g)
    if det == 0:
        return None
    id_ = F(1) / det
    return [[(e*i_-f*h)*id_, (c*h-b*i_)*id_, (b*f-c*e)*id_],
            [(f*g-d*i_)*id_, (a*i_-c*g)*id_, (c*d-a*f)*id_],
            [(d*h-e*g)*id_, (b*g-a*h)*id_, (a*e-b*d)*id_]]

def mat_mul_F(A, B):
    n = len(A); m = len(B[0]); k = len(B)
    return [[sum(A[i][l]*B[l][j] for l in range(k)) for j in range(m)] for i in range(n)]

print("=== Fully-twisted Ore intertwiner T(n) = Φ_tw(n) · Φ_S(n)⁻¹ ===\n")

t_vals = []
for n in range(NMAX):
    Phi_S = [[DD[n], DE[n], EE[n]],
             [DD[n+1], DE[n+1], EE[n+1]],
             [DD[n+2], DE[n+2], EE[n+2]]]

    Phi_tw = [[q_tw(n,0), q_tw(n,1), q_tw(n,2)],
              [q_tw(n+1,0), q_tw(n+1,1), q_tw(n+1,2)],
              [q_tw(n+2,0), q_tw(n+2,1), q_tw(n+2,2)]]

    Phi_S_inv = inv3(Phi_S)
    if Phi_S_inv is None:
        t_vals.append(None)
        continue

    T = mat_mul_F(Phi_tw, Phi_S_inv)
    t0, t1, t2 = T[0]
    t_vals.append((t0, t1, t2))

    if n < 15:
        print(f"n={n}:")
        print(f"  t₀ = {float(t0):.8e}  (num: {len(str(abs(t0.numerator)))} digits, den: {len(str(t0.denominator))} digits)")
        print(f"  t₁ = {float(t1):.8e}  (num: {len(str(abs(t1.numerator)))} digits, den: {len(str(t1.denominator))} digits)")
        print(f"  t₂ = {float(t2):.8e}  (num: {len(str(abs(t2.numerator)))} digits, den: {len(str(t2.denominator))} digits)")

# Check: do t_i grow polynomially or still exponentially?
print("\n=== Growth check: |t₀(n)| ===")
for n in [0, 5, 10, 15, 20, 25, 30, 35]:
    if n < NMAX and t_vals[n] is not None:
        t0 = t_vals[n][0]
        print(f"  n={n}: {float(t0):.6e}  (log10 = {len(str(abs(t0.numerator)))-1:.0f} digits)")

# Consecutive ratios
print("\n=== t₀(n+1)/t₀(n) — should approach a constant if rational ===")
for n in range(min(NMAX-1, 20)):
    if t_vals[n] is None or t_vals[n+1] is None: continue
    t0_n = t_vals[n][0]
    t0_np1 = t_vals[n+1][0]
    if t0_n != 0:
        r = t0_np1 / t0_n
        print(f"  n={n}: {float(r):.12f}")

# Try rational function fitting
print("\n=== Rational function fit for t₀(n) ===")
t0_data = [(n, t_vals[n][0]) for n in range(NMAX) if t_vals[n] is not None]

def try_rational_fit(data, dp, dq):
    num_unknowns = (dp + 1) + dq
    if len(data) < num_unknowns + 2:
        return None
    rows = []; rhs = []
    for n_val, r_val in data[:num_unknowns + 2]:
        row = [F(n_val)**j for j in range(dp+1)] + [-r_val * F(n_val)**j for j in range(1, dq+1)]
        rows.append(row)
        rhs.append(r_val)
    m = len(rows); ncols = num_unknowns
    mat = [list(rows[i]) + [rhs[i]] for i in range(m)]
    pivots = []; col = 0
    for row_idx in range(min(m, ncols)):
        if col >= ncols: break
        pivot_row = None
        for r in range(row_idx, m):
            if mat[r][col] != 0: pivot_row = r; break
        if pivot_row is None: col += 1; continue
        mat[row_idx], mat[pivot_row] = mat[pivot_row], mat[row_idx]
        pv = mat[row_idx][col]
        for c in range(ncols+1): mat[row_idx][c] /= pv
        for r in range(m):
            if r != row_idx and mat[r][col] != 0:
                fac = mat[r][col]
                for c in range(ncols+1): mat[r][c] -= fac * mat[row_idx][c]
        pivots.append(col); col += 1
    if len(pivots) < ncols: return None
    solution = [mat[i][ncols] for i in range(ncols)]
    p_coeffs = solution[:dp+1]
    q_coeffs = [F(1)] + solution[dp+1:]
    for n_val, r_val in data[num_unknowns+2:]:
        p_val = sum(c * F(n_val)**j for j,c in enumerate(p_coeffs))
        q_val = sum(c * F(n_val)**j for j,c in enumerate(q_coeffs))
        if q_val == 0: continue
        if p_val / q_val != r_val: return None
    return p_coeffs, q_coeffs

for dp in range(20):
    for dq in range(max(0, dp-5), dp+3):
        if dp + 1 + dq + 4 > len(t0_data): break
        result = try_rational_fit(t0_data, dp, dq)
        if result is not None:
            p_coeffs, q_coeffs = result
            print(f"  FOUND t₀: deg(p)={dp}, deg(q)={dq}")
            # Clear denominators and simplify
            all_c = p_coeffs + q_coeffs
            denoms = [c.denominator for c in all_c if c != 0]
            L = reduce(lambda a,b: a*b//gcd(a,b), denoms) if denoms else 1
            p_int = [int(c * L) for c in p_coeffs]
            q_int = [int(c * L) for c in q_coeffs]
            g_val = reduce(gcd, [abs(x) for x in p_int + q_int if x != 0])
            p_int = [x // g_val for x in p_int]
            q_int = [x // g_val for x in q_int]
            from sympy import Symbol, factor as sfact
            N = Symbol('N')
            p_sym = sum(c * N**j for j,c in enumerate(p_int))
            q_sym = sum(c * N**j for j,c in enumerate(q_int))
            print(f"  Factored p(N): {sfact(p_sym)}")
            print(f"  Factored q(N): {sfact(q_sym)}")
            break
    else:
        continue
    break
else:
    print("  No rational fit found up to degree 19")

# Same for t₁
print("\n=== Rational function fit for t₁(n) ===")
t1_data = [(n, t_vals[n][1]) for n in range(NMAX) if t_vals[n] is not None]
for dp in range(20):
    for dq in range(max(0, dp-5), dp+3):
        if dp + 1 + dq + 4 > len(t1_data): break
        result = try_rational_fit(t1_data, dp, dq)
        if result is not None:
            p_coeffs, q_coeffs = result
            print(f"  FOUND t₁: deg(p)={dp}, deg(q)={dq}")
            all_c = p_coeffs + q_coeffs
            denoms = [c.denominator for c in all_c if c != 0]
            L = reduce(lambda a,b: a*b//gcd(a,b), denoms) if denoms else 1
            p_int = [int(c * L) for c in p_coeffs]
            q_int = [int(c * L) for c in q_coeffs]
            g_val = reduce(gcd, [abs(x) for x in p_int + q_int if x != 0])
            p_int = [x // g_val for x in p_int]
            q_int = [x // g_val for x in q_int]
            from sympy import Symbol, factor as sfact
            N = Symbol('N')
            p_sym = sum(c * N**j for j,c in enumerate(p_int))
            q_sym = sum(c * N**j for j,c in enumerate(q_int))
            print(f"  Factored p(N): {sfact(p_sym)}")
            print(f"  Factored q(N): {sfact(q_sym)}")
            break
    else:
        continue
    break
else:
    print("  No rational fit found up to degree 19")

# Same for t₂
print("\n=== Rational function fit for t₂(n) ===")
t2_data = [(n, t_vals[n][2]) for n in range(NMAX) if t_vals[n] is not None]
for dp in range(20):
    for dq in range(max(0, dp-5), dp+3):
        if dp + 1 + dq + 4 > len(t2_data): break
        result = try_rational_fit(t2_data, dp, dq)
        if result is not None:
            p_coeffs, q_coeffs = result
            print(f"  FOUND t₂: deg(p)={dp}, deg(q)={dq}")
            all_c = p_coeffs + q_coeffs
            denoms = [c.denominator for c in all_c if c != 0]
            L = reduce(lambda a,b: a*b//gcd(a,b), denoms) if denoms else 1
            p_int = [int(c * L) for c in p_coeffs]
            q_int = [int(c * L) for c in q_coeffs]
            g_val = reduce(gcd, [abs(x) for x in p_int + q_int if x != 0])
            p_int = [x // g_val for x in p_int]
            q_int = [x // g_val for x in q_int]
            from sympy import Symbol, factor as sfact
            N = Symbol('N')
            p_sym = sum(c * N**j for j,c in enumerate(p_int))
            q_sym = sum(c * N**j for j,c in enumerate(q_int))
            print(f"  Factored p(N): {sfact(p_sym)}")
            print(f"  Factored q(N): {sfact(q_sym)}")
            break
    else:
        continue
    break
else:
    print("  No rational fit found up to degree 19")

print("\nDone.")
