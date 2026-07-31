#!/usr/bin/env python3
"""Check if u(n) = (-16)^n is a solution, and explore r(n) structure."""
from fractions import Fraction
import sys

# Rebuild the exact recurrence (fast, ~30s)
def M_mat_qq(n):
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

A_row = [Fraction(33750), Fraction(-36000), Fraction(9000)]
T = [[Fraction(1 if i==j else 0) for j in range(3)] for i in range(3)]
q_vals = []
for N in range(120):
    q = sum(A_row[k] * T[k][0] for k in range(3))
    q_vals.append(q)
    T_new = [[Fraction(0)]*3 for _ in range(3)]
    M = M_mat_qq(N)
    for i in range(3):
        for j in range(3):
            for k in range(3):
                T_new[i][j] += T[i][k] * M[k][j]
    T = T_new

# Gaussian elimination for recurrence
order = 3; degs = [28, 21, 14, 7]
n_unknowns = sum(d+1 for d in degs)
rows = []
for N in range(n_unknowns + 10):
    row = []
    for j in range(order + 1):
        for k in range(degs[j] + 1):
            row.append(Fraction(N)**k * q_vals[N + j])
    rows.append(row)

mat = [list(row) for row in rows]
n_rows = len(mat); n_cols = len(mat[0])
pivot_cols = []
for col in range(n_cols):
    pivot_row = None
    for row in range(len(pivot_cols), n_rows):
        if mat[row][col] != 0:
            pivot_row = row; break
    if pivot_row is None: continue
    mat[pivot_row], mat[len(pivot_cols)] = mat[len(pivot_cols)], mat[pivot_row]
    pivot_row = len(pivot_cols); pivot_cols.append(col)
    piv = mat[pivot_row][col]
    for row in range(n_rows):
        if row == pivot_row: continue
        if mat[row][col] != 0:
            factor = mat[row][col] / piv
            for c in range(n_cols):
                mat[row][c] -= factor * mat[pivot_row][c]

free_cols = [c for c in range(n_cols) if c not in pivot_cols]
null_vec = [Fraction(0)] * n_cols
null_vec[free_cols[0]] = Fraction(1)
for i in range(len(pivot_cols) - 1, -1, -1):
    pc = pivot_cols[i]
    val = -sum(mat[i][c] * null_vec[c] for c in range(n_cols) if c != pc) / mat[i][pc]
    null_vec[pc] = val

idx = 0; polys = []
for j in range(order + 1):
    coeffs = null_vec[idx:idx + degs[j] + 1]
    idx += degs[j] + 1
    polys.append(coeffs)

def eval_cj(j, n):
    return sum(polys[j][k] * Fraction(n)**k for k in range(len(polys[j])))

# === CHECK 1: Does u(n) = (-16)^n satisfy the recurrence? ===
print("=== Check: is (-16)^n a solution? ===")
for n in range(10):
    # c3(n)*(-16)^(n+3) + c2(n)*(-16)^(n+2) + c1(n)*(-16)^(n+1) + c0(n)*(-16)^n = 0?
    # Divide by (-16)^n: c3(n)*(-4096) + c2(n)*(256) + c1(n)*(-16) + c0(n) = 0?
    val = eval_cj(3,n)*(-4096) + eval_cj(2,n)*256 + eval_cj(1,n)*(-16) + eval_cj(0,n)
    print(f"  n={n}: {val}")
    if val != 0:
        print("  ← NOT zero! (-16)^n is NOT a solution.")
        break

# === CHECK 2: Gauged recurrence d_j(n) = (-16)^j * c_j(n) ===
# d3(n) = -4096*c3(n), d2(n) = 256*c2(n), d1(n) = -16*c1(n), d0(n) = c0(n)
# Check if d0+d1+d2+d3 = 0 for all n (constants killing)
print("\n=== Check: gauged (S-1) condition ===")
for n in range(5):
    val = eval_cj(0,n) - 16*eval_cj(1,n) + 256*eval_cj(2,n) - 4096*eval_cj(3,n)
    print(f"  n={n}: d0+d1+d2+d3 = {val}")

# === APPROACH: Factor c3 and c0 to find the Petkovšek roots ===
print("\n=== Factoring c3 (degree 7) ===")
# c3(n) has degree 7. Check rational roots using rational root theorem.
# c3 has LC=1, constant term = polys[3][0]
c3_const = polys[3][0]
print(f"  c3 constant term: {c3_const}")
print(f"  c3 coefficients: {polys[3]}")

# Check half-integer roots: c3(-k/2) = 0?
for num in range(-20, 20):
    for den in [1, 2]:
        n = Fraction(num, den)
        val = eval_cj(3, n)
        if val == 0:
            print(f"  Root: n = {n}")

# === Factor c0 (degree 28) ===
print("\n=== Factoring c0 (degree 28): rational/half-integer roots ===")
for num in range(-20, 5):
    for den in [1, 2]:
        n = Fraction(num, den)
        val = eval_cj(0, n)
        if val == 0:
            print(f"  Root: n = {n}")

# === Compute the gauged v-recurrence and try to find hyper solution ===
# v(n) = u(n)/(-16)^n, recurrence: d3v(n+3) + d2v(n+2) + d1v(n+1) + d0v(n) = 0
# where dj = (-16)^j * cj
# Look for solution v(n) with v(n+1)/v(n) rational → find r_v(n) = v(n+1)/v(n)
# Then the original gauge is r(n) = -16 * r_v(n)

# Check if v(n) = (n!)^7 is close to a solution
# v(n) = (n!)^7, v(n+1)/v(n) = (n+1)^7
# Check: d3(n)*(n+3)^7*(n+2)^7*(n+1)^7 + d2(n)*(n+2)^7*(n+1)^7 + d1(n)*(n+1)^7 + d0(n) = 0?
print("\n=== Check: is v(n) = (n!)^7 a solution of gauged recurrence? ===")
for n in range(5):
    d3 = -4096*eval_cj(3,n)
    d2 = 256*eval_cj(2,n)
    d1 = -16*eval_cj(1,n)
    d0 = eval_cj(0,n)
    val = d3*(n+3)**7*(n+2)**7*(n+1)**7 + d2*(n+2)**7*(n+1)**7 + d1*(n+1)**7 + d0
    print(f"  n={n}: {val}")
    
# Check more general (a)_n^7 with a = 1
# v(n) = ((1)_n)^7 = (n!)^7, same as above

# Try v(n) = Γ(n+1/2)^a * Γ(n+1)^b * ...
# Actually, let me just compute exact v(n) = q(n)/(-16)^n and see if v(n+1)/v(n) is recognizable
print("\n=== Computing v(n) = q(n)/(-16)^n and ratios ===")
for n in range(15):
    v_n = q_vals[n] / Fraction(-16)**n
    v_n1 = q_vals[n+1] / Fraction(-16)**(n+1)
    if v_n != 0:
        rv = v_n1 / v_n
        print(f"  n={n}: v(n+1)/v(n) = {float(rv):.10f}")
