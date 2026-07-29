#!/usr/bin/env python3
"""
Factor the L_25 recurrence coefficients to identify Pochhammer structure.
Extracts the recurrence from CMF, then factors each polynomial over Q.
"""
from fractions import Fraction as F
from sympy import Symbol, Rational, factor, Poly, ZZ, QQ
from sympy import factorint

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

# Compute scalar sequence
print("Computing CMF sequence...")
N_MAX = 90
I3 = [[1,0,0],[0,1,0],[0,0,1]]
prod = [row[:] for row in I3]
q = []
for N in range(N_MAX + 5):
    q.append(prod[0][0])  # first row, first column
    if N < N_MAX + 4:
        prod = mat_mul(prod, M_int(N))

# Extract recurrence with degree pattern (28,21,14,7)
deg = [28, 21, 14, 7]
n_coeffs = [d+1 for d in deg]
total = sum(n_coeffs)

print(f"Extracting recurrence (total {total} unknowns)...")

# Build linear system
n_eqns = total + 10
A_sys = []
for eq_idx in range(n_eqns):
    N = eq_idx
    row = []
    for i in range(4):
        for m in range(deg[i]+1):
            row.append(F(N**m * q[N+i]))
    A_sys.append(row)

# Fix last unknown = 1
A_mat = [row[:-1] for row in A_sys]
b_vec = [-row[-1] for row in A_sys]
n_vars = total - 1

# Gaussian elimination
aug = [A_mat[i][:] + [b_vec[i]] for i in range(len(A_mat))]
n_r = len(aug)
n_c = n_vars

piv_row = 0
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
            factor_val = aug[r][col] / pivot_val
            for c2 in range(n_c + 1):
                aug[r][c2] -= factor_val * aug[piv_row][c2]
    piv_row += 1

solution = [F(0)] * n_c
for r in range(piv_row):
    pc = -1
    for c in range(n_c):
        if aug[r][c] != 0:
            pc = c
            break
    if pc >= 0:
        solution[pc] = aug[r][n_c] / aug[r][pc]
solution.append(F(1))

# Build sympy polynomials
N_sym = Symbol('N')
polys = []
idx = 0
for i in range(4):
    coeffs = solution[idx:idx+deg[i]+1]
    idx += deg[i]+1
    # Convert to sympy Rational
    p = sum(Rational(c.numerator, c.denominator) * N_sym**m for m, c in enumerate(coeffs))
    polys.append(p)

# Factor each polynomial
print("\n=== Factored recurrence coefficients ===")
for i in range(4):
    print(f"\nc_{i}(N) [degree {deg[i]}]:")
    p = polys[i]
    # First clear denominators
    poly_obj = Poly(p, N_sym, domain='QQ')
    # Get integer polynomial by multiplying by LCM of denominators
    coeffs_list = poly_obj.all_coeffs()
    denoms = [Rational(c).denominator for c in coeffs_list]
    from math import gcd
    from functools import reduce
    lcm_d = reduce(lambda a, b: a * b // gcd(a, b), denoms)
    p_int = p * lcm_d
    print(f"  Denominator factor: 1/{lcm_d}")

    # Factor over QQ
    f = factor(p_int)
    print(f"  Factored: {f}")

    # Also try Poly.factor_list for more detail
    poly_int = Poly(p_int, N_sym, domain='ZZ')
    try:
        fl = poly_int.factor_list()
        print(f"  Factor list: {fl}")
    except Exception as e:
        print(f"  (factor_list failed: {e})")

# Also compute the Pochhammer gauge
# We need the ratio c_3(N+k) products to match c_0(N) structure
print("\n\n=== Gauge analysis ===")
print("For the gauge g_N = (-16)^N * prod(shifted Pochhammer):")
print(f"c_3 leading: {polys[3].coeff(N_sym, deg[3])}")
print(f"c_0 leading: {polys[0].coeff(N_sym, deg[0])}")
print(f"Ratio: {polys[0].coeff(N_sym, deg[0]) / polys[3].coeff(N_sym, deg[3])}")
print(f"= 4096 (expected, since Poincaré product = -4096)")
