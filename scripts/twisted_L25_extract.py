#!/usr/bin/env python3
"""
Problem 2.5: Extract the twisted operator L₂₅♯ from the CMF sequence.

Twist: H_n = (-16)^n · (2)_n² · (3)_n² · (5/2)_n · (7/2)_n²
Twisted sequence: q♯(n) = q(n) / H_n
Expected degree pattern: (13, 13, 13, 13) (uniform)

After extraction, factor each coefficient to identify P₆, P₉, P₁₀, B(n).
"""
from fractions import Fraction as F
from sympy import Symbol, Rational, factor, Poly, expand, simplify, prod as sprod
from functools import reduce
from math import gcd

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

# Pochhammer (rising factorial)
def pochhammer(a_num, a_den, n):
    """Compute (a)_n = a(a+1)...(a+n-1) where a = a_num/a_den, return as Fraction."""
    result = F(1)
    for k in range(n):
        result *= F(a_num + k * a_den, a_den)
    return result

def H(n):
    """Twist factor H_n = (-16)^n · (2)_n² · (3)_n² · (5/2)_n · (7/2)_n²"""
    if n == 0:
        return F(1)
    neg16_n = F(-16)**n
    poch_2 = pochhammer(2, 1, n)     # (2)_n = 2·3·...·(n+1) = (n+1)!/1!
    poch_3 = pochhammer(3, 1, n)     # (3)_n = 3·4·...·(n+2) = (n+2)!/2!
    poch_5_2 = pochhammer(5, 2, n)   # (5/2)_n
    poch_7_2 = pochhammer(7, 2, n)   # (7/2)_n
    return neg16_n * poch_2**2 * poch_3**2 * poch_5_2 * poch_7_2**2

# Compute CMF sequence
print("Computing CMF sequence...")
N_MAX = 90
I3 = [[1,0,0],[0,1,0],[0,0,1]]
prod_mat = [row[:] for row in I3]
q_raw = []
for N in range(N_MAX + 5):
    q_raw.append(prod_mat[0][0])
    if N < N_MAX + 4:
        prod_mat = mat_mul(prod_mat, M_int(N))

# Compute twisted sequence
print("Computing twist factors and twisted sequence...")
H_vals = [H(n) for n in range(N_MAX + 5)]
q_twist = []
for n in range(N_MAX + 5):
    if H_vals[n] == 0:
        q_twist.append(F(0))
    else:
        q_twist.append(F(q_raw[n]) / H_vals[n])

print(f"q_raw[0..4] = {q_raw[:5]}")
print(f"H[0..4] = {[H_vals[i] for i in range(5)]}")
print(f"q♯[0..4] = {[float(q_twist[i]) for i in range(5)]}")

# Verify that δ(n) = H(n+1)/H(n) matches the expected formula
print("\nVerifying twist ratio δ(n) = -2(n+2)²(n+3)²(2n+5)(2n+7)²:")
for n in range(5):
    delta_actual = H_vals[n+1] / H_vals[n] if H_vals[n] != 0 else None
    delta_expected = F(-2) * (n+2)**2 * (n+3)**2 * (2*n+5) * (2*n+7)**2
    if delta_actual == delta_expected:
        print(f"  n={n}: ✓")
    else:
        print(f"  n={n}: ✗ actual={delta_actual}, expected={delta_expected}")

# Extract recurrence with degree pattern (13, 13, 13, 13)
print("\n=== Extracting twisted recurrence (degree 13,13,13,13) ===")
deg = [13, 13, 13, 13]
n_coeffs = [d+1 for d in deg]
total = sum(n_coeffs)  # 56
print(f"Total unknowns: {total}")

n_eqns = total + 15  # overdetermined
A_sys = []
for eq_idx in range(n_eqns):
    N = eq_idx
    row = []
    for i in range(4):
        for m in range(deg[i]+1):
            row.append(F(N)**m * q_twist[N+i])
    A_sys.append(row)

# Fix last unknown = 1 (leading coeff of c_3)
A_mat = [row[:-1] for row in A_sys]
b_vec = [-row[-1] for row in A_sys]
n_vars = total - 1

print(f"System: {len(A_mat)} equations × {n_vars} unknowns")
print("Solving with Gaussian elimination over Q...")

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
            fac = aug[r][col] / pivot_val
            for c2 in range(n_c + 1):
                aug[r][c2] -= fac * aug[piv_row][c2]
    piv_row += 1
    if col % 10 == 0:
        print(f"  pivot {col}/{n_c}...")

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
    p = sum(Rational(c.numerator, c.denominator) * N_sym**m for m, c in enumerate(coeffs))
    polys.append(p)

# Verify recurrence on twisted sequence
print("\n=== Verifying twisted recurrence ===")
def eval_poly_F(coeffs, N):
    return sum(c * F(N)**m for m, c in enumerate(coeffs))

max_residual = F(0)
all_ok = True
poly_coeffs = []
idx = 0
for i in range(4):
    poly_coeffs.append(solution[idx:idx+deg[i]+1])
    idx += deg[i]+1

for N in range(min(70, len(q_twist)-3)):
    val = sum(eval_poly_F(poly_coeffs[i], N) * q_twist[N+i] for i in range(4))
    if val != 0:
        print(f"  NONZERO residual at N={N}: {val}")
        all_ok = False
        break

if all_ok:
    print("  All residuals exactly ZERO. Twisted recurrence verified!")

# Factor each coefficient
print("\n=== Factored twisted recurrence coefficients ===")
for i in range(4):
    p = polys[i]
    poly_obj = Poly(p, N_sym, domain='QQ')
    coeffs_list = poly_obj.all_coeffs()
    denoms = [Rational(c).denominator for c in coeffs_list]
    lcm_d = reduce(lambda a, b: a * b // gcd(a, b), denoms)
    p_int = expand(p * lcm_d)
    print(f"\nℓ_{i}(N) [degree {deg[i]}]:")
    print(f"  Denom factor: 1/{lcm_d}")
    f = factor(p_int)
    print(f"  Factored: {f}")

# Also print raw sympy factor_list for detailed structure
print("\n=== Detailed factor lists ===")
for i in range(4):
    p = polys[i]
    poly_obj = Poly(p, N_sym, domain='QQ')
    coeffs_list = poly_obj.all_coeffs()
    denoms = [Rational(c).denominator for c in coeffs_list]
    lcm_d = reduce(lambda a, b: a * b // gcd(a, b), denoms)
    p_int = Poly(expand(p * lcm_d), N_sym, domain='ZZ')
    try:
        fl = p_int.factor_list()
        print(f"\nℓ_{i}: content={fl[0]}, factors:")
        for fac, mult in fl[1]:
            print(f"  ({fac.as_expr()})^{mult}  [degree {fac.degree()}]")
    except Exception as e:
        print(f"\nℓ_{i}: factor_list failed: {e}")

print("\nDone.")
