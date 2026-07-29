#!/usr/bin/env python3
"""
P2.5: Check if the CMF differential operator factors as L_K ∘ D.

The integrated-K module has ODE: k(1-k²)Y''' + (1-3k²)Y'' - kY' = 0
which factors as [k(1-k²)D² + (1-3k²)D - k] ∘ D.

After pullback through k(z) = 4√(2z)/(1-z), the ODE becomes an order-3
differential equation in z. We need to verify this matches the CMF.

Strategy: compute the ODE from the CMF recurrence by computing many terms
of f(z) = Σ Q_n z^n and using numerical differentiation to extract the
operator coefficients.
"""
from fractions import Fraction
from mpmath import mp, mpf, matrix as mpmat, power, sqrt, polyroots, log10
import sys

mp.dps = 100

# Step 1: Compute Q_{N,1} terms using exact Fraction arithmetic
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

# Compute Q_{N,1} for N=0..NTERMS
NTERMS = 150
print(f"Computing {NTERMS} terms of Q_{{N,1}}...", flush=True)

q_row = [Fraction(33750), Fraction(-36000), Fraction(9000)]
Q = [q_row[0]]

for N in range(NTERMS - 1):
    M = M_entries(N)
    q_new = [sum(q_row[i]*M[i][j] for i in range(3)) for j in range(3)]
    q_row = q_new
    Q.append(q_row[0])

print(f"  Done. Q[0]={Q[0]}, Q[1]={Q[1]}", flush=True)

# Step 2: Find the scalar recurrence of order 3
# Q[n] satisfies c₃(n)Q[n+3] + c₂(n)Q[n+2] + c₁(n)Q[n+1] + c₀(n)Q[n] = 0
# with deg(c_k) = 7(3-k), so deg(c₃)=0·7=0... wait no.
# degree pattern (28,21,14,7) means deg(c₃)=7, deg(c₂)=14, deg(c₁)=21, deg(c₀)=28
# Wait, the proof says degree pattern (28,21,14,7) which means:
# c₀ has degree 28, c₁ has degree 21, c₂ has degree 14, c₃ has degree 7.
# Actually no - it says "degree pattern (28,21,14,7)" for the order-3 recurrence
# with 4 coefficients. Leading coefficient c₃ has degree 28, trailing c₀ has degree 7.
# Wait, let me re-read. The proof says:
# "The scalar recurrence extracted via Casorati minors has order 3,
#  degree pattern (28,21,14,7)..."
# This means c₃(n) has degree 7, c₂(n) has degree 14, c₁(n) has degree 21, c₀(n) has degree 28?
# Or c₃ has 28, c₂ has 21, c₁ has 14, c₀ has 7?
#
# The proof also says "Newton polygon has constant slope -7". With 4 terms at shifts
# n+3, n+2, n+1, n, the Newton polygon has points (0, deg c₃), (1, deg c₂), (2, deg c₁), (3, deg c₀).
# Constant slope -7 means deg c_k = 28 - 7k: c₃ at k=0 is 28, c₂ at k=1 is 21, c₁ at k=2 is 14, c₀ at k=3 is 7.
# So c₃ has degree 28, c₀ has degree 7.

# Total unknowns: (28+1) + (21+1) + (14+1) + (7+1) = 29+22+15+8 = 74 coefficients
# minus 1 for overall scaling = 73 free parameters
# We need at least 73 equations: n = 0, 1, ..., 72 gives 73 equations.

print("\nFinding scalar recurrence...", flush=True)
deg3, deg2, deg1, deg0 = 28, 21, 14, 7  # degrees of c₃, c₂, c₁, c₀

# Build the system: for each n, c₃(n)Q[n+3] + c₂(n)Q[n+2] + c₁(n)Q[n+1] + c₀(n)Q[n] = 0
# Unknowns: coefficients of c_k(n) = Σ a_{k,j} n^j for j=0..deg_k
# Total unknowns
total_unk = (deg3+1) + (deg2+1) + (deg1+1) + (deg0+1)
print(f"  Total unknowns: {total_unk}")

neqs = total_unk + 5  # extra equations for overdetermined check
assert NTERMS >= neqs + 3, f"Need {neqs+3} terms but only have {NTERMS}"

# Build matrix: row for each n, columns are the unknown coefficients
# For n: c₃(n)Q[n+3] + ... = Σⱼ a₃ⱼ n^j Q[n+3] + Σⱼ a₂ⱼ n^j Q[n+2] + ...
from fractions import Fraction

print(f"  Building {neqs} × {total_unk} matrix...", flush=True)
rows = []
for n in range(neqs):
    row = []
    for j in range(deg3+1):
        row.append(Fraction(n)**j * Q[n+3])
    for j in range(deg2+1):
        row.append(Fraction(n)**j * Q[n+2])
    for j in range(deg1+1):
        row.append(Fraction(n)**j * Q[n+1])
    for j in range(deg0+1):
        row.append(Fraction(n)**j * Q[n])
    rows.append(row)

# Gaussian elimination to find the kernel
print(f"  Gaussian elimination (exact fractions)...", flush=True)
# Use numpy-like approach but with exact fractions
m = len(rows)
nc = total_unk
mat = [list(row) for row in rows]

# Forward elimination
pivot_cols = []
for col in range(nc):
    # Find non-zero entry in this column
    found = -1
    for r in range(len(pivot_cols), m):
        if mat[r][col] != 0:
            found = r
            break
    if found == -1:
        continue

    # Swap
    mat[len(pivot_cols)], mat[found] = mat[found], mat[len(pivot_cols)]
    pivot_row = len(pivot_cols)
    pivot_cols.append(col)

    # Eliminate
    pivot_val = mat[pivot_row][col]
    for r in range(m):
        if r != pivot_row and mat[r][col] != 0:
            factor = mat[r][col] / pivot_val
            for c2 in range(nc):
                mat[r][c2] -= factor * mat[pivot_row][c2]

    if len(pivot_cols) % 10 == 0:
        print(f"    {len(pivot_cols)}/{nc} pivots...", flush=True)

rank = len(pivot_cols)
nullity = nc - rank
print(f"  Rank = {rank}, Nullity = {nullity}")

if nullity != 1:
    print(f"  ERROR: Expected nullity 1 but got {nullity}")
    sys.exit(1)

# Extract the kernel vector
free_col = set(range(nc)) - set(pivot_cols)
free_col = list(free_col)[0]
print(f"  Free column: {free_col}")

# Back-substitute
kernel = [Fraction(0)] * nc
kernel[free_col] = Fraction(1)
for i in range(rank-1, -1, -1):
    pc = pivot_cols[i]
    val = Fraction(0)
    for c2 in range(nc):
        if c2 != pc:
            val += mat[i][c2] * kernel[c2]
    kernel[pc] = -val / mat[i][pc]

# Extract the polynomial coefficients
c3_coeffs = kernel[:deg3+1]
c2_coeffs = kernel[deg3+1:deg3+1+deg2+1]
c1_coeffs = kernel[deg3+1+deg2+1:deg3+1+deg2+1+deg1+1]
c0_coeffs = kernel[deg3+1+deg2+1+deg1+1:]

def poly_eval(coeffs, n):
    return sum(c * Fraction(n)**k for k, c in enumerate(coeffs))

# Verify the recurrence
print("\nVerifying recurrence...", flush=True)
for n in range(neqs, neqs + 10):
    val = poly_eval(c3_coeffs, n)*Q[n+3] + poly_eval(c2_coeffs, n)*Q[n+2] + \
          poly_eval(c1_coeffs, n)*Q[n+1] + poly_eval(c0_coeffs, n)*Q[n]
    if val != 0:
        print(f"  FAIL at n={n}: residual = {val}")
        break
else:
    print(f"  All 10 extra equations satisfied. Recurrence verified!")

# Print Poincaré analysis
print("\n=== Poincaré polynomial check ===")
# Leading coefficient of c₃(n) is at degree 28
lc3 = c3_coeffs[deg3]  # coefficient of n^28
lc2 = c2_coeffs[deg2]  # coefficient of n^21
lc1 = c1_coeffs[deg1]  # coefficient of n^14
lc0 = c0_coeffs[deg0]  # coefficient of n^7

# Poincaré polynomial: lc3·ξ³ + lc2·ξ² + lc1·ξ + lc0 = 0
# Normalize by lc3
if lc3 != 0:
    a = lc2/lc3
    b = lc1/lc3
    c = lc0/lc3
    print(f"  ξ³ + {float(a):.6f}ξ² + {float(b):.6f}ξ + {float(c):.6f}")
    print(f"  Expected: ξ³ + 560ξ² + 8960ξ + 4096")
    # Check
    print(f"  Ratios: {float(a)/560:.10f}, {float(b)/8960:.10f}, {float(c)/4096:.10f}")

# Also check if c₃ factors: the leading coeff of c₃(n) at n^28
# and the trailing coeff of c₀(n) at n^7
# should give the Poincaré roots
print(f"\n  lc₃ (n^28 coeff) = {float(lc3):.6e}")
print(f"  lc₀ (n^7 coeff)  = {float(lc0):.6e}")

# Print degree structure
for name, coeffs, deg in [("c₃", c3_coeffs, deg3), ("c₂", c2_coeffs, deg2),
                           ("c₁", c1_coeffs, deg1), ("c₀", c0_coeffs, deg0)]:
    # Find actual degree
    actual_deg = deg
    while actual_deg > 0 and coeffs[actual_deg] == 0:
        actual_deg -= 1
    print(f"  {name}: expected deg {deg}, actual deg {actual_deg}")

print("\n=== ODE factorization check ===")
print("The integrated-K ODE factors as L_K ∘ D where")
print("  L_K = k(1-k²)D² + (1-3k²)D - k")
print("After pullback k(z) = 4√(2z)/(1-z):")
print("  k=0 at z=0")
print("  k=1 at z=17-12√2")
print("  k=-1 at z=17+12√2")
print("  k=∞ at z=1")

# The scalar recurrence ←→ differential equation
# n·Q_n ←→ z·f'(z)  (under g.f. Σ Q_n z^n = f(z))
# The recurrence c₃(n)Q_{n+3}+...+c₀(n)Q_n = 0
# translates to a differential equation involving z, f, f', f'', f'''
# with polynomial coefficients in z.

# For the factorization check, we need the differential equation.
# This requires converting the recurrence to an ODE using the map:
# n^k Q_{n+s} ←→ polynomial in (z, zD_z) applied to f(z)

# This is a standard computation. Let me do it numerically.
# At z₀ near 0, compute f, f', f'', f''' and check the ODE.

print("\nChecking ODE factorization numerically...", flush=True)

# Compute f(z) = Σ Q_n z^n at z = 10^{-4} (well inside convergence radius)
z0 = mpf('0.0001')
f0 = sum(mpf(int(Q[n].numerator))/mpf(int(Q[n].denominator)) * z0**n for n in range(NTERMS))
f1 = sum(mpf(int(Q[n].numerator))/mpf(int(Q[n].denominator)) * n * z0**(n-1) for n in range(1, NTERMS))
f2 = sum(mpf(int(Q[n].numerator))/mpf(int(Q[n].denominator)) * n*(n-1) * z0**(n-2) for n in range(2, NTERMS))
f3 = sum(mpf(int(Q[n].numerator))/mpf(int(Q[n].denominator)) * n*(n-1)*(n-2) * z0**(n-3) for n in range(3, NTERMS))

print(f"  At z={z0}: f = {mp.nstr(f0, 20)}")

# Check: does f(z) satisfy the pullback of k(1-k²)Y''' + (1-3k²)Y'' - kY' = 0?
# With k(z) = 4√(2z)/(1-z), compute Y = f(z), Y' = df/dk = (df/dz)/(dk/dz)
k_z = 4*sqrt(2*z0)/(1-z0)
# dk/dz = 2√2(1+z)/(√z·(1-z)²)
dk_dz = 2*sqrt(2)*(1+z0)/(sqrt(z0)*(1-z0)**2)
# d²k/dz²
# k = 4√2 · z^{1/2} · (1-z)^{-1}
# dk/dz = 4√2 · [(1/2)z^{-1/2}(1-z)^{-1} + z^{1/2}(1-z)^{-2}]
#        = 4√2 · z^{-1/2}(1-z)^{-2} · [(1-z)/2 + z]
#        = 4√2 · z^{-1/2}(1-z)^{-2} · (1+z)/2
#        = 2√2 · (1+z) / (√z · (1-z)²)  ✓

# Y = f, dY/dk = f'/dk' = f₁/dk_dz
Y_k1 = f1 / dk_dz

# For the ODE check, I need dY/dk, d²Y/dk², d³Y/dk³
# dY/dk = (dY/dz)/(dk/dz)
# d²Y/dk² = d/dk[dY/dk] = (1/(dk/dz)) · d/dz[f'/dk']
# This gets complex. Let me use the chain rule properly.

# Let θ = z·d/dz and Θ = k·d/dk
# Y satisfies in k-variable: k(1-k²)Y''' + (1-3k²)Y'' - kY' = 0
# We can also write this as: [(1-k²)Θ(Θ-1) - (Θ+1)k²] ∘ (1/k)Θ Y = 0

# Alternative: evaluate the ODE directly at many z-points
# Use the substitution Y(k(z)) = f(z)
# Chain rule: Y'(k) = f'(z)/k'(z)
# Y''(k) = [f''(z)k'(z) - f'(z)k''(z)] / [k'(z)]³
# Y'''(k) = {...} / [k'(z)]⁵ + ...

# This is messy. Let me instead check numerically:
# Does the generating function of Q_n satisfy the pullback ODE?

# I'll compute Y'', Y''' via the chain rule
# Let me use multiple z-values and finite differences for simplicity

def k_of_z(z):
    return 4*sqrt(2*z)/(1-z)

def dk_of_z(z):
    return 2*sqrt(2)*(1+z)/(sqrt(z)*(1-z)**2)

def d2k_of_z(z):
    # d²k/dz² by finite difference from dk/dz
    h = z * mpf('1e-15')
    return (dk_of_z(z+h) - dk_of_z(z-h)) / (2*h)

# Compute Y and its k-derivatives at several z-values
z_vals = [mpf('0.001'), mpf('0.005'), mpf('0.01')]
print(f"\nODE residual check at several z-values:")

for z0 in z_vals:
    # Compute f(z₀) and derivatives
    f_val = sum(mpf(int(Q[n].numerator))/mpf(int(Q[n].denominator)) * z0**n for n in range(NTERMS))
    f1_val = sum(mpf(int(Q[n].numerator))/mpf(int(Q[n].denominator)) * n * z0**(n-1) for n in range(1, NTERMS))
    f2_val = sum(mpf(int(Q[n].numerator))/mpf(int(Q[n].denominator)) * n*(n-1) * z0**(n-2) for n in range(2, NTERMS))
    f3_val = sum(mpf(int(Q[n].numerator))/mpf(int(Q[n].denominator)) * n*(n-1)*(n-2) * z0**(n-3) for n in range(3, NTERMS))

    k0 = k_of_z(z0)
    dk = dk_of_z(z0)
    d2k = d2k_of_z(z0)

    # Y' = f'/dk
    Y1 = f1_val / dk
    # Y'' = (f'' - Y'·d2k) / dk²
    Y2 = (f2_val - Y1 * d2k) / dk**2
    # Y''' by finite difference of Y''
    h = z0 * mpf('1e-12')
    z_p = z0 + h
    z_m = z0 - h

    f1_p = sum(mpf(int(Q[n].numerator))/mpf(int(Q[n].denominator)) * n * z_p**(n-1) for n in range(1, NTERMS))
    f2_p = sum(mpf(int(Q[n].numerator))/mpf(int(Q[n].denominator)) * n*(n-1) * z_p**(n-2) for n in range(2, NTERMS))
    dk_p = dk_of_z(z_p)
    d2k_p = d2k_of_z(z_p)
    Y1_p = f1_p / dk_p
    Y2_p = (f2_p - Y1_p * d2k_p) / dk_p**2

    f1_m = sum(mpf(int(Q[n].numerator))/mpf(int(Q[n].denominator)) * n * z_m**(n-1) for n in range(1, NTERMS))
    f2_m = sum(mpf(int(Q[n].numerator))/mpf(int(Q[n].denominator)) * n*(n-1) * z_m**(n-2) for n in range(2, NTERMS))
    dk_m = dk_of_z(z_m)
    d2k_m = d2k_of_z(z_m)
    Y1_m = f1_m / dk_m
    Y2_m = (f2_m - Y1_m * d2k_m) / dk_m**2

    Y3 = (Y2_p - Y2_m) / (2*h) / dk

    # ODE: k(1-k²)Y''' + (1-3k²)Y'' - kY' = 0
    residual = k0*(1-k0**2)*Y3 + (1-3*k0**2)*Y2 - k0*Y1
    relative = residual / (abs(k0*Y1) + 1)

    print(f"  z={mp.nstr(z0,4)}: k={mp.nstr(k0,8)}, "
          f"|residual| = {mp.nstr(abs(residual), 8)}, "
          f"|relative| = {mp.nstr(abs(relative), 8)}")

print("\nIf residuals are small, f(z) = Σ Q_n z^n satisfies the pullback of the integrated-K ODE.")
print("This would prove the module identification and hence L = G.")

print("\nDone.")
