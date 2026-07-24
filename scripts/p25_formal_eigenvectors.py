#!/usr/bin/env python3
"""P2.5: Compute formal adjoint eigenvectors w± as 1/n series.

From Q4855 §6: expand B̄(n) in powers of 1/n, compute w± to high order,
then take the cross product v = w₊ × w₋ to get the neutral invariant line.
The cross product has coefficients fixed by √2 → -√2, so it's rational.
Then reconstruct v(n) as a rational function of n.
"""
from mpmath import mp, mpf, matrix, sqrt, nstr
import sys

mp.dps = 100

# The balanced matrix B̄(n) has entries that are rational functions of n.
# We need to expand each entry as a Laurent series in 1/n.
# B̄(n) = B₀ + B₁/n + B₂/n² + ...

# First, let me compute B̄(n) symbolically using exact integer arithmetic.
# Then expand each entry as a power series in 1/n.

# M(n) entries are polynomials in n of degree ≤ 7.
# δ(n) is a polynomial of degree 7.
# D_n = diag(1, n+1, (n+1)²)
# B̄(n) = D_n⁻¹ M(n) D_{n+1} / δ(n)

# B̄[i,j](n) = M[i,j](n) * (n+2)^j_scale / ((n+1)^i_scale * δ(n))
# where i_scale = [0,1,2][i], j_scale = [0,1,2][j]

# For the Laurent expansion: we need M[i,j](n) / (n+1)^i_scale * (n+2)^j_scale / δ(n)

# Let me work with mpmath for now and compute B̄(n) for integer n,
# then extract the leading matrix B₀ and correction terms.

sqrt2 = sqrt(2)

def M_poly(n):
    """M(n) entries as exact mpf values."""
    m11 = (-2*n-5)*(n+3)**2*(136*n**4+1424*n**3+5548*n**2+9551*n+6141)
    m12 = 384*n**6+6384*n**5+44168*n**4+162698*n**3+336377*n**2+369933*n+169011
    m13 = -480*n**4-4980*n**3-19210*n**2-32690*n-20730
    m21 = (n+2)**2*(n+3)**2*(4*n+10)*(48*n**3+386*n**2+1017*n+879)
    m22 = (n+2)**2*(-272*n**5-3848*n**4-21732*n**3-61184*n**2-85761*n-47808)
    m23 = (n+2)**2*(320*n**3+2540*n**2+6610*n+5640)
    m31 = (-4*n-10)*(n+2)**2*(n+3)**2*(32*n**4+302*n**3+1037*n**2+1530*n+813)
    m32 = (n+2)**2*(192*n**6+2984*n**5+19116*n**4+64452*n**3+120256*n**2+117279*n+46476)
    m33 = (n+2)**2*(-16*n**5-408*n**4-2912*n**3-8884*n**2-12254*n-6240)
    return [[m11,m12,m13],[m21,m22,m23],[m31,m32,m33]]

def delta_poly(n):
    return -2*(n+2)**2*(n+3)**2*(2*n+5)*(2*n+7)**2

# Compute B̄ at a large n to extract the leading matrix B₀
# B₀ = lim_{n→∞} B̄(n)
# The entries stabilize as n → ∞

N_test = 10000
M_test = M_poly(mpf(N_test))
d_test = delta_poly(mpf(N_test))
n1 = mpf(N_test + 1)
n2 = mpf(N_test + 2)

B0 = matrix(3, 3)
for i in range(3):
    for j in range(3):
        scale_inv = [1, 1/n1, 1/n1**2][i]
        scale_next = [1, n2, n2**2][j]
        B0[i,j] = M_test[i][j] * scale_inv * scale_next / d_test

print("Leading matrix B₀ = lim B̄(n):")
for i in range(3):
    print(f"  [{nstr(B0[i,0], 15)}, {nstr(B0[i,1], 15)}, {nstr(B0[i,2], 15)}]")

# Eigenvalues of B₀
tr = B0[0,0] + B0[1,1] + B0[2,2]
cofsum = (B0[0,0]*B0[1,1] - B0[0,1]*B0[1,0]
        + B0[0,0]*B0[2,2] - B0[0,2]*B0[2,0]
        + B0[1,1]*B0[2,2] - B0[1,2]*B0[2,1])
det = (B0[0,0]*(B0[1,1]*B0[2,2]-B0[1,2]*B0[2,1])
     - B0[0,1]*(B0[1,0]*B0[2,2]-B0[1,2]*B0[2,0])
     + B0[0,2]*(B0[1,0]*B0[2,1]-B0[1,1]*B0[2,0]))

from mpmath import polyroots
roots = polyroots([1, -tr, cofsum, -det])
print(f"\nEigenvalues of B₀:")
for i, r in enumerate(roots):
    print(f"  λ_{i} = {nstr(r, 20)}")
    # Check: are these 1, 17+12√2, 17-12√2?
    for val, name in [(1, "1"), (17+12*sqrt2, "17+12√2"), (17-12*sqrt2, "17-12√2")]:
        if abs(r - val) < 1e-10:
            print(f"       = {name}")

# Right eigenvectors of B₀
print(f"\nRight eigenvectors of B₀:")
for lam in roots:
    C = matrix(3, 3)
    for i in range(3):
        for j in range(3):
            C[i,j] = B0[i,j] - (lam if i==j else 0)
    row0 = [C[0,0], C[0,1], C[0,2]]
    row1 = [C[1,0], C[1,1], C[1,2]]
    v1 = row0[1]*row1[2] - row0[2]*row1[1]
    v2 = row0[2]*row1[0] - row0[0]*row1[2]
    v3 = row0[0]*row1[1] - row0[1]*row1[0]
    norm = abs(v1) + abs(v2) + abs(v3)
    if norm > 0:
        v1, v2, v3 = v1/v3, v2/v3, mpf(1)
    print(f"  λ={nstr(lam,8)}: v = ({nstr(v1,15)}, {nstr(v2,15)}, {nstr(v3,15)})")

# Left eigenvectors of B₀
print(f"\nLeft eigenvectors of B₀ᵀ:")
B0T = B0.T
for lam in roots:
    C = matrix(3, 3)
    for i in range(3):
        for j in range(3):
            C[i,j] = B0T[i,j] - (lam if i==j else 0)
    row0 = [C[0,0], C[0,1], C[0,2]]
    row1 = [C[1,0], C[1,1], C[1,2]]
    v1 = row0[1]*row1[2] - row0[2]*row1[1]
    v2 = row0[2]*row1[0] - row0[0]*row1[2]
    v3 = row0[0]*row1[1] - row0[1]*row1[0]
    norm = abs(v1) + abs(v2) + abs(v3)
    if norm > 0:
        v1, v2, v3 = v1/v3, v2/v3, mpf(1)
    print(f"  λ={nstr(lam,8)}: w = ({nstr(v1,15)}, {nstr(v2,15)}, {nstr(v3,15)})")

# Now the key computation: formal expansion of B̄(n) in 1/n
# B̄(n) = B₀ + B₁/n + B₂/n² + ...
# Extract B_k by computing B̄(n)·n^k at several large n and interpolating

print("\n" + "="*80)
print("Computing Laurent expansion coefficients B_k")
print("="*80)

# Method: compute B̄(n) at n = N, N+1, ..., N+K for large N,
# then use the Vandermonde system to extract B₀, B₁, ...

K_ORDER = 20  # Number of 1/n terms to extract

def balanced_entry(n, i, j):
    M = M_poly(mpf(n))
    d = delta_poly(mpf(n))
    n1 = mpf(n+1)
    n2 = mpf(n+2)
    scale_inv = [1, 1/n1, 1/n1**2][i]
    scale_next = [1, n2, n2**2][j]
    return M[i][j] * scale_inv * scale_next / d

# Alternative: direct Laurent extraction
# B̄[i,j](n) is a ratio of polynomials in n.
# Let me compute it for specific (i,j) and expand.
# The numerator of B̄[i,j] is M[i,j](n) * (n+2)^(j_scale)
# The denominator is δ(n) * (n+1)^(i_scale)
# where i_scale = [0,1,2][i], j_scale = [0,1,2][j]

# To expand P(n)/Q(n) in 1/n, compute P and Q as polynomials,
# then do polynomial long division.

# For now, let me just evaluate B̄ at many large n and do Richardson extrapolation.
# Use n = 1000, 2000, ..., 1000*K to get K terms.

BASE = 2000
NPTS = K_ORDER + 5

B_data = []
for idx in range(NPTS):
    n = BASE + idx * 100
    B = matrix(3, 3)
    for i in range(3):
        for j in range(3):
            B[i,j] = balanced_entry(n, i, j)
    B_data.append((n, B))

# For each entry (i,j), fit B̄[i,j](n) = Σ_{k=0}^{K} c_k / n^k
# using the data points. This is a Vandermonde system in 1/n.

B_coeffs = [[None]*3 for _ in range(3)]
for i in range(3):
    for j in range(3):
        # Fit c_0 + c_1/n + ... + c_K/n^K to the data
        A_mat = matrix(K_ORDER+1, K_ORDER+1)
        b_vec = matrix(K_ORDER+1, 1)
        for idx in range(K_ORDER+1):
            n = B_data[idx][0]
            inv_n = mpf(1) / n
            for k in range(K_ORDER+1):
                A_mat[idx, k] = inv_n**k
            b_vec[idx, 0] = B_data[idx][1][i,j]
        try:
            x = mp.lu_solve(A_mat, b_vec)
            B_coeffs[i][j] = [x[k] for k in range(K_ORDER+1)]
        except:
            B_coeffs[i][j] = None

# Print B₀ (from expansion)
print("\nB₀ (from Laurent expansion):")
for i in range(3):
    row = [nstr(B_coeffs[i][j][0], 15) if B_coeffs[i][j] else "?" for j in range(3)]
    print(f"  [{', '.join(row)}]")

# Check: B₀ should have eigenvalues 1, 17±12√2
print("\nB₁ (correction at order 1/n):")
for i in range(3):
    row = [nstr(B_coeffs[i][j][1], 15) if B_coeffs[i][j] else "?" for j in range(3)]
    print(f"  [{', '.join(row)}]")

# The right eigenvectors of B₀ at the two outer eigenvalues are w₊ and w₋
# (after conjugation by √2 → -√2).
# To find the neutral line, compute the projector onto the neutral eigenspace.

# The leading neutral right eigenvector is e₃ = (0, 0, 1).
# The correction to the neutral eigenvector at order 1/n:
# (B₀ - I) v₁ = -B₁ v₀ (projected onto complement of v₀)
# where v₀ = e₃, so B₁ v₀ = column 2 of B₁

# Actually, let me use the recursive procedure directly.
# v(n) = v₀ + v₁/n + v₂/n² + ...
# r₀(n) = r₀₀ + r₀₁/n + r₀₂/n² + ...
# B̄(n) v(n) = r₀(n) v(n+1)
# where v(n+1) = v₀ + v₁/(n+1) + ... and 1/(n+1) = 1/n - 1/n² + 1/n³ - ...

# Leading: B₀ v₀ = r₀₀ v₀ → r₀₀ = 1, v₀ = e₃ ✓

# Order 1/n: B₁ v₀ + B₀ v₁ = r₀₁ v₀ + r₀₀ (v₁ - v₁)
# Wait, v(n+1) = v₀ + v₁/(n+1) + v₂/(n+1)² + ...
# 1/(n+1) = (1/n)(1/(1+1/n)) = (1/n)(1 - 1/n + 1/n² - ...)
# So v₁/(n+1) = v₁/n - v₁/n² + v₁/n³ - ...
# v₂/(n+1)² = v₂/n² - 2v₂/n³ + ...

# B̄(n) v(n) = r₀(n) v(n+1)
# (B₀ + B₁/n + B₂/n² + ...)(v₀ + v₁/n + v₂/n² + ...)
# = (r₀₀ + r₀₁/n + r₀₂/n² + ...)(v₀ + v₁/n(1 - 1/n + ...) + v₂/n²(1 - 2/n + ...) + ...)

# Order 0: B₀ v₀ = r₀₀ v₀  →  r₀₀ = 1 (neutral eigenvalue)

# Order 1/n: B₁ v₀ + B₀ v₁ = r₀₁ v₀ + r₀₀ v₁ = r₀₁ v₀ + v₁
# → (B₀ - I) v₁ = r₀₁ v₀ - B₁ v₀
# → (B₀ - I) v₁ = (r₀₁ I - B₁) v₀

# Left eigenvector l₀ of B₀ for eigenvalue 1: l₀ (B₀ - I) = 0
# Solvability: l₀ (r₀₁ I - B₁) v₀ = 0
# → r₀₁ l₀·v₀ = l₀ B₁ v₀
# → r₀₁ = l₀ B₁ v₀ / (l₀ · v₀)

# Let me compute this.
print("\n" + "="*80)
print("Formal expansion of neutral eigenvector")
print("="*80)

# B₀, B₁, B₂, ... as 3×3 matrices
def B_matrix(k):
    B = matrix(3, 3)
    for i in range(3):
        for j in range(3):
            if B_coeffs[i][j] is not None:
                B[i,j] = B_coeffs[i][j][k]
    return B

B0_mat = B_matrix(0)
B1_mat = B_matrix(1)
I3 = matrix([[1,0,0],[0,1,0],[0,0,1]])

# Verify B₀ eigenvalue 1 for e₃
Be3 = matrix([[B0_mat[i,2]] for i in range(3)])
print(f"B₀ · e₃ = ({nstr(Be3[0,0],15)}, {nstr(Be3[1,0],15)}, {nstr(Be3[2,0],15)})")
print(f"Expected (0, 0, 1)")

# Left eigenvector for eigenvalue 1
# l₀ · B₀ = l₀  →  l₀ (B₀ - I) = 0
# Compute from B₀ᵀ
B0mI = B0_mat - I3
# Null vector of (B₀-I)ᵀ
C = B0mI.T
row0 = [C[0,0], C[0,1], C[0,2]]
row1 = [C[1,0], C[1,1], C[1,2]]
l1 = row0[1]*row1[2] - row0[2]*row1[1]
l2 = row0[2]*row1[0] - row0[0]*row1[2]
l3 = row0[0]*row1[1] - row0[1]*row1[0]
# Normalize so l·e₃ = 1 (i.e., l3 = 1)
if abs(l3) > 1e-50:
    l1, l2, l3 = l1/l3, l2/l3, mpf(1)
l0_vec = [l1, l2, l3]
print(f"\nLeft eigenvector l₀ for λ=1: ({nstr(l1,15)}, {nstr(l2,15)}, {nstr(l3,15)})")
print(f"l₀ · e₃ = {nstr(l3,15)} (should be 1)")

# r₀₁ = l₀ B₁ e₃ / (l₀ · e₃)
B1e3 = [sum(B1_mat[i,j] * (1 if j==2 else 0) for j in range(3)) for i in range(3)]
l0_B1_e3 = sum(l0_vec[i] * B1e3[i] for i in range(3))
l0_e3 = l3
r01 = l0_B1_e3 / l0_e3
print(f"\nr₀₁ = l₀ · B₁ · e₃ / (l₀ · e₃) = {nstr(r01, 20)}")
print(f"Expected: ≈ -3 (formal index)")

# Solve for v₁ from (B₀ - I) v₁ = (r₀₁ I - B₁) e₃
rhs = matrix(3, 1)
for i in range(3):
    rhs[i,0] = (r01 if i==2 else 0) - B1_mat[i,2]

print(f"\nRHS = (r₀₁ I - B₁) e₃ = ({nstr(rhs[0,0],15)}, {nstr(rhs[1,0],15)}, {nstr(rhs[2,0],15)})")

# (B₀ - I) v₁ = RHS
# B₀ - I is singular (rank 2), so we need the pseudoinverse
# Constrain: l₀ · v₁ = 0 (normalization: keep v₃ = 1 at leading order)
# Actually, we can set v₁₃ = 0 (don't perturb the leading component)

# Augmented system: (B₀ - I | l₀ᵀ) (v₁ | 0)ᵀ = (RHS | 0)
# Wait, let me just replace the last row with l₀ and set last RHS to 0
A_aug = matrix(3, 3)
for i in range(2):
    for j in range(3):
        A_aug[i,j] = B0mI[i,j]
for j in range(3):
    A_aug[2,j] = l0_vec[j]
b_aug = matrix(3, 1)
b_aug[0,0] = rhs[0,0]
b_aug[1,0] = rhs[1,0]
b_aug[2,0] = mpf(0)

v1 = mp.lu_solve(A_aug, b_aug)
print(f"\nv₁ = ({nstr(v1[0],20)}, {nstr(v1[1],20)}, {nstr(v1[2],20)})")

# Verify: (B₀-I) v₁ should equal RHS
check = B0mI * v1
print(f"Check (B₀-I)v₁ = ({nstr(check[0,0],15)}, {nstr(check[1,0],15)}, {nstr(check[2,0],15)})")
print(f"RHS            = ({nstr(rhs[0,0],15)}, {nstr(rhs[1,0],15)}, {nstr(rhs[2,0],15)})")

# Now iterate to higher orders
# General recursion at order k:
# Σ_{j=0}^{k} B_j v_{k-j} = Σ_{j=0}^{k} r₀_j [shift-corrected v_{k-j}]
# The shift correction: v_m/(n+1)^m expands as v_m/n^m * (1 - m/n + ...)
# More precisely, if v(n+1) = Σ v_k/(n+1)^k, then
# v_k/(n+1)^k = v_k · n^{-k} · (1+1/n)^{-k}
# = v_k · n^{-k} · Σ_{j≥0} (-k choose j) n^{-j}
# = v_k · Σ_{j≥0} C(-k,j) n^{-(k+j)}

# So at order 1/n^p, the RHS contributes:
# Σ_{k+j=p} r₀_k · C(-m, j) · v_m  where m = p-k-j... this is getting complicated.
# Let me use the direct approach: compute v(n) numerically and do Padé.

print("\n" + "="*80)
print("Numerical neutral eigenvector for Padé reconstruction")
print("="*80)

from mpmath import polyroots as pr

def neutral_eigen(n):
    """Returns (λ_neutral, v1/v3, v2/v3) for the balanced matrix at n."""
    n = mpf(n)
    M = M_poly(n)
    d = delta_poly(n)
    n1 = n + 1
    n2 = n + 2
    B = matrix(3, 3)
    for i in range(3):
        for j in range(3):
            scale_inv = [1, 1/n1, 1/n1**2][i]
            scale_next = [1, n2, n2**2][j]
            B[i,j] = M[i][j] * scale_inv * scale_next / d

    tr = B[0,0] + B[1,1] + B[2,2]
    cof = (B[0,0]*B[1,1] - B[0,1]*B[1,0]
         + B[0,0]*B[2,2] - B[0,2]*B[2,0]
         + B[1,1]*B[2,2] - B[1,2]*B[2,1])
    det = (B[0,0]*(B[1,1]*B[2,2]-B[1,2]*B[2,1])
         - B[0,1]*(B[1,0]*B[2,2]-B[1,2]*B[2,0])
         + B[0,2]*(B[1,0]*B[2,1]-B[1,1]*B[2,0]))
    roots = pr([1, -tr, cof, -det])
    idx = min(range(3), key=lambda i: abs(roots[i] - 1))
    lam = roots[idx]

    C = matrix(3, 3)
    for i in range(3):
        for j in range(3):
            C[i,j] = B[i,j] - (lam if i==j else 0)
    r0 = [C[0,0], C[0,1], C[0,2]]
    r1 = [C[1,0], C[1,1], C[1,2]]
    v1 = r0[1]*r1[2] - r0[2]*r1[1]
    v2 = r0[2]*r1[0] - r0[0]*r1[2]
    v3 = r0[0]*r1[1] - r0[1]*r1[0]
    if abs(v3) > 1e-50:
        return lam, v1/v3, v2/v3
    return lam, v1, v2

# Compute v1(n)·n and v2(n)·n for large n to see if they stabilize
print(f"\n{'n':>4} {'v1·n':>20} {'v2·n':>20} {'v1·n²':>20} {'v2·n²':>20}")
for n in [10, 20, 30, 40, 50, 60, 80, 100, 150, 200, 300, 500]:
    lam, v1, v2 = neutral_eigen(n)
    print(f"{n:4d} {nstr(v1*n,15):>20} {nstr(v2*n,15):>20} {nstr(v1*n*n,15):>20} {nstr(v2*n*n,15):>20}")
