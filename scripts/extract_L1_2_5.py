#!/usr/bin/env python3
"""Problem 2.5: Extract the order-2 recurrence L₁ from the difference sequence v_k = q_{k+1} - q_k.

If L = L₁·(S-1) and q_k satisfies L, then v_k = q_{k+1} - q_k satisfies L₁.
L₁ is order 2 with polynomial coefficients. We find these coefficients by solving
a linear system over Q.

Expected: degree pattern dropping by 7, so L₁ has coefficients of degree ≤ 21, 14, 7.
"""
from mpmath import mp, mpf, nstr, matrix, catalan
from fractions import Fraction

mp.dps = 500  # very high precision for rational reconstruction

def M(n):
    """3x3 matrix M(n) with exact rational entries."""
    n = mpf(n)
    m11 = (-2*n-5)*(n+3)**2 * (136*n**4 + 1424*n**3 + 5548*n**2 + 9551*n + 6141)
    m12 = 384*n**6 + 6384*n**5 + 44168*n**4 + 162698*n**3 + 336377*n**2 + 369933*n + 169011
    m13 = -480*n**4 - 4980*n**3 - 19210*n**2 - 32690*n - 20730

    m21 = (n+2)**2*(n+3)**2*(4*n+10)*(48*n**3 + 386*n**2 + 1017*n + 879)
    m22 = (n+2)**2*(-272*n**5 - 3848*n**4 - 21732*n**3 - 61184*n**2 - 85761*n - 47808)
    m23 = (n+2)**2*(320*n**3 + 2540*n**2 + 6610*n + 5640)

    m31 = (-4*n-10)*(n+2)**2*(n+3)**2*(32*n**4 + 302*n**3 + 1037*n**2 + 1530*n + 813)
    m32 = (n+2)**2*(192*n**6 + 2984*n**5 + 19116*n**4 + 64452*n**3 + 120256*n**2 + 117279*n + 46476)
    m33 = (n+2)**2*(-16*n**5 - 408*n**4 - 2912*n**3 - 8884*n**2 - 12254*n - 6240)

    return matrix([[m11, m12, m13],
                   [m21, m22, m23],
                   [m31, m32, m33]])

A_mat = matrix([[mpf(30921), mpf(-32972), mpf(8240)],
                [mpf(33750), mpf(-36000), mpf(9000)]])

# Compute q_N = Q_{N,1} for N = 0, ..., N_max
N_max = 100
print(f"Computing q_N for N=0..{N_max} at {mp.dps} digit precision...")

q_vals = []
T = matrix([[1,0,0],[0,1,0],[0,0,1]])
for N in range(N_max + 1):
    AT = A_mat * T
    q_val = AT[1, 0]
    q_vals.append(q_val)
    T = T * M(N)

# Compute difference sequence v_k = q_{k+1} - q_k
v_vals = [q_vals[k+1] - q_vals[k] for k in range(N_max)]
print(f"Computed {len(v_vals)} values of v_k")

# Now find the order-2 recurrence for v_k:
# α₂(k) v_{k+2} + α₁(k) v_{k+1} + α₀(k) v_k = 0
# where α_i(k) are polynomials in k.

# Strategy: if deg(α₂) = d₂, deg(α₁) = d₁, deg(α₀) = d₀,
# total unknowns = (d₂+1) + (d₁+1) + (d₀+1).
# Each k value gives one linear equation.
# We need total unknowns + 1 equations (to find up to scalar).

# Expected degree pattern: (21, 14, 7). Total unknowns = 22+15+8 = 45.
# But that's a LOT. Let me start smaller and work up.

# First: try SMALLER degree patterns and see which one works.
# Poincaré normalization: for large k, the leading terms of α_i determine the
# Poincaré characteristic equation. The roots should be 17±12√2.

# For order 2 with roots 17±12√2 (product = 1, sum = 34):
# The Poincaré equation is t² - 34t + 1 = 0.
# This means: leading(α₂)/leading(α₂) * t² + leading(α₁)/leading(α₂) * t + leading(α₀)/leading(α₂) → t² - 34t + 1
# In other words, for the normalized recurrence:
# v_{k+2} + (α₁/α₂)v_{k+1} + (α₀/α₂)v_k = 0
# As k→∞: α₁/α₂ → -34 · k^{d₁-d₂} and α₀/α₂ → 1 · k^{d₀-d₂}

# Check from the data: v_{k+2}/(v_k) for large k should approach prod of roots = 1 after normalization.
print("\n=== Poincaré analysis of v_k ===")
for k in [20, 30, 40, 50, 60, 70]:
    if k+2 < len(v_vals) and v_vals[k] != 0:
        r = v_vals[k+2] / v_vals[k]
        print(f"  v[{k+2}]/v[{k}] = {nstr(r, 15)}")
        # This should behave like k^{2*(d₂-d₀)} * (17+12√2)² (dominant squared)
        # or rather, v_{k+2} ≈ α₁/α₂ · v_{k+1} (dominant mode)

# Better: consecutive ratio
print("\nConsecutive ratio v[k+1]/v[k]:")
for k in [30, 40, 50, 60, 70, 80]:
    if k+1 < len(v_vals) and v_vals[k] != 0:
        r = v_vals[k+1] / v_vals[k]
        # After dividing by k^(d₂-d₀), this should approach the dominant Poincaré root
        print(f"  v[{k+1}]/v[{k}] = {nstr(r, 15)}, /k^7 = {nstr(r/k**7, 10)}")

# So we need d₂ - d₀ = 7 to absorb the k⁷ growth.
# Let's try several degree patterns:
print("\n=== Trying degree patterns ===")

def try_recurrence(d2, d1, d0, num_data=None):
    """Try to find order-2 recurrence with given degree pattern.
    Returns True if found, prints coefficients."""
    total_unknowns = (d2+1) + (d1+1) + (d0+1)
    if num_data is None:
        num_data = total_unknowns + 5  # some extra for overdetermined check

    # Build linear system: for each k, the equation
    # Σ_{j=0}^{d2} c2_j k^j · v_{k+2} + Σ_{j=0}^{d1} c1_j k^j · v_{k+1} + Σ_{j=0}^{d0} c0_j k^j · v_k = 0
    # Unknowns: c2_0,...,c2_{d2}, c1_0,...,c1_{d1}, c0_0,...,c0_{d0}

    # Use k = 2, 3, ..., num_data+1 (avoid k=0,1 which might be boundary)
    k_start = 5  # start away from boundary
    mat_rows = []
    for k_idx in range(num_data):
        k = k_start + k_idx
        if k+2 >= len(v_vals):
            break
        row = []
        for j in range(d2+1):
            row.append(mpf(k)**j * v_vals[k+2])
        for j in range(d1+1):
            row.append(mpf(k)**j * v_vals[k+1])
        for j in range(d0+1):
            row.append(mpf(k)**j * v_vals[k])
        mat_rows.append(row)

    n_rows = len(mat_rows)
    n_cols = total_unknowns

    if n_rows < n_cols:
        print(f"  ({d2},{d1},{d0}): Not enough data ({n_rows} < {n_cols})")
        return False

    # Solve homogeneous system: find null space
    # Method: use SVD-like approach via QR, or just solve Ax=0 by fixing one variable
    # Fix the last coefficient = 1 and solve the (n_cols-1) × (n_cols-1) system.

    # Actually, let's use a different approach: compute the right nullvector
    # by solving the rectangular system.

    # Form the matrix
    M_sys = matrix(n_rows, n_cols)
    for i in range(n_rows):
        for j in range(n_cols):
            M_sys[i, j] = mat_rows[i][j]

    # Try: fix c0_{d0} = 1, solve for the rest
    # Move the last column to the RHS
    rhs = matrix(n_rows, 1)
    M_red = matrix(n_rows, n_cols - 1)
    for i in range(n_rows):
        rhs[i, 0] = -M_sys[i, n_cols-1]
        for j in range(n_cols - 1):
            M_red[i, j] = M_sys[i, j]

    try:
        from mpmath import lu_solve as mplu
        # Use least squares for overdetermined system
        # M_red^T M_red x = M_red^T rhs
        MTM = M_red.T * M_red
        MTr = M_red.T * rhs
        x = mp.lu_solve(MTM, MTr)

        # Reconstruct full solution
        coeffs = [x[i,0] for i in range(n_cols-1)] + [mpf(1)]

        # Check residuals
        residuals = []
        for k_idx in range(n_rows):
            k = k_start + k_idx
            res = sum(coeffs[j] * mat_rows[k_idx][j] for j in range(n_cols))
            residuals.append(abs(res))

        max_res = max(residuals)
        avg_res = sum(residuals) / len(residuals)

        if max_res < mpf(10)**(-mp.dps//2):
            print(f"  ({d2},{d1},{d0}): FOUND! max_residual = {float(max_res):.3e}")

            # Extract polynomial coefficients
            c2_coeffs = coeffs[:d2+1]
            c1_coeffs = coeffs[d2+1:d2+1+d1+1]
            c0_coeffs = coeffs[d2+1+d1+1:]

            print(f"    α₂(k) coeffs (const..k^{d2}): {[nstr(c,10) for c in c2_coeffs]}")
            print(f"    α₁(k) coeffs (const..k^{d1}): {[nstr(c,10) for c in c1_coeffs]}")
            print(f"    α₀(k) coeffs (const..k^{d0}): {[nstr(c,10) for c in c0_coeffs]}")

            # Verify on held-out data
            print("    Held-out verification:")
            for k in range(70, min(80, len(v_vals)-2)):
                alpha2 = sum(c2_coeffs[j] * k**j for j in range(d2+1))
                alpha1 = sum(c1_coeffs[j] * k**j for j in range(d1+1))
                alpha0 = sum(c0_coeffs[j] * k**j for j in range(d0+1))
                res = alpha2*v_vals[k+2] + alpha1*v_vals[k+1] + alpha0*v_vals[k]
                print(f"      k={k}: residual = {float(abs(res)):.3e}")

            return True
        else:
            print(f"  ({d2},{d1},{d0}): Failed, max_residual = {float(max_res):.3e}")
            return False
    except Exception as e:
        print(f"  ({d2},{d1},{d0}): Solver failed: {e}")
        return False

# Try degree patterns
for d2, d1, d0 in [(7, 7, 7), (14, 7, 7), (7, 14, 7), (21, 14, 7),
                    (14, 14, 7), (14, 14, 14), (21, 21, 14), (21, 14, 14),
                    (28, 21, 14)]:
    if try_recurrence(d2, d1, d0):
        break
