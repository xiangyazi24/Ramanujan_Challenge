#!/usr/bin/env python3
"""P2.7: Derive the 4-term recurrence for step-line MOP Q_n(1) values,
then compare the Poincaré polynomial with the P2.7 recurrence.

The step-line MOP polynomials Q_{(n2,n3)} for measures
μ₂ = (-log t)dt, μ₃ = ½(log²t)dt satisfy a 4-term recurrence
on the step line. If this recurrence matches the P2.7 operator
(after gauge/pullback), the proof is complete.
"""
from fractions import Fraction
import sys

def mu2_moment(k):
    return Fraction(1, (k+1)**2)

def mu3_moment(k):
    return Fraction(1, (k+1)**3)

def det_exact(M):
    n = len(M)
    M = [row[:] for row in M]
    det = Fraction(1)
    for col in range(n):
        pivot = None
        for row in range(col, n):
            if M[row][col] != 0:
                pivot = row
                break
        if pivot is None:
            return Fraction(0)
        if pivot != col:
            M[col], M[pivot] = M[pivot], M[col]
            det = -det
        det *= M[col][col]
        inv = Fraction(1, M[col][col])
        for row in range(col+1, n):
            factor = M[row][col] * inv
            for j in range(col+1, n):
                M[row][j] -= factor * M[col][j]
            M[row][col] = Fraction(0)
    return det

def mop_Q_at_1(n2, n3):
    """Compute Q_{(n2,n3)}(1) exactly."""
    N = n2 + n3
    if N == 0:
        return Fraction(1)

    moment_rows = []
    for r in range(n2):
        row = [mu2_moment(r + c) for c in range(N)]
        moment_rows.append(row)
    for r in range(n3):
        row = [mu3_moment(r + c) for c in range(N)]
        moment_rows.append(row)

    Delta = det_exact(moment_rows)
    if Delta == 0:
        return None

    full_rows = []
    for r in range(n2):
        row = [mu2_moment(r + c) for c in range(N + 1)]
        full_rows.append(row)
    for r in range(n3):
        row = [mu3_moment(r + c) for c in range(N + 1)]
        full_rows.append(row)
    full_rows.append([Fraction(1)] * (N + 1))  # t=1: [1,1,...,1]

    full_det = det_exact(full_rows)
    return full_det / Delta

# Step line
def step_line(max_N):
    indices = [(0, 0)]
    for m in range(1, max_N + 1):
        if m % 2 == 1:
            indices.append(((m+1)//2, (m-1)//2))
        else:
            indices.append((m//2, m//2))
    return indices

# Compute Q_N(1) on step line
MAX_STEP = 20
print(f"Computing Q_N(1) on step line for N=0..{MAX_STEP}...", flush=True)
step = step_line(MAX_STEP)
q_vals = []
for i, (n2, n3) in enumerate(step):
    N = n2 + n3
    val = mop_Q_at_1(n2, n3)
    q_vals.append(val)
    if val is not None:
        print(f"  N={N}: ({n2},{n3}), Q(1) = {float(val):.15e}", flush=True)
    else:
        print(f"  N={N}: ({n2},{n3}), SINGULAR", flush=True)
        break

# Derive the 4-term recurrence coefficients
# Ansatz: a_n q[n+3] + b_n q[n+2] + c_n q[n+1] + d_n q[n] = 0
# For each n, solve for (a,b,c,d) up to scaling
print(f"\n=== Deriving 4-term recurrence coefficients ===", flush=True)

def solve_4term(q, n):
    """Given q[n], q[n+1], q[n+2], q[n+3], find (a,b,c,d) s.t. a*q3+b*q2+c*q1+d*q0=0."""
    # Set a=1, solve b*q2 + c*q1 + d*q0 = -q3
    # 3x3 system: no, it's 1 equation with 3 unknowns, underdetermined for single n
    # Need to collect recurrence from MULTIPLE q values
    pass

# For a 4-term recurrence with polynomial coefficients a(n), b(n), c(n), d(n),
# we need to find these polynomials.
# First, let's check the CONSTANT-coefficient case (leading Poincaré behavior).
# The recurrence q[n+3] = α q[n+2] + β q[n+1] + γ q[n]
# gives characteristic equation λ³ - α λ² - β λ - γ = 0

# Compute successive ratios
print(f"\nRatios q[n+1]/q[n]:")
for i in range(1, len(q_vals)):
    if q_vals[i-1] and q_vals[i-1] != 0 and q_vals[i]:
        r = q_vals[i] / q_vals[i-1]
        print(f"  q[{i}]/q[{i-1}] = {float(r):.12f}")

# Estimate the Poincaré polynomial from the asymptotic ratios
# The dominant root should converge
print(f"\n=== Asymptotic Poincaré analysis ===")
print("Ratios should converge to the dominant root of the Poincaré polynomial.")
if len(q_vals) >= 15:
    # Use the last few ratios to estimate dominant root
    ratios = []
    for i in range(max(len(q_vals)-5, 10), len(q_vals)):
        if q_vals[i-1] and q_vals[i-1] != 0 and q_vals[i]:
            r = q_vals[i] / q_vals[i-1]
            ratios.append(float(r))
    if ratios:
        print(f"  Last ratios: {[f'{r:.10f}' for r in ratios]}")
        print(f"  Approximate dominant root: {ratios[-1]:.10f}")

# P2.7 dominant root: one real root of 4μ³-220μ²+8μ-1
# μ₀ ≈ 54.96
# But we're looking at q_n/q_{n-1}, which should approach μ₀/64 ≈ 0.8588
# The MOP ratios converge to ~0.33, which is DIFFERENT.
# So the MOP needs a pullback to match P2.7.

# Let's check: what is the MOP Poincaré polynomial?
# For 4-term recurrence: a*λ³ + b*λ² + c*λ + d = 0
# From the ratios ~ 0.33, the dominant root is about 1/3.

# Check if the MOP system has rational recurrence coefficients
# by verifying a linear dependence among 4 consecutive values
# with polynomial n-dependent coefficients.

# Test: for each n, compute the vector (q[n], q[n+1], q[n+2], q[n+3])
# and check for a 4-term relation.
# With constant coefficients: a q[n+3] + b q[n+2] + c q[n+1] + d q[n] = 0
# Set d=1, solve the 3x3 system from n=0,1,2

if len(q_vals) >= 7:
    print(f"\n=== Testing constant-coefficient 4-term recurrence ===")
    # From n=0: a q3 + b q2 + c q1 + q0 = 0
    # From n=1: a q4 + b q3 + c q2 + q1 = 0
    # From n=2: a q5 + b q4 + c q3 + q2 = 0
    M_sys = []
    rhs = []
    for n in range(3):
        M_sys.append([q_vals[n+3], q_vals[n+2], q_vals[n+1]])
        rhs.append(-q_vals[n])

    # Solve 3x3 system
    # Using Cramer's rule with exact Fraction arithmetic
    det_M = det_exact(M_sys)
    print(f"  System determinant: {det_M}")
    if det_M != 0:
        a = det_exact([[rhs[i]] + M_sys[i][1:] for i in range(3)]) / det_M
        b = det_exact([M_sys[i][:1] + [rhs[i]] + M_sys[i][2:] for i in range(3)]) / det_M
        c = det_exact([M_sys[i][:2] + [rhs[i]] for i in range(3)]) / det_M
        print(f"  a = {a} = {float(a):.10f}")
        print(f"  b = {b} = {float(b):.10f}")
        print(f"  c = {c} = {float(c):.10f}")
        print(f"  d = 1")

        # Verify
        print(f"\n  Verification:")
        max_err = Fraction(0)
        for n in range(min(len(q_vals) - 3, 15)):
            res = a * q_vals[n+3] + b * q_vals[n+2] + c * q_vals[n+1] + q_vals[n]
            if res != 0:
                print(f"    n={n}: residual = {float(res):.6e} (FAIL)")
                max_err = max(max_err, abs(res))
            else:
                if n < 5:
                    print(f"    n={n}: OK")

        if max_err == 0:
            print(f"  Constant-coefficient 4-term recurrence FOUND!")
            print(f"  Characteristic polynomial: {float(a):.10f} λ³ + {float(b):.10f} λ² + {float(c):.10f} λ + 1 = 0")
        else:
            print(f"  Not constant-coefficient (expected for polynomial coefficients).")

    # Now try variable-coefficient: a(n) q[n+3] + b(n) q[n+2] + c(n) q[n+1] + d(n) q[n] = 0
    # where a,b,c,d are polynomials in n.
    # Try degree 0 (constant) first, then degree 1, etc.
    print(f"\n=== Testing variable-coefficient 4-term recurrence ===")
    for deg in range(0, 5):
        nparams = 4 * (deg + 1)
        nequations = len(q_vals) - 3
        if nequations < nparams:
            print(f"  degree {deg}: need {nparams} params, only {nequations} equations. Skip.")
            continue

        # Build the system: for each n=0,...,nequations-1:
        # sum_{j=0}^3 sum_{k=0}^deg  alpha_{j,k} * n^k * q[n+j] = 0
        # Unknowns: alpha_{j,k} for j=0..3, k=0..deg
        # Normalize by setting alpha_{3,deg} = 1 (leading coeff of a(n))
        # Then we have nparams-1 unknowns and nequations equations

        from fractions import Fraction
        nunk = nparams - 1
        neq = min(nequations, nunk + 5)  # use a few extra for verification

        A_mat = []
        b_vec = []
        for n in range(neq):
            row = []
            for j in range(4):
                for k in range(deg + 1):
                    if j == 3 and k == deg:
                        continue  # this is the normalized one
                    row.append(Fraction(n)**k * q_vals[n + j])
            A_mat.append(row)
            b_vec.append(-Fraction(n)**deg * q_vals[n + 3])

        if len(A_mat[0]) == 0:
            continue

        # Solve using least squares (but with exact arithmetic, use first nunk equations)
        A_solve = [row[:nunk] for row in A_mat[:nunk]]
        b_solve = b_vec[:nunk]

        det_sys = det_exact(A_solve)
        if det_sys == 0:
            print(f"  degree {deg}: singular system")
            continue

        # Check residual on extra equations
        # Actually, just solve and verify on all
        # Use Gaussian elimination
        aug = [A_solve[i][:] + [b_solve[i]] for i in range(nunk)]
        n_aug = len(aug)
        for col in range(n_aug):
            pivot = None
            for row in range(col, n_aug):
                if aug[row][col] != 0:
                    pivot = row
                    break
            if pivot is None:
                break
            aug[col], aug[pivot] = aug[pivot], aug[col]
            for row in range(n_aug):
                if row != col and aug[row][col] != 0:
                    factor = aug[row][col] / aug[col][col]
                    for j in range(n_aug + 1):
                        aug[row][j] -= factor * aug[col][j]

        sol = [aug[i][n_aug] / aug[i][i] for i in range(nunk)]

        # Reconstruct coefficients
        coeffs = {}
        idx = 0
        for j in range(4):
            for k in range(deg + 1):
                if j == 3 and k == deg:
                    coeffs[(j, k)] = Fraction(1)
                else:
                    coeffs[(j, k)] = sol[idx]
                    idx += 1

        # Verify
        max_res = Fraction(0)
        ok_count = 0
        for n in range(neq):
            res = Fraction(0)
            for j in range(4):
                poly_val = sum(coeffs[(j, k)] * Fraction(n)**k for k in range(deg + 1))
                res += poly_val * q_vals[n + j]
            if res == 0:
                ok_count += 1
            max_res = max(max_res, abs(res))

        if max_res == 0:
            print(f"  degree {deg}: EXACT recurrence found! ({ok_count}/{neq} verified)")
            for j in range(4):
                poly = [coeffs[(j, k)] for k in range(deg + 1)]
                print(f"    coeff_{j}(n) = {[str(c) for c in poly]}")
            break
        else:
            print(f"  degree {deg}: max residual = {float(max_res):.6e}")

print("\nDone.")
