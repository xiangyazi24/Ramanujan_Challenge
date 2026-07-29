#!/usr/bin/env python3
"""P2.7: Compute step-line type-II multiple-orthogonal polynomials
for the pair (μ₂, μ₃) where μ₂ = (-log t)dt, μ₃ = ½(log²t)dt.

Then compare Q_N(1) values and recurrence with P2.7.
"""
from fractions import Fraction
import sys

# Moments: ∫₀¹ t^k (-log t) dt = 1/(k+1)²
#           ∫₀¹ t^k ½(log²t) dt = 1/(k+1)³
def mu2_moment(k):
    return Fraction(1, (k+1)**2)

def mu3_moment(k):
    return Fraction(1, (k+1)**3)

def det_exact(M):
    """Compute determinant of a matrix of Fractions using Gaussian elimination."""
    n = len(M)
    M = [row[:] for row in M]  # copy
    det = Fraction(1)
    for col in range(n):
        # Find pivot
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

def mop_Q_at_point(n2, n3, t_val):
    """Evaluate Q_(n2,n3)(t_val) using the determinant formula."""
    N = n2 + n3
    if N == 0:
        return Fraction(1)

    # Build moment matrix (N rows, N columns)
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

    # Build (N+1) x (N+1) matrix with last row = [1, t, t², ..., t^N]
    full_rows = []
    for r in range(n2):
        row = [mu2_moment(r + c) for c in range(N + 1)]
        full_rows.append(row)
    for r in range(n3):
        row = [mu3_moment(r + c) for c in range(N + 1)]
        full_rows.append(row)

    # Last row: [t^0, t^1, ..., t^N]
    t_pow = Fraction(1)
    last_row = []
    for c in range(N + 1):
        last_row.append(t_pow)
        t_pow *= t_val
    full_rows.append(last_row)

    full_det = det_exact(full_rows)
    return full_det / Delta

# Step line: (0,0), (1,0), (1,1), (2,1), (2,2), (3,2), (3,3), ...
def step_line_indices(max_steps):
    indices = [(0, 0)]
    n = 0
    while len(indices) <= max_steps:
        n += 1
        if n % 2 == 1:
            indices.append(((n+1)//2, (n-1)//2))
        else:
            indices.append((n//2, n//2))
    return indices

# Compute Q_N(1) for step-line indices
print("=== Step-line MOP polynomials: Q(1) values ===", flush=True)
step = step_line_indices(15)

q_at_1 = []
for i, (n2, n3) in enumerate(step):
    N = n2 + n3
    val = mop_Q_at_point(n2, n3, Fraction(1))
    q_at_1.append(val)
    print(f"  Step {i}: ({n2},{n3}), N={N}, Q(1) = {val}", flush=True)
    if val is None:
        break

# P2.7 initial values
q_p27 = [
    Fraction(-215040420000),
    Fraction(-167282265043404, 905),
    Fraction(-964185327658080, 6071),
]

print(f"\n=== P2.7 q_n values ===")
for n, q in enumerate(q_p27):
    print(f"  q_{n} = {q} = {float(q):.6e}")

# Check ratios q_p27[n] / Q_step[n](1)
print(f"\n=== Ratios q_n(P2.7) / Q_step(1) ===")
for n in range(min(len(q_p27), len(q_at_1))):
    if q_at_1[n] and q_at_1[n] != 0:
        ratio = q_p27[n] / q_at_1[n]
        print(f"  n={n}: q_{n}/Q_{n}(1) = {ratio} = {float(ratio):.6e}")

# Compute the 4-term recurrence on the step line
# If Q_N(1) satisfies a_n Q_{N}(1) + b_n Q_{N-1}(1) + c_n Q_{N-2}(1) + d_n Q_{N-3}(1) = 0
# Find a_n, b_n, c_n by solving:
# a_n Q_N(1) + b_n Q_{N-1}(1) + c_n Q_{N-2}(1) + d_n Q_{N-3}(1) = 0
print(f"\n=== Testing 4-term recurrence on Q_step(1) ===", flush=True)
if len(q_at_1) >= 7:
    for n in range(3, min(len(q_at_1), 12)):
        v = [q_at_1[n-3], q_at_1[n-2], q_at_1[n-1], q_at_1[n]]
        if all(x is not None for x in v):
            # Check if they satisfy a linear relation
            # We can compute the "recurrence ratio" as -v[3]*(v[1]*v[0]) / ...
            # Actually let's just print the ratios v[n]/v[n-1]
            if v[2] != 0:
                r1 = v[3] / v[2]
                print(f"  N={n}: Q_{n}(1)/Q_{n-1}(1) = {float(r1):.10f}")

# Also compute mixed moment ∫₀¹ ∫₀¹ Q_N(xy) (1-½log(xy))/(1-xy) dx dy
# = Q(1)L - (P_{2,Q}(1) + P_{3,Q}(1))
# For this we need the full polynomial, not just Q(1).
# Let me compute Q as a polynomial (list of coefficients)
print(f"\n=== Computing MOP polynomial coefficients ===", flush=True)

def mop_Q_coeffs(n2, n3):
    """Return coefficients [c_0, c_1, ..., c_N] of Q_(n2,n3)(t)."""
    N = n2 + n3
    if N == 0:
        return [Fraction(1)]

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

    # Cofactor expansion along the last row
    coeffs = []
    for col in range(N + 1):
        # Build minor matrix: all rows, all columns except col
        minor_rows = []
        for r in range(n2):
            row = [mu2_moment(r + c) for c in range(N + 1) if c != col]
            minor_rows.append(row)
        for r in range(n3):
            row = [mu3_moment(r + c) for c in range(N + 1) if c != col]
            minor_rows.append(row)
        sign = Fraction((-1)**(N + col))
        minor_det = det_exact(minor_rows)
        coeffs.append(sign * minor_det / Delta)

    return coeffs

# Compute mixed moment pair for a polynomial R(t) = sum c_k t^k
# W[R] = R(1)*L - sum c_k * h_k where h_k = H_k^(2) + H_k^(3)
def mixed_moment_pair(coeffs):
    """Return (q, p) such that W[R] = q*L - p."""
    q = sum(coeffs)
    h_vals = []
    for k in range(len(coeffs)):
        hk = sum(Fraction(1, j**2) + Fraction(1, j**3) for j in range(1, k+1))
        h_vals.append(hk)
    p = sum(coeffs[k] * h_vals[k] for k in range(len(coeffs)))
    return q, p

print("Step-line MOP: Q(1) and mixed moment W[Q]:")
for i in range(min(8, len(step))):
    n2, n3 = step[i]
    N = n2 + n3
    coeffs = mop_Q_coeffs(n2, n3)
    if coeffs is None:
        print(f"  Step {i}: ({n2},{n3}) — singular")
        continue
    q, p = mixed_moment_pair(coeffs)
    q_1 = sum(coeffs)
    print(f"  Step {i}: ({n2},{n3}), N={N}")
    print(f"    Q(1) = {q_1}")
    print(f"    W[Q] = {q_1}·L - {p}")
    print(f"    coeffs = {[str(c) for c in coeffs[:5]]}{'...' if len(coeffs)>5 else ''}")

print("\nDone.")
