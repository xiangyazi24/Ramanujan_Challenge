#!/usr/bin/env python3
"""Problem 2.5: Sym²(Delannoy) analysis.

Computes:
1. The invariant conic J₀ for the leading matrix C
2. Numerical verification of the Ore intertwiner
3. The gauge between challenge and Sym²(Delannoy) sequences
"""
from mpmath import mp, mpf, matrix, catalan, nstr, sqrt, log, pi, fac
import numpy as np
from fractions import Fraction as F

mp.dps = 100

# ---- Leading matrix C from Q4772 (eq 4.3) ----
C = np.array([
    [17, -24, 0],
    [-12, 17, 0],
    [8, -12, 1]
], dtype=float)

print("=== Leading matrix C ===")
eigvals = np.linalg.eigvals(C)
print(f"Eigenvalues: {eigvals}")
# Should be 1, 17+12√2, 17-12√2

# ---- Find J₀ satisfying C^T J₀ C = J₀ ----
# J₀ is a symmetric 3×3 matrix. The equation C^T J₀ C = J₀ is linear in the
# 6 independent entries of J₀. We solve it.

# Let J = [[a,b,c],[b,d,e],[c,e,f]]
# C^T J C = J gives 6 equations for 6 unknowns
# This is equivalent to: (C⊗C)^T vec(J) = vec(J)
# where ⊗ is the Kronecker product

# More precisely, for symmetric J, we use the vectorization of the upper triangle.
# But let's just solve it directly with the full 9x9 system.

CkC = np.kron(C.T, C.T)  # (C⊗C)^T acts on vec(J)
# The equation is C^T J C = J, i.e., vec(C^T J C) = vec(J)
# vec(C^T J C) = (C^T ⊗ C^T) vec(J)
# So we need the eigenspace of (C^T ⊗ C^T) with eigenvalue 1

I9 = np.eye(9)
M_sys = CkC - I9
# Null space of M_sys
U, S, Vt = np.linalg.svd(M_sys)
null_mask = S < 1e-10
null_vecs = Vt[null_mask]  # rows of Vt corresponding to zero singular values

print(f"\n=== Null space dimension of C^T⊗C^T - I: {null_mask.sum()} ===")
for i, v in enumerate(null_vecs):
    J = v.reshape(3,3)
    Jsym = (J + J.T) / 2
    print(f"\nNull vector {i}: {v}")
    print(f"As symmetric matrix:\n{Jsym}")
    print(f"det = {np.linalg.det(Jsym):.6f}")

# ---- Sym²(Delannoy) recurrence ----
# From Q4772 eq 2.3:
# (n+3)²(2n+3) U_{n+3} - (2n+5)(35n²+140n+131) U_{n+2}
# + (2n+3)(35n²+140n+131) U_{n+1} - (2n+5)(n+1)² U_n = 0

def sym2_del_coeffs(n):
    """Return (c0, c1, c2, c3) for U recurrence at index n."""
    c3 = (n+3)**2 * (2*n+3)
    c2 = -(2*n+5) * (35*n**2 + 140*n + 131)
    c1 = (2*n+3) * (35*n**2 + 140*n + 131)
    c0 = -(2*n+5) * (n+1)**2
    return (c0, c1, c2, c3)

# Central Delannoy numbers
def delannoy(N):
    D = [mpf(0)] * (N+1)
    D[0] = 1
    D[1] = 3
    for n in range(1, N):
        D[n+1] = (3*(2*n+1)*D[n] - n*D[n-1]) / (n+1)
    return D

N = 80
D = delannoy(N)
print("\n=== Central Delannoy numbers ===")
print(f"D[0..5] = {[int(D[i]) for i in range(6)]}")

# D_n² sequence
Dsq = [D[n]**2 for n in range(N+1)]
print(f"D²[0..5] = {[int(Dsq[i]) for i in range(6)]}")

# Verify Sym²(Delannoy) recurrence on D_n²
print("\n=== Verify Sym²(Delannoy) recurrence on D_n² ===")
max_err = 0
for n in range(N-3):
    c0, c1, c2, c3 = sym2_del_coeffs(n)
    resid = c3*Dsq[n+3] + c2*Dsq[n+2] + c1*Dsq[n+1] + c0*Dsq[n]
    if abs(resid) > max_err:
        max_err = float(abs(resid))
print(f"Max residual: {max_err:.3e}")

# ---- Problem 2.5 CMF matrix ----
def M25(n):
    """3x3 matrix M(n) for Problem 2.5."""
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

# Initial matrix A (2×3)
A25 = matrix([[mpf(30921), mpf(-32972), mpf(8240)],
              [mpf(33750), mpf(-36000), mpf(9000)]])

# ---- Compute challenge scalar sequence q_N ----
print("\n=== Challenge scalar sequence q_N (first column) ===")
q_vals = []
T = matrix([[1,0,0],[0,1,0],[0,0,1]])
for N_idx in range(60):
    AT = A25 * T
    q_val = AT[1, 0]
    q_vals.append(q_val)
    T = T * M25(N_idx)

print(f"q[0] = {nstr(q_vals[0], 20)}")
print(f"q[1] = {nstr(q_vals[1], 20)}")
print(f"q[2] = {nstr(q_vals[2], 20)}")

# ---- Compute gauge ratio q_N / ((-16)^N * Dsq[N]) ----
print("\n=== Gauge analysis: q_N / ((-16)^N * D_N²) ===")
for N_idx in range(min(15, len(q_vals))):
    if Dsq[N_idx] != 0:
        ratio = q_vals[N_idx] / (mpf(-16)**N_idx * Dsq[N_idx])
        print(f"  N={N_idx}: ratio = {nstr(ratio, 20)}")

# ---- Instead try: find q_N / (g_N * U_N) where U_N is a LINEAR COMBINATION
#      of D_n², D_n * E_n, E_n² (where E_n is the second Delannoy solution) ----
# The second solution E_n of the Delannoy recurrence has E_n ~ C * (3-2√2)^n / √n
# For practical computation, we need E_n via the recurrence with different initial conditions.

# Two independent solutions of (n+1)D_{n+1} = 3(2n+1)D_n - nD_{n-1}:
# y1 = D_n (Delannoy), y2 = something with y2_0 = 0, y2_1 = 1

def delannoy_pair(N):
    """Compute two independent solutions of the Delannoy recurrence."""
    y1 = [mpf(0)] * (N+1)
    y2 = [mpf(0)] * (N+1)
    y1[0] = 1; y1[1] = 3
    y2[0] = 0; y2[1] = 1
    for n in range(1, N):
        y1[n+1] = (3*(2*n+1)*y1[n] - n*y1[n-1]) / (n+1)
        y2[n+1] = (3*(2*n+1)*y2[n] - n*y2[n-1]) / (n+1)
    return y1, y2

y1, y2 = delannoy_pair(N)
print(f"\ny1[0..5] = {[nstr(y1[i],10) for i in range(6)]}")
print(f"y2[0..5] = {[nstr(y2[i],10) for i in range(6)]}")

# Three Sym² solutions: y1², y1*y2, y2²
sym2_sols = {
    'y1^2': [y1[n]**2 for n in range(N+1)],
    'y1*y2': [y1[n]*y2[n] for n in range(N+1)],
    'y2^2': [y2[n]**2 for n in range(N+1)],
}

# Try to find A, B, C such that q_N / g_N = A*y1² + B*y1*y2 + C*y2²
# where g_N = (-16)^N * (N!)^7 * correction
# First try g_N = (-16)^N * fac(N)^7

print("\n=== Trying gauge g_N = (-16)^N * (N!)^7 ===")
gvals = [mpf(-16)**n * fac(n)**7 for n in range(20)]
r_vals = [q_vals[n] / gvals[n] if gvals[n] != 0 else None for n in range(20)]

# r_vals should be a linear combination of sym2 solutions
# Solve: r[n] = A * y1[n]² + B * y1[n]*y2[n] + C * y2[n]²
# Using n = 0, 1, 2 to find A, B, C

from mpmath import lu_solve

if r_vals[0] is not None:
    mat = matrix([
        [y1[0]**2, y1[0]*y2[0], y2[0]**2],
        [y1[1]**2, y1[1]*y2[1], y2[1]**2],
        [y1[2]**2, y1[2]*y2[2], y2[2]**2],
    ])
    rhs = matrix([r_vals[0], r_vals[1], r_vals[2]])
    try:
        coeffs = lu_solve(mat, rhs)
        print(f"Coefficients: A={nstr(coeffs[0],20)}, B={nstr(coeffs[1],20)}, C={nstr(coeffs[2],20)}")

        # Verify at higher indices
        for n in range(3, 15):
            predicted = coeffs[0]*y1[n]**2 + coeffs[1]*y1[n]*y2[n] + coeffs[2]*y2[n]**2
            actual = r_vals[n]
            if actual is not None and predicted != 0:
                rel_err = float(abs((actual - predicted)/predicted))
                print(f"  n={n}: rel_err = {rel_err:.3e}")
    except Exception as e:
        print(f"Solve failed: {e}")

# ---- Try different gauges ----
# g_N = (-16)^N * product of Pochhammer terms
# The formal exponent for c₀=-16 mode is 33/2, so the power is n^{33/2}
# For Sym², the power is n^{-1} (from D_n ~ n^{-1/2}, so D_n² ~ n^{-1})
# So gauge must contribute n^{33/2 - (-1)} = n^{35/2}
# Factor n^7 from (n!)^7 ~ n^{7n+7/2} e^{-7n} (2π)^{7/2}
# So gauge growth ~ (-16)^n n^{7n} (Stirling) which accounts for the exponential
# The power correction from (n!)^7 gives n^{7/2}, but we need n^{35/2},
# so we need additional n^{14} worth of Pochhammer.

# Let's try: g_N = (-16)^N * Γ(N+1)^7 * Γ(N+3/2)^a * ...
# More systematically: just look at q_N / (sym2_sol) ratio

print("\n=== Direct ratio q_N / Dsq_N ===")
for n in range(min(10, len(q_vals))):
    if Dsq[n] != 0:
        r = q_vals[n] / Dsq[n]
        print(f"  N={n}: q/D² = {nstr(r, 30)}")

print("\n=== Ratio of consecutive q_N / ((-16)^N * D_N²) ===")
for n in range(1, min(15, len(q_vals))):
    if Dsq[n] != 0 and Dsq[n-1] != 0:
        r_n = q_vals[n] / (mpf(-16)**n * Dsq[n])
        r_nm1 = q_vals[n-1] / (mpf(-16)**(n-1) * Dsq[n-1])
        if r_nm1 != 0:
            ratio = r_n / r_nm1
            print(f"  N={n}: ratio[N]/ratio[N-1] = {nstr(ratio, 30)}")
