#!/usr/bin/env python3
"""Problem 2.5: Test gauge transforms to find (S-1) factor.

Key insight from Q4735: The limiting matrix C_∞ has eigenvalues
17+12√2, 17-12√2, and 1. The eigenvalue-1 mode corresponds to a
solution with only polynomial (not exponential) growth.

We need to find h(N) such that q̃_N = q_N / h(N) satisfies an order-3
recurrence where (S-1) is a right factor.

Candidate gauges:
  h₁(N) = A052795(N+1) = (6(N+1))! / (5(N+1)+1)! = (6N+6)!/(5N+6)!
  h₂(N) = product of the GCD factors: (N+2)^5 (N+3)^5 (N+4)^2 (2N+3)(2N+5)^2(2N+7)^2
"""
from mpmath import mp, mpf, fac, nstr, matrix
mp.dps = 200

def M_mat(n):
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
    return matrix([[m11, m12, m13], [m21, m22, m23], [m31, m32, m33]])

A_mat = matrix([[mpf(30921), mpf(-32972), mpf(8240)],
                [mpf(33750), mpf(-36000), mpf(9000)]])

# Compute q_N values
N_max = 60
T = matrix([[1,0,0],[0,1,0],[0,0,1]])
q_vals = []
for N in range(N_max + 1):
    AT = A_mat * T
    q_vals.append(AT[1,0])  # q_N = second row, first column
    T = T * M_mat(N)

# 1. Ratios q_{N+1}/q_N
print("=== q_{N+1}/q_N ratios ===")
for N in [5, 10, 15, 20, 25, 30, 40, 50]:
    if N+1 < len(q_vals) and q_vals[N] != 0:
        ratio = q_vals[N+1] / q_vals[N]
        print(f"  N={N:2d}: q_{N+1}/q_N = {nstr(ratio, 15)}")

# 2. Ratio / N^7 (should approach eigenvalue of C_∞)
print("\n=== q_{N+1}/(q_N * N^7) — should approach a Poincaré root ===")
for N in [10, 15, 20, 25, 30, 40, 50]:
    if N+1 < len(q_vals) and q_vals[N] != 0:
        ratio = q_vals[N+1] / (q_vals[N] * mpf(N)**7)
        print(f"  N={N:2d}: ratio/N^7 = {nstr(ratio, 15)}")

# 3. Test A052795 gauge
print("\n=== A052795 gauge: h(N) = (6N)!/(5N+1)! ===")
def A052795(n):
    return fac(6*n) / fac(5*n + 1)

for N in [3, 5, 10, 15, 20]:
    if N < len(q_vals):
        h = A052795(N+1)  # shifted
        gauged = q_vals[N] / h
        print(f"  N={N:2d}: q_N/A052795(N+1) = {nstr(gauged, 15)}")

# 4. Let's look at det(M(N)) / eigenvalue product
# det = -8*(n+1)*(n+2)^6*(n+3)^5*(2n+3)^2*(2n+5)^3*(2n+7)^4
# Product of eigenvalues = -det (for 3x3)
# If eigenvalues are λ₁ ~ (17+12√2)n^7, λ₂ ~ (17-12√2)n^7, λ₃ ~ n^7
# then det ~ (17+12√2)(17-12√2)(1) * n^21 = 1 * n^21 = n^21
# But actual det LC = -4096. So the eigenvalue scaling is off.
# Let me check: det degree = 21, LC = -4096.
# (17+12√2)(17-12√2) = 289-288 = 1.
# So eigenvalue product ~ a³ n^21 where a is from scaling, a³ = -4096 → a = -16.
# So eigenvalues are ~ -16(17+12√2)n^d₁, -16(17-12√2)n^d₂, -16*1*n^d₃
# with d₁+d₂+d₃ = 21.

# 5. Try to find the eigenvalue-1 mode by looking at the matrix eigenvector.
print("\n=== Eigenvalues of M(N)/N^7 for large N ===")
for N in [20, 30, 50]:
    MN = M_mat(N)
    # Scale: M[i,j] has degree 7+i-j (roughly). So D^{-1} M D where D = diag(1, N, N^2)?
    # Actually M[0,0] ~ n^7, M[0,1] ~ n^6, M[0,2] ~ n^4
    # M[1,0] ~ n^8, M[1,1] ~ n^7, M[1,2] ~ n^5
    # M[2,0] ~ n^9, M[2,1] ~ n^8, M[2,2] ~ n^7
    # So D = diag(1, n, n^2) works: D^{-1} M D ~ n^7 * C_∞
    D = matrix([[1, 0, 0], [0, mpf(N), 0], [0, 0, mpf(N)**2]])
    Dinv = matrix([[1, 0, 0], [0, 1/mpf(N), 0], [0, 0, 1/mpf(N)**2]])
    C = Dinv * MN * D / mpf(N)**7
    # Eigenvalues
    import numpy as np
    C_np = np.array([[float(C[i,j]) for j in range(3)] for i in range(3)])
    eigvals = np.linalg.eigvals(C_np)
    print(f"  N={N}: eigenvalues of D^-1 M D / N^7 = {sorted([float(e.real) for e in eigvals], reverse=True)}")

# The scaling D=diag(1,n,n²) normalizes the matrix.
# This means the solution vector has components scaling as (1, n, n²) * eigenvalue^N.
# For the q_N component (1st component, scaled by 1):
#   q_N ~ C₁ * λ₁^N + C₂ * λ₂^N + C₃ * λ₃^N
# where λ₁ ~ -16(17+12√2)N^7, λ₂ ~ -16(17-12√2)N^7, λ₃ ~ -16*N^7
# (using the det LC = -4096 = (-16)^3 implies each eigenvalue ~ -16 * n^7 * poincare_root)

# 6. The "eigenvalue 1 mode" solution has growth q_N ~ C * (-16)^N * N^{7N} * ...
# which is superexponential. The gauge h(N) should absorb this.
# Natural candidate: h(N) = (-16)^N * prod_{k=0}^{N-1} k^7 = (-16)^N * (N-1)!^7
# But this is too crude. Let's look at exact eigenvalue at finite N.

# 7. Direct check: is 1 a root of the normalized char poly at each N?
print("\n=== Is 1 a root of the normalized Poincaré poly at finite N? ===")
for N in [5, 10, 20, 30]:
    MN = M_mat(N)
    D = matrix([[1, 0, 0], [0, mpf(N), 0], [0, 0, mpf(N)**2]])
    Dinv = matrix([[1, 0, 0], [0, 1/mpf(N), 0], [0, 0, 1/mpf(N)**2]])
    C = Dinv * MN * D / mpf(N)**7
    # Evaluate char poly at lambda=1
    I3 = matrix([[1,0,0],[0,1,0],[0,0,1]])
    det_C_minus_I = (C - I3).det()
    print(f"  N={N}: det(C/n^7 - I) = {nstr(det_C_minus_I, 10)}")
