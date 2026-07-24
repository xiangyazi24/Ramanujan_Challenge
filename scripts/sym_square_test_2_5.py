#!/usr/bin/env python3
"""Test the symmetric-square condition for Problem 2.5 matrix M(n).

A 3x3 matrix M(n) is a symmetric square (up to coboundary) iff there
exists a symmetric 3x3 matrix J(n) and scalar lambda(n) such that:
  M(n)^T · J(n+1) · M(n) = lambda(n) · J(n)

This is equivalent to M preserving a quadratic form up to scaling.

Approach: propagate J forward using M(n) and check if it stabilizes
to a rational function of n.
"""
from mpmath import mp, mpf, matrix, nstr
import numpy as np

mp.dps = 100

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

def mat_to_np(M):
    return np.array([[float(M[i,j]) for j in range(3)] for i in range(3)])

# Method 1: Find J by solving the eigenvalue problem at a specific n.
# M^T J M = lambda J can be rewritten as (M^T ⊗ M) vec(J) = lambda vec(J)
# restricted to symmetric matrices.

# The 6-dimensional symmetric representation:
# J = [[J11, J12, J13], [J12, J22, J23], [J13, J23, J33]]
# We index: (11, 12, 13, 22, 23, 33)

def sym_kronecker(M_np):
    """Compute the 6x6 matrix representing M^T ⊗ M on symmetric matrices."""
    # (M^T J M)_{ij} = sum_{k,l} M_{ki} J_{kl} M_{lj}
    # For symmetric J, we use basis {E11, E12+E21, E13+E31, E22, E23+E32, E33}
    # Indices: 0=(1,1), 1=(1,2), 2=(1,3), 3=(2,2), 4=(2,3), 5=(3,3)

    idx_to_ij = [(0,0), (0,1), (0,2), (1,1), (1,2), (2,2)]

    K = np.zeros((6, 6))
    for a, (i, j) in enumerate(idx_to_ij):
        for b, (k, l) in enumerate(idx_to_ij):
            # (M^T J M)_{ij} with J = e_k e_l^T + e_l e_k^T (if k!=l) or e_k e_k^T (if k==l)
            # = sum_r sum_s M_{ri} (delta_{rk}delta_{sl} + delta_{rl}delta_{sk}) M_{sj}
            # = M_{ki} M_{lj} + M_{li} M_{kj}  (if k!=l)
            # = M_{ki} M_{kj}  (if k==l)
            if k == l:
                val = M_np[k, i] * M_np[k, j]
            else:
                val = M_np[k, i] * M_np[l, j] + M_np[l, i] * M_np[k, j]

            # Adjust for the normalization of the off-diagonal basis elements
            # We're using J_{kl} as coefficient (not 2*J_{kl}), so:
            K[a, b] = val

    return K

# Test for several values of n
print("=== Symmetric-square eigenvalue test ===")
print("For each n, compute eigenvalues of sym-Kronecker (M^T ⊗ M)|_sym")
print("If M is symmetric square, one eigenvalue should be det(N)^2 = (det M)^{2/3}")
print()

for n_val in [0, 1, 2, 5, 10, 20]:
    M_n = M_mat(n_val)
    M_np = mat_to_np(M_n)
    K = sym_kronecker(M_np)
    eigvals = np.linalg.eigvals(K)
    eigvals_sorted = sorted(eigvals, key=lambda x: abs(x), reverse=True)
    det_M = np.linalg.det(M_np)
    print(f"n={n_val}: det(M)={det_M:.6e}")
    print(f"  K eigenvalues: {[f'{e:.6e}' for e in eigvals_sorted]}")

    # For a symmetric square M = Sym²(N), the eigenvalues of K are:
    # lambda_i * lambda_j for i <= j where lambda_i are eigenvalues of N⊗N|_sym
    # Actually, for Sym²(N), the eigenvalue of the Kronecker product on symmetric
    # matrices is det(N)^2, and the eigenvalue ratios should be specific.

    # More useful: check if there's a common eigenvector J across different n values.
    # The INVARIANT CONIC test: if J is n-independent, all n-values share the same eigenvector.

# Method 2: Direct propagation.
# Start with J(0) = identity (arbitrary), compute J(1) = M(0)^T J(0) M(0) / lambda(0)
# where lambda(0) is chosen to normalize J.
print("\n=== Direct propagation of J(n) ===")
print("Start J(0) = I, propagate J(n+1) via M(n)^T J(n) M(n) / lambda(n)")

J = np.eye(3)
for n_val in range(20):
    M_n = M_mat(n_val)
    M_np = mat_to_np(M_n)
    J_new = M_np.T @ J @ M_np
    # Normalize by J_new[0,0] to keep it bounded
    lam = J_new[0, 0]
    J_new = J_new / lam
    if n_val in [0, 1, 2, 5, 10, 15, 19]:
        print(f"n={n_val}: lambda={lam:.6e}")
        print(f"  J(n+1) normalized:")
        for i in range(3):
            print(f"    [{J_new[i,0]:.10f}, {J_new[i,1]:.10f}, {J_new[i,2]:.10f}]")
    J = J_new

# If J stabilizes, the matrix is (asymptotically) preserving a quadratic form.
# The EXACT test requires J(n) to be rational in n, not just asymptotically stable.

# Method 3: Check det(M)^{2/3} at integer n
print("\n\n=== Check if det(M(n))^{2/3} is rational ===")
for n_val in [0, 1, 2, 3, 5, 10]:
    M_n = M_mat(n_val)
    M_np = mat_to_np(M_n)
    det_M = np.linalg.det(M_np)
    if det_M < 0:
        det_23 = -((-det_M)**(2/3))
    else:
        det_23 = det_M**(2/3)
    print(f"  n={n_val}: det(M) = {det_M:.6e}, det^(2/3) = {det_23:.6e}")
