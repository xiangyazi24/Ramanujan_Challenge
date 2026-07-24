#!/usr/bin/env python3
"""Problem 2.5: Extract the scalar recurrence from the CMF product.

Compute Q_{N,1} = e_1 · M(0)·M(1)·...·M(N-1) · e_1^T  (first column of product)
and find the minimal recurrence satisfied by this sequence.

Use Berlekamp-Massey over Q: find c_0,...,c_r such that
c_r(n) Q_{n+r} + ... + c_0(n) Q_n = 0  for all n,
where c_j are polynomial functions of n.
"""
from fractions import Fraction as F

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

def mat_mul_int(A, B):
    return [[sum(A[i][k]*B[k][j] for k in range(3)) for j in range(3)] for i in range(3)]

# Compute Q_{N,j} for j=0,1,2 (three columns of the product)
print("Computing CMF product sequences...")
N_MAX = 60
I3 = [[1,0,0],[0,1,0],[0,0,1]]
prod = [row[:] for row in I3]

# Q[j][N] = (e_1 · M(0)···M(N-1))_j = prod[0][j]
Q = [[], [], []]

for N in range(N_MAX + 1):
    for j in range(3):
        Q[j].append(prod[0][j])
    if N < N_MAX:
        prod = mat_mul_int(prod, M_int(N))
    if N % 10 == 0:
        print(f"  N={N}: Q[0]={Q[0][-1]}, len(digits)~{len(str(abs(Q[0][-1])))}")

print(f"\nComputed {N_MAX+1} terms.")

# Now find the recurrence: assume order r with polynomial coefficients of degree d.
# c_r(n) Q_{n+r} + ... + c_0(n) Q_n = 0
# where c_j(n) = sum_{m=0}^d a_{j,m} n^m
#
# This gives a linear system in the a_{j,m} coefficients.
# We need (r+1)(d+1) unknowns, and N_MAX - r equations.

seq = Q[0]  # first column

print("\n=== Finding minimal recurrence for Q[0] ===")

for r in range(3, 8):
    for d in range(5, 30):
        n_unknowns = (r + 1) * (d + 1)
        n_equations = N_MAX - r

        if n_equations < n_unknowns + 5:
            continue

        # Build the system: for each n, sum_{j=0}^r sum_{m=0}^d a_{j,m} n^m Q[n+j] = 0
        rows = []
        for n in range(n_equations):
            row = []
            for j in range(r + 1):
                for m in range(d + 1):
                    row.append(F(n**m * seq[n + j]))
            rows.append(row)

        # Gaussian elimination to find the null space
        mat = [row[:] for row in rows]
        n_r = len(mat)
        n_c = n_unknowns

        pivots = []
        piv_row = 0
        for col in range(n_c):
            found = -1
            for rr in range(piv_row, n_r):
                if mat[rr][col] != 0:
                    found = rr
                    break
            if found == -1:
                continue
            mat[found], mat[piv_row] = mat[piv_row], mat[found]
            pv = mat[piv_row][col]
            for rr in range(n_r):
                if rr != piv_row and mat[rr][col] != 0:
                    f = F(mat[rr][col]) / F(pv)
                    for cc in range(n_c):
                        mat[rr][cc] = F(mat[rr][cc]) - f * F(mat[piv_row][cc])
            pivots.append((piv_row, col))
            piv_row += 1

        rank = len(pivots)
        null_dim = n_c - rank

        if null_dim == 0:
            continue

        # Found a recurrence!
        print(f"  r={r}, d={d}: rank={rank}, null_dim={null_dim}")

        if null_dim >= 1:
            # Extract a null vector
            free_cols = [c for c in range(n_c) if c not in [p[1] for p in pivots]]

            # Set the first free variable to 1, rest to 0
            sol = [F(0)] * n_c
            sol[free_cols[0]] = F(1)

            # Back-substitute
            for pr, pc in reversed(pivots):
                val = F(0)
                for cc in range(pc + 1, n_c):
                    val += F(mat[pr][cc]) * sol[cc]
                sol[pc] = -val / F(mat[pr][pc])

            # Extract c_j(n) polynomials
            print(f"\n  Recurrence found: order {r}, coeff degree {d}")
            for j in range(r + 1):
                coeffs = [sol[j * (d + 1) + m] for m in range(d + 1)]
                nonzero = [m for m in range(d + 1) if coeffs[m] != 0]
                if nonzero:
                    max_deg = max(nonzero)
                    print(f"    c_{j}(n): degree {max_deg}, leading = {coeffs[max_deg]}")

            # Verify on all available data
            n_verify = 0
            for n in range(N_MAX - r):
                val = sum(sum(sol[j*(d+1)+m] * F(n**m) for m in range(d+1)) * F(seq[n+j]) for j in range(r+1))
                if val != 0:
                    print(f"    VERIFICATION FAILED at n={n}")
                    break
                n_verify += 1

            print(f"    Verified on {n_verify} points (all available)")

            # Check: is c_r(n) nonzero for all n >= 0?
            # Evaluate c_r at n=0,...,10
            print(f"\n    c_r(n) at n=0,...,10:")
            for n in range(11):
                cr = sum(sol[r*(d+1)+m] * F(n**m) for m in range(d+1))
                print(f"      c_{r}({n}) = {cr}")

            break
    else:
        continue
    break

print("\nDone.")
