#!/usr/bin/env python3
"""P2.5: Extract the scalar recurrence from the CMF matrix.

Compute Q_{N,0} (first component of Q_row · M(0)·M(1)·...·M(N-1))
for many N values, then use Berlekamp-Massey to find the scalar recurrence.
After finding it, compare with the Delannoy-square recurrence after H_n twist.
"""
from mpmath import mp, mpf, matrix, nstr, rf

mp.dps = 200

def M_exact(n):
    n = int(n)
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

# Compute using exact integer arithmetic (no floating point!)
# Q_row starts as [33750, -36000, 9000]
# Q_{N+1,j} = Σ_k Q_{N,k} · M(N)[k,j]

# Actually, we want the sequence Q_{N,0} for the scalar recurrence.
# But we need all 3 components to iterate.

print("Computing Q_{N,j} for N = 0..120 using exact integers...")

q = [33750, -36000, 9000]
Q0_seq = [33750]  # Q_{0,0}

for N in range(120):
    M = M_exact(N)
    new_q = [0, 0, 0]
    for j in range(3):
        for k in range(3):
            new_q[j] += q[k] * M[k][j]
    q = new_q
    Q0_seq.append(q[0])
    if N < 5 or N % 20 == 19:
        print(f"  N={N+1}: Q0 has {len(str(abs(q[0])))} digits")

print(f"  Total terms: {len(Q0_seq)}")

# Now find the scalar recurrence for Q0_seq.
# The recurrence is: c3(n) Q0(n+3) + c2(n) Q0(n+2) + c1(n) Q0(n+1) + c0(n) Q0(n) = 0
# where c0, c1, c2, c3 are polynomials in n.
# From the summary: degree staircase (28, 21, 14, 7) means
# deg c3 = 28, deg c2 = 21, deg c1 = 14, deg c0 = 7.

# Total coefficients: 29 + 22 + 15 + 8 = 74 unknowns.
# Each value of n gives one linear equation.
# We need at least 74 equations, so need Q0(0), ..., Q0(76).

# Actually, let me first verify the order and degree by trying smaller degrees.
# Start with the expected (28, 21, 14, 7) staircase.

# But first, let me try normalizing by H_n to get the "nice" recurrence.
def H_int(n):
    """H_n = (-16)^n (2)_n^2 (3)_n^2 (5/2)_n (7/2)_n^2, computed exactly as a rational number."""
    # Use exact integer/rational arithmetic
    from fractions import Fraction
    result = Fraction(1)
    for k in range(n):
        # (-16) * (k+2)^2 * (k+3)^2 * (k+5/2) * (k+7/2)^2
        result *= Fraction(-16)
        result *= Fraction(k+2)**2
        result *= Fraction(k+3)**2
        result *= Fraction(2*k+5, 2)
        result *= Fraction(2*k+7, 2)**2
    return result

# Compute Q̂_n = Q0_n / H_n for small n
print("\nQ̂_n = Q0_n / H_n for n = 0..10:")
from fractions import Fraction

Q0_rat = [Fraction(x) for x in Q0_seq[:30]]

for n in range(15):
    Hn = H_int(n)
    if Hn != 0:
        q_hat = Q0_rat[n] / Hn
    else:
        q_hat = Q0_rat[n]
    print(f"  n={n:2d}: Q̂ = {float(q_hat):.6f}")

# Check if Q̂_n satisfies the Delannoy-square recurrence!
# (2n+3)(n+3)² d_{n+3} - (2n+5)(35n²+140n+131) d_{n+2}
#   + (2n+3)(35n²+140n+131) d_{n+1} - (2n+5)(n+1)² d_n = 0
print("\nChecking if Q̂_n satisfies Delannoy-square recurrence:")
Q_hat = []
for n in range(30):
    Hn = H_int(n)
    Q_hat.append(Q0_rat[n] / Hn if Hn else Q0_rat[n])

for n in range(25):
    P_n = 35*n**2 + 140*n + 131
    residual = ((2*n+3)*(n+3)**2 * Q_hat[n+3]
              - (2*n+5)*P_n * Q_hat[n+2]
              + (2*n+3)*P_n * Q_hat[n+1]
              - (2*n+5)*(n+1)**2 * Q_hat[n])
    print(f"  n={n:2d}: residual = {float(residual):.6e}" + (" ✓" if abs(float(residual)) < 1e-10 else " ✗"))

# If Q̂_n does NOT satisfy the Delannoy-square recurrence,
# let's find what recurrence Q̂_n does satisfy.
print("\n" + "="*60)
print("Finding the actual scalar recurrence for Q̂_n")
print("="*60)

# Try order 3 with polynomial coefficients of increasing degree
from fractions import Fraction

def find_recurrence(seq, order, max_deg, verbose=True):
    """Find recurrence Σ_{i=0}^{order} c_i(n) seq[n+i] = 0
    with deg c_i ≤ max_deg."""
    num_coeffs = (order + 1) * (max_deg + 1)
    num_eqs = len(seq) - order
    if num_eqs < num_coeffs + 5:
        return None

    # Use exact rational arithmetic
    # Build the system: for each n, Σ_{i,j} a_{i,j} n^j seq[n+i] = 0
    # Row for each n, columns for (i, j)

    # Use overdetermined system, solve via nullspace
    # Just check: for each candidate degree, build the matrix and find null vectors

    A = []
    for n in range(min(num_eqs, num_coeffs + 10)):
        row = []
        for i in range(order + 1):
            for j in range(max_deg + 1):
                row.append(Fraction(n)**j * seq[n+i])
        A.append(row)

    # Find the nullspace using Gaussian elimination
    m = len(A)
    nc = len(A[0])
    # Convert to matrix
    import copy
    mat = [copy.copy(row) for row in A]

    # Gaussian elimination
    pivot_cols = []
    r = 0
    for c in range(nc):
        # Find pivot
        found = -1
        for i in range(r, m):
            if mat[i][c] != Fraction(0):
                found = i
                break
        if found < 0:
            continue
        mat[r], mat[found] = mat[found], mat[r]
        pivot_cols.append(c)
        for i in range(m):
            if i != r and mat[i][c] != Fraction(0):
                factor = mat[i][c] / mat[r][c]
                for j in range(nc):
                    mat[i][j] -= factor * mat[r][j]
        r += 1

    rank = r
    null_dim = nc - rank
    if verbose:
        print(f"  Order {order}, max_deg {max_deg}: rank={rank}/{nc}, null_dim={null_dim}")

    if null_dim > 0:
        # Extract null vector: find a free column
        free_cols = [c for c in range(nc) if c not in pivot_cols]
        # For the first free column, set it to 1, others to 0
        x = [Fraction(0)] * nc
        fc = free_cols[0]
        x[fc] = Fraction(1)
        for pr_idx in range(len(pivot_cols)-1, -1, -1):
            pc = pivot_cols[pr_idx]
            s = Fraction(0)
            for j in range(pc+1, nc):
                s += mat[pr_idx][j] * x[j]
            x[pc] = -s / mat[pr_idx][pc]

        # Verify
        max_res = max(abs(sum(A[n][c]*x[c] for c in range(nc)))
                     for n in range(m))
        if verbose:
            print(f"    Max residual: {float(max_res):.2e}")

        # Also verify on held-out data
        holdout_res = Fraction(0)
        for n in range(num_coeffs + 10, min(num_eqs, num_coeffs + 20)):
            val = Fraction(0)
            for i in range(order + 1):
                c_i = sum(x[i*(max_deg+1)+j] * Fraction(n)**j for j in range(max_deg+1))
                val += c_i * seq[n+i]
            holdout_res = max(holdout_res, abs(val))
        if verbose:
            print(f"    Holdout residual: {float(holdout_res):.2e}")

        if holdout_res == 0:
            # Extract polynomial coefficients
            polys = []
            for i in range(order + 1):
                coeffs = [x[i*(max_deg+1)+j] for j in range(max_deg+1)]
                polys.append(coeffs)
            return polys

    return None

# Try order 3 with degrees 3, 4, 5, ...
for deg in range(2, 10):
    result = find_recurrence(Q_hat[:30], 3, deg)
    if result is not None:
        print(f"\n  *** FOUND recurrence of order 3, degree {deg} ***")
        for i, poly in enumerate(result):
            nonzero = [(j, poly[j]) for j in range(len(poly)) if poly[j] != 0]
            print(f"    c_{i}(n) = " + " + ".join(f"({c})*n^{j}" for j, c in nonzero))
        break
