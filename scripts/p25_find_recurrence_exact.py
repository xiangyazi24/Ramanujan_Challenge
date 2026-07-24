#!/usr/bin/env python3
"""Re-derive the scalar recurrence for Q̂_n = Q_{N,0}/H_n using EXACT arithmetic."""
from fractions import Fraction

def M_entries(n):
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

# H_{n+1}/H_n = -2(n+2)^2(n+3)^2(2n+5)(2n+7)^2
def delta(n):
    return Fraction(-2*(n+2)**2*(n+3)**2*(2*n+5)*(2*n+7)**2)

# Compute Q̂_n for n = 0..NMAX
NMAX = 80
q = [Fraction(33750), Fraction(-36000), Fraction(9000)]
Q0_raw = [Fraction(33750)]
H_vals = [Fraction(1)]

for N in range(NMAX):
    M = M_entries(N)
    new_q = [Fraction(0)]*3
    for j in range(3):
        for k in range(3):
            new_q[j] += q[k] * M[k][j]
    q = new_q
    Q0_raw.append(q[0])
    H_vals.append(H_vals[-1] * delta(N))

Q_hat = [Q0_raw[n] / H_vals[n] for n in range(NMAX+1)]

print(f"Computed {NMAX+1} exact Q̂_n values")
print(f"Q̂_0 = {Q_hat[0]}")
print(f"Q̂_1 = {Q_hat[1]}")
print(f"Q̂_2 = {Q_hat[2]}")
print(f"Q̂_10 = {float(Q_hat[10]):.15e}")

# Find minimal recurrence: try order r=1,2,3 with increasing degree d
import sys

def find_recurrence(order, max_deg, Q, verbose=True):
    """Find recurrence of given order and polynomial degree max_deg."""
    num_unknowns = (order + 1) * (max_deg + 1)
    num_equations = len(Q) - order
    if num_equations < num_unknowns:
        if verbose:
            print(f"  Not enough data: {num_equations} equations < {num_unknowns} unknowns")
        return None

    # Build matrix: row for each n, columns for coefficients a_{s,k} where c_s(n) = Σ a_{s,k} n^k
    # Equation: Σ_{s=0}^{order} c_s(n) Q[n+s] = 0
    # i.e., Σ_{s=0}^{order} Σ_{k=0}^{max_deg} a_{s,k} n^k Q[n+s] = 0

    rows = []
    for n in range(min(num_equations, num_unknowns + 10)):
        row = []
        for s in range(order + 1):
            for k in range(max_deg + 1):
                row.append(Fraction(n)**k * Q[n+s])
            # end k
        # end s
        rows.append(row)

    actual_eqs = len(rows)
    nc = num_unknowns

    # Gaussian elimination over Q
    mat = [row[:] for row in rows]
    nr = len(mat)
    pivots = []

    for c in range(nc):
        # Find nonzero pivot
        piv = -1
        for i in range(len(pivots), nr):
            if mat[i][c] != 0:
                piv = i
                break
        if piv == -1:
            continue
        mat[len(pivots)], mat[piv] = mat[piv], mat[len(pivots)]
        piv_row = len(pivots)
        piv_val = mat[piv_row][c]
        # Eliminate
        for i in range(nr):
            if i != piv_row and mat[i][c] != 0:
                f = Fraction(mat[i][c], piv_val)
                for j in range(nc):
                    mat[i][j] -= f * mat[piv_row][j]
        pivots.append(c)

    rank = len(pivots)
    null_dim = nc - rank

    if verbose:
        print(f"  Order {order}, degree {max_deg}: {nc} unknowns, {actual_eqs} equations, rank={rank}, null_dim={null_dim}")

    if null_dim == 0:
        return None

    # Extract null vector
    free_cols = [c for c in range(nc) if c not in pivots]
    x = [Fraction(0)] * nc
    fc = free_cols[0]
    x[fc] = Fraction(1)

    # Back-substitute
    for pr_idx in range(len(pivots)-1, -1, -1):
        pc = pivots[pr_idx]
        s = sum(mat[pr_idx][j] * x[j] for j in range(pc+1, nc))
        x[pc] = -s / mat[pr_idx][pc]

    # Verify on ALL available data
    max_n_verify = len(Q) - order
    bad = False
    for n in range(max_n_verify):
        residual = Fraction(0)
        for s in range(order + 1):
            c_s = Fraction(0)
            for k in range(max_deg + 1):
                c_s += x[s*(max_deg+1) + k] * Fraction(n)**k
            residual += c_s * Q[n+s]
        if residual != 0:
            if verbose:
                print(f"  FAILED verification at n={n}: residual != 0")
            bad = True
            break

    if not bad:
        if verbose:
            print(f"  VERIFIED on n=0..{max_n_verify-1}")
        return x, max_deg, order
    return None

# Search
for order in [1, 2, 3]:
    for deg in range(1, 16):
        print(f"\nTrying order={order}, degree={deg}...")
        result = find_recurrence(order, deg, Q_hat)
        if result is not None:
            x, max_deg, ord_ = result
            print(f"\n*** FOUND: order={ord_}, degree={max_deg} ***")

            # Extract and print coefficients
            for s in range(ord_ + 1):
                coeffs = []
                for k in range(max_deg + 1):
                    coeffs.append(x[s*(max_deg+1) + k])
                # Clear denominators
                from math import gcd
                denoms = [c.denominator for c in coeffs if c != 0]
                if denoms:
                    lcm = denoms[0]
                    for d in denoms[1:]:
                        lcm = lcm * d // gcd(lcm, d)
                    int_coeffs = [int(c * lcm) for c in coeffs]
                    # Divide by GCD
                    g = 0
                    for ic in int_coeffs:
                        g = gcd(g, abs(ic))
                    if g > 0:
                        int_coeffs = [ic // g for ic in int_coeffs]
                    print(f"  c_{s}(n) = {int_coeffs} (coefficients of 1, n, n^2, ...)")

            # Print leading coefficients
            print("\nLeading coefficients (degree-{} terms):".format(max_deg))
            for s in range(ord_ + 1):
                lc = x[s*(max_deg+1) + max_deg]
                print(f"  lc(c_{s}) = {lc}")

            # Poincaré polynomial
            lcs = [x[s*(max_deg+1) + max_deg] for s in range(ord_ + 1)]
            print(f"\nPoincaré polynomial: {lcs[0]} + {lcs[1]}ξ + {lcs[2]}ξ² + {lcs[3] if ord_>=3 else ''}ξ³")

            sys.exit(0)

        if (order+1)*(deg+1) > 60:
            print(f"  Too many unknowns ({(order+1)*(deg+1)}), stopping degree search for order={order}")
            break

print("\nNo recurrence found!")
