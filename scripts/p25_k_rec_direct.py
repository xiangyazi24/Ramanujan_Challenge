#!/usr/bin/env python3
"""
P2.5: Find k-recurrence for f(k) directly.

Instead of computing f(k) from triangular inversion of Q̂_N,
use a different approach: the f(k) are determined by the CMF
matrix and the initial condition. The e₁-trajectory column 1
values u₁(N) satisfy:

u₁(N) = Σ_{k=0}^N α₁(k) B(N,k)

where α₁(k) = f(k) (since u₁ is the e₁-trajectory of column 1).

For the q-row: f_q(k) = 33750 α₁(k) - 36000 α₂(k) + 9000 α₃(k)

So first compute α_j(k) for j=1,2,3 and many k values,
then search for recurrence.

Strategy: compute 60 terms, search for rec of order 3-5 with
polynomial coefficients of degree up to 12.
"""
from fractions import Fraction as Fr
from math import comb
import sys

def M_entries(n):
    m11 = (-2*n-5)*(n+3)**2*(136*n**4+1424*n**3+5548*n**2+9551*n+6141)
    m12 = 384*n**6+6384*n**5+44168*n**4+162698*n**3+336377*n**2+369933*n+169011
    m13 = -(480*n**4+4980*n**3+19210*n**2+32690*n+20730)
    m21 = (n+2)**2*(n+3)**2*(4*n+10)*(48*n**3+386*n**2+1017*n+879)
    m22 = (n+2)**2*(-272*n**5-3848*n**4-21732*n**3-61184*n**2-85761*n-47808)
    m23 = (n+2)**2*(320*n**3+2540*n**2+6610*n+5640)
    m31 = (-4*n-10)*(n+2)**2*(n+3)**2*(32*n**4+302*n**3+1037*n**2+1530*n+813)
    m32 = (n+2)**2*(192*n**6+2984*n**5+19116*n**4+64452*n**3+120256*n**2+117279*n+46476)
    m33 = (n+2)**2*(-16*n**5-408*n**4-2912*n**3-8884*n**2-12254*n-6240)
    return [[m11,m12,m13],[m21,m22,m23],[m31,m32,m33]]

def delta_H(n):
    return -2*(n+2)**2*(n+3)**2*(2*n+5)*(2*n+7)**2

def B(N, k):
    if k < 0 or k > N:
        return Fr(0)
    return Fr(2**k * comb(2*k, k) * comb(N, k) * comb(N+k, k))

KMAX = 60
print(f"Computing CMF trajectories for N=0..{KMAX}...", flush=True)

# Compute all three e_j trajectories (column 1 only)
rows = [[Fr(1), Fr(0), Fr(0)],
        [Fr(0), Fr(1), Fr(0)],
        [Fr(0), Fr(0), Fr(1)]]

# u_j(N) = e_j · Π M_H(m) · e₁
u = {j: [rows[j][0]] for j in range(3)}

for N in range(KMAX):
    M = M_entries(N)
    d = Fr(delta_H(N))
    MH = [[Fr(M[i][j]) / d for j in range(3)] for i in range(3)]
    for j in range(3):
        r = rows[j]
        new_r = [sum(r[i]*MH[i][k] for i in range(3)) for k in range(3)]
        rows[j] = new_r
        u[j].append(new_r[0])
    if (N+1) % 10 == 0:
        print(f"  N={N+1} done", flush=True)

# Triangular inversion: α_j(k) from u_j(N)
print(f"\nTriangular inversion for α_j(k)...", flush=True)
alpha = {j: [] for j in range(3)}
for j in range(3):
    for K in range(KMAX + 1):
        rhs = u[j][K]
        for k in range(K):
            rhs -= alpha[j][k] * B(K, k)
        bKK = B(K, K)
        alpha[j].append(rhs / bKK)

# Verify
for j in range(3):
    for N in [0, 5, 10, 20]:
        if N <= KMAX:
            val = sum(alpha[j][k] * B(N, k) for k in range(N+1))
            if val != u[j][N]:
                print(f"  VERIFY FAIL: j={j}, N={N}")
                break
    else:
        print(f"  α_{j} verified", flush=True)

# Now search for k-recurrence for α₁(k)
# Try: c_r(k)α(k+r) + ... + c_0(k)α(k) = 0
# where c_i(k) are polynomials of degree ≤ d

def search_recurrence(seq, name, max_order=5, max_deg=12):
    """Search for linear recurrence with polynomial coefficients."""
    L = len(seq)
    print(f"\nSearching for k-recurrence of {name} (L={L})...", flush=True)

    for order in range(2, max_order + 1):
        for deg in range(1, max_deg + 1):
            # Number of unknowns per coefficient: deg+1
            # Total unknowns: (order+1)*(deg+1)
            n_unknowns = (order + 1) * (deg + 1)
            # Number of equations: L - order (one per valid k)
            n_eqs = L - order
            if n_eqs < n_unknowns + 3:
                continue

            # Build matrix: for each k in [0, L-order-1],
            # Σ_{i=0}^{order} c_i(k) · seq[k+i] = 0
            # c_i(k) = Σ_{j=0}^{deg} a_{i,j} k^j
            # So: Σ_i Σ_j a_{i,j} k^j seq[k+i] = 0
            # Unknowns: a_{i,j} for i=0..order, j=0..deg

            # Use modular arithmetic first for speed (detect nullity)
            P = 2**61 - 1  # Mersenne prime

            import numpy as np

            # Build matrix mod P
            mat = []
            for k in range(min(n_eqs, n_unknowns + 5)):
                row = []
                for i in range(order + 1):
                    s = seq[k + i]
                    # Convert fraction to int mod P
                    s_mod = (s.numerator * pow(s.denominator, P-2, P)) % P
                    for j in range(deg + 1):
                        row.append((pow(k, j, P) * s_mod) % P)
                mat.append(row)

            # Gaussian elimination mod P
            mat_np = [[int(x) for x in row] for row in mat]
            nrows = len(mat_np)
            ncols = n_unknowns
            pivots = []
            for col in range(ncols):
                # Find pivot
                found = False
                for r in range(len(pivots), nrows):
                    if mat_np[r][col] % P != 0:
                        # Swap
                        mat_np[len(pivots)], mat_np[r] = mat_np[r], mat_np[len(pivots)]
                        # Eliminate
                        inv = pow(mat_np[len(pivots)][col], P-2, P)
                        for r2 in range(nrows):
                            if r2 != len(pivots) and mat_np[r2][col] % P != 0:
                                factor = (mat_np[r2][col] * inv) % P
                                for c2 in range(ncols):
                                    mat_np[r2][c2] = (mat_np[r2][c2] - factor * mat_np[len(pivots)][c2]) % P
                        pivots.append(col)
                        found = True
                        break
                if not found:
                    pass

            rank = len(pivots)
            nullity = ncols - rank
            if nullity > 0:
                print(f"  order={order}, deg={deg}: nullity={nullity} "
                      f"({n_unknowns} unknowns, rank={rank})")
                if nullity == 1:
                    print(f"  *** FOUND: unique recurrence of order {order}, "
                          f"degree {deg} ***")
                    return order, deg
                elif nullity >= 2:
                    print(f"  (nullity > 1, try lower degree)")

    print(f"  No recurrence found up to order {max_order}, degree {max_deg}")
    return None, None

# Search for α₁(k) recurrence
o, d = search_recurrence(alpha[0], "α₁(k)", max_order=4, max_deg=10)

if o is not None:
    # Also check α₂ and α₃
    print("\nChecking if α₂ and α₃ satisfy the same recurrence...")
    o2, d2 = search_recurrence(alpha[1], "α₂(k)", max_order=o, max_deg=d)
    o3, d3 = search_recurrence(alpha[2], "α₃(k)", max_order=o, max_deg=d)

# Also search for f(k) = q-row combination
q = [Fr(33750), Fr(-36000), Fr(9000)]
p = [Fr(30921), Fr(-32972), Fr(8240)]

f_vals = [sum(q[j] * alpha[j][k] for j in range(3)) for k in range(KMAX + 1)]
g_vals = [sum(p[j] * alpha[j][k] for j in range(3)) for k in range(KMAX + 1)]

print(f"\nf(k) positivity: all positive? {all(f > 0 for f in f_vals)}")
print(f"g(k)/f(k) for last 5 terms:")
from mpmath import mp, mpf, catalan
mp.dps = 50
for k in range(KMAX - 4, KMAX + 1):
    r = mpf(g_vals[k].numerator) / mpf(g_vals[k].denominator) / (mpf(f_vals[k].numerator) / mpf(f_vals[k].denominator))
    print(f"  k={k}: g/f = {mp.nstr(r, 40)}")
print(f"  G   = {mp.nstr(catalan, 40)}")

o_f, d_f = search_recurrence(f_vals, "f(k)", max_order=4, max_deg=10)

print("\nDone.")
