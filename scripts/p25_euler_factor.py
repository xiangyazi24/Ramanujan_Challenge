#!/usr/bin/env python3
"""P2.5: Reconstruct the Euler operator L_5 of F(z) = Σ Q̂_n z^n,
then test whether D_n² (Delannoy squares) satisfies L_5.
If yes → L_{3,D} is a right factor of L_5."""
from fractions import Fraction
import sys

# CMF matrix entries
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

def H(n):
    if n == 0:
        return Fraction(1)
    val = Fraction(1)
    for k in range(n):
        f = Fraction(-16) * Fraction(k+2)**2 * Fraction(k+3)**2 * Fraction(2*k+5, 2) * Fraction(2*k+7, 2)**2
        val *= f
    return val

# Compute Q̂_n for n = 0..NMAX
NMAX = 170
print(f"Computing Q̂_n for n = 0..{NMAX}...", flush=True)
qhat = [None] * (NMAX + 1)
row = [Fraction(33750), Fraction(-36000), Fraction(9000)]
qhat[0] = Fraction(33750)
for n in range(NMAX):
    M = M_entries(n)
    new_row = [Fraction(0)] * 3
    for j in range(3):
        for k in range(3):
            new_row[j] += row[k] * Fraction(M[k][j])
    row = new_row
    h = H(n + 1)
    qhat[n + 1] = Fraction(row[0], h) if h != 0 else None
    if (n + 1) % 20 == 0:
        print(f"  n={n+1} done", flush=True)

# Central Delannoy numbers D_n = P_n(3)
print("Computing D_n²...", flush=True)
D = [Fraction(0)] * (NMAX + 1)
D[0] = Fraction(1)
D[1] = Fraction(3)
for n in range(1, NMAX):
    D[n+1] = (Fraction(6*n+3) * D[n] - Fraction(n) * D[n-1]) / Fraction(n+1)
D2 = [d*d for d in D]

# ================================================================
# Part 1: Reconstruct Euler operator of order ORD with z-degree ZDEG
# ================================================================
# Try order 5, z-degree 11 first (72 unknowns as in Q4865)
# If kernel dim != 1, try other parameters

def find_euler_operator(seq, order, zdeg, training_end=None):
    """Find Euler operator Σ z^j P_j(θ) annihilating seq.
    Returns list of (j,m,c) tuples and the coefficient matrix."""
    labels = [(j, m) for j in range(zdeg + 1) for m in range(order + 1)]
    num_unknowns = len(labels)

    if training_end is None:
        training_end = len(seq) - 20

    # Build matrix: for each n, row is [S_j(n-j) * a_{n-j}] = [(n-j)^m * a_{n-j}]
    rows = []
    for n in range(zdeg, training_end):
        row = []
        for j, m in labels:
            if n - j >= 0 and n - j < len(seq) and seq[n-j] is not None:
                row.append(Fraction(n-j)**m * seq[n-j])
            else:
                row.append(Fraction(0))
        rows.append(row)

    nrows = len(rows)
    ncols = num_unknowns
    print(f"  Matrix: {nrows} x {ncols}", flush=True)

    # Gaussian elimination to find kernel (exact over Q)
    # Augment: work on the matrix directly
    mat = [list(r) for r in rows]

    # Row echelon form
    pivot_cols = []
    r = 0
    for c in range(ncols):
        # Find pivot
        found = -1
        for i in range(r, nrows):
            if mat[i][c] != 0:
                found = i
                break
        if found == -1:
            continue
        # Swap
        mat[r], mat[found] = mat[found], mat[r]
        pivot_cols.append(c)
        # Eliminate
        pivot = mat[r][c]
        for i in range(nrows):
            if i != r and mat[i][c] != 0:
                factor = Fraction(mat[i][c], pivot)
                for cc in range(ncols):
                    mat[i][cc] -= factor * mat[r][cc]
        r += 1
        if r % 10 == 0:
            print(f"    Pivot {r}/{ncols}", flush=True)

    rank = len(pivot_cols)
    kernel_dim = ncols - rank
    print(f"  Rank = {rank}, kernel dim = {kernel_dim}", flush=True)

    if kernel_dim != 1:
        return None, kernel_dim

    # Extract kernel vector
    free_cols = [c for c in range(ncols) if c not in pivot_cols]
    assert len(free_cols) == 1
    fc = free_cols[0]

    # Set free variable = 1, solve for pivot variables
    kern = [Fraction(0)] * ncols
    kern[fc] = Fraction(1)

    # Back-substitution
    for i in range(rank - 1, -1, -1):
        pc = pivot_cols[i]
        val = Fraction(0)
        for cc in range(ncols):
            if cc != pc:
                val += mat[i][cc] * kern[cc]
        kern[pc] = -Fraction(val, mat[i][pc])

    # Clear denominators
    from math import gcd
    denoms = [k.denominator for k in kern if k != 0]
    from functools import reduce
    def lcm(a, b):
        return a * b // gcd(a, b)
    L = reduce(lcm, denoms, 1)
    kern_int = [int(k * L) for k in kern]
    # Make primitive
    nz = [abs(x) for x in kern_int if x != 0]
    g = reduce(gcd, nz)
    kern_int = [x // g for x in kern_int]
    # Make first nonzero positive
    for x in kern_int:
        if x != 0:
            if x < 0:
                kern_int = [-x for x in kern_int]
            break

    # Verify on held-out terms
    print(f"  Verifying on held-out terms {training_end}..{len(seq)-1}...", flush=True)
    bad = []
    for n in range(training_end, len(seq)):
        res = Fraction(0)
        for c_val, (j, m) in zip(kern_int, labels):
            if c_val != 0 and n - j >= 0 and n - j < len(seq) and seq[n-j] is not None:
                res += Fraction(c_val) * Fraction(n-j)**m * seq[n-j]
        if res != 0:
            bad.append(n)
    if bad:
        print(f"  FAILED on terms: {bad[:5]}")
        return None, -1
    print(f"  All held-out terms verified!", flush=True)

    return list(zip(labels, kern_int)), kernel_dim


# Try order 5, z-degree 11
print("\n=== Euler operator reconstruction (order=5, z-degree=11) ===", flush=True)
result, kdim = find_euler_operator(qhat, 5, 11, training_end=140)

if result is None:
    print(f"Order 5, z-degree 11 failed (kdim={kdim}). Trying alternatives...")
    for order, zdeg in [(3, 13), (4, 12), (5, 12), (6, 11), (3, 11), (4, 11)]:
        print(f"\nTrying order={order}, z-degree={zdeg}...", flush=True)
        result, kdim = find_euler_operator(qhat, order, zdeg, training_end=140)
        if result is not None:
            print(f"SUCCESS with order={order}, z-degree={zdeg}")
            break

if result is None:
    print("Could not find Euler operator! Exiting.")
    sys.exit(1)

# Extract the Euler polynomials P_j(T)
euler_coeffs = {}
for (j, m), c in result:
    if c != 0:
        if j not in euler_coeffs:
            euler_coeffs[j] = {}
        euler_coeffs[j][m] = c

print("\nEuler operator: L = Σ z^j P_j(θ)")
max_order = max(m for (j, m), c in result if c != 0)
max_zdeg = max(j for (j, m), c in result if c != 0)
print(f"  Actual order (max m with nonzero coeff): {max_order}")
print(f"  Actual z-degree: {max_zdeg}")

for j in sorted(euler_coeffs.keys()):
    terms = []
    for m in sorted(euler_coeffs[j].keys()):
        terms.append(f"{euler_coeffs[j][m]}*T^{m}")
    print(f"  P_{j}(T) = {' + '.join(terms)}")

# Indicial polynomial (z=0 part)
print("\nIndicial polynomial P_0(T):")
if 0 in euler_coeffs:
    for m in sorted(euler_coeffs[0].keys()):
        print(f"  T^{m}: {euler_coeffs[0][m]}")

# ================================================================
# Part 2: Test whether D_n² satisfies L_5
# ================================================================
print("\n=== Testing D_n² against L_5 ===", flush=True)

def eval_euler(seq, coeffs_list, n):
    """Evaluate Σ c_{j,m} (n-j)^m seq[n-j]."""
    res = Fraction(0)
    for (j, m), c in coeffs_list:
        if c != 0 and 0 <= n - j < len(seq) and seq[n-j] is not None:
            res += Fraction(c) * Fraction(n-j)**m * seq[n-j]
    return res

bad_dn2 = []
for n in range(max_zdeg, min(NMAX - 5, 150)):
    res = eval_euler(D2, result, n)
    if res != 0:
        bad_dn2.append(n)
        if len(bad_dn2) <= 3:
            print(f"  FAIL at n={n}: residual = {res}")

if not bad_dn2:
    print(f"  D_n² satisfies L_5 for all tested n (up to {min(NMAX-5,150)})!")
    print("  => L_{3,D} (Brafman/Delannoy operator) is a RIGHT FACTOR of L_5")
    print("  => L_5 = Q_2 · L_{3,D}")
else:
    print(f"  D_n² does NOT satisfy L_5 ({len(bad_dn2)} failures)")

# ================================================================
# Part 3: Test whether D_n · Q_n(3) satisfies L_5
# ================================================================
print("\n=== Testing D_n · Q_n(3) (Legendre second-kind) against L_5 ===", flush=True)

# Q_n(3) = Legendre function of second kind at x=3
# Recurrence: (n+1)Q_{n+1}(3) = (2n+1)·3·Q_n(3) - n·Q_{n-1}(3)
# Q_0(3) = (1/2)log((3+1)/(3-1)) = (1/2)log(2) ... but this needs mpmath
# For exact test, compute as mpmath with very high precision
from mpmath import mp, mpf, log, nstr
mp.dps = 200

Qleg = [mpf(0)] * (NMAX + 1)
Qleg[0] = log(2) / 2
Qleg[1] = mpf(3) * Qleg[0] - 1
for n in range(1, NMAX):
    Qleg[n+1] = (mpf(6*n+3) * Qleg[n] - mpf(n) * Qleg[n-1]) / mpf(n+1)

DQ = [mpf(int(D[n])) * Qleg[n] for n in range(NMAX + 1)]

print("Testing D_n·Q_n(3) against L_5 (numerical, 200 dps)...")
for n in [20, 50, 80, 100, 130]:
    res = mpf(0)
    for (j, m), c in result:
        if c != 0 and 0 <= n - j < len(DQ):
            res += mpf(c) * mpf(n-j)**m * DQ[n-j]
    if abs(res) > 0:
        digits = -float(mp.log10(abs(res))) if abs(res) > mpf(10)**(-200) else 200
    else:
        digits = 200
    print(f"  n={n}: |residual| ≈ 10^{-digits:.0f}")

# ================================================================
# Part 4: If L_{3,D} is a right factor, apply L_{3,D} to Q̂_n
# to get g_n = [L_{3,D}(F)](z^n), check if g_n satisfies order-2 recurrence
# ================================================================
if not bad_dn2:
    print("\n=== Applying L_{3,D} to Q̂_n to get g_n ===", flush=True)
    # L_{3,D} = (2θ-3)θ² - z(2θ+1)(35θ²-9) + z²(2θ+1)(35θ²+70θ+26) - z³(2θ+5)(θ+1)²
    # S_0(T) = 2T³ - 3T²
    # S_1(T) = -(2T+1)(35T²-9) = -70T³ - 35T² + 18T + 9
    # S_2(T) = (2T+1)(35T²+70T+26) = 70T³ + 175T² + 122T + 26
    # S_3(T) = -(2T+5)(T+1)² = -2T³ - 9T² - 12T - 5

    def S0(T):
        return 2*T**3 - 3*T**2
    def S1(T):
        return -70*T**3 - 35*T**2 + 18*T + 9
    def S2(T):
        return 70*T**3 + 175*T**2 + 122*T + 26
    def S3(T):
        return -2*T**3 - 9*T**2 - 12*T - 5

    g = [None] * (NMAX + 1)
    for n in range(3, NMAX + 1):
        g[n] = S0(Fraction(n)) * qhat[n] + S1(Fraction(n-1)) * qhat[n-1] + \
               S2(Fraction(n-2)) * qhat[n-2] + S3(Fraction(n-3)) * qhat[n-3]

    # Check if g_n is identically zero
    all_zero = all(g[n] == 0 for n in range(3, NMAX + 1) if g[n] is not None)
    if all_zero:
        print("  g_n = L_{3,D}(Q̂_n) is IDENTICALLY ZERO!")
        print("  => Q̂_n ∈ ker(L_{3,D}), meaning Q̂_n satisfies the Delannoy ODE")
        print("  This would mean Q̂_n is a linear combination of D_n², D_n·Q_n(3), and the log² branch")
    else:
        print(f"  g_n is NOT identically zero. Sample values:")
        for n in [5, 10, 20, 30]:
            if g[n] is not None:
                print(f"    g_{n} = {float(g[n]):.6e}")

        # Check if g_n satisfies a second-order Euler recurrence
        print("\n  Looking for order-2 Euler recurrence for g_n...")
        # Try order 2, z-degree D for various D
        for zdeg2 in range(3, 12):
            res2, kd2 = find_euler_operator(g, 2, zdeg2, training_end=min(NMAX-20, 120))
            if res2 is not None:
                print(f"  Found order-2 Euler operator with z-degree {zdeg2}!")
                max_m2 = max(m for (j,m), c in res2 if c != 0)
                max_j2 = max(j for (j,m), c in res2 if c != 0)
                print(f"  Actual order: {max_m2}, actual z-degree: {max_j2}")
                break
        else:
            print("  No order-2 Euler recurrence found up to z-degree 11")

print("\nDone.")
