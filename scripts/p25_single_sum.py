#!/usr/bin/env python3
"""P2.5: Search for Q̂_n = Σ_k R(n,k) F_D(n,k) with harmonic terms.

F_D(n,k) = 2^k C(2k,k) C(n,k) C(n+k,k) is the Delannoy summand.
R(n,k) may involve harmonic numbers, odd partial sums, etc.
"""
from fractions import Fraction
import math

# CMF column-0 Q values
def M_entries(n):
    n = Fraction(n)
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
    n = Fraction(n)
    return Fraction(-2)*(n+2)**2*(n+3)**2*(2*n+5)*(2*n+7)**2

def MH_at(n):
    M = M_entries(n)
    d = delta_H(n)
    return [[M[i][j] / d for j in range(3)] for i in range(3)]

NMAX = 30

print("Computing CMF Q̂_n values...", flush=True)
q_row = [Fraction(33750), Fraction(-36000), Fraction(9000)]
cmf_q = []
for N in range(NMAX):
    cmf_q.append(q_row[0])
    MH = MH_at(N)
    q_new = [sum(q_row[i] * MH[i][j] for i in range(3)) for j in range(3)]
    q_row = q_new
print(f"  {len(cmf_q)} values", flush=True)

# Delannoy summand
def F_D(n, k):
    return Fraction(2)**k * Fraction(math.comb(2*k, k)) * Fraction(math.comb(n, k)) * Fraction(math.comb(n+k, k))

# D_n^2 check
for n in range(5):
    dn2 = sum(F_D(n, k) for k in range(n+1))
    P_n = [1, 3, 13, 63, 321]
    print(f"  D_{n}^2 = {dn2}, P_{n}(3)^2 = {P_n[n]**2}, match = {dn2 == P_n[n]**2}")

# Harmonic functions
def H(n):
    """H_n = Σ_{j=1}^n 1/j"""
    return sum(Fraction(1, j) for j in range(1, n+1))

def H2(n):
    """H_n^{(2)} = Σ_{j=1}^n 1/j²"""
    return sum(Fraction(1, j**2) for j in range(1, n+1))

def H_nk(n, k):
    """H_n - H_k"""
    return sum(Fraction(1, j) for j in range(k+1, n+1))

def odd_partial(k):
    """Σ_{j=0}^{k-1} (-1)^j/(2j+1)²"""
    return sum(Fraction((-1)**j, (2*j+1)**2) for j in range(k))

def odd_H(k):
    """Σ_{j=0}^{k-1} 1/(2j+1) = H_{2k}/2 - H_k/2 + ... actually just sum"""
    return sum(Fraction(1, 2*j+1) for j in range(k))

def central_H(n, k):
    """Σ_{j=1}^k 1/(j(2j-1)) or similar"""
    if k == 0:
        return Fraction(0)
    return sum(Fraction(1, j * (2*j-1)) for j in range(1, k+1))

# Precompute basis functions
print("\nPrecomputing basis functions at (n,k) points...", flush=True)

# Basis functions for R(n,k): each is a function (n,k) → value
# We try: R(n,k) = Σ_i α_i(n) · f_i(n,k) where f_i are known functions
# and α_i(n) are polynomial in n

basis_funcs = {
    '1': lambda n, k: Fraction(1),
    'H_2k': lambda n, k: sum(Fraction(1, j) for j in range(1, 2*k+1)) if k > 0 else Fraction(0),
    'H_k': lambda n, k: H(k),
    'H_n': lambda n, k: H(n),
    'H_nk': lambda n, k: H_nk(n, k),
    'H2_k': lambda n, k: H2(k),
    'H2_n': lambda n, k: H2(n),
    'S_k': lambda n, k: odd_partial(k),
    'oddH_k': lambda n, k: odd_H(k),
    'k': lambda n, k: Fraction(k),
    '1/(2k+1)': lambda n, k: Fraction(1, 2*k+1),
    'H_n+k': lambda n, k: H(n + k),
}

# First, just compute Q̂_n / D_n² to understand the scale
print("\n=== Q̂_n / D_n² ===", flush=True)
for n in range(min(10, NMAX)):
    dn2 = sum(F_D(n, k) for k in range(n+1))
    ratio = cmf_q[n] / dn2 if dn2 != 0 else None
    print(f"  n={n}: Q̂={float(cmf_q[n]):.6e}, D²={dn2}, ratio={float(ratio):.6f}" if ratio else f"  n={n}: D²=0")

# Compute residuals: Q̂_n - D_n² · polynomial_in_n
# Maybe Q̂_n ≈ D_n² · (a + b·n) → check
print("\n=== Fitting Q̂_n / D_n² to polynomial ===", flush=True)
ratios = []
for n in range(NMAX):
    dn2 = sum(F_D(n, k) for k in range(n+1))
    if dn2 != 0:
        ratios.append(cmf_q[n] / dn2)
    else:
        ratios.append(None)

# Check if ratios are polynomial in n
# Compute finite differences
diffs = list(ratios[:15])
print("  Ratios:", [float(r) for r in diffs[:10] if r])
for order in range(1, 5):
    new_diffs = [diffs[i+1] - diffs[i] for i in range(len(diffs)-1) if diffs[i] and diffs[i+1]]
    diffs = new_diffs
    print(f"  Δ^{order}: {[float(d) for d in diffs[:8]]}")

# Search: Q̂_n = Σ_k F_D(n,k) · R(n,k)
# where R(n,k) = Σ c_i · n^a · f_i(n,k)
# Start with simple: R(n,k) = c0 + c1·H_k + c2·S_k + c3·k/(2k+1) + c4·n + ...
# with n-polynomial coefficients

print("\n=== Single-sum search ===", flush=True)

# Build basis at each (n, k) and compute expected sums
# For each basis function f_i and n-power p, compute Σ_k F_D(n,k)·n^p·f_i(n,k)
# Then solve for coefficients

n_basis = ['1', 'H_k', 'S_k', 'H_n', 'H2_k', 'oddH_k', 'k', 'H_nk', 'H_n+k']

for n_poly_deg in range(4):
    for nb in range(2, len(n_basis)+1):
        selected = n_basis[:nb]
        n_unknowns = nb * (n_poly_deg + 1)
        n_train = n_unknowns + 3
        n_holdout = 3

        if n_train + n_holdout > NMAX:
            break

        # Build system
        A_rows = []
        b_vec = []

        for n in range(n_train + n_holdout):
            row = []
            for fi, fname in enumerate(selected):
                f = basis_funcs[fname]
                for p in range(n_poly_deg + 1):
                    # Σ_k F_D(n,k) · n^p · f(n,k)
                    val = Fraction(0)
                    for k in range(n + 1):
                        val += F_D(n, k) * Fraction(n)**p * f(n, k)
                    row.append(val)
            A_rows.append(row)
            b_vec.append(cmf_q[n])

        # Solve
        m = n_train
        aug = [list(A_rows[i]) + [b_vec[i]] for i in range(m)]
        n_cols = n_unknowns

        pivot_cols = []
        row_idx = 0
        for col in range(n_cols):
            found = -1
            for rr in range(row_idx, m):
                if aug[rr][col] != 0:
                    found = rr
                    break
            if found == -1:
                continue
            aug[row_idx], aug[found] = aug[found], aug[row_idx]
            piv = aug[row_idx][col]
            for j2 in range(n_cols + 1):
                aug[row_idx][j2] /= piv
            for rr in range(m):
                if rr == row_idx:
                    continue
                if aug[rr][col] == 0:
                    continue
                factor = aug[rr][col]
                for j2 in range(n_cols + 1):
                    aug[rr][j2] -= factor * aug[row_idx][j2]
            pivot_cols.append(col)
            row_idx += 1

        rank = len(pivot_cols)
        if rank < n_unknowns:
            continue

        x = [Fraction(0)] * n_unknowns
        for pi, pc in enumerate(pivot_cols):
            x[pc] = aug[pi][n_cols]

        # Holdout
        all_ok = True
        for n in range(n_train, n_train + n_holdout):
            pred = Fraction(0)
            for fi, fname in enumerate(selected):
                f = basis_funcs[fname]
                for p in range(n_poly_deg + 1):
                    coeff = x[fi * (n_poly_deg + 1) + p]
                    val = Fraction(0)
                    for k in range(n + 1):
                        val += F_D(n, k) * Fraction(n)**p * f(n, k)
                    pred += coeff * val
            if pred != cmf_q[n]:
                all_ok = False
                break

        if all_ok:
            print(f"\n*** MATCH: basis={selected}, n-deg={n_poly_deg} ***", flush=True)
            for fi, fname in enumerate(selected):
                coeffs = [x[fi*(n_poly_deg+1)+p] for p in range(n_poly_deg+1)]
                nonzero = [(p, c) for p, c in enumerate(coeffs) if c != 0]
                if nonzero:
                    terms = [f"({c})·n^{p}" if p > 0 else f"{c}" for p, c in nonzero]
                    print(f"  coeff of {fname}: {' + '.join(terms)}")
            break
    else:
        continue
    break

print("\nDone.")
