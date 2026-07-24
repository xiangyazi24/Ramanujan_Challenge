#!/usr/bin/env python3
"""P2.5 AZ certificate search — pure Python, no Sage.
Implements the universal 3×3 matrix AZ equation from Q4875.

The equation:
  E(n)(uv-Q)M_H(n) = Q[uQ A_u + vQ B_v + (1+t²)C_t + (n+1+(n-1)u²v²)(A+B)]

where Q = 1-u²v², and A,B,C are 3×3 polynomial matrices in (u,v,t)
with polynomial-in-n coefficients.

All arithmetic is done over F_p, p = 2^61 - 1.
"""
import sys
from collections import defaultdict

P_MOD = (1 << 61) - 1  # 2305843009213693951

def modinv(a, p=P_MOD):
    return pow(a % p, p - 2, p)

# --- Multivariate polynomial over F_p ---
# Represented as dict: {(a,b,c): coeff} for u^a v^b t^c

def poly_zero():
    return {}

def poly_const(c):
    c = c % P_MOD
    if c == 0:
        return {}
    return {(0, 0, 0): c}

def poly_var(idx):
    e = [0, 0, 0]
    e[idx] = 1
    return {tuple(e): 1}

U_VAR, V_VAR, T_VAR = poly_var(0), poly_var(1), poly_var(2)

def poly_add(a, b):
    r = dict(a)
    for k, v in b.items():
        r[k] = (r.get(k, 0) + v) % P_MOD
        if r[k] == 0:
            del r[k]
    return r

def poly_sub(a, b):
    r = dict(a)
    for k, v in b.items():
        r[k] = (r.get(k, 0) - v) % P_MOD
        if r[k] == 0:
            del r[k]
    return r

def poly_scale(a, c):
    c = c % P_MOD
    if c == 0:
        return {}
    return {k: (v * c) % P_MOD for k, v in a.items() if (v * c) % P_MOD != 0}

def poly_mul(a, b):
    r = {}
    for ka, va in a.items():
        for kb, vb in b.items():
            k = (ka[0]+kb[0], ka[1]+kb[1], ka[2]+kb[2])
            r[k] = (r.get(k, 0) + va * vb) % P_MOD
    # Clean zeros
    return {k: v for k, v in r.items() if v != 0}

def poly_deriv(a, var_idx):
    r = {}
    for k, v in a.items():
        if k[var_idx] == 0:
            continue
        new_k = list(k)
        coeff = (v * k[var_idx]) % P_MOD
        new_k[var_idx] -= 1
        new_k = tuple(new_k)
        r[new_k] = (r.get(new_k, 0) + coeff) % P_MOD
    return {k: v for k, v in r.items() if v != 0}

# Q = 1 - u²v²
Q_POLY = poly_sub(poly_const(1), {(2, 2, 0): 1})
# suv = u²v²
SUV = {(2, 2, 0): 1}
# u*v
UV = {(1, 1, 0): 1}
# 1+t²
ONE_PLUS_T2 = poly_add(poly_const(1), {(0, 0, 2): 1})

# --- CMF matrix data ---
def p_data(z):
    """Polynomial building blocks at z (integer mod p)."""
    z = z % P_MOD
    p11 = (136*pow(z,4,P_MOD) + 1424*pow(z,3,P_MOD) + 5548*pow(z,2,P_MOD) + 9551*z + 6141) % P_MOD
    p12 = (384*pow(z,6,P_MOD) + 6384*pow(z,5,P_MOD) + 44168*pow(z,4,P_MOD) + 162698*pow(z,3,P_MOD)
           + 336377*pow(z,2,P_MOD) + 369933*z + 169011) % P_MOD
    p13 = (480*pow(z,4,P_MOD) + 4980*pow(z,3,P_MOD) + 19210*pow(z,2,P_MOD) + 32690*z + 20730) % P_MOD
    p21 = (48*pow(z,3,P_MOD) + 386*pow(z,2,P_MOD) + 1017*z + 879) % P_MOD
    p22 = (272*pow(z,5,P_MOD) + 3848*pow(z,4,P_MOD) + 21732*pow(z,3,P_MOD) + 61184*pow(z,2,P_MOD)
           + 85761*z + 47808) % P_MOD
    p23 = (320*pow(z,3,P_MOD) + 2540*pow(z,2,P_MOD) + 6610*z + 5640) % P_MOD
    p31 = (32*pow(z,4,P_MOD) + 302*pow(z,3,P_MOD) + 1037*pow(z,2,P_MOD) + 1530*z + 813) % P_MOD
    p32 = (192*pow(z,6,P_MOD) + 2984*pow(z,5,P_MOD) + 19116*pow(z,4,P_MOD) + 64452*pow(z,3,P_MOD)
           + 120256*pow(z,2,P_MOD) + 117279*z + 46476) % P_MOD
    p33 = (16*pow(z,5,P_MOD) + 408*pow(z,4,P_MOD) + 2912*pow(z,3,P_MOD) + 8884*pow(z,2,P_MOD)
           + 12254*z + 6240) % P_MOD
    return p11, p12, p13, p21, p22, p23, p31, p32, p33


def D_M(z):
    z = z % P_MOD
    return (2 * pow(z+2, 2, P_MOD) % P_MOD * pow(z+3, 2, P_MOD) % P_MOD
            * ((2*z+5) % P_MOD) % P_MOD * pow(2*z+7, 2, P_MOD) % P_MOD) % P_MOD


def normalized_matrix_mod(nn):
    """M_H(n) = M(n)/delta(n) where delta = -D_M, computed mod p."""
    z = nn % P_MOD
    p11, p12, p13, p21, p22, p23, p31, p32, p33 = p_data(z)

    z2 = (z + 2) % P_MOD
    z3 = (z + 3) % P_MOD
    w5 = (2*z + 5) % P_MOD
    w7 = (2*z + 7) % P_MOD
    z2sq = pow(z2, 2, P_MOD)
    z3sq = pow(z3, 2, P_MOD)
    w7sq = pow(w7, 2, P_MOD)

    d1 = modinv((2 * z2sq % P_MOD * w7sq % P_MOD) % P_MOD)
    d2 = modinv((2 * z3sq % P_MOD * w5 % P_MOD * w7sq % P_MOD) % P_MOD)
    d3 = modinv(w7sq)
    d4 = modinv((2 * z3sq % P_MOD * w5 % P_MOD * w7sq % P_MOD) % P_MOD)

    # From the Sage code (normalized_matrix):
    # Row 0: [p11/(2(z+2)²(2z+7)²), -p12/(2(z+2)²(z+3)²(2z+5)(2z+7)²), p13/(2(z+2)²(z+3)²(2z+5)(2z+7)²)]
    # Row 1: [-p21/(2z+7)², p22/(2(z+3)²(2z+5)(2z+7)²), -p23/(2(z+3)²(2z+5)(2z+7)²)]
    # Row 2: [p31/(2z+7)², -p32/(2(z+3)²(2z+5)(2z+7)²), p33/(2(z+3)²(2z+5)(2z+7)²)]

    inv_d_full = modinv((2 * z2sq % P_MOD * z3sq % P_MOD * w5 % P_MOD * w7sq % P_MOD) % P_MOD)

    M = [
        [(p11 * d1) % P_MOD, (P_MOD - p12 * inv_d_full % P_MOD) % P_MOD, (p13 * inv_d_full) % P_MOD],
        [(P_MOD - p21 * d3 % P_MOD) % P_MOD, (p22 * d4) % P_MOD, (P_MOD - p23 * d4 % P_MOD) % P_MOD],
        [(p31 * d3) % P_MOD, (P_MOD - p32 * d4 % P_MOD) % P_MOD, (p33 * d4) % P_MOD],
    ]
    return M


# --- Certificate denominator ---
def cert_den_value(nn, mode):
    z = nn % P_MOD
    if mode == "one":
        return 1
    z1 = (z + 1) % P_MOD
    z2 = (z + 2) % P_MOD
    z3 = (z + 3) % P_MOD
    w3 = (2*z + 3) % P_MOD
    w5 = (2*z + 5) % P_MOD
    w7 = (2*z + 7) % P_MOD
    if mode == "squarefree":
        return (z1 * z2 % P_MOD * z3 % P_MOD * w3 % P_MOD * w5 % P_MOD * w7 % P_MOD) % P_MOD
    if mode == "matrix":
        return (pow(z2, 2, P_MOD) * pow(z3, 2, P_MOD) % P_MOD * w5 % P_MOD * pow(w7, 2, P_MOD) % P_MOD) % P_MOD
    if mode == "orbit":
        return (z1 * w3 % P_MOD * pow(z2, 2, P_MOD) % P_MOD * pow(z3, 2, P_MOD) % P_MOD
                * w5 % P_MOD * pow(w7, 2, P_MOD) % P_MOD) % P_MOD
    raise ValueError(mode)


def cert_den_degree(mode):
    return {"one": 0, "squarefree": 6, "matrix": 7, "orbit": 9}[mode]


# --- Basis ---
# Minimal basis (parity-reduced):
# A_basis: [1, uv, u²+v², u²-v², t²]
# A_swap:  [1, uv, u²+v², -(u²-v²), t²]
# C_basis: [t]

A_BASIS = [
    poly_const(1),
    UV,
    poly_add({(2, 0, 0): 1}, {(0, 2, 0): 1}),  # u²+v²
    poly_sub({(2, 0, 0): 1}, {(0, 2, 0): 1}),   # u²-v²
    {(0, 0, 2): 1},  # t²
]

A_SWAP = [
    poly_const(1),
    UV,
    poly_add({(2, 0, 0): 1}, {(0, 2, 0): 1}),  # u²+v²
    poly_sub({(0, 2, 0): 1}, {(2, 0, 0): 1}),   # -(u²-v²) = v²-u²
    {(0, 0, 2): 1},  # t²
]

C_BASIS = [
    {(0, 0, 1): 1},  # t
]


def labels_for_degree(d):
    labels = []
    for r in range(d + 1):
        for j in range(len(A_BASIS)):
            labels.append(('A', j, r))
        for j in range(len(C_BASIS)):
            labels.append(('C', j, r))
    return labels


def contribution(label, nn):
    """Q × RHS contribution of one numerator coefficient."""
    kind, j, r = label
    z = nn % P_MOD
    scalar = pow(z, r, P_MOD)
    # Kpoly = (n+1) + (n-1)u²v²
    Kpoly = poly_add(poly_const((z + 1) % P_MOD), poly_scale(SUV, (z - 1) % P_MOD))

    if kind == 'A':
        m = A_BASIS[j]
        ms = A_SWAP[j]
        # inner = u*Q*m_u + v*Q*ms_v + Kpoly*(m+ms)
        m_u = poly_deriv(m, 0)
        ms_v = poly_deriv(ms, 1)
        term1 = poly_mul(poly_mul(U_VAR, Q_POLY), m_u)
        term2 = poly_mul(poly_mul(V_VAR, Q_POLY), ms_v)
        term3 = poly_mul(Kpoly, poly_add(m, ms))
        inner = poly_add(poly_add(term1, term2), term3)
        return poly_scale(poly_mul(Q_POLY, inner), scalar)

    if kind == 'C':
        m = C_BASIS[j]
        m_t = poly_deriv(m, 2)
        inner = poly_mul(ONE_PLUS_T2, m_t)
        return poly_scale(poly_mul(Q_POLY, inner), scalar)

    raise ValueError(kind)


def targets_at_n(nn, den_mode):
    """E(n)·(uv - Q)·M_H(n), returned as 9 polynomials (flattened 3×3)."""
    E = cert_den_value(nn, den_mode)
    MH = normalized_matrix_mod(nn)
    # uv - Q = uv - (1 - u²v²) = uv + u²v² - 1
    uv_minus_Q = poly_sub(UV, Q_POLY)  # = uv - 1 + u²v²

    targets = []
    for i in range(3):
        for j in range(3):
            targets.append(poly_scale(poly_mul(uv_minus_Q, poly_const(MH[i][j])), E))
    return targets


def good_n(nn, den_mode):
    z = nn % P_MOD
    return D_M(z) != 0 and cert_den_value(nn, den_mode) != 0


def first_good_ns(count, den_mode, start=0):
    ans = []
    k = start
    while len(ans) < count:
        if good_n(k, den_mode):
            ans.append(k)
        k += 1
    return ans


# --- Build and solve ---
def build_system(d, train_ns, den_mode):
    labels = labels_for_degree(d)
    ncols = len(labels)
    rows = []
    rhs_rows = []

    for ni, nn in enumerate(train_ns):
        cols = [contribution(label, nn) for label in labels]
        targets = targets_at_n(nn, den_mode)

        support = set()
        for f in cols:
            support.update(f.keys())
        for f in targets:
            support.update(f.keys())

        for exponent in sorted(support):
            row = [(f.get(exponent, 0)) % P_MOD for f in cols]
            b = [(f.get(exponent, 0)) % P_MOD for f in targets]
            if any(c != 0 for c in row) or any(c != 0 for c in b):
                rows.append(row)
                rhs_rows.append(b)

    return rows, rhs_rows, labels


def gauss_rank_and_solve(rows, rhs_rows):
    """Gaussian elimination over F_p. Returns (rank, aug_rank, solution_or_None).
    Solution is a list of 9 column vectors (one per RHS).
    """
    if not rows:
        return 0, 0, None

    m = len(rows)
    ncols = len(rows[0])
    nrhs = len(rhs_rows[0]) if rhs_rows else 0

    # Augmented matrix: [L | B]
    aug = []
    for i in range(m):
        aug.append(list(rows[i]) + list(rhs_rows[i]))

    total_cols = ncols + nrhs
    pivot_cols = []
    pivot_rows = []

    row_idx = 0
    for col in range(ncols):
        # Find pivot
        found = -1
        for r in range(row_idx, m):
            if aug[r][col] != 0:
                found = r
                break
        if found == -1:
            continue
        aug[row_idx], aug[found] = aug[found], aug[row_idx]
        pivot_cols.append(col)
        pivot_rows.append(row_idx)

        inv_piv = modinv(aug[row_idx][col])
        for j in range(total_cols):
            aug[row_idx][j] = (aug[row_idx][j] * inv_piv) % P_MOD

        for r in range(m):
            if r == row_idx:
                continue
            if aug[r][col] == 0:
                continue
            factor = aug[r][col]
            for j in range(total_cols):
                aug[r][j] = (aug[r][j] - factor * aug[row_idx][j]) % P_MOD

        row_idx += 1

    rank_L = len(pivot_cols)

    # Check augmented rank: look for rows where L part is zero but B part is not
    aug_rank = rank_L
    for r in range(row_idx, m):
        if any(aug[r][ncols + j] != 0 for j in range(nrhs)):
            aug_rank += 1

    consistent = (rank_L == aug_rank)

    if not consistent:
        return rank_L, aug_rank, None

    # Extract solution: for each RHS column
    solutions = []
    for rhs_col in range(nrhs):
        x = [0] * ncols
        for pi, (pr, pc) in enumerate(zip(pivot_rows, pivot_cols)):
            x[pc] = aug[pr][ncols + rhs_col]
        solutions.append(x)

    return rank_L, aug_rank, solutions


def degree_bound_after_clearing(d, den_mode):
    return max(cert_den_degree(den_mode) + 7, d + 8)


def solve_system(d, den_mode):
    proof_count = degree_bound_after_clearing(d, den_mode) + 1
    holdout_count = 8
    all_ns = first_good_ns(proof_count + holdout_count, den_mode)
    train_ns = all_ns[:proof_count]
    holdout_ns = all_ns[proof_count:]

    print(f"\n{'='*60}", flush=True)
    print(f"  den_mode={den_mode}, n-degree={d}", flush=True)
    print(f"  train n-values: {train_ns[:5]}...  (count={len(train_ns)})", flush=True)
    print(f"  holdout: {holdout_ns}", flush=True)

    labels = labels_for_degree(d)
    print(f"  unknowns/entry: {len(labels)}", flush=True)

    rows, rhs_rows, labels = build_system(d, train_ns, den_mode)
    print(f"  equations: {len(rows)}", flush=True)

    rank_L, aug_rank, solutions = gauss_rank_and_solve(rows, rhs_rows)

    consistent = (rank_L == aug_rank)
    affine_dim = len(labels) - rank_L if consistent else None

    print(f"  rank={rank_L}, aug_rank={aug_rank}, consistent={consistent}", flush=True)
    if consistent:
        print(f"  affine_dim={affine_dim}", flush=True)

    if not consistent:
        print(f"  INCONSISTENT — no certificate at this degree/denominator", flush=True)
        return None

    # Holdout verification
    print(f"  Verifying on holdout n-values...", flush=True)
    all_pass = True
    for nn in holdout_ns:
        cols = [contribution(label, nn) for label in labels]
        targets = targets_at_n(nn, den_mode)
        for e in range(9):
            got = {}
            for k in range(len(labels)):
                if solutions[e][k] != 0:
                    got = poly_add(got, poly_scale(cols[k], solutions[e][k]))
            # Compare got with targets[e]
            diff = poly_sub(got, targets[e])
            if diff:
                print(f"  HOLDOUT FAIL: n={nn}, entry={divmod(e,3)}", flush=True)
                all_pass = False
                break
        if not all_pass:
            break

    if all_pass:
        print(f"  HOLDOUT PASS ✓", flush=True)
        return {
            'd': d,
            'den_mode': den_mode,
            'labels': labels,
            'solutions': solutions,
            'rank': rank_L,
            'affine_dim': affine_dim,
        }
    else:
        print(f"  HOLDOUT FAIL — spurious solution", flush=True)
        return None


# --- Rational reconstruction ---
def rational_recon(a, p=P_MOD):
    """Rational reconstruction: find r/s with |r|,s < sqrt(p/2) such that r/s ≡ a (mod p)."""
    a = a % p
    if a == 0:
        return (0, 1)
    old_r, r = p, a
    old_s, s = 0, 1
    bound = int(p**0.5)
    while r > bound:
        q = old_r // r
        old_r, r = r, old_r - q * r
        old_s, s = s, old_s - q * s
    if s < 0:
        r, s = -r, -s
    return (r, s)


def display_solution(candidate):
    labels = candidate['labels']
    solutions = candidate['solutions']
    print(f"\n  Certificate coefficients (rational reconstruction):")
    for e in range(9):
        i, j = divmod(e, 3)
        nonzero = [(k, solutions[e][k]) for k in range(len(labels)) if solutions[e][k] != 0]
        if nonzero:
            print(f"  Entry ({i},{j}):")
            for k, val in nonzero:
                kind, idx, r = labels[k]
                rr, ss = rational_recon(val)
                basis_name = f"A[{idx}]" if kind == 'A' else f"C[{idx}]"
                if ss == 1:
                    print(f"    n^{r} * {basis_name}: {rr}")
                else:
                    print(f"    n^{r} * {basis_name}: {rr}/{ss}")


# --- Main ---
def main():
    print(f"P2.5 AZ certificate search")
    print(f"prime = 2^61-1 = {P_MOD}")
    print(f"A_basis size = {len(A_BASIS)}, C_basis size = {len(C_BASIS)}")
    print(f"slots at n-degree 0: {len(A_BASIS) + len(C_BASIS)} = 6")

    den_modes = ["one", "squarefree", "matrix", "orbit"]
    degrees = [0, 1, 2, 3, 4]

    for den_mode in den_modes:
        for d in degrees:
            candidate = solve_system(d, den_mode)
            if candidate is not None:
                display_solution(candidate)
                print(f"\n  *** CERTIFICATE FOUND: den_mode={den_mode}, n-degree={d} ***")
                return

    print(f"\nNo certificate found in the ansatz ladder.")
    print(f"Consider: larger n-degree, bigger basis, or different denominator.")


if __name__ == '__main__':
    main()
