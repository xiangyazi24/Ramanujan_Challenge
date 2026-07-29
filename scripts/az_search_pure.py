#!/usr/bin/env python3
"""Pure Python AZ certificate search for P2.5 (no Sage dependency).

Works over GF(p) using Python's native modular arithmetic.
Polynomials in (u,v,t) represented as dict: (i,j,k) -> coefficient.
"""
import sys
from collections import defaultdict

P = (1 << 61) - 1  # Mersenne prime

def mod(x):
    return x % P

def modinv(a):
    return pow(a, P - 2, P)

# Polynomial operations: poly is dict (exp_u, exp_v, exp_t) -> coeff mod P
def poly_zero():
    return {}

def poly_const(c):
    c = mod(c)
    return {(0,0,0): c} if c else {}

def poly_var(var):
    if var == 'u': return {(1,0,0): 1}
    if var == 'v': return {(0,1,0): 1}
    if var == 't': return {(0,0,1): 1}

def poly_add(a, b):
    r = dict(a)
    for k, v in b.items():
        r[k] = mod(r.get(k, 0) + v)
        if r[k] == 0: del r[k]
    return r

def poly_sub(a, b):
    r = dict(a)
    for k, v in b.items():
        r[k] = mod(r.get(k, 0) - v)
        if r[k] == 0: del r[k]
    return r

def poly_scale(a, c):
    c = mod(c)
    if c == 0: return {}
    return {k: mod(v * c) for k, v in a.items() if mod(v * c)}

def poly_mul(a, b):
    r = {}
    for (i1,j1,k1), c1 in a.items():
        for (i2,j2,k2), c2 in b.items():
            key = (i1+i2, j1+j2, k1+k2)
            r[key] = mod(r.get(key, 0) + c1 * c2)
            if r[key] == 0 and key in r: del r[key]
    return r

def poly_deriv(a, var):
    idx = {'u': 0, 'v': 1, 't': 2}[var]
    r = {}
    for exp, c in a.items():
        if exp[idx] == 0: continue
        new_exp = list(exp)
        coeff = mod(c * exp[idx])
        new_exp[idx] -= 1
        key = tuple(new_exp)
        r[key] = mod(r.get(key, 0) + coeff)
        if r[key] == 0 and key in r: del r[key]
    return r

def poly_neg(a):
    return {k: mod(-v) for k, v in a.items()}

# Build basis monomials
u = poly_var('u')
v = poly_var('v')
t = poly_var('t')
one = poly_const(1)

# u+v, u*v, u²+v², t(u+v), t², t
upv = poly_add(u, v)
uv = poly_mul(u, v)
u2pv2 = poly_add(poly_mul(u, u), poly_mul(v, v))
tupv = poly_mul(t, upv)
t2 = poly_mul(t, t)

# u-v, t(u-v), (u-v)(u+v)
umv = poly_sub(u, v)
tumv = poly_mul(t, umv)
umv_upv = poly_mul(umv, upv)

SYM_BASIS = [one, upv, uv, u2pv2, tupv, t2, t]
ANTI_BASIS = [umv, tumv, umv_upv]
A_BASIS = SYM_BASIS + ANTI_BASIS
A_SWAP_BASIS = SYM_BASIS + [poly_neg(m) for m in ANTI_BASIS]

# C basis: only t-dependent monomials
C_BASIS = [tupv, t2, t]

# Q = 1 - u²v²
u2v2 = poly_mul(uv, uv)
Q = poly_sub(one, u2v2)

# Matrix data
def p_data(n):
    n = mod(n)
    p11 = mod(136*pow(n,4,P) + 1424*pow(n,3,P) + 5548*pow(n,2,P) + 9551*n + 6141)
    p12 = mod(384*pow(n,6,P) + 6384*pow(n,5,P) + 44168*pow(n,4,P) + 162698*pow(n,3,P) + 336377*pow(n,2,P) + 369933*n + 169011)
    p13 = mod(480*pow(n,4,P) + 4980*pow(n,3,P) + 19210*pow(n,2,P) + 32690*n + 20730)
    p21 = mod(48*pow(n,3,P) + 386*pow(n,2,P) + 1017*n + 879)
    p22 = mod(272*pow(n,5,P) + 3848*pow(n,4,P) + 21732*pow(n,3,P) + 61184*pow(n,2,P) + 85761*n + 47808)
    p23 = mod(320*pow(n,3,P) + 2540*pow(n,2,P) + 6610*n + 5640)
    p31 = mod(32*pow(n,4,P) + 302*pow(n,3,P) + 1037*pow(n,2,P) + 1530*n + 813)
    p32 = mod(192*pow(n,6,P) + 2984*pow(n,5,P) + 19116*pow(n,4,P) + 64452*pow(n,3,P) + 120256*pow(n,2,P) + 117279*n + 46476)
    p33 = mod(16*pow(n,5,P) + 408*pow(n,4,P) + 2912*pow(n,3,P) + 8884*pow(n,2,P) + 12254*n + 6240)
    return p11, p12, p13, p21, p22, p23, p31, p32, p33

def M_H_entry(n, i, j):
    """Normalized matrix M_H(n)[i,j] over GF(P)."""
    n = mod(n)
    pp = p_data(n)
    # Signs: row 0 has (+, -, +), row 1 has (-, +, -), row 2 has (+, -, +)
    signs = [[1, P-1, 1], [P-1, 1, P-1], [1, P-1, 1]]
    # Denominators
    n2 = mod(n + 2)
    n3 = mod(n + 3)
    n27 = mod(2*n + 7)
    n25 = mod(2*n + 5)

    den_map = {
        (0,0): mod(2 * n2 * n2 % P * (n27 * n27 % P) % P),
        (0,1): mod(2 * n2 * n2 % P * (n3 * n3 % P) % P * n25 % P * (n27 * n27 % P) % P),
        (0,2): mod(2 * n2 * n2 % P * (n3 * n3 % P) % P * n25 % P * (n27 * n27 % P) % P),
        (1,0): mod(n27 * n27 % P),
        (1,1): mod(2 * (n3 * n3 % P) % P * n25 % P * (n27 * n27 % P) % P),
        (1,2): mod(2 * (n3 * n3 % P) % P * n25 % P * (n27 * n27 % P) % P),
        (2,0): mod(n27 * n27 % P),
        (2,1): mod(2 * (n3 * n3 % P) % P * n25 % P * (n27 * n27 % P) % P),
        (2,2): mod(2 * (n3 * n3 % P) % P * n25 % P * (n27 * n27 % P) % P),
    }

    idx = i * 3 + j
    num = mod(signs[i][j] * pp[idx])
    den = den_map[(i,j)]
    return mod(num * modinv(den))

def good_n(nn):
    n = mod(nn)
    n2 = mod(n+2); n3 = mod(n+3); n25 = mod(2*n+5); n27 = mod(2*n+7)
    return mod(n2 * n3 % P * n25 % P * n27 % P) != 0

def contribution_poly(kind, basis_idx, n_power, nn):
    """Compute the RHS contribution polynomial for one unknown."""
    z = mod(nn)
    scalar = pow(z, n_power, P)
    # Kpoly = (n+1) + (n-1)*u²v²
    Kpoly = poly_add(poly_const(mod(z + 1)), poly_scale(u2v2, mod(z - 1)))

    if kind == 'A':
        m = A_BASIS[basis_idx]
        ms = A_SWAP_BASIS[basis_idx]
        # inner = u*Q*∂_u(m) + v*Q*∂_v(ms) + Kpoly*(m+ms)
        du_m = poly_deriv(m, 'u')
        dv_ms = poly_deriv(ms, 'v')
        term1 = poly_mul(poly_mul(u, Q), du_m)
        term2 = poly_mul(poly_mul(v, Q), dv_ms)
        term3 = poly_mul(Kpoly, poly_add(m, ms))
        inner = poly_add(poly_add(term1, term2), term3)
        return poly_scale(poly_mul(Q, inner), scalar)

    if kind == 'C':
        m = C_BASIS[basis_idx]
        # inner = (1+t²)*∂_t(m)
        dt_m = poly_deriv(m, 't')
        t2p1 = poly_add(one, t2)
        inner = poly_mul(t2p1, dt_m)
        return poly_scale(poly_mul(Q, inner), scalar)

def cert_den(nn, mode):
    """Certificate denominator E(n)."""
    n = mod(nn)
    if mode == "one":
        return 1
    n1 = mod(n+1); n2 = mod(n+2); n3 = mod(n+3)
    n23 = mod(2*n+3); n25 = mod(2*n+5); n27 = mod(2*n+7)
    return mod(n1 * n2 % P * n3 % P * n23 % P * n25 % P * n27 % P)

def target_poly(nn, i, j, den_mode="one"):
    """Target: E(n) * (uv - Q) * M_H(n)[i,j]."""
    uv_minus_Q = poly_sub(uv, Q)
    scalar = mod(cert_den(nn, den_mode) * M_H_entry(nn, i, j) % P)
    return poly_scale(uv_minus_Q, scalar)

def build_and_solve(i, j, d, train_ns, den_mode="one"):
    """Build and solve the linear system for entry (i,j) at n-degree d."""
    n_A = len(A_BASIS)
    n_C = len(C_BASIS)
    labels = []
    for r in range(d + 1):
        for idx in range(n_A):
            labels.append(('A', idx, r))
        for idx in range(n_C):
            labels.append(('C', idx, r))

    n_unknowns = len(labels)

    # Collect all equations
    rows = []
    rhs_vals = []

    for nn in train_ns:
        # Contributions stay as-is; only target gets E(n)
        cols = [contribution_poly(k, idx, r, nn) for k, idx, r in labels]
        tgt = target_poly(nn, i, j, den_mode)

        # Collect all monomials
        all_monoms = set(tgt.keys())
        for col in cols:
            all_monoms.update(col.keys())

        for exp in sorted(all_monoms):
            row = [col.get(exp, 0) for col in cols]
            b = tgt.get(exp, 0)
            if any(c != 0 for c in row) or b != 0:
                rows.append(row)
                rhs_vals.append(b)

    n_eqs = len(rows)
    print(f"  entry ({i},{j}): {n_eqs} equations, {n_unknowns} unknowns", flush=True)

    # Gaussian elimination over GF(P)
    # Augmented matrix
    aug = [row + [b] for row, b in zip(rows, rhs_vals)]

    pivots = []
    pivot_col = 0
    for row_idx in range(len(aug)):
        if pivot_col >= n_unknowns:
            break
        # Find pivot
        found = -1
        for k in range(row_idx, len(aug)):
            if aug[k][pivot_col] != 0:
                found = k
                break
        if found == -1:
            pivot_col += 1
            continue
        # Swap
        aug[row_idx], aug[found] = aug[found], aug[row_idx]
        # Scale pivot row
        inv_piv = modinv(aug[row_idx][pivot_col])
        aug[row_idx] = [mod(x * inv_piv) for x in aug[row_idx]]
        # Eliminate
        for k in range(len(aug)):
            if k == row_idx: continue
            if aug[k][pivot_col] == 0: continue
            factor = aug[k][pivot_col]
            aug[k] = [mod(aug[k][c] - factor * aug[row_idx][c]) for c in range(len(aug[k]))]
        pivots.append((row_idx, pivot_col))
        pivot_col += 1

    rank = len(pivots)

    # Check consistency: any row with all-zero LHS but nonzero RHS?
    consistent = True
    for row in aug:
        if all(row[c] == 0 for c in range(n_unknowns)) and row[n_unknowns] != 0:
            consistent = False
            break

    affine_dim = n_unknowns - rank if consistent else None
    print(f"    rank={rank}, consistent={consistent}, affine_dim={affine_dim}", flush=True)

    if not consistent:
        return None

    # Extract solution (set free variables to 0)
    solution = [0] * n_unknowns
    for ri, ci in pivots:
        solution[ci] = aug[ri][n_unknowns]

    return solution

def verify(solution, labels, i, j, test_ns, den_mode="one"):
    for nn in test_ns:
        lhs = poly_zero()
        for coeff, label in zip(solution, labels):
            if coeff == 0: continue
            k, idx, r = label
            cp = contribution_poly(k, idx, r, nn)
            lhs = poly_add(lhs, poly_scale(cp, coeff))
        tgt = target_poly(nn, i, j, den_mode)
        diff = poly_sub(lhs, tgt)
        if diff:
            return False
    return True

def main():
    print(f"Prime: 2^61-1 = {P}", flush=True)
    print(f"A monomials: {len(A_BASIS)}, C monomials: {len(C_BASIS)}", flush=True)

    for den_mode in ["one", "squarefree"]:
        for d in range(6):
            n_labels = (d + 1) * (len(A_BASIS) + len(C_BASIS))
            den_deg = 0 if den_mode == "one" else 6
            n_train = max(den_deg + 8, d + 9) + 1
            n_holdout = 4

            train_ns = []
            k = 0
            while len(train_ns) < n_train + n_holdout:
                if good_n(k):
                    if den_mode == "squarefree":
                        if cert_den(k, den_mode) == 0:
                            k += 1
                            continue
                    train_ns.append(k)
                k += 1
            holdout = train_ns[n_train:]
            train_ns = train_ns[:n_train]

            print(f"\n{'='*60}", flush=True)
            print(f"den={den_mode}, d={d}, unknowns={n_labels}, train={len(train_ns)}", flush=True)

            all_ok = True
            for i in range(3):
                for j in range(3):
                    sol = build_and_solve(i, j, d, train_ns, den_mode)
                    if sol is None:
                        print(f"    INCONSISTENT at ({i},{j})", flush=True)
                        all_ok = False
                        break
                    labels = []
                    for r in range(d + 1):
                        for idx in range(len(A_BASIS)):
                            labels.append(('A', idx, r))
                        for idx in range(len(C_BASIS)):
                            labels.append(('C', idx, r))
                    ok = verify(sol, labels, i, j, holdout, den_mode)
                    print(f"    holdout: {'PASS' if ok else 'FAIL'}", flush=True)
                    if not ok:
                        all_ok = False
                        break
                if not all_ok:
                    break

            if all_ok:
                print(f"\n*** CERTIFICATE FOUND: den={den_mode}, d={d} ***", flush=True)
                return

    print("\nNo certificate found.", flush=True)

if __name__ == '__main__':
    main()
