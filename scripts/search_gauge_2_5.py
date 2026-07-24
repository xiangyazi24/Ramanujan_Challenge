#!/usr/bin/env python3
"""Systematic search for gauge r(n) = -16*P(n)/Q(n) with various Q denominators.

Q can only have roots at n = -7/2, -9/2, -11/2 (from rational root of c3 and shifts).
"""
from fractions import Fraction
from mpmath import mp, mpf, nstr, matrix as mp_matrix, lu_solve
mp.dps = 150

# Rebuild exact recurrence
def M_mat_qq(n):
    n = Fraction(n)
    m11 = (-2*n-5)*(n+3)**2 * (136*n**4 + 1424*n**3 + 5548*n**2 + 9551*n + 6141)
    m12 = 384*n**6 + 6384*n**5 + 44168*n**4 + 162698*n**3 + 336377*n**2 + 369933*n + 169011
    m13 = -480*n**4 - 4980*n**3 - 19210*n**2 - 32690*n - 20730
    m21 = (n+2)**2*(n+3)**2*(4*n+10)*(48*n**3 + 386*n**2 + 1017*n + 879)
    m22 = (n+2)**2*(-272*n**5 - 3848*n**4 - 21732*n**3 - 61184*n**2 - 85761*n - 47808)
    m23 = (n+2)**2*(320*n**3 + 2540*n**2 + 6610*n + 5640)
    m31 = (-4*n-10)*(n+2)**2*(n+3)**2*(32*n**4 + 302*n**3 + 1037*n**2 + 1530*n + 813)
    m32 = (n+2)**2*(192*n**6 + 2984*n**5 + 19116*n**4 + 64452*n**3 + 120256*n**2 + 117279*n + 46476)
    m33 = (n+2)**2*(-16*n**5 - 408*n**4 - 2912*n**3 - 8884*n**2 - 12254*n - 6240)
    return [[m11,m12,m13],[m21,m22,m23],[m31,m32,m33]]

A_row = [Fraction(33750), Fraction(-36000), Fraction(9000)]
T = [[Fraction(1 if i==j else 0) for j in range(3)] for i in range(3)]
q_vals = []
for N in range(120):
    q = sum(A_row[k] * T[k][0] for k in range(3))
    q_vals.append(q)
    T_new = [[Fraction(0)]*3 for _ in range(3)]
    M = M_mat_qq(N)
    for i in range(3):
        for j in range(3):
            for k in range(3):
                T_new[i][j] += T[i][k] * M[k][j]
    T = T_new

# Gaussian elimination for recurrence
order = 3; degs = [28, 21, 14, 7]
n_unknowns = sum(d+1 for d in degs)
rows = []
for N in range(n_unknowns + 10):
    row = []
    for j in range(order + 1):
        for k in range(degs[j] + 1):
            row.append(Fraction(N)**k * q_vals[N + j])
    rows.append(row)
mat_g = [list(row) for row in rows]
n_rows = len(mat_g); n_cols = len(mat_g[0])
pivot_cols = []
for col in range(n_cols):
    pivot_row = None
    for row in range(len(pivot_cols), n_rows):
        if mat_g[row][col] != 0:
            pivot_row = row; break
    if pivot_row is None: continue
    mat_g[pivot_row], mat_g[len(pivot_cols)] = mat_g[len(pivot_cols)], mat_g[pivot_row]
    pr = len(pivot_cols); pivot_cols.append(col)
    piv = mat_g[pr][col]
    for row in range(n_rows):
        if row == pr: continue
        if mat_g[row][col] != 0:
            f = mat_g[row][col] / piv
            for c in range(n_cols):
                mat_g[row][c] -= f * mat_g[pr][c]
free_cols = [c for c in range(n_cols) if c not in pivot_cols]
null_vec = [Fraction(0)] * n_cols
null_vec[free_cols[0]] = Fraction(1)
for i in range(len(pivot_cols) - 1, -1, -1):
    pc = pivot_cols[i]
    val = -sum(mat_g[i][c] * null_vec[c] for c in range(n_cols) if c != pc) / mat_g[i][pc]
    null_vec[pc] = val
idx = 0; polys_qq = []
for j in range(order + 1):
    coeffs = null_vec[idx:idx + degs[j] + 1]
    idx += degs[j] + 1
    polys_qq.append(coeffs)

# Convert to mpf
polys_mp = [[mpf(c.numerator)/mpf(c.denominator) for c in p] for p in polys_qq]

def eval_poly_mp(coeffs, x):
    val = mpf(0)
    for c in reversed(coeffs):
        val = val * x + c
    return val

def c_j(j, n):
    return eval_poly_mp(polys_mp[j], n)

# Functional equation:
# F(n) = c3(n)*r(n)*r(n+1)*r(n+2) + c2(n)*r(n)*r(n+1) + c1(n)*r(n) + c0(n) = 0
# where r(n) = -16 * P(n) / Q(n)
# Q(n) is fixed (enumerated), P(n) = n^{dQ+7} + s1*n^{dQ+6} + ... + s_{dQ+7}

# Try each Q candidate
Q_candidates = [
    ("1", []),
    ("(2n+7)", [Fraction(-7,2)]),
    ("(2n+9)", [Fraction(-9,2)]),
    ("(2n+11)", [Fraction(-11,2)]),
    ("(2n+7)(2n+9)", [Fraction(-7,2), Fraction(-9,2)]),
    ("(2n+7)(2n+11)", [Fraction(-7,2), Fraction(-11,2)]),
    ("(2n+9)(2n+11)", [Fraction(-9,2), Fraction(-11,2)]),
    ("(2n+7)(2n+9)(2n+11)", [Fraction(-7,2), Fraction(-9,2), Fraction(-11,2)]),
    ("(2n+7)^2", [Fraction(-7,2), Fraction(-7,2)]),
    ("(2n+7)^2(2n+9)", [Fraction(-7,2), Fraction(-7,2), Fraction(-9,2)]),
    ("(2n+7)^2(2n+11)", [Fraction(-7,2), Fraction(-7,2), Fraction(-11,2)]),
    ("(2n+7)^3", [Fraction(-7,2), Fraction(-7,2), Fraction(-7,2)]),
]

def eval_Q(roots, n):
    val = mpf(1)
    for r in roots:
        val *= (n - mpf(r.numerator)/mpf(r.denominator))
    return val

for Q_name, Q_roots in Q_candidates:
    dQ = len(Q_roots)
    dP = dQ + 7
    n_params = dP  # P is monic of degree dP, so dP free coefficients

    # Use n_params evaluation points
    eval_points = list(range(n_params + 5))  # overdetermined

    def F_val(n_val, sigma):
        P_val = mpf(1)
        nv = mpf(n_val)
        for i in range(dP):
            P_val = P_val * nv + sigma[i]

        Q_val = eval_Q(Q_roots, nv)
        Q_val1 = eval_Q(Q_roots, nv+1)
        Q_val2 = eval_Q(Q_roots, nv+2)

        r0 = mpf(-16) * P_val / Q_val if Q_val != 0 else mpf('inf')

        P_val1 = mpf(1)
        for i in range(dP):
            P_val1 = P_val1 * (nv+1) + sigma[i]
        P_val2 = mpf(1)
        for i in range(dP):
            P_val2 = P_val2 * (nv+2) + sigma[i]

        r1 = mpf(-16) * P_val1 / Q_val1 if Q_val1 != 0 else mpf('inf')
        r2 = mpf(-16) * P_val2 / Q_val2 if Q_val2 != 0 else mpf('inf')

        return c_j(3, nv)*r0*r1*r2 + c_j(2, nv)*r0*r1 + c_j(1, nv)*r0 + c_j(0, nv)

    # Use half-integer initial guess: roots spread around 1,...,dP/2
    sigma0 = [mpf(0)] * dP
    # Start with simple guess: P(n) = (n+1)(n+3/2)(n+2)(n+5/2)(n+3)(n+7/2)(n+4)... extended
    from itertools import combinations
    from functools import reduce
    import operator
    guess_roots = [mpf(k)/2 for k in range(2, 2+2*dP, 2)][:dP]  # 1, 2, 3, ...
    # Adjust to half-integers for variety
    if dP <= 7:
        guess_roots = [1, 1.5, 2, 2.5, 3, 3.5, 4][:dP]
    else:
        guess_roots = [0.5 + k*0.5 for k in range(dP)]

    # Compute sigma from roots
    def roots_to_sigma(roots):
        n = len(roots)
        sigma = []
        for k in range(1, n+1):
            s = mpf(0)
            for combo in combinations(range(n), k):
                prod = mpf(1)
                for idx in combo:
                    prod *= roots[idx]
                s += prod
            sigma.append(s)
        return sigma

    sigma = roots_to_sigma(guess_roots)

    # Skip Q candidates that cause division by zero in eval points
    skip = False
    for n_val in eval_points[:n_params]:
        Q_val = eval_Q(Q_roots, mpf(n_val))
        if abs(Q_val) < mpf('1e-50'):
            skip = True
            break
    if skip:
        # Shift evaluation points
        eval_points = [n + 10 for n in eval_points]

    # Newton iteration
    converged = False
    for iteration in range(50):
        # Evaluate F at n_params points
        F_vals = []
        for n_val in eval_points[:n_params]:
            try:
                fv = F_val(n_val, sigma)
                F_vals.append(fv)
            except:
                break
        if len(F_vals) < n_params:
            break

        max_res = max(abs(fv) for fv in F_vals)
        if max_res < mpf('1e-100'):
            converged = True
            break

        # Jacobian
        eps = mpf('1e-50')
        J = mp_matrix(n_params, n_params)
        for j_col in range(n_params):
            sp = list(sigma); sm = list(sigma)
            sp[j_col] += eps; sm[j_col] -= eps
            for i_row in range(n_params):
                fp = F_val(eval_points[i_row], sp)
                fm = F_val(eval_points[i_row], sm)
                J[i_row, j_col] = (fp - fm) / (2*eps)

        b = mp_matrix(n_params, 1)
        for i in range(n_params):
            b[i, 0] = -F_vals[i]

        try:
            delta = lu_solve(J, b)
            sigma = [sigma[j] + delta[j, 0] for j in range(n_params)]
        except:
            break

    if converged:
        # Verify at MORE points
        max_verify = mpf(0)
        for n_val in range(30):
            Q_val = eval_Q(Q_roots, mpf(n_val))
            if abs(Q_val) > mpf('1e-50'):
                fv = abs(F_val(n_val, sigma))
                if fv > max_verify:
                    max_verify = fv

        if max_verify < mpf('1e-50'):
            print(f"\n✓✓✓ FOUND GAUGE with Q = {Q_name} ✓✓✓")
            print(f"  P coefficients (below leading n^{dP}):")
            for i, s in enumerate(sigma):
                print(f"    σ_{i+1} = {nstr(s, 30)}")

            # Find roots of P
            import numpy as np
            coeffs_np = [1.0] + [float(s) for s in sigma]
            roots_P = np.roots(coeffs_np)
            print(f"  P roots:")
            for i, r in enumerate(sorted(roots_P, key=lambda x: x.real)):
                if abs(r.imag) < 1e-8:
                    # Check half-integer
                    x = r.real
                    best_frac = None; best_err = 1e10
                    for num in range(-15, 15):
                        for den in [1, 2, 3, 4, 5, 6]:
                            err = abs(x + num/den)
                            if err < best_err:
                                best_frac = f"-{num}/{den}" if den > 1 else str(-num)
                                best_err = err
                    print(f"    a_{i+1} = {x:.12f}  ≈ {best_frac} (err {best_err:.2e})")
                else:
                    print(f"    a_{i+1} = {r.real:.12f} + {r.imag:.12f}i  (COMPLEX)")

            # Print r(n) formula
            print(f"\n  r(n) = -16 * P(n) / ({Q_name})")
            break
        else:
            print(f"  Q={Q_name}: Newton converged but verification FAILED (max F = {nstr(max_verify, 5)})")
    else:
        status = f"max|F|={nstr(max(abs(F_val(n_val, sigma)) for n_val in eval_points[:3]), 5)}" if n_params > 0 else "skip"
        print(f"  Q={Q_name}: Newton did NOT converge ({status})")
else:
    print("\nNo gauge found with tested Q candidates.")
