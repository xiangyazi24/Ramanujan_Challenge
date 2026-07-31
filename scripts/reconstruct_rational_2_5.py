#!/usr/bin/env python3
"""Rational function reconstruction from integer evaluations of r_int(n).

Uses the numerical values of r_int(n) computed via dominant-mode cancellation.
Tries Padé-type reconstruction for various degree pairs.
"""
from mpmath import mp, mpf, nstr, matrix as mp_matrix, lu_solve, log10
mp.dps = 600

# Recompute r_int values (same as extract_rint_2_5.py but stored)
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
    return [[m11,m12,m13],[m21,m22,m23],[m31,m32,m33]]

N_max = 55
T = [[mpf(1) if i==j else mpf(0) for j in range(3)] for i in range(3)]
u = [[] for _ in range(3)]

for N in range(N_max + 3):
    for i in range(3):
        u[i].append(T[i][0])
    if N < N_max + 2:
        M = M_mat(N)
        T_new = [[mpf(0)]*3 for _ in range(3)]
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    T_new[i][j] += T[i][k] * M[k][j]
        T = T_new

R = u[0][50] / u[1][50]
w = [u[0][n] - R * u[1][n] for n in range(len(u[0]))]

# Compute r_int values
r_int_vals = {}
for n in range(45):
    if abs(w[n]) > mpf('1e-500'):
        r_int_vals[n] = w[n+1] / w[n]

print(f"Computed r_int at {len(r_int_vals)} points")

# Precision check: at n=0, how many digits are reliable?
# w[0] = u0[0] - R*u1[0] = 1 - R*0 = 1 (if u1[0] = 0)
# w[1] = u0[1] - R*u1[1]
print(f"w[0] = {nstr(w[0], 20)}")
print(f"w[1] = {nstr(w[1], 20)}")
print(f"r_int(0) = {nstr(r_int_vals[0], 40)}")

# Rational function reconstruction: r_int(n) * Q(n) + 16*P(n) = 0
# i.e., r_int(n) * Q(n) = -16*P(n)
# P monic of degree dP, Q monic of degree dQ = dP - 7

# For a given (dP, dQ), set up: Q(n)*r_int(n) + 16*P(n) = 0
# Q(n) = n^dQ + q_{dQ-1}*n^{dQ-1} + ... + q_0
# P(n) = n^dP + p_{dP-1}*n^{dP-1} + ... + p_0
# Unknowns: p_0,...,p_{dP-1}, q_0,...,q_{dQ-1} = dP + dQ unknowns
# From the leading terms: -16 * n^dP = r_int(n) * n^dQ, which gives -16*n^7 ≈ r_int(n)/n^dQ
# This means the leading coefficient of P is 1 and of Q is 1

# At each data point n_i:
# r_int(n_i) * (n_i^dQ + q_{dQ-1}*n_i^{dQ-1} + ... + q_0) + 16*(n_i^dP + p_{dP-1}*n_i^{dP-1} + ... + p_0) = 0
# i.e., sum_j q_j * r_int(n_i) * n_i^j + sum_k p_k * 16 * n_i^k = -r_int(n_i)*n_i^dQ - 16*n_i^dP

for dQ in range(0, 8):
    dP = dQ + 7
    n_unknowns = dP + dQ  # (dP-1 + 1) + (dQ-1 + 1) = dP + dQ free coefficients
    if n_unknowns == 0:
        # dP=7, dQ=0: r_int(n) = -16*P(n) where P = n^7 + p6*n^6 + ... + p0
        # 7 unknowns
        n_unknowns = 7

    # Use overdetermined system with data points n=0,...,n_unknowns+5
    data_points = list(range(min(n_unknowns + 10, 45)))
    n_eq = len(data_points)

    # Build matrix
    mat = mp_matrix(n_eq, n_unknowns)
    rhs = mp_matrix(n_eq, 1)

    for idx, n_val in enumerate(data_points):
        if n_val not in r_int_vals:
            continue
        rv = r_int_vals[n_val]
        nv = mpf(n_val)

        # Columns: p_0, p_1, ..., p_{dP-1}, q_0, q_1, ..., q_{dQ-1}
        col = 0
        for k in range(dP):
            mat[idx, col] = mpf(16) * nv**k
            col += 1
        for j in range(dQ):
            mat[idx, col] = rv * nv**j
            col += 1

        rhs[idx, 0] = -rv * nv**dQ - mpf(16) * nv**dP

    # Solve via least squares (pseudoinverse)
    # A^T A x = A^T b
    ATA = mat.T * mat
    ATb = mat.T * rhs

    try:
        sol = lu_solve(ATA, ATb)

        # Compute residual
        residual = mat * sol - rhs
        max_res = max(abs(residual[i, 0]) for i in range(n_eq))
        max_rhs = max(abs(rhs[i, 0]) for i in range(n_eq))
        rel_res = max_res / max_rhs if max_rhs > 0 else max_res

        if rel_res < mpf('1e-10'):
            status = "✓✓✓ MATCH"
        elif rel_res < mpf('1e-3'):
            status = "close"
        else:
            status = "no fit"

        print(f"dP={dP}, dQ={dQ}: max|residual|/max|rhs| = {nstr(rel_res, 5)}  [{status}]")

        if rel_res < mpf('1e-10'):
            # Print the coefficients
            P_coeffs = [sol[k, 0] for k in range(dP)] + [mpf(1)]  # p_0,...,p_{dP-1}, 1
            Q_coeffs = [sol[dP + j, 0] for j in range(dQ)] + [mpf(1)] if dQ > 0 else [mpf(1)]
            
            print(f"  P(n) = n^{dP} + {' + '.join(nstr(P_coeffs[k],15) + '*n^' + str(k) for k in range(dP-1,-1,-1))}")
            if dQ > 0:
                print(f"  Q(n) = n^{dQ} + {' + '.join(nstr(Q_coeffs[j],15) + '*n^' + str(j) for j in range(dQ-1,-1,-1))}")
            
            # Check if coefficients are close to simple rationals
            print(f"  P coefficients (checking for half-integers):")
            for k in range(dP):
                val = P_coeffs[k]
                # Check if 2*val is close to an integer
                val2 = val * 2
                nearest_int = round(float(val2))
                err = abs(val2 - nearest_int)
                print(f"    p_{k} = {nstr(val, 20)}  2*p_{k} ≈ {nearest_int} (err {nstr(err, 5)})")
            
            break
    except Exception as e:
        print(f"dP={dP}, dQ={dQ}: solve failed ({e})")

