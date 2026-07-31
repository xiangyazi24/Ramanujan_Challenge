#!/usr/bin/env python3
"""Check if r_int is rational by extending degree search and checking residual pattern."""
from mpmath import mp, mpf, nstr, matrix as mp_matrix, lu_solve
mp.dps = 400

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
r_int_vals = {}
for n in range(45):
    if abs(w[n]) > mpf('1e-300'):
        r_int_vals[n] = w[n+1] / w[n]

# Systematic degree scan
print("Degree scan: log10(relative residual) vs (dP, dQ)")
print(f"{'dQ':>4} {'dP':>4} {'n_unk':>6} {'log10(resid)':>14} {'verdict':>10}")
print("-" * 50)

for dQ in range(0, 16):
    dP = dQ + 7
    n_unk = dP + dQ
    if n_unk == 0:
        n_unk = 7
    
    data_points = list(range(min(n_unk + 10, 44)))
    n_eq = len(data_points)
    if n_unk >= n_eq:
        break
    
    mat = mp_matrix(n_eq, n_unk)
    rhs = mp_matrix(n_eq, 1)
    
    for idx, n_val in enumerate(data_points):
        if n_val not in r_int_vals:
            continue
        rv = r_int_vals[n_val]
        nv = mpf(n_val)
        col = 0
        for k in range(dP):
            mat[idx, col] = mpf(16) * nv**k
            col += 1
        for j in range(dQ):
            mat[idx, col] = rv * nv**j
            col += 1
        rhs[idx, 0] = -rv * nv**dQ - mpf(16) * nv**dP
    
    try:
        ATA = mat.T * mat
        ATb = mat.T * rhs
        sol = lu_solve(ATA, ATb)
        residual = mat * sol - rhs
        max_res = max(abs(residual[i, 0]) for i in range(n_eq))
        max_rhs = max(abs(rhs[i, 0]) for i in range(n_eq))
        rel_res = max_res / max_rhs if max_rhs > 0 else max_res
        logr = float(mp.log10(rel_res)) if rel_res > 0 else -999
        
        if logr < -50:
            verdict = "RATIONAL?"
        elif logr < -10:
            verdict = "close"
        else:
            verdict = "no"
        
        print(f"{dQ:>4} {dP:>4} {n_unk:>6} {logr:>14.1f} {verdict:>10}")
    except:
        print(f"{dQ:>4} {dP:>4} {n_unk:>6} {'FAIL':>14}")

# Also test: is r_int(n)/(-16*n^7 - 16*33/2*n^6) closer to a simpler rational function?
print("\n\nTesting: r_int(n) / (-16*(n+a1)...(n+a_k)) for Pochhammer guesses...")
# Try: does r_int(n) / (-16*(n+1)^2*(n+2)^2*(n+3)^3/n^0) look rational?

# Actually, let me instead check: take r_int values at n=15,...,40 (where recessive is negligible)
# and try reconstruction there only (higher precision)
print("\n\nReconstruction using n=10..40 only (recessive contribution < 10^{-15}):")
for dQ in range(0, 12):
    dP = dQ + 7
    n_unk = dP + dQ
    data_points = list(range(10, min(10 + n_unk + 10, 44)))
    n_eq = len(data_points)
    if n_unk >= n_eq:
        break
    
    mat = mp_matrix(n_eq, n_unk)
    rhs = mp_matrix(n_eq, 1)
    
    for idx, n_val in enumerate(data_points):
        if n_val not in r_int_vals:
            continue
        rv = r_int_vals[n_val]
        nv = mpf(n_val)
        col = 0
        for k in range(dP):
            mat[idx, col] = mpf(16) * nv**k
            col += 1
        for j in range(dQ):
            mat[idx, col] = rv * nv**j
            col += 1
        rhs[idx, 0] = -rv * nv**dQ - mpf(16) * nv**dP
    
    try:
        ATA = mat.T * mat
        ATb = mat.T * rhs
        sol = lu_solve(ATA, ATb)
        residual = mat * sol - rhs
        max_res = max(abs(residual[i, 0]) for i in range(n_eq))
        max_rhs = max(abs(rhs[i, 0]) for i in range(n_eq))
        rel_res = max_res / max_rhs if max_rhs > 0 else max_res
        logr = float(mp.log10(rel_res)) if rel_res > 0 else -999
        
        if logr < -50:
            verdict = "RATIONAL?"
        elif logr < -10:
            verdict = "close"
        else:
            verdict = "no"
        
        print(f"  dQ={dQ}, dP={dP}: log10(resid) = {logr:.1f}  [{verdict}]")
    except:
        print(f"  dQ={dQ}, dP={dP}: FAIL")

