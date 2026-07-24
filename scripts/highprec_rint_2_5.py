#!/usr/bin/env python3
"""High-precision extraction of r_int(n) to settle rationality.

Strategy: compute at mp.dps=2000 up to N=160.
R converges to ~200 digits at n=130 (error ~(16/543)^130 ≈ 10^{-200}).
Then w = u1 - R*u2 has dominant mode cancelled to 200 digits.
Forward-iterate w to n=250 for r_int at n=130..240.
Recessive contamination at n=130: (0.47/16)^130 ≈ 10^{-200}. So 200-digit accuracy.
"""
from mpmath import mp, mpf, nstr
import time
mp.dps = 2000

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

# Instead of storing the full matrix, compute the first column only
# V(N+1) = M(N) * V(N), where V = first column of T

# We need TWO independent solutions: u1 and u2
# u1: V(0) = [1, 0, 0]
# u2: V(0) = [0, 1, 0]

N_max = 165

print(f"Computing two forward solutions with mp.dps={mp.dps} up to N={N_max}...")
t0 = time.time()

v1 = [mpf(1), mpf(0), mpf(0)]
v2 = [mpf(0), mpf(1), mpf(0)]
u1_vals = [v1[0]]  # u1(0)
u2_vals = [v2[0]]  # u2(0)

for N in range(N_max + 5):
    M = M_mat(N)
    # v_new = M * v
    def mat_vec(M, v):
        return [sum(M[i][j]*v[j] for j in range(3)) for i in range(3)]

    v1 = mat_vec(M, v1)
    v2 = mat_vec(M, v2)
    u1_vals.append(v1[0])
    u2_vals.append(v2[0])

    if N % 20 == 0:
        elapsed = time.time() - t0
        ndig = len(str(abs(int(v1[0])))) if abs(v1[0]) > 1 else 0
        print(f"  N={N}: |u1| ~ 10^{ndig}, elapsed {elapsed:.1f}s")

elapsed = time.time() - t0
print(f"Done in {elapsed:.1f}s")

# Compute R = u1(N_R)/u2(N_R) for various N_R to check convergence
print("\nDominant ratio R = u1(n)/u2(n) convergence:")
for n in [80, 90, 100, 110, 120, 130, 140, 150, 160]:
    if n < len(u2_vals) and abs(u2_vals[n]) > 0:
        R_n = u1_vals[n] / u2_vals[n]
        print(f"  n={n}: R = {nstr(R_n, 60)}")

# Use R at n=160 (converged to ~240+ digits)
R = u1_vals[160] / u2_vals[160]
print(f"\nUsing R = u1(160)/u2(160)")

# Compute w(n) = u1(n) - R*u2(n) for n=0..N_max
w = [u1_vals[n] - R * u2_vals[n] for n in range(len(u1_vals))]

# Compute r_int(n) = w(n+1)/w(n)
r_int = {}
for n in range(N_max + 2):
    if abs(w[n]) > mpf('1e-1500'):
        r_int[n] = w[n+1] / w[n]

# Print a few key values
print("\nr_int(n) for small n:")
for n in [0, 1, 2, 5, 10, 20]:
    if n in r_int:
        print(f"  r_int({n}) = {nstr(r_int[n], 50)}")

print("\nr_int(n)/(-16*n^7) for moderate n:")
for n in [20, 50, 80, 100, 130, 150]:
    if n in r_int:
        ratio = r_int[n] / (-16 * mpf(n)**7)
        print(f"  n={n}: {nstr(ratio, 30)}")

# Now do rational reconstruction with HIGH PRECISION data
# Use n=130..160 where accuracy is best
from mpmath import matrix as mp_matrix, lu_solve

print("\n\nRational reconstruction from n=130..160 data:")
for dQ in range(0, 15):
    dP = dQ + 7
    n_unk = dP + dQ

    data_points = list(range(130, min(130 + n_unk + 5, 162)))
    n_eq = len(data_points)
    if n_unk >= n_eq:
        continue

    mat = mp_matrix(n_eq, n_unk)
    rhs = mp_matrix(n_eq, 1)

    for idx, n_val in enumerate(data_points):
        if n_val not in r_int:
            continue
        rv = r_int[n_val]
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

        if logr < -100:
            verdict = "★★★ RATIONAL"
        elif logr < -50:
            verdict = "★ likely"
        elif logr < -20:
            verdict = "close"
        else:
            verdict = "no"

        print(f"  dP={dP:>3}, dQ={dQ:>3}: log10(resid) = {logr:>8.1f}  [{verdict}]")
    except Exception as e:
        print(f"  dP={dP:>3}, dQ={dQ:>3}: FAIL ({e})")
