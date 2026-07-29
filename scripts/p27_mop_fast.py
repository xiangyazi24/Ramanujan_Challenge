#!/usr/bin/env python3
"""P2.7: Fast MOP step-line Q(1) computation using mpmath.
Compute enough values to identify the polynomial degree of the
4-term recurrence coefficients.
"""
import mpmath as mp
mp.mp.dps = 300

def mu2_moment(k):
    return mp.mpf(1) / mp.mpf(k+1)**2

def mu3_moment(k):
    return mp.mpf(1) / mp.mpf(k+1)**3

def mop_Q_at_1(n2, n3):
    N = n2 + n3
    if N == 0:
        return mp.mpf(1)

    moment_rows = []
    for r in range(n2):
        row = [mu2_moment(r + c) for c in range(N)]
        moment_rows.append(row)
    for r in range(n3):
        row = [mu3_moment(r + c) for c in range(N)]
        moment_rows.append(row)

    Delta = mp.det(mp.matrix(moment_rows))
    if abs(Delta) < mp.mpf(10)**(-200):
        return None

    full_rows = []
    for r in range(n2):
        row = [mu2_moment(r + c) for c in range(N + 1)]
        full_rows.append(row)
    for r in range(n3):
        row = [mu3_moment(r + c) for c in range(N + 1)]
        full_rows.append(row)
    full_rows.append([mp.mpf(1)] * (N + 1))

    full_det = mp.det(mp.matrix(full_rows))
    return full_det / Delta

def step_line(max_N):
    indices = [(0, 0)]
    for m in range(1, max_N + 1):
        if m % 2 == 1:
            indices.append(((m+1)//2, (m-1)//2))
        else:
            indices.append((m//2, m//2))
    return indices

MAX_STEP = 40
print(f"Computing Q_N(1) for N=0..{MAX_STEP}...", flush=True)
step = step_line(MAX_STEP)
q_vals = []
for i, (n2, n3) in enumerate(step):
    N = n2 + n3
    val = mop_Q_at_1(n2, n3)
    q_vals.append(val)
    if i % 5 == 0 or i < 5:
        print(f"  N={N}: Q(1) = {mp.nstr(val, 15)}", flush=True)

# Filter out None values — find the longest contiguous non-None prefix
valid_count = 0
for v in q_vals:
    if v is None:
        break
    valid_count += 1
q_vals = q_vals[:valid_count]
print(f"\nComputed {len(q_vals)} valid values.", flush=True)

if len(q_vals) < 7:
    print("Not enough valid values for recurrence search.")
    import sys; sys.exit(0)

# Now find the recurrence degree
# 4-term recurrence: Σ_{j=0}^3 P_j(n) q[n+j] = 0
# where P_j has degree d in n.
# Total unknowns: 4(d+1), normalize one → 4(d+1)-1
# Equations: len(q_vals) - 3

for deg in range(4, 15):
    nparams = 4 * (deg + 1) - 1
    neq = len(q_vals) - 3
    if neq < nparams + 3:
        print(f"  degree {deg}: need {nparams} params, only {neq} equations. Not enough.")
        continue

    # Build system
    A_mat = []
    b_vec = []
    for n in range(neq):
        row = []
        for j in range(4):
            for k in range(deg + 1):
                if j == 3 and k == deg:
                    continue
                row.append(mp.mpf(n)**k * q_vals[n + j])
        A_mat.append(row)
        b_vec.append(-mp.mpf(n)**deg * q_vals[n + 3])

    # Solve first nparams equations
    A_solve = mp.matrix([r[:nparams] for r in A_mat[:nparams]])
    b_solve = mp.matrix([b_vec[i] for i in range(nparams)])

    try:
        sol = mp.lu_solve(A_solve, b_solve)
    except:
        print(f"  degree {deg}: singular system")
        continue

    # Check residual on remaining equations
    max_res = mp.mpf(0)
    for n in range(nparams, min(neq, nparams + 5)):
        res = b_vec[n]
        for i in range(nparams):
            res -= A_mat[n][i] * sol[i]
        max_res = max(max_res, abs(res))

    print(f"  degree {deg}: max verification residual = {mp.nstr(max_res, 6)}", flush=True)

    if max_res < mp.mpf(10)**(-100):
        print(f"  *** FOUND: polynomial degree {deg} ***")

        # Extract Poincaré polynomial
        # Leading coefficient of each P_j
        coeffs = {}
        idx = 0
        for j in range(4):
            for k in range(deg + 1):
                if j == 3 and k == deg:
                    coeffs[(j, k)] = mp.mpf(1)
                else:
                    coeffs[(j, k)] = sol[idx]
                    idx += 1

        # Poincaré: lc(P_3) λ³ + lc(P_2) λ² + lc(P_1) λ + lc(P_0) = 0
        lc = [coeffs[(j, deg)] for j in range(4)]
        print(f"  Leading coefficients: {[mp.nstr(c, 10) for c in lc]}")

        # Characteristic polynomial: lc[3]λ³ + lc[2]λ² + lc[1]λ + lc[0]
        # Divide by lc[3]:
        char_poly = [lc[j]/lc[3] for j in range(4)]
        print(f"  Poincaré: λ³ + {mp.nstr(char_poly[2], 10)}λ² + {mp.nstr(char_poly[1], 10)}λ + {mp.nstr(char_poly[0], 10)} = 0")

        roots = mp.polyroots([char_poly[3], char_poly[2], char_poly[1], char_poly[0]])
        print(f"  Roots: {[mp.nstr(r, 10) for r in roots]}")
        break

print("\nDone.")
