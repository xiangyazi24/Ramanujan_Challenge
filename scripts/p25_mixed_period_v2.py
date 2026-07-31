#!/usr/bin/env python3
"""P2.5: Mixed-period ansatz with CORRECT globally-normalized recurrence coefficients.
Poincaré polynomial: ξ³-35ξ²+35ξ-1 = (ξ-1)(ξ²-34ξ+1).
"""
from mpmath import mp, mpf, log, nstr
mp.dps = 200

# CORRECT globally-normalized coefficients
c_coeffs = [
    [-170972650800, -826494925500, -1792449886332, -2317972607944, -2000297648936,
     -1219354055500, -541255279788, -177419351856, -43002662976, -7620091136,
     -960400960, -81589760, -4190208, -98304],
    [8781630505200, 38850314624124, 78557994908508, 96136040496551, 79442239242197,
     46814452218572, 20241514501104, 6502490145168, 1552168938336, 271943188864,
     33995217088, 2871763456, 146952192, 3440640],
    [-21132458248680, -87529225645944, -165451256319618, -189073879129764, -145809619841418,
     -80164318460172, -32338316008004, -9694892892592, -2160716677664, -353683596544,
     -41340724928, -3268370944, -156684288, -3440640],
    [587448626688, 2442715444224, 4635428285664, 5317694979920, 4116150568664,
     2270943978716, 919036676572, 276298241680, 61721801728, 10120470656,
     1184128064, 93632000, 4485120, 98304],
]

def eval_c(i, n):
    val = mpf(0)
    nk = mpf(1)
    for coeff in c_coeffs[i]:
        val += coeff * nk
        nk *= n
    return val

NMAX = 200
print(f"Computing P_n(3) and Q_n(3) for n = 0..{NMAX}...")

P = [mpf(0)] * (NMAX + 1)
Q = [mpf(0)] * (NMAX + 1)
P[0] = mpf(1); P[1] = mpf(3)
Q[0] = log(2) / 2; Q[1] = 3 * Q[0] - 1

for n in range(1, NMAX):
    P[n+1] = ((6*n+3)*P[n] - n*P[n-1]) / (n+1)
    Q[n+1] = ((6*n+3)*Q[n] - n*Q[n-1]) / (n+1)

def L_basis(i, j, n):
    """Apply recurrence to P_{n+i}(3)·Q_{n+j}(3)."""
    val = mpf(0)
    for s in range(4):
        ni = n + s + i
        nj = n + s + j
        if 0 <= ni <= NMAX and 0 <= nj <= NMAX:
            val += eval_c(s, n) * P[ni] * Q[nj]
        else:
            return None
    return val

# Also check: does the correct recurrence annihilate D_n^2, PQ, Q^2?
print("\nChecking if correct recurrence annihilates standard symmetric-square basis:")
for name, fi, fj in [("D_n^2", lambda n: P[n], lambda n: P[n]),
                       ("P_nQ_n", lambda n: P[n], lambda n: Q[n]),
                       ("Q_n^2", lambda n: Q[n], lambda n: Q[n])]:
    max_res = mpf(0)
    for n in range(50):
        val = mpf(0)
        for s in range(4):
            val += eval_c(s, n) * fi(n+s) * fj(n+s)
        max_res = max(max_res, abs(val))
    print(f"  L[{name}]: max residual = {nstr(max_res, 6)}")

# Search for polynomial r_{ij}(n) such that J_n = Σ r_{ij}(n) P_{n+i}Q_{n+j} satisfies L[J_n]=0
print("\n" + "="*60)
print("Mixed-period ansatz search with CORRECT recurrence")
print("="*60)

for deg in range(10):
    total_unknowns = 9 * (deg + 1)
    num_eq = total_unknowns + 20
    
    rows = []
    for eq_n in range(num_eq):
        n = eq_n
        row = [mpf(0)] * total_unknowns
        ok = True
        for i in range(3):
            for j in range(3):
                val = L_basis(i, j, n)
                if val is None:
                    ok = False
                    break
                for d in range(deg + 1):
                    col = d * 9 + i * 3 + j
                    row[col] = (mpf(n) ** d) * val
            if not ok:
                break
        if ok:
            rows.append(row)
    
    actual_eqs = len(rows)
    nc = total_unknowns
    
    # Gaussian elimination
    mat = [row[:] for row in rows[:nc + 10]]
    nr = len(mat)
    rank = 0
    for c in range(nc):
        best = -1
        best_val = mpf(0)
        for i in range(rank, nr):
            if abs(mat[i][c]) > best_val:
                best_val = abs(mat[i][c])
                best = i
        if best_val < mpf('1e-100'):
            continue
        mat[rank], mat[best] = mat[best], mat[rank]
        pv = mat[rank][c]
        for j2 in range(nc):
            mat[rank][j2] /= pv
        for i in range(nr):
            if i != rank and abs(mat[i][c]) > mpf('1e-150'):
                f = mat[i][c]
                for j2 in range(nc):
                    mat[i][j2] -= f * mat[rank][j2]
        rank += 1
    
    null_dim = nc - rank
    trivial_dim = deg + 1  # Wronskian
    effective = null_dim - trivial_dim
    
    print(f"\nDegree {deg}: {nc} unknowns, {actual_eqs} eqs, rank={rank}, null_dim={null_dim}, "
          f"trivial={trivial_dim}, effective={effective}")
    
    if effective > 0:
        print(f"\n*** NONTRIVIAL SOLUTION AT DEGREE {deg}! ***")
        # Extract null vectors
        mat2 = [row[:] for row in rows[:nc + 10]]
        r2 = 0
        pivot_cols = []
        for c in range(nc):
            best = -1
            best_val = mpf(0)
            for i in range(r2, len(mat2)):
                if abs(mat2[i][c]) > best_val:
                    best_val = abs(mat2[i][c])
                    best = i
            if best_val < mpf('1e-100'):
                continue
            mat2[r2], mat2[best] = mat2[best], mat2[r2]
            pivot_cols.append(c)
            pv = mat2[r2][c]
            for j2 in range(nc):
                mat2[r2][j2] /= pv
            for i in range(len(mat2)):
                if i != r2 and abs(mat2[i][c]) > mpf('1e-150'):
                    f = mat2[i][c]
                    for j2 in range(nc):
                        mat2[i][j2] -= f * mat2[r2][j2]
            r2 += 1
        
        free_cols = [c for c in range(nc) if c not in pivot_cols]
        print(f"  Free columns: {free_cols}")
        
        # Extract ALL null vectors
        for fidx, fc in enumerate(free_cols):
            x = [mpf(0)] * nc
            x[fc] = mpf(1)
            for pr_idx in range(len(pivot_cols)-1, -1, -1):
                pc = pivot_cols[pr_idx]
                s = sum(mat2[pr_idx][j2] * x[j2] for j2 in range(pc+1, nc))
                x[pc] = -s
            
            # Check if this is a Wronskian-trivial solution
            # Wronskian gives r_{ij} = f(n)(δ_{i<j} - δ_{i>j})/(n+j-n-i) type
            is_trivial = True
            for i in range(3):
                if abs(x[0*9 + i*3 + i]) > mpf('1e-100'):  # diagonal must be 0
                    is_trivial = False
            
            # Verify residual
            max_res = mpf(0)
            for eq_n in range(min(actual_eqs, 100)):
                n = eq_n
                val = mpf(0)
                for d in range(deg + 1):
                    for i in range(3):
                        for j in range(3):
                            col = d * 9 + i * 3 + j
                            lb = L_basis(i, j, n)
                            if lb is not None:
                                val += x[col] * (mpf(n) ** d) * lb
                max_res = max(max_res, abs(val))
            
            label = "TRIVIAL" if is_trivial else "NONTRIVIAL"
            print(f"\n  Null vector {fidx} ({label}), residual: {nstr(max_res, 6)}")
            
            if not is_trivial:
                print("  r_{ij} coefficients:")
                for i in range(3):
                    for j in range(3):
                        poly = [x[d*9 + i*3 + j] for d in range(deg+1)]
                        if any(abs(c) > mpf('1e-50') for c in poly):
                            terms = []
                            for d, c in enumerate(poly):
                                if abs(c) > mpf('1e-50'):
                                    cf = float(c)
                                    if d == 0:
                                        terms.append(f"{cf:.10g}")
                                    elif d == 1:
                                        terms.append(f"{cf:.10g}·n")
                                    else:
                                        terms.append(f"{cf:.10g}·n^{d}")
                            print(f"    r_{{{i},{j}}}(n) = {' + '.join(terms)}")
        break
    
    if null_dim == 0 and nc > 100:
        print("  Too many unknowns, stopping")
        break

