#!/usr/bin/env python3
"""P2.5: Mixed-period ansatz — find polynomial r_{ij}(n) such that
J_n = Σ_{0≤i,j≤2} r_{ij}(n) P_{n+i}(3) Q_{n+j}(3)
satisfies the degree-13 scalar recurrence for Q̂_n.

Uses modular arithmetic for fast rank computation, then exact for extraction.
"""
from mpmath import mp, mpf, log
mp.dps = 120

c_coeffs = [
    [-42743162700, -206623731375, -448112471583, -579493151986, -500074412234,
     -304838513875, -135313819947, -44354837964, -10750665744, -1905022784,
     -240100240, -20397440, -1047552, -24576],
    [8781630505200, 38850314624124, 78557994908508, 96136040496551, 79442239242197,
     46814452218572, 20241514501104, 6502490145168, 1552168938336, 271943188864,
     33995217088, 2871763456, 146952192, 3440640],
    [-10566229124340, -43764612822972, -82725628159809, -94536939564882, -72904809920709,
     -40082159230086, -16169158004002, -4847446446296, -1080358338832, -176841798272,
     -20670362464, -1634185472, -78342144, -1720320],
    [146862156672, 610678861056, 1158857071416, 1329423744980, 1029037642166,
     567735994679, 229759169143, 69074560420, 15430450432, 2530117664,
     296032016, 23408000, 1121280, 24576],
]

def eval_c(i, n):
    val = mpf(0)
    nk = mpf(1)
    for coeff in c_coeffs[i]:
        val += coeff * nk
        nk *= n
    return val

NMAX = 120
print(f"Computing P_n(3) and Q_n(3) for n = 0..{NMAX}...")

P = [mpf(0)] * (NMAX + 1)
Q = [mpf(0)] * (NMAX + 1)
P[0] = mpf(1)
P[1] = mpf(3)
Q[0] = log(2) / 2
Q[1] = 3 * Q[0] - 1

for n in range(1, NMAX):
    P[n + 1] = ((6 * n + 3) * P[n] - n * P[n - 1]) / (n + 1)
    Q[n + 1] = ((6 * n + 3) * Q[n] - n * Q[n - 1]) / (n + 1)

print(f"  P_10 = {P[10]}, Q_10 = {Q[10]}")
from mpmath import nstr
print(f"  P_20*Q_20 = {nstr(P[20]*Q[20], 15)}")
print(f"  1/(4√2·20) = {nstr(1/(4*mpf(2)**mpf('0.5')*20), 15)}")

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

def check_rank_numerical(num_unknowns, max_deg, num_eqs):
    """Check if polynomial r_{ij}(n) of degree max_deg gives null space."""
    total_unknowns = 9 * (max_deg + 1)
    assert total_unknowns == num_unknowns

    # Build matrix rows: for each n, equation is
    # Σ_{deg=0}^{max_deg} Σ_{i,j=0}^{2} r_{ij,deg} n^deg L[B_{ij}](n) = 0
    rows = []
    for eq_n in range(num_eqs):
        n = eq_n
        row = [mpf(0)] * total_unknowns
        ok = True
        for i in range(3):
            for j in range(3):
                val = L_basis(i, j, n)
                if val is None:
                    ok = False
                    break
                for d in range(max_deg + 1):
                    col = d * 9 + i * 3 + j
                    row[col] = (mpf(n) ** d) * val
            if not ok:
                break
        if ok:
            rows.append(row)

    actual_eqs = len(rows)
    print(f"  {total_unknowns} unknowns, {actual_eqs} equations (needed ≥ {total_unknowns})")
    if actual_eqs < total_unknowns:
        print("  Insufficient equations!")
        return -1

    # Gaussian elimination with pivoting (numerical)
    mat = [row[:] for row in rows[:total_unknowns + 5]]
    nr = len(mat)
    nc = total_unknowns
    rank = 0
    for c in range(nc):
        # Find pivot
        best = -1
        best_val = mpf(0)
        for i in range(rank, nr):
            if abs(mat[i][c]) > best_val:
                best_val = abs(mat[i][c])
                best = i
        if best_val < mpf('1e-50'):
            continue
        mat[rank], mat[best] = mat[best], mat[rank]
        pv = mat[rank][c]
        for j in range(nc):
            mat[rank][j] /= pv
        for i in range(nr):
            if i != rank and abs(mat[i][c]) > mpf('1e-100'):
                f = mat[i][c]
                for j in range(nc):
                    mat[i][j] -= f * mat[rank][j]
        rank += 1

    null_dim = total_unknowns - rank
    print(f"  rank = {rank}, null_dim = {null_dim}")

    if null_dim > 0:
        # Extract null vector
        pivot_cols = []
        r2 = 0
        mat2 = [row[:] for row in rows[:total_unknowns + 5]]
        for c in range(nc):
            best = -1
            best_val = mpf(0)
            for i in range(r2, len(mat2)):
                if abs(mat2[i][c]) > best_val:
                    best_val = abs(mat2[i][c])
                    best = i
            if best_val < mpf('1e-50'):
                continue
            mat2[r2], mat2[best] = mat2[best], mat2[r2]
            pivot_cols.append(c)
            pv = mat2[r2][c]
            for j in range(nc):
                mat2[r2][j] /= pv
            for i in range(len(mat2)):
                if i != r2 and abs(mat2[i][c]) > mpf('1e-100'):
                    f = mat2[i][c]
                    for j in range(nc):
                        mat2[i][j] -= f * mat2[r2][j]
            r2 += 1

        free_cols = [c for c in range(nc) if c not in pivot_cols]
        print(f"  Free columns: {free_cols}")

        # Extract one null vector
        x = [mpf(0)] * nc
        fc = free_cols[0]
        x[fc] = mpf(1)
        for pr_idx in range(len(pivot_cols) - 1, -1, -1):
            pc = pivot_cols[pr_idx]
            s = sum(mat2[pr_idx][j] * x[j] for j in range(pc + 1, nc))
            x[pc] = -s

        # Verify on ALL equations
        max_res = mpf(0)
        for eq_n in range(min(actual_eqs, 80)):
            n = eq_n
            val = mpf(0)
            for d in range(max_deg + 1):
                for i in range(3):
                    for j in range(3):
                        col = d * 9 + i * 3 + j
                        lb = L_basis(i, j, n)
                        if lb is not None:
                            val += x[col] * (mpf(n) ** d) * lb
            max_res = max(max_res, abs(val))

        print(f"  Verification max residual: {nstr(max_res, 6)}")

        # Display the r_{ij} polynomials
        print("\n  r_{ij}(n) coefficients:")
        for i in range(3):
            for j in range(3):
                poly_coeffs = []
                for d in range(max_deg + 1):
                    col = d * 9 + i * 3 + j
                    poly_coeffs.append(x[col])
                nonzero = any(abs(c) > mpf('1e-50') for c in poly_coeffs)
                if nonzero:
                    terms = []
                    for d, c in enumerate(poly_coeffs):
                        if abs(c) > mpf('1e-50'):
                            if d == 0:
                                terms.append(f"{float(c):.8g}")
                            elif d == 1:
                                terms.append(f"{float(c):.8g}·n")
                            else:
                                terms.append(f"{float(c):.8g}·n^{d}")
                    print(f"    r_{{{i},{j}}}(n) = {' + '.join(terms)}")

    return null_dim

# Search with increasing polynomial degree
for deg in range(8):
    num_unk = 9 * (deg + 1)
    num_eq = num_unk + 15
    print(f"\n{'='*60}")
    print(f"Degree {deg}: r_{{ij}}(n) polynomial of degree ≤ {deg} ({num_unk} unknowns)")
    print(f"{'='*60}")
    nd = check_rank_numerical(num_unk, deg, num_eq)
    trivial_dim = deg + 1  # Wronskian identity gives d+1 trivial null vectors
    effective = nd - trivial_dim
    print(f"  Trivial dim (Wronskian): {trivial_dim}, effective null dim: {effective}")
    if effective > 0:
        print(f"\n*** NONTRIVIAL solution at degree {deg}! ***")
        break
    if nd < 0:
        print("  Stopping: insufficient data")
        break
