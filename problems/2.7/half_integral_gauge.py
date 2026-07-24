#!/usr/bin/env python3
"""
Verify the half-integral gauge factorization from Q5197.

Key claim: G(n) = R(n) · D_h(n) where:
  h_n = (4)_n^3 / [(5/2)_n (7/2)_n (9/2)_n]
  D_h(n) = diag(h_{n+2}, h_{n+1}, h_n)
  R(n) ∈ GL_3(Q(n))  — a RATIONAL matrix

If R(n) is rational, then the Zudilin error bounds transfer through
the gauge, proving c₀(e) = 0 unconditionally.

Method: compute R(n) = state_P(n) · D_h(n)^{-1} · state_Z(n)^{-1}
at many values of n and check if entries are rational functions.
"""
from fractions import Fraction as F
from mpmath import mp, mpf, matrix as mpmatrix, lu_solve, det as mpdet, fabs

mp.dps = 120

def pochhammer_exact(x, n):
    """(x)_n = x(x+1)...(x+n-1), exact rational."""
    result = F(1)
    for i in range(n):
        result *= (x + i)
    return result

def h_exact(n):
    """h_n = (4)_n^3 / [(5/2)_n (7/2)_n (9/2)_n]"""
    num = pochhammer_exact(F(4), n) ** 3
    den = (pochhammer_exact(F(5,2), n) *
           pochhammer_exact(F(7,2), n) *
           pochhammer_exact(F(9,2), n))
    return num / den

# === P2.7 coefficients ===
def A_c(n): return 1024*(2*n+5)**4*(2*n+7)**3*(2*n+9)**3*(946*n*n+6407*n+10860)
def B_c(n): return 128*(2*n+7)**3*(2*n+9)**3*(104060*n**6+1745370*n**5+12145238*n**4+44886481*n**3+92943995*n**2+102256019*n+46709052)
def C_c(n): return 16*(n+3)**4*(2*n+9)**3*(3784*n**5+57792*n**4+351019*n**3+1059230*n**2+1587211*n+944620)
def D_c(n): return (n+3)**4*(n+4)**6*(946*n*n+4515*n+5399)

# === Zudilin coefficients ===
def QZ(n): return 946*n**2 - 731*n + 153
def MZ(n): return 104060*n**6+127710*n**5+12788*n**4-34525*n**3-8482*n**2+3298*n+1071
def NZ(n): return 3784*n**5-1032*n**4-1925*n**3+853*n**2+328*n-184
def RZ(n): return 946*n**2+1161*n+368

N = 30

def zudilin_terms(init, N):
    u = list(init)
    for n in range(2, N):
        m = n
        d = 2 * QZ(m) * (2*m+1) * (m+1)**3
        nxt = F(2*MZ(m)) * u[n] + F(-2*m*NZ(m)) * u[n-1] + F(RZ(m)*m*(m-1)**3) * u[n-2]
        u.append(nxt / F(d))
    return u

def p27_terms(init, N):
    u = list(init)
    for n in range(2, N):
        nxt = F(B_c(n), A_c(n)) * u[n] + F(-C_c(n-1), A_c(n-1)) * u[n-1] + F(D_c(n-2), A_c(n-2)) * u[n-2]
        u.append(nxt)
    return u

# Zudilin solutions
b  = zudilin_terms([F(1), F(7), F(163)], N)
bt = zudilin_terms([F(0), F(23,2), F(2145,8)], N)
btt = zudilin_terms([F(0), F(17,2), F(3135,16)], N)

# P2.7 solutions (UNSCALED q_n, p_n, s_n)
q = p27_terms([F(-215040420000), F(-167282265043404, 905), F(-964185327658080, 6071)], N)
p = p27_terms([F(-612218384750), F(-9525021973931919, 18100), F(-29561828382772029, 65380)], N)

# Scaled P2.7 solutions: q̂_n = 64^n q_n, etc.
qhat = [F(64)**n * q[n] for n in range(N)]
phat = [F(64)**n * p[n] for n in range(N)]

# Third P2.7 solution with s_0=1, s_1=0, s_2=0
s = p27_terms([F(1), F(0), F(0)], N)
shat = [F(64)**n * s[n] for n in range(N)]

# Zudilin "mixed" solution m = b̃ + b̃̃
m_z = [bt[n] + btt[n] for n in range(N)]

# Precompute h values
h = [h_exact(n) for n in range(N)]

print("=== h_n check ===")
for n in range(6):
    hf = float(h[n])
    print(f"  h_{n} = {h[n]} ≈ {hf:.6f}")

print("\n=== Verify h_n growth ===")
for n in range(5, 15):
    ratio = float(h[n]) / float(h[n-1]) if h[n-1] != 0 else float('inf')
    expected = float(F(n+3)**3 / (F(n+F(3,2)) * F(n+F(5,2)) * F(n+F(7,2))))
    print(f"  h_{n}/h_{n-1} = {ratio:.10f}, expected r({n-1}) = {expected:.10f}")

def frac_to_mpf(x):
    if isinstance(x, F):
        return mpf(x.numerator) / mpf(x.denominator)
    return mpf(x)

# Compute R(n) = state_P(n) · D_h(n)^{-1} · state_Z(n)^{-1}
print("\n=== R(n) = state_P(n) · D_h(n)^{-1} · state_Z(n)^{-1} ===")

R_entries_exact = {(i,j): [] for i in range(3) for j in range(3)}
R_ns = []

for n in range(20):
    # Zudilin state: columns [b, bt+btt, third]
    # Use b, m_z = bt+btt, and bt (or btt) as third
    # Actually: we should use b, bt, btt as three independent solutions
    Z = mpmatrix(3, 3)
    for col, seq in enumerate([b, bt, btt]):
        Z[0, col] = frac_to_mpf(seq[n+2])
        Z[1, col] = frac_to_mpf(seq[n+1])
        Z[2, col] = frac_to_mpf(seq[n])

    # P2.7 state: columns [q̂, p̂, ŝ]
    P = mpmatrix(3, 3)
    for col, seq in enumerate([qhat, phat, shat]):
        P[0, col] = frac_to_mpf(seq[n+2])
        P[1, col] = frac_to_mpf(seq[n+1])
        P[2, col] = frac_to_mpf(seq[n])

    # D_h(n)^{-1} = diag(1/h_{n+2}, 1/h_{n+1}, 1/h_n)
    Dh_inv = mpmatrix(3, 3)
    Dh_inv[0, 0] = frac_to_mpf(F(1) / h[n+2])
    Dh_inv[1, 1] = frac_to_mpf(F(1) / h[n+1])
    Dh_inv[2, 2] = frac_to_mpf(F(1) / h[n])

    # R(n) = P · (D_h(n) · Z)^{-1} = P · Z^{-1} · D_h(n)^{-1}
    # Wait: G(n) = R(n) · D_h(n), so state_P(n) = G(n) · state_Z(n) = R(n) · D_h(n) · state_Z(n)
    # Therefore R(n) = state_P(n) · state_Z(n)^{-1} · D_h(n)^{-1}

    d = mpdet(Z)
    if fabs(d) < mpf(10)**(-40):
        print(f"  n={n}: Zudilin matrix singular")
        continue

    # Compute G(n) = state_P(n) · state_Z(n)^{-1}
    G = mpmatrix(3, 3)
    for row in range(3):
        rhs = mpmatrix(3, 1)
        for j in range(3):
            rhs[j, 0] = P[row, j]
        x = lu_solve(Z.T, rhs)
        for j in range(3):
            G[row, j] = x[j, 0]

    # R(n) = G(n) · D_h(n)^{-1}
    R = G * Dh_inv

    R_ns.append(n)
    for i in range(3):
        for j in range(3):
            R_entries_exact[(i,j)].append(float(R[i,j]))

    if n <= 5 or n >= 18:
        print(f"\nR({n}):")
        for i in range(3):
            vals = [f"{float(R[i,j]):20.6f}" for j in range(3)]
            print(f"  [{', '.join(vals)}]")

# Check if R(n) is approximately constant
print("\n\n=== Is R(n) approximately constant? ===")
for i in range(3):
    for j in range(3):
        vals = R_entries_exact[(i,j)]
        if len(vals) < 3:
            continue
        mn, mx = min(vals), max(vals)
        spread = mx - mn
        avg = sum(vals) / len(vals)
        rel = spread / (abs(avg) + 1e-300)
        print(f"  R[{i},{j}]: min={mn:.6e}, max={mx:.6e}, spread={spread:.2e}, rel_spread={rel:.2e}")

# Try rational function identification
print("\n\n=== Rational function identification for R(n) entries ===")
import numpy as np

for i in range(3):
    for j in range(3):
        vals = R_entries_exact[(i,j)]
        ns = R_ns
        if len(vals) < 6:
            continue

        for deg in range(8):
            ns_arr = np.array(ns, dtype=float)
            gs_arr = np.array(vals, dtype=float)
            if len(ns_arr) < deg + 2:
                break
            coeffs = np.polyfit(ns_arr, gs_arr, deg)
            pred_all = np.polyval(coeffs, ns_arr)
            resid = np.max(np.abs(pred_all - gs_arr))
            rel = resid / (np.max(np.abs(gs_arr)) + 1e-300)
            if rel < 1e-8:
                print(f"  R[{i},{j}]: POLYNOMIAL deg {deg}, rel resid = {rel:.2e}")
                if deg <= 3:
                    print(f"    coeffs = {coeffs}")
                break
        else:
            # Try Padé [p,q] with p+q <= 6
            best_rel = 1.0
            best_desc = ""
            from numpy.linalg import lstsq
            for p_deg in range(5):
                for q_deg in range(5):
                    if p_deg + q_deg > 6:
                        continue
                    n_unknowns = (p_deg + 1) + q_deg  # q_0 = 1
                    if len(vals) < n_unknowns + 2:
                        continue
                    # f(n)(1 + b1*n + ... + bq*n^q) = a0 + a1*n + ... + ap*n^p
                    A_mat = np.zeros((len(vals), n_unknowns))
                    b_vec = np.array(vals)
                    for idx, (nn, gv) in enumerate(zip(ns, vals)):
                        for d in range(p_deg + 1):
                            A_mat[idx, d] = nn**d
                        for d in range(1, q_deg + 1):
                            A_mat[idx, p_deg + 1 + d - 1] = -gv * nn**d
                    sol, _, _, _ = lstsq(A_mat, b_vec, rcond=None)
                    pred = []
                    for nn in ns:
                        denom = 1 + sum(sol[p_deg + 1 + d - 1] * nn**d for d in range(1, q_deg + 1))
                        numer = sum(sol[d] * nn**d for d in range(p_deg + 1))
                        pred.append(numer/denom if abs(denom) > 1e-20 else 1e30)
                    resid = max(abs(pp - v) for pp, v in zip(pred, vals))
                    rel = resid / (max(abs(v) for v in vals) + 1e-300)
                    if rel < best_rel:
                        best_rel = rel
                        best_desc = f"Padé [{p_deg},{q_deg}]"
            if best_rel < 1e-8:
                print(f"  R[{i},{j}]: {best_desc}, rel resid = {best_rel:.2e}")
            else:
                print(f"  R[{i},{j}]: no good fit, best rel = {best_rel:.2e}, vals[0:4] = {vals[:4]}")

print("\n=== det R(n) ===")
for idx, n in enumerate(R_ns[:20]):
    if idx >= len(R_entries_exact[(0,0)]):
        break
    R_mat = [[R_entries_exact[(i,j)][idx] for j in range(3)] for i in range(3)]
    d = (R_mat[0][0]*(R_mat[1][1]*R_mat[2][2]-R_mat[1][2]*R_mat[2][1])
        -R_mat[0][1]*(R_mat[1][0]*R_mat[2][2]-R_mat[1][2]*R_mat[2][0])
        +R_mat[0][2]*(R_mat[1][0]*R_mat[2][1]-R_mat[1][1]*R_mat[2][0]))
    print(f"  det R({n}) = {d:.10e}")
