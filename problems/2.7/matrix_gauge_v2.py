#!/usr/bin/env python3
"""
Matrix gauge G(n) between P2.7 and Zudilin companion systems.
Using UNSCALED sequences (same Poincaré polynomial).

q_n (P2.7) and b_n (Zudilin) have the SAME Poincaré polynomial:
4ν³ - 220ν² + 8ν - 1 = 0

So a rational matrix gauge G(n) should exist if the modules are isomorphic.
G(n) = state_P(n) · state_Z(n)^{-1}
"""
from fractions import Fraction as F
from mpmath import mp, mpf, matrix as mpmatrix, lu_solve, det as mpdet, log10, fabs

mp.dps = 80

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

# Compute sequences with exact rational arithmetic
N = 25

def zudilin_terms(init, N):
    u = list(init)
    for n in range(2, N):
        m = n
        d = 2 * QZ(m) * (2*m+1) * (m+1)**3
        nxt = F(2*MZ(m)) * u[n] + F(-2*m*NZ(m)) * u[n-1] + F(RZ(m)*m*(m-1)**3) * u[n-2]
        u.append(nxt / d)
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

# P2.7 solutions (UNSCALED)
q = p27_terms([F(-215040420000), F(-167282265043404, 905), F(-964185327658080, 6071)], N)
p = p27_terms([F(-612218384750), F(-9525021973931919, 18100), F(-29561828382772029, 65380)], N)
s = p27_terms([F(1), F(0), F(0)], N)

print("=== Verify: dominant multiplier comparison ===")
for n in range(5, 12):
    rho_q = float(q[n+1]) / float(q[n])
    rho_b = float(b[n+1]) / float(b[n])
    print(f"  n={n}: q_{n+1}/q_n = {rho_q:.10f}, b_{n+1}/b_n = {rho_b:.10f}")

def frac_to_mpf(x):
    if isinstance(x, F):
        return mpf(x.numerator) / mpf(x.denominator)
    return mpf(x)

# Compute G(n) = state_P(n) · state_Z(n)^{-1}
print("\n=== G(n) at n=0,...,15 ===")
G_entries = {(i,j): [] for i in range(3) for j in range(3)}
G_ns = []

for n in range(16):
    # Zudilin state: columns are [b, bt, btt] at indices n+2, n+1, n
    Z = mpmatrix(3, 3)
    for col, seq in enumerate([b, bt, btt]):
        Z[0, col] = frac_to_mpf(seq[n+2])
        Z[1, col] = frac_to_mpf(seq[n+1])
        Z[2, col] = frac_to_mpf(seq[n])

    # P2.7 state: columns are [q, p, s]
    P = mpmatrix(3, 3)
    for col, seq in enumerate([q, p, s]):
        P[0, col] = frac_to_mpf(seq[n+2])
        P[1, col] = frac_to_mpf(seq[n+1])
        P[2, col] = frac_to_mpf(seq[n])

    # G(n) = P · Z^{-1}
    d = mpdet(Z)
    if fabs(d) < mpf(10)**(-40):
        print(f"  n={n}: Zudilin matrix nearly singular, det = {d}")
        continue

    # Solve Z^T · G^T = P^T for each row of G
    G = mpmatrix(3, 3)
    for row in range(3):
        rhs = mpmatrix(3, 1)
        for j in range(3):
            rhs[j, 0] = P[row, j]
        # Solve Z^T x = rhs
        ZT = Z.T
        x = lu_solve(ZT, rhs)
        for j in range(3):
            G[row, j] = x[j, 0]

    G_ns.append(n)
    for i in range(3):
        for j in range(3):
            G_entries[(i,j)].append(float(G[i,j]))

    # Print G(n) compactly
    if n <= 5 or n >= 14:
        print(f"\nG({n}):")
        for i in range(3):
            vals = [f"{float(G[i,j]):15.6f}" for j in range(3)]
            print(f"  [{', '.join(vals)}]")

# Check if G(n) is approximately constant
print("\n\n=== Is G(n) approximately constant? ===")
for i in range(3):
    for j in range(3):
        vals = G_entries[(i,j)]
        if len(vals) < 3:
            continue
        mn, mx = min(vals), max(vals)
        spread = mx - mn
        avg = sum(vals) / len(vals)
        rel = spread / (abs(avg) + 1e-300)
        print(f"  G[{i},{j}]: min={mn:.6e}, max={mx:.6e}, spread={spread:.2e}, rel_spread={rel:.2e}")

# Try to identify rational functions: G[i,j](n) = P(n)/Q(n)
print("\n\n=== Rational function identification ===")
for i in range(3):
    for j in range(3):
        vals = G_entries[(i,j)]
        ns = G_ns
        if len(vals) < 6:
            continue

        # Try polynomial fit deg 0 (constant)
        for deg in range(6):
            import numpy as np
            ns_arr = np.array(ns[:deg+2], dtype=float)
            gs_arr = np.array(vals[:deg+2], dtype=float)
            if len(ns_arr) < deg + 1:
                break
            coeffs = np.polyfit(ns_arr, gs_arr, deg)
            pred_all = np.polyval(coeffs, np.array(ns, dtype=float))
            resid = np.max(np.abs(pred_all - np.array(vals, dtype=float)))
            rel = resid / (np.max(np.abs(np.array(vals, dtype=float))) + 1e-300)
            if rel < 1e-8:
                print(f"  G[{i},{j}]: POLYNOMIAL deg {deg}, rel resid = {rel:.2e}")
                print(f"    coeffs = {coeffs}")
                break
        else:
            # Try Padé approximant: f(n) = (a0 + a1*n) / (1 + b1*n)
            # Cross-multiply: f(n)(1 + b1*n) = a0 + a1*n
            # f(n) + f(n)*b1*n = a0 + a1*n
            # Three unknowns, use three points
            if len(vals) >= 5:
                from numpy.linalg import lstsq
                # f(n) = (a0 + a1*n + a2*n^2) / (1 + b1*n + b2*n^2)
                # f(n) + f(n)*b1*n + f(n)*b2*n^2 = a0 + a1*n + a2*n^2
                A_mat = np.zeros((len(vals), 5))
                b_vec = np.array(vals)
                for idx, (nn, gv) in enumerate(zip(ns, vals)):
                    A_mat[idx, 0] = 1
                    A_mat[idx, 1] = nn
                    A_mat[idx, 2] = nn**2
                    A_mat[idx, 3] = -gv * nn
                    A_mat[idx, 4] = -gv * nn**2
                sol, _, _, _ = lstsq(A_mat, b_vec, rcond=None)
                pred = []
                for nn in ns:
                    denom = 1 + sol[3]*nn + sol[4]*nn**2
                    numer = sol[0] + sol[1]*nn + sol[2]*nn**2
                    pred.append(numer/denom if abs(denom) > 1e-20 else 1e30)
                resid = max(abs(p - v) for p, v in zip(pred, vals))
                rel = resid / (max(abs(v) for v in vals) + 1e-300)
                if rel < 1e-8:
                    print(f"  G[{i},{j}]: RATIONAL [2,2], rel resid = {rel:.2e}")
                    print(f"    num = {sol[0]:.6e} + {sol[1]:.6e}*n + {sol[2]:.6e}*n^2")
                    print(f"    den = 1 + {sol[3]:.6e}*n + {sol[4]:.6e}*n^2")
                else:
                    print(f"  G[{i},{j}]: no good fit, rel resid = {rel:.2e}, first values: {vals[:4]}")

print("\n=== Check: det G(n) ===")
for n in range(min(16, len(G_entries[(0,0)]))):
    G_mat = [[G_entries[(i,j)][n] for j in range(3)] for i in range(3)]
    d = (G_mat[0][0]*(G_mat[1][1]*G_mat[2][2]-G_mat[1][2]*G_mat[2][1])
        -G_mat[0][1]*(G_mat[1][0]*G_mat[2][2]-G_mat[1][2]*G_mat[2][0])
        +G_mat[0][2]*(G_mat[1][0]*G_mat[2][1]-G_mat[1][1]*G_mat[2][0]))
    print(f"  det G({n}) = {d:.10e}")
