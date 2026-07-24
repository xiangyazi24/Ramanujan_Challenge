#!/usr/bin/env sage
"""
P2.7: Test Hermite-Padé approach from Q4904.
Check if type-II HP approximants for (Li_2, Li_3) measures match P2.7 recurrence.
"""
from sage.all import *

R.<t> = PolynomialRing(QQ)

def Q_rs(r, s):
    """Type-II Hermite-Padé polynomial of multi-index (r,s)"""
    N = r + s
    if N == 0:
        return R(1)

    M = matrix(QQ, N, N + 1)
    row = 0
    for k in range(r):
        for m in range(N + 1):
            M[row, m] = QQ(1) / QQ((k + m + 1)^2)
        row += 1
    for k in range(s):
        for m in range(N + 1):
            M[row, m] = QQ(1) / QQ((k + m + 1)^3)
        row += 1

    ker = M.right_kernel().basis()
    if len(ker) != 1:
        print(f"  WARNING: kernel dim = {len(ker)} at (r,s)=({r},{s})")
        return None
    v = ker[0]
    v = v / v[N]  # monic normalization
    return sum(v[m] * t^m for m in range(N + 1))

def Q_step(N):
    """Step-line HP polynomial"""
    return Q_rs((N + 1)//2, N//2)

# Compute HP polynomials and their values at t=1
print("Computing step-line Hermite-Padé polynomials...")
q_values = []
for N in range(15):
    P_poly = Q_step(N)
    if P_poly is None:
        q_values.append(None)
        continue
    q_val = P_poly(1)
    q_values.append(q_val)
    print(f"  N={N}: Q_step = {P_poly}, q={q_val}")

# Check: do q_values satisfy a 4-term recurrence?
print("\nChecking for 4-term recurrence in q_n = Q_step(n)(1)...")
# A_n q_{n+3} + B_n q_{n+2} + C_n q_{n+1} + D_n q_n = 0

# For comparison: P2.7 recurrence coefficients
# A_n = 1024(2n+5)^4(2n+7)^3(2n+9)^3(946n²+6407n+10860)
# etc. - very large coefficients

# Let's check ratios instead
# The recurrence for P2.7 has Poincaré roots that are
# squared Cooper level-11 roots

# First check if q_n has a 3-term pattern (order 2) or 4-term (order 3)
print("\nRatios q_{n+1}/q_n:")
for n in range(min(12, len(q_values)-1)):
    if q_values[n] != 0 and q_values[n+1] is not None:
        print(f"  n={n}: {float(q_values[n+1]/q_values[n]):.10f}")

# Try to find a recurrence of order 3 (4-term)
print("\nSearching for order-3 recurrence with polynomial coefficients...")
# For each trial degree d, try to find c_3(n), c_2(n), c_1(n), c_0(n)
# of degrees d, d, d, d such that c_3(n)q[n+3] + ... = 0

for max_deg in range(1, 8):
    # total unknowns = 4*(max_deg+1) - 1 = 4*max_deg + 3
    n_unk = 4*(max_deg + 1)
    n_eqs = min(len(q_values) - 3, n_unk + 2)

    if n_eqs < n_unk:
        continue

    mat = matrix(QQ, n_eqs, n_unk)
    for n_idx in range(n_eqs):
        for j in range(max_deg + 1):
            if q_values[n_idx + 3] is not None:
                mat[n_idx, j] = QQ(n_idx)^j * q_values[n_idx + 3]
            if q_values[n_idx + 2] is not None:
                mat[n_idx, max_deg + 1 + j] = QQ(n_idx)^j * q_values[n_idx + 2]
            if q_values[n_idx + 1] is not None:
                mat[n_idx, 2*(max_deg + 1) + j] = QQ(n_idx)^j * q_values[n_idx + 1]
            if q_values[n_idx] is not None:
                mat[n_idx, 3*(max_deg + 1) + j] = QQ(n_idx)^j * q_values[n_idx]

    ker = mat.right_kernel()
    if ker.dimension() > 0:
        print(f"  deg={max_deg}: kernel dim = {ker.dimension()}")
        break
    else:
        print(f"  deg={max_deg}: no recurrence")

# Now test the P2.7 challenge recurrence directly with the q_values
print("\n=== Testing P2.7 recurrence directly ===")
# P2.7: u_{n+1} = (B_n/A_n) u_n - (C_{n-1}/A_{n-1}) u_{n-1} + (D_{n-2}/A_{n-2}) u_{n-2}
# With:
# A_n = 1024(2n+5)^4(2n+7)^3(2n+9)^3(946n^2+6407n+10860)
# B_n = 128(2n+7)^3(2n+9)^3(104060n^6+1745370n^5+12145238n^4+44886481n^3+92943995n^2+102256019n+46709052)
# C_n = 16(n+3)^4(2n+9)^3(3784n^5+57792n^4+351019n^3+1059230n^2+1587211n+944620)
# D_n = (n+3)^4(n+4)^6(946n^2+4515n+5399)

def A_p27(n):
    return QQ(1024)*(2*n+5)^4*(2*n+7)^3*(2*n+9)^3*(946*n^2+6407*n+10860)

def B_p27(n):
    return QQ(128)*(2*n+7)^3*(2*n+9)^3*(104060*n^6+1745370*n^5+12145238*n^4+44886481*n^3+92943995*n^2+102256019*n+46709052)

def C_p27(n):
    return QQ(16)*(n+3)^4*(2*n+9)^3*(3784*n^5+57792*n^4+351019*n^3+1059230*n^2+1587211*n+944620)

def D_p27(n):
    return (n+3)^4*(n+4)^6*(946*n^2+4515*n+5399)

# The recurrence is: A_n u_{n+1} - B_n u_n + C_{n-1} u_{n-1} - D_{n-2} u_{n-2} = 0
# Wait, let me re-read the challenge. It says:
# u_{n+1} = (B_n/A_n) u_n - (C_{n-1}/A_{n-1}) u_{n-1} + (D_{n-2}/A_{n-2}) u_{n-2}
# So: A_n u_{n+1} = B_n u_n - C_{n-1} u_{n-1} + D_{n-2} u_{n-2}

# Test with challenge initial values:
# q_0 = -215040420000
# q_1 = -167282265043404/905
# q_2 = -96418532765808/6071

q0_challenge = QQ(-215040420000)
q1_challenge = QQ(-167282265043404, 905)
q2_challenge = QQ(-964185327658080, 6071)

# Compute q3:
# A_2 q_3 = B_2 q_2 - C_1 q_1 + D_0 q_0
q3_challenge = (B_p27(2)*q2_challenge - C_p27(1)*q1_challenge + D_p27(0)*q0_challenge) / A_p27(2)
print(f"q₃ (challenge) = {q3_challenge}")

# Compare with HP values
print(f"\nHP q-values vs challenge (need gauge/normalization):")
print(f"  HP q_0 = {q_values[0]}")
print(f"  Challenge q_0 = {q0_challenge}")
if q_values[0] != 0:
    ratio = q0_challenge / q_values[0]
    print(f"  Ratio challenge/HP = {ratio}")

# The HP values start at 1 (monic normalization), challenge values are huge
# They won't match directly without a gauge transformation g(n)
# q_challenge(n) = g(n) * q_HP(n)

# Check if the Poincaré roots match
# P2.7 Poincaré roots come from the ratio of leading coefficients
lc_A = 1024  # leading term ~n^12
lc_D = 1     # leading term ~n^12
# Poincaré polynomial from recurrence A*u_{n+1} = B*u_n - C*u_{n-1} + D*u_{n-2}
# is determined by leading coefficients at each shift
# Actually the Poincaré polynomial for the 4-term recurrence uses the
# leading coefficients of A, B, C, D evaluated at n→∞

# Leading terms:
# A_n ~ 1024 * 2^4 * 2^3 * 2^3 * 946 * n^(2+4+3+3) = 1024*16*8*8*946 * n^12
# B_n ~ 128 * 2^3 * 2^3 * 104060 * n^(3+3+6) = 128*8*8*104060 * n^12
# C_n ~ 16 * 1 * 2^3 * 3784 * n^(4+3+5) = 16*8*3784 * n^12
# D_n ~ 1 * 1 * 946 * n^(4+6+2) = 946 * n^12

# Poincaré polynomial: A*ξ³ - B*ξ² + C*ξ - D = 0 (from the shifted form)
# Wait, need to be more careful with the Poincaré analysis

# At leading order, A(n)u_{n+1} = B(n)u_n - C(n-1)u_{n-1} + D(n-2)u_{n-2}
# The shift form: A(n)u_{n+1} - B(n)u_n + C(n-1)u_{n-1} - D(n-2)u_{n-2} = 0
# Poincaré poly: A_∞ ξ - B_∞ + C_∞/ξ - D_∞/ξ² = 0
# Multiply by ξ²: A_∞ ξ³ - B_∞ ξ² + C_∞ ξ - D_∞ = 0

lc_A_full = 1024 * 16 * 8 * 8 * 946
lc_B_full = 128 * 8 * 8 * 104060
lc_C_full = 16 * 8 * 3784
lc_D_full = 946

print(f"\nPoincaré coefficients (leading):")
print(f"  A_∞ = {lc_A_full}")
print(f"  B_∞ = {lc_B_full}")
print(f"  C_∞ = {lc_C_full}")
print(f"  D_∞ = {lc_D_full}")

# Poincaré polynomial: A ξ³ - B ξ² + C ξ - D = 0
x = var('x')
poincare = lc_A_full * x^3 - lc_B_full * x^2 + lc_C_full * x - lc_D_full
print(f"\nPoincaré poly: {poincare}")
roots = solve(poincare == 0, x)
print(f"Roots: {[r.rhs().n() for r in roots]}")

print("\nDone.")
