#!/usr/bin/env python3
"""
Find the EXACT rational s values for the third P2.7 solution correspondence.
Uses exact Fraction arithmetic throughout — no floating point.
"""
from fractions import Fraction as F

# === Recurrence coefficients (exact) ===
def A_c(n): return 1024*(2*n+5)**4*(2*n+7)**3*(2*n+9)**3*(946*n*n+6407*n+10860)
def B_c(n): return 128*(2*n+7)**3*(2*n+9)**3*(104060*n**6+1745370*n**5+12145238*n**4+44886481*n**3+92943995*n**2+102256019*n+46709052)
def C_c(n): return 16*(n+3)**4*(2*n+9)**3*(3784*n**5+57792*n**4+351019*n**3+1059230*n**2+1587211*n+944620)
def D_c(n): return (n+3)**4*(n+4)**6*(946*n*n+4515*n+5399)

def QZ(n): return 946*n**2 - 731*n + 153
def MZ(n): return 104060*n**6+127710*n**5+12788*n**4-34525*n**3-8482*n**2+3298*n+1071
def NZ(n): return 3784*n**5-1032*n**4-1925*n**3+853*n**2+328*n-184
def RZ(n): return 946*n**2+1161*n+368

def r_val(n):
    return F(8*(n+4)**3, (2*n+5)*(2*n+7)*(2*n+9))

# 3x3 matrix operations over Fraction
def mat_mul(A, B):
    return [[sum(A[i][k]*B[k][j] for k in range(3)) for j in range(3)] for i in range(3)]

def mat_inv(M):
    a,b,c = M[0]; d,e,f = M[1]; g,h,i = M[2]
    det = a*(e*i-f*h) - b*(d*i-f*g) + c*(d*h-e*g)
    inv_det = F(1) / det
    return [
        [(e*i-f*h)*inv_det, (c*h-b*i)*inv_det, (b*f-c*e)*inv_det],
        [(f*g-d*i)*inv_det, (a*i-c*g)*inv_det, (c*d-a*f)*inv_det],
        [(d*h-e*g)*inv_det, (b*g-a*h)*inv_det, (a*e-b*d)*inv_det],
    ]

def mat_vec(M, v):
    return [sum(M[i][j]*v[j] for j in range(3)) for i in range(3)]

def mat_scale(s, M):
    return [[s*M[i][j] for j in range(3)] for i in range(3)]

def mat_add(A, B):
    return [[A[i][j]+B[i][j] for j in range(3)] for i in range(3)]

def mat_sub(A, B):
    return [[A[i][j]-B[i][j] for j in range(3)] for i in range(3)]

# Companion matrices
def CP_at(n):
    alpha = F(64) * F(B_c(n+2), A_c(n+2))
    beta = F(-64**2) * F(C_c(n+1), A_c(n+1))
    gamma = F(64**3) * F(D_c(n), A_c(n))
    return [[alpha, beta, gamma], [F(1), F(0), F(0)], [F(0), F(1), F(0)]]

def CZ_at(n):
    m = n + 2
    den = F(QZ(m) * (2*m+1) * (m+1)**3)
    alpha = F(MZ(m)) / den
    beta = F(-m * NZ(m)) / den
    gamma = F(RZ(m) * m * (m-1)**3, 2) / den
    return [[alpha, beta, gamma], [F(1), F(0), F(0)], [F(0), F(1), F(0)]]

# Source matrix Z0 = [z_b, z_2, z_3] (columns)
Z0 = [
    [F(163), F(2145,8), F(3135,16)],
    [F(7),   F(23,2),   F(17,2)],
    [F(1),   F(0),      F(0)],
]
Z0_inv = mat_inv(Z0)

# Target vectors
xq = [F(64)**2 * F(-964185327658080, 6071), F(64) * F(-167282265043404, 905), F(-215040420000)]
xp = [F(64)**2 * F(-29561828382772029, 65380), F(64) * F(-9525021973931919, 18100), F(-612218384750)]

print("xq exact:")
for k in range(3):
    print(f"  xq[{k}] = {xq[k]} = {float(xq[k]):.6e}")
print("xp exact:")
for k in range(3):
    print(f"  xp[{k}] = {xp[k]} = {float(xp[k]):.6e}")

# P2.7 forward iteration (exact)
def p27_fwd(init, N):
    u = list(init)
    for n in range(2, N):
        nxt = F(B_c(n), A_c(n)) * u[n] - F(C_c(n-1), A_c(n-1)) * u[n-1] + F(D_c(n-2), A_c(n-2)) * u[n-2]
        u.append(nxt)
    return u

# Zudilin forward iteration (exact)
def zud_fwd(init, N):
    u = list(init)
    for n in range(2, N):
        d = F(2 * QZ(n) * (2*n+1) * (n+1)**3)
        nxt = (F(2*MZ(n)) * u[n] + F(-2*n*NZ(n)) * u[n-1] + F(RZ(n)*n*(n-1)**3) * u[n-2]) / d
        u.append(nxt)
    return u

N = 15

print("\nComputing sequences (exact)...")
b  = zud_fwd([F(1), F(7), F(163)], N+5)
bt = zud_fwd([F(0), F(23,2), F(2145,8)], N+5)
btt= zud_fwd([F(0), F(17,2), F(3135,16)], N+5)

q_raw = p27_fwd([F(-215040420000), F(-167282265043404, 905), F(-964185327658080, 6071)], N+5)
p_raw = p27_fwd([F(-612218384750), F(-9525021973931919, 18100), F(-29561828382772029, 65380)], N+5)

qhat = [F(64)**n * q_raw[n] for n in range(N+5)]
phat = [F(64)**n * p_raw[n] for n in range(N+5)]

# Gauge propagation: R(n+1) = (1/r(n)) * C_P(n) * R(n) * C_Z(n)^{-1}
def propagate_R(R0, n_max):
    R_vals = [R0]
    for n in range(n_max):
        CZ_n = CZ_at(n)
        CP_n = CP_at(n)
        CZ_inv = mat_inv(CZ_n)
        R_CZinv = mat_mul(R_vals[-1], CZ_inv)
        R_next = mat_scale(F(1) / r_val(n), mat_mul(CP_n, R_CZinv))
        R_vals.append(R_next)
    return R_vals

# R0(s) = [xq, s, xp-s] · Z0^{-1}
# We parameterize s = (s2, s1, s0) and solve for the gauge equation to hold.
# The gauge equation R(n+1)·r(n)·C_Z(n) = C_P(n)·R(n) must hold as an identity.
# Since R(n) depends linearly on s, and the gauge equation is linear in R,
# the residual at each n is AFFINE in s.

# Compute R(0; s) = P0(s) · Z0^{-1}
# P0(s) = [xq, s, xp-s]
# P0(s)[i][0] = xq[i], P0(s)[i][1] = s[i], P0(s)[i][2] = xp[i] - s[i]

# R0(s) = P0(s) · Z0^{-1} is affine in (s0, s1, s2):
# R0 = R0_base + s0 * Delta0 + s1 * Delta1 + s2 * Delta2

# R0_base: s = (0, 0, 0)
P0_base = [[xq[i], F(0), xp[i]] for i in range(3)]
R0_base = mat_mul(P0_base, Z0_inv)

# Delta_k: derivative w.r.t. s_k
# d(P0)/d(s_k): column 1 gets e_k, column 2 gets -e_k
R0_deltas = []
for k in range(3):
    P0_d = [[F(0)]*3 for _ in range(3)]
    P0_d[k][1] = F(1)
    P0_d[k][2] = F(-1)
    R0_deltas.append(mat_mul(P0_d, Z0_inv))

# Propagate all components
print("Propagating gauge equations...")
R_base = propagate_R(R0_base, N)
R_deltas = [propagate_R(R0_deltas[k], N) for k in range(3)]

# The gauge equation should vanish: R(n+1)·r(n)·C_Z(n) - C_P(n)·R(n) = 0
# For the full R = R_base + Σ s_k R_delta_k, the residual at n is:
# res(n) = [R_base(n+1) + Σ s_k R_dk(n+1)] · r(n) · C_Z(n) - C_P(n) · [R_base(n) + Σ s_k R_dk(n)]
#
# Each R component ALREADY satisfies the propagation R(n+1) = (1/r)·C_P·R(n)·C_Z^{-1}
# by construction! So the gauge equation is automatically satisfied.
# The propagation IS the gauge equation.

# This means: for ANY s, the propagated R satisfies the gauge equation exactly.
# The question is whether R(n) is a RATIONAL FUNCTION of n (entries are p(n)/q(n) for polynomials p, q).

# Since R(n) satisfies the gauge equation for ALL n by construction,
# and the gauge equation has rational coefficients,
# if R(n) is rational for sufficiently many n values, it extends to a rational function.

# But every R with rational R(0) has rational entries at every integer n!
# (Because C_P(n) and C_Z(n) have rational entries for integer n.)

# The question is: among the 3-parameter family of R's, is there one where
# the entries are rational FUNCTIONS (polynomial ratios), not just rational numbers at each integer.

# For R to be a rational function of n, the entries R_ij(n) must satisfy:
# R_ij(n) = p_ij(n) / d(n) for some polynomials p_ij, d.
# This means R_ij(n) · d(n) is a polynomial in n for all i,j.

# Let's try: check if the entries of R(n; s=s_opt) look like rational functions
# by computing many values and trying to fit.

# First: use the exact propagation to compute R(n) values for the optimal s.
# The optimal s from the least-squares (exact Fraction result):

print("\nComputing optimal s via least squares...")

# Build normal equations
G = [[F(0)]*3 for _ in range(3)]
bvec = [F(0)]*3

for n in range(2, N):
    for k in range(3):
        for l in range(3):
            for i in range(3):
                for j in range(3):
                    G[k][l] += R_deltas[k][n][i][j] * R_deltas[l][n][i][j]
        for i in range(3):
            for j in range(3):
                bvec[k] -= R_base[n][i][j] * R_deltas[k][n][i][j]

G_inv = mat_inv(G)
c_opt = mat_vec(G_inv, bvec)

print(f"\nExact optimal s (Fraction):")
for k in range(3):
    print(f"  s[{k}] = {c_opt[k]}")
    print(f"        = {float(c_opt[k]):.15e}")
    print(f"        numerator digits: {len(str(abs(c_opt[k].numerator)))}")
    print(f"        denominator:      {c_opt[k].denominator}")

# Compute R(n) with optimal s
print("\nPropagating R with optimal s...")
R0_opt = [[R0_base[i][j] + sum(c_opt[k]*R0_deltas[k][i][j] for k in range(3)) for j in range(3)] for i in range(3)]
R_opt = propagate_R(R0_opt, N)

print("R(n) entries (exact rational):")
for n in [0, 1, 2, 3]:
    print(f"\n  R({n}):")
    for i in range(3):
        for j in range(3):
            val = R_opt[n][i][j]
            print(f"    R[{i},{j}] = {float(val):.6e}  (denom: {val.denominator})")

# Check: do ALL R(n) entries share a common denominator pattern?
print("\n\nDenominator analysis:")
for n in range(N):
    denoms = set()
    for i in range(3):
        for j in range(3):
            denoms.add(R_opt[n][i][j].denominator)
    lcm_d = 1
    for d in denoms:
        from math import gcd
        lcm_d = lcm_d * d // gcd(lcm_d, d)
    print(f"  n={n}: lcm of denominators = {lcm_d}, #digits = {len(str(lcm_d))}")

# The s values from least squares are NOT the true minimum for a polynomial gauge.
# The true s should make R(n) exactly a rational function p(n)/q(n).
# This means: for the correct s, the NUMERATORS D(n)*R_ij(n) should be polynomial in n.

# Let's try a different approach: direct fitting.
# Assume R_ij(n) = p_ij(n) / D(n) where D(n) is a product of known factors.
# The Sage code from Q5202 tries D = Delta denominator, etc.

# But first, let's check: what do the R(n) denominators look like?
print("\n\nFactorization of R(0) denominators:")
R0 = R_opt[0]
for i in range(3):
    for j in range(3):
        d = R0[i][j].denominator
        n = R0[i][j].numerator
        print(f"  R0[{i},{j}]: num={n}, den={d}")
