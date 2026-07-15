#!/usr/bin/env python3
"""
Rational gauge search using the CORRECT rank-one h-twist.

Key correction from Q5202: the h-twist is rank-one (h_n · I_3), giving
  C_Z^{(h)}(n) = r(n) · C_Z(n)
where r(n) = h_{n+1}/h_n = (n+4)^3 / [(n+5/2)(n+7/2)(n+9/2)]

NOT the diagonal state-vector twist diag(h_{n+2}, h_{n+1}, h_n).

The gauge equation is: R(n+1) · r(n) · C_Z(n) = C_P(n) · R(n)
"""
from fractions import Fraction as F
from mpmath import mp, mpf, nstr, log10, fabs

mp.dps = 80

# === P2.7 recurrence coefficients ===
def A_c(n): return 1024*(2*n+5)**4*(2*n+7)**3*(2*n+9)**3*(946*n*n+6407*n+10860)
def B_c(n): return 128*(2*n+7)**3*(2*n+9)**3*(104060*n**6+1745370*n**5+12145238*n**4+44886481*n**3+92943995*n**2+102256019*n+46709052)
def C_c(n): return 16*(n+3)**4*(2*n+9)**3*(3784*n**5+57792*n**4+351019*n**3+1059230*n**2+1587211*n+944620)
def D_c(n): return (n+3)**4*(n+4)**6*(946*n*n+4515*n+5399)

# Zudilin recurrence
def QZ(n): return 946*n**2 - 731*n + 153
def MZ(n): return 104060*n**6+127710*n**5+12788*n**4-34525*n**3-8482*n**2+3298*n+1071
def NZ(n): return 3784*n**5-1032*n**4-1925*n**3+853*n**2+328*n-184
def RZ(n): return 946*n**2+1161*n+368

# r(n) = (n+4)^3 / [(n+5/2)(n+7/2)(n+9/2)]
#       = (n+4)^3 / [(2n+5)(2n+7)(2n+9)/8]
#       = 8(n+4)^3 / [(2n+5)(2n+7)(2n+9)]
def r_val(n):
    return F(8*(n+4)**3, (2*n+5)*(2*n+7)*(2*n+9))

# Companion matrices at integer n, using exact Fractions
def CP_at(n):
    """P2.7 scaled companion matrix at step n."""
    alpha = F(64) * F(B_c(n+2), A_c(n+2))
    beta = F(-64**2) * F(C_c(n+1), A_c(n+1))
    gamma = F(64**3) * F(D_c(n), A_c(n))
    return [[alpha, beta, gamma], [F(1), F(0), F(0)], [F(0), F(1), F(0)]]

def CZ_at(n):
    """Zudilin companion matrix at step n (m = n+2)."""
    m = n + 2
    den = F(QZ(m) * (2*m+1) * (m+1)**3)
    alpha = F(MZ(m)) / den
    beta = F(-m * NZ(m)) / den
    gamma = F(RZ(m) * m * (m-1)**3, 2) / den
    return [[alpha, beta, gamma], [F(1), F(0), F(0)], [F(0), F(1), F(0)]]

def mat_mul(A, B):
    """3x3 matrix multiply."""
    return [[sum(A[i][k]*B[k][j] for k in range(3)) for j in range(3)] for i in range(3)]

def mat_inv(M):
    """3x3 matrix inverse using exact Fractions."""
    a,b,c = M[0]
    d,e,f = M[1]
    g,h,i = M[2]
    det = a*(e*i-f*h) - b*(d*i-f*g) + c*(d*h-e*g)
    if det == 0:
        raise ValueError("Singular matrix")
    inv_det = F(1) / det
    return [
        [(e*i-f*h)*inv_det, (c*h-b*i)*inv_det, (b*f-c*e)*inv_det],
        [(f*g-d*i)*inv_det, (a*i-c*g)*inv_det, (c*d-a*f)*inv_det],
        [(d*h-e*g)*inv_det, (b*g-a*h)*inv_det, (a*e-b*d)*inv_det],
    ]

def mat_scale(s, M):
    """Scalar * matrix."""
    return [[s*M[i][j] for j in range(3)] for i in range(3)]

def mat_vec(M, v):
    """Matrix * vector."""
    return [sum(M[i][j]*v[j] for j in range(3)) for i in range(3)]

# Forward iteration (exact Fractions)
def iterate_fwd(init, N, recurrence='P27'):
    u = list(init)
    for n_step in range(2, N):
        if recurrence == 'P27':
            n = n_step
            nxt = F(B_c(n), A_c(n)) * u[n] - F(C_c(n-1), A_c(n-1)) * u[n-1] + F(D_c(n-2), A_c(n-2)) * u[n-2]
        else:  # Zudilin
            n = n_step
            d = F(2 * QZ(n) * (2*n+1) * (n+1)**3)
            nxt = F(2*MZ(n)) * u[n] / d + F(-2*n*NZ(n)) * u[n-1] / d + F(RZ(n)*n*(n-1)**3) * u[n-2] / d
        u.append(nxt)
    return u

N_MAX = 25

# Zudilin solutions
print("Computing Zudilin solutions...")
b  = iterate_fwd([F(1), F(7), F(163)], N_MAX+5, 'Zudilin')
bt = iterate_fwd([F(0), F(23,2), F(2145,8)], N_MAX+5, 'Zudilin')
btt= iterate_fwd([F(0), F(17,2), F(3135,16)], N_MAX+5, 'Zudilin')

# P2.7 solutions
print("Computing P2.7 solutions...")
q_raw = iterate_fwd([F(-215040420000), F(-167282265043404, 905), F(-964185327658080, 6071)], N_MAX+5, 'P27')
p_raw = iterate_fwd([F(-612218384750), F(-9525021973931919, 18100), F(-29561828382772029, 65380)], N_MAX+5, 'P27')

# Scaled: q̂_n = 64^n q_n
qhat = [F(64)**n * q_raw[n] for n in range(N_MAX+5)]
phat = [F(64)**n * p_raw[n] for n in range(N_MAX+5)]

# Three basis third solutions for P2.7
print("Computing basis third solutions...")
basis_inits = [[F(1), F(0), F(0)], [F(0), F(1), F(0)], [F(0), F(0), F(1)]]
basis_seqs = []
for init in basis_inits:
    seq = iterate_fwd(init, N_MAX+5, 'P27')
    seq_hat = [F(64)**n * seq[n] for n in range(N_MAX+5)]
    basis_seqs.append(seq_hat)

# Source matrix Z0 = [z_b, z_2, z_3] (columns)
# z_b = (b_2, b_1, b_0) = (163, 7, 1)
# z_2 = (b̃_2, b̃_1, b̃_0) = (2145/8, 23/2, 0)
# z_3 = (b̃̃_2, b̃̃_1, b̃̃_0) = (3135/16, 17/2, 0)
Z0 = [
    [F(163), F(2145,8), F(3135,16)],
    [F(7),   F(23,2),   F(17,2)],
    [F(1),   F(0),      F(0)],
]
Z0_inv = mat_inv(Z0)

# det Z0 should be 825/32
det_Z0 = (Z0[0][0]*(Z0[1][1]*Z0[2][2]-Z0[1][2]*Z0[2][1])
         - Z0[0][1]*(Z0[1][0]*Z0[2][2]-Z0[1][2]*Z0[2][0])
         + Z0[0][2]*(Z0[1][0]*Z0[2][1]-Z0[1][1]*Z0[2][0]))
print(f"det Z0 = {det_Z0} (should be 825/32 = {F(825,32)})")

# Target vectors
xq = [F(64)**2 * q_raw[2], F(64) * q_raw[1], q_raw[0]]
xp = [F(64)**2 * p_raw[2], F(64) * p_raw[1], p_raw[0]]

print(f"\nxq = {[float(x) for x in xq]}")
print(f"xp = {[float(x) for x in xp]}")

# Zudilin state vectors at step n
def state_Z(n):
    """3x3 matrix: columns are state vectors of b, b̃, b̃̃ at step n."""
    return [
        [b[n+2], bt[n+2], btt[n+2]],
        [b[n+1], bt[n+1], btt[n+1]],
        [b[n],   bt[n],   btt[n]],
    ]

# Gauge propagation: R(n+1) = (1/r(n)) * C_P(n) * R(n) * C_Z(n)^{-1}
def propagate_R(R0, n_max):
    """Propagate gauge from R(0) to R(n_max)."""
    R_vals = [R0]
    for n in range(n_max):
        CZ_n = CZ_at(n)
        CP_n = CP_at(n)
        CZ_inv = mat_inv(CZ_n)
        R_CZinv = mat_mul(R_vals[-1], CZ_inv)
        R_next = mat_scale(F(1) / r_val(n), mat_mul(CP_n, R_CZinv))
        R_vals.append(R_next)
    return R_vals

# Compute R(0; s) = [xq, s, xp-s] · Z0^{-1}
# where s = (s2, s1, s0) parameterizes the third P2.7 solution
# and Z0^{-1} = inv(source basis matrix)

# With the rank-one twist, R(0) maps Zudilin states to P2.7 states:
#   R(0) · h_0 · X_Z(u, 0) = X_P(v, 0)
# But h_0 = (4)_0^3 / [(5/2)_0(7/2)_0(9/2)_0] = 1/1 = 1
# Actually h_0 = product over k=0,...,-1 = empty product = 1
# So R(0) · Z0 = P0 (target fundamental matrix at n=0)
# => R(0) = P0 · Z0^{-1}

# h_0 = 1 (empty product)
print(f"\nh_0 = 1 (empty Pochhammer product)")

# For each choice of s, R(0) = [xq, s, xp-s] · Z0^{-1}
# Since P0 is linear in s, R(0) is affine in (s0, s1, s2).

# R_fixed = [xq, 0, xp] · Z0^{-1}
P0_fixed = [
    [xq[0], F(0), xp[0]],
    [xq[1], F(0), xp[1]],
    [xq[2], F(0), xp[2]],
]
R0_fixed = mat_mul(P0_fixed, Z0_inv)

# R_delta_k = change from adding e_k to the second column and subtracting from third
# delta_k contribution to P0: column 1 gets e_k, column 2 gets -e_k
# P0_delta_k[i][1] = e_k[i], P0_delta_k[i][2] = -e_k[i]
R0_deltas = []
for k in range(3):
    P0_d = [[F(0)]*3 for _ in range(3)]
    P0_d[k][1] = F(1)
    P0_d[k][2] = F(-1)
    R0_deltas.append(mat_mul(P0_d, Z0_inv))

print("\n=== Propagating gauge with naive s = (0, 0, 0) [i.e., second col = 0] ===")
R_naive = propagate_R(R0_fixed, 15)
for n in range(0, 16, 3):
    maxR = max(abs(float(R_naive[n][i][j])) for i in range(3) for j in range(3))
    print(f"  n={n:3d}: max|R| = {maxR:.4e}")

# Now propagate with each delta separately to understand the linear structure
print("\n=== Propagating each delta separately ===")
for dk in range(3):
    R_dk = propagate_R(R0_deltas[dk], 15)
    for n in [0, 5, 10, 15]:
        maxR = max(abs(float(R_dk[n][i][j])) for i in range(3) for j in range(3))
        print(f"  delta_{dk}, n={n:3d}: max|R| = {maxR:.4e}")

# Optimize: minimize ||R_fixed(n) + s0*R_d0(n) + s1*R_d1(n) + s2*R_d2(n)||
# across multiple n values
print("\n=== Multi-point Frobenius norm optimization ===")

# Propagate all components
R_f = propagate_R(R0_fixed, 20)
R_d = [propagate_R(R0_deltas[k], 20) for k in range(3)]

# Build normal equations for least-squares
# G_kl = Σ_n Σ_ij R_dk(n)_ij * R_dl(n)_ij
# b_k = -Σ_n Σ_ij R_f(n)_ij * R_dk(n)_ij

G = [[F(0)]*3 for _ in range(3)]
bvec = [F(0)]*3

n_range = range(5, 20)
for n in n_range:
    for k in range(3):
        for l in range(3):
            for i in range(3):
                for j in range(3):
                    G[k][l] += R_d[k][n][i][j] * R_d[l][n][i][j]
        for i in range(3):
            for j in range(3):
                bvec[k] -= R_f[n][i][j] * R_d[k][n][i][j]

# Solve G c = b
G_inv = mat_inv(G)
c_opt = mat_vec(G_inv, bvec)
print(f"Optimal s = ({float(c_opt[0]):.15e}, {float(c_opt[1]):.15e}, {float(c_opt[2]):.15e})")

# Check if c_opt are recognizable rationals
for k in range(3):
    val = c_opt[k]
    found = False
    for d in range(1, 10000):
        frac = val * d
        numer_approx = round(float(frac))
        if abs(float(frac - numer_approx)) < 1e-10:
            print(f"  c[{k}] ≈ {numer_approx}/{d}")
            found = True
            break
    if not found:
        print(f"  c[{k}] = {float(val):.20e} (not a simple rational)")

# Compute R with optimal s
R0_opt = [[R0_fixed[i][j] + sum(c_opt[k]*R0_deltas[k][i][j] for k in range(3)) for j in range(3)] for i in range(3)]
R_opt = propagate_R(R0_opt, 20)

print("\n=== Growth with optimal s ===")
for n in range(0, 20):
    maxR = max(abs(float(R_opt[n][i][j])) for i in range(3) for j in range(3))
    logR = log10(mpf(maxR)) if maxR > 0 else mpf(-999)
    print(f"  n={n:3d}: max|R| = {maxR:.6e}  log10 = {float(logR):.2f}")

# Check if R(n) entries for the optimal s are recognizable rational functions
# For each entry R_ij(n), compute the values at n=0,...,19 and try to fit a
# rational function p(n)/q(n)
print("\n=== Checking rationality of R(n) entries ===")
from mpmath import pslq

for i in range(3):
    for j in range(3):
        vals = [float(R_opt[n][i][j]) for n in range(20)]
        # Check if constant
        if all(abs(vals[k] - vals[0]) < 1e-10 * max(1, abs(vals[0])) for k in range(len(vals))):
            print(f"  R[{i},{j}]: CONSTANT ≈ {vals[0]:.6e}")
            continue
        # Check if polynomial (successive differences should vanish)
        diffs = vals[:]
        for deg in range(10):
            new_diffs = [diffs[k+1] - diffs[k] for k in range(len(diffs)-1)]
            if all(abs(d) < 1e-10 * max(1, max(abs(v) for v in vals)) for d in new_diffs):
                print(f"  R[{i},{j}]: degree-{deg} polynomial")
                break
            diffs = new_diffs
        else:
            # Check growth rate
            if all(abs(v) > 0 for v in vals[1:]):
                ratios = [abs(vals[k+1]/vals[k]) for k in range(5, min(15, len(vals)-1)) if abs(vals[k]) > 1e-50]
                if ratios:
                    avg_ratio = sum(ratios)/len(ratios)
                    print(f"  R[{i},{j}]: growth ratio ≈ {avg_ratio:.6f}, vals[0]={vals[0]:.4e}, vals[10]={vals[10]:.4e}")
                else:
                    print(f"  R[{i},{j}]: vals[0]={vals[0]:.4e}")
