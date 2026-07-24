#!/usr/bin/env python3
"""
Search for a rational matrix gauge R(n) between h-twisted Zudilin and P2.7.

Approach: parameterize the third P2.7 solution with 3 free parameters (s0,s1,s2),
compute R(n) = Phi_P(n) · [D_h(n) · Phi_Z(n)]^{-1} for n=0,...,N,
and search for (s0,s1,s2) that makes R(n) polynomially bounded (not exponentially growing).
"""
from mpmath import mp, mpf, matrix as mpm, eye, det, lu_solve, fabs, log10, nstr
from fractions import Fraction as F

mp.dps = 100

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

# h_n = (4)_n^3 / [(5/2)_n (7/2)_n (9/2)_n]
def h(n):
    val = mpf(1)
    for k in range(n):
        val *= mpf(4+k)**3 / (mpf(5,2)+k) / (mpf(7,2)+k) / (mpf(9,2)+k)
    return val

N_MAX = 30

# Forward iteration
def zudilin_fwd(init, N):
    u = list(init)
    for n in range(2, N):
        d = mpf(2 * QZ(n) * (2*n+1) * (n+1)**3)
        nxt = mpf(2*MZ(n)) * u[n] + mpf(-2*n*NZ(n)) * u[n-1] + mpf(RZ(n)*n*(n-1)**3) * u[n-2]
        u.append(nxt / d)
    return u

def p27_fwd(init, N):
    u = list(init)
    for n in range(2, N):
        nxt = mpf(B_c(n))/mpf(A_c(n)) * u[n] - mpf(C_c(n-1))/mpf(A_c(n-1)) * u[n-1] + mpf(D_c(n-2))/mpf(A_c(n-2)) * u[n-2]
        u.append(nxt)
    return u

# Compute all sequences
b  = zudilin_fwd([mpf(1), mpf(7), mpf(163)], N_MAX+5)
bt = zudilin_fwd([mpf(0), mpf(23)/2, mpf(2145)/8], N_MAX+5)
btt= zudilin_fwd([mpf(0), mpf(17)/2, mpf(3135)/16], N_MAX+5)

q0 = mpf(-215040420000)
q1 = mpf(-167282265043404) / mpf(905)
q2 = mpf(-964185327658080) / mpf(6071)
p0 = mpf(-612218384750)
p1 = mpf(-9525021973931919) / mpf(18100)
p2 = mpf(-29561828382772029) / mpf(65380)

q_raw = p27_fwd([q0, q1, q2], N_MAX+5)
p_raw = p27_fwd([p0, p1, p2], N_MAX+5)

# Scaled P2.7 sequences: q̂_n = 64^n q_n, p̂_n = 64^n p_n
qhat = [mpf(64)**n * q_raw[n] for n in range(N_MAX+5)]
phat = [mpf(64)**n * p_raw[n] for n in range(N_MAX+5)]

# h values
h_vals = [h(n) for n in range(N_MAX+5)]

# D_h(n) = diag(h(n+2), h(n+1), h(n))
def Dh(n):
    return mpm([[h_vals[n+2], 0, 0], [0, h_vals[n+1], 0], [0, 0, h_vals[n]]])

# State vectors: X(u, n) = (u_{n+2}, u_{n+1}, u_n)^T
def state(u, n):
    return mpm([u[n+2], u[n+1], u[n]])

# Zudilin fundamental matrix (columns = state vectors of b, b̃, b̃̃)
def Phi_Z(n):
    return mpm([
        [b[n+2], bt[n+2], btt[n+2]],
        [b[n+1], bt[n+1], btt[n+1]],
        [b[n],   bt[n],   btt[n]]
    ])

# P2.7 fundamental matrix (columns = state vectors of q̂, p̂, s)
# s is the third solution parameterized by (s0, s1, s2)
def Phi_P(n, s):
    return mpm([
        [qhat[n+2], phat[n+2], s[n+2]],
        [qhat[n+1], phat[n+1], s[n+1]],
        [qhat[n],   phat[n],   s[n]]
    ])

# Compute R(n) = Phi_P(n) · [D_h(n) · Phi_Z(n)]^{-1}
def R_matrix(n, s):
    M = Dh(n) * Phi_Z(n)
    return Phi_P(n, s) * M**(-1)

# === Test 1: s = (1, 0, 0) — the naive third solution ===
print("=== Test 1: s = (1, 0, 0) ===")
s_naive = p27_fwd([mpf(1), mpf(0), mpf(0)], N_MAX+5)
s_naive_hat = [mpf(64)**n * s_naive[n] for n in range(N_MAX+5)]

for n in [0, 5, 10, 15, 20]:
    R = R_matrix(n, s_naive_hat)
    maxR = max(fabs(R[i,j]) for i in range(3) for j in range(3))
    print(f"  n={n}: max|R(n)| = {nstr(maxR, 10)}, log10 = {nstr(log10(maxR), 5)}")

# === Test 2: try s = b̃̃ mapped to P2.7 ===
# Hypothesis: if b̃̃ maps to some specific P2.7 solution, that solution
# should make R(n) rational.
# But we don't know what b̃̃ maps to!

# === Approach: search over (s0, s1, s2) to minimize max|R(20)| ===
# R(n) is LINEAR in s. So R(20) = A + s0*B0 + s1*B1 + s2*B2
# where A is the contribution from the known (q̂, p̂) columns and
# B0, B1, B2 are the contributions from basis third solutions.

print("\n=== Computing R(n) dependence on (s0, s1, s2) ===")

# Compute R for basis third solutions
basis_inits = [
    [mpf(1), mpf(0), mpf(0)],
    [mpf(0), mpf(1), mpf(0)],
    [mpf(0), mpf(0), mpf(1)],
]
basis_seqs = []
for init in basis_inits:
    seq = p27_fwd(init, N_MAX+5)
    seq_hat = [mpf(64)**n * seq[n] for n in range(N_MAX+5)]
    basis_seqs.append(seq_hat)

# R_basis[k] = R(n) when third column uses basis_seqs[k]
# R(n, s) = R(n, 0_third) + Σ s_k * R(n, basis_k_third_only)

# Actually, it's easier: Phi_P has the third column as the s sequence.
# R = Phi_P · M^{-1} where M = D_h · Phi_Z
# R_ij = Σ_l Phi_P_il · M^{-1}_lj
# The third column of Phi_P is s, so R depends linearly on (s_{n+2}, s_{n+1}, s_n).
# And s is linear in (s0, s1, s2).

# More precisely: R(n) = R_fixed(n) + s0 * R_e0(n) + s1 * R_e1(n) + s2 * R_e2(n)
# where R_fixed uses s=0 (impossible since Phi_P becomes rank 2), so let's use
# R with the i-th basis as the third column and subtract.

# R(n, s=c0*e0+c1*e1+c2*e2) = c0*R(n,e0) + c1*R(n,e1) + c2*R(n,e2)
# But this isn't right because the first two columns are FIXED.

# Let me rewrite: Phi_P(n, s) = [state(q̂), state(p̂), state(s)]
# R(n, s) = Phi_P(n, s) · M(n)^{-1}
# The dependence on s is through the third column of Phi_P.

# R_ij(n, s) = Φ_P(n)_i1 * Minv(n)_1j + Φ_P(n)_i2 * Minv(n)_2j + s_{n+2-i} * Minv(n)_3j
# Wait, let me index properly.
# Φ_P = [[q̂_{n+2}, p̂_{n+2}, s_{n+2}],
#         [q̂_{n+1}, p̂_{n+1}, s_{n+1}],
#         [q̂_n, p̂_n, s_n]]
# So Φ_P[0,2] = s_{n+2}, Φ_P[1,2] = s_{n+1}, Φ_P[2,2] = s_n
# R = Φ_P · Minv
# R[i,j] = Σ_l Φ_P[i,l] * Minv[l,j]
#         = Φ_P[i,0]*Minv[0,j] + Φ_P[i,1]*Minv[1,j] + Φ_P[i,2]*Minv[2,j]
# The first two terms are fixed (independent of s).
# R[i,j] = R_fixed[i,j] + s_{n+2-i} * Minv[2,j]
#
# So R is AFFINE in the state vector of s at step n.
# And the state vector of s at step n is a LINEAR combination of (s_0, s_1, s_2)
# via the recurrence propagation matrix.

# Let me just compute R for the three basis cases and the zero case.

# R with zero third column (Phi_P singular, but we can still compute the contribution)
# Actually, R_fixed = [q̂ state, p̂ state, 0] · Minv = just the first two columns of R contribute
# Let me compute it properly.

n_target = 20

# Compute M inverse at n_target
M_n = Dh(n_target) * Phi_Z(n_target)
Minv = M_n**(-1)

print(f"\nMinv(n={n_target}) third row:")
for j in range(3):
    print(f"  Minv[2,{j}] = {nstr(Minv[2,j], 15)}")

# For each basis third solution, compute R(n_target):
for k in range(3):
    R_k = R_matrix(n_target, basis_seqs[k])
    maxR = max(fabs(R_k[i,j]) for i in range(3) for j in range(3))
    print(f"\nBasis e{k}: max|R({n_target})| = {nstr(maxR, 10)}, log10 = {nstr(log10(maxR), 5)}")
    print(f"  R entries:")
    for i in range(3):
        for j in range(3):
            print(f"    R[{i},{j}] = {nstr(R_k[i,j], 12)}")

# Now find the linear combination s = c0*e0 + c1*e1 + c2*e2 that minimizes max|R(n)|
# for several values of n simultaneously.

# Strategy: for R(n) to be bounded polynomially, the exponential growth must cancel.
# The dominant solution grows like ρ₀^n ≈ 55^n, so after 64^n scaling, it grows like (55*64)^n.
# This contamination must vanish.

# Better approach: compute R(n) for n = n_target using s parameterized by (c0, c1, c2),
# and find the (c0, c1, c2) that minimizes the Frobenius norm of R.

# R(n, c) = R_fixed(n) + c0*R_delta0(n) + c1*R_delta1(n) + c2*R_delta2(n)
# where R_fixed uses s=0 (but that makes Phi_P singular!).

# Alternative: just use R(n, e_k) for k=0,1,2 as the three "basis R matrices".
# Then R(n, s=c0*e0+c1*e1+c2*e2) = c0*R(n,e0) + c1*R(n,e1) + c2*R(n,e2)
# This is because R depends linearly on the third column of Phi_P,
# and the third column of Phi_P is the state vector of s, which is linear in (c0,c1,c2).
# Wait, but Phi_P also has the FIRST two columns as q̂ and p̂, which don't depend on s.
# So R = (qhat_state, phat_state, s_state) · Minv
#      = qhat_state ⊗ Minv[0,:] + phat_state ⊗ Minv[1,:] + s_state ⊗ Minv[2,:]
# The first two terms are fixed. The third term is linear in s_state,
# which is linear in (c0, c1, c2).

# So: R(n, c) = R_0(n) + c0*Delta_0(n) + c1*Delta_1(n) + c2*Delta_2(n)
# where R_0 is R with s=0, and Delta_k is the change when adding e_k.

# But R with s=0 makes Phi_P singular, so det(Phi_P)=0 and R may not be well-defined.
# Actually R = Phi_P · Minv, and if Phi_P is singular that just means R is singular.
# The product is still well-defined (just not invertible).

# R_0 = [qhat_state, phat_state, 0] · Minv
R0 = mpm(3, 3)
for i in range(3):
    qh_i = [qhat[n_target+2], qhat[n_target+1], qhat[n_target]][i]
    ph_i = [phat[n_target+2], phat[n_target+1], phat[n_target]][i]
    for j in range(3):
        R0[i,j] = qh_i * Minv[0,j] + ph_i * Minv[1,j]

# Delta_k: contribution from e_k to the third column
# Delta_k[i,j] = (e_k state at n_target)_i * Minv[2,j]
Deltas = []
for k in range(3):
    Dk = mpm(3, 3)
    s_k = basis_seqs[k]
    for i in range(3):
        s_i = [s_k[n_target+2], s_k[n_target+1], s_k[n_target]][i]
        for j in range(3):
            Dk[i,j] = s_i * Minv[2,j]
    Deltas.append(Dk)

# Now minimize ||R_0 + c0*D0 + c1*D1 + c2*D2||_F^2 over (c0, c1, c2)
# This is a linear least squares problem.

# Vectorize: let x = [R0_ij], delta_k = [Dk_ij], c = [c0, c1, c2]
# ||x + Σ c_k delta_k||^2 = ||x||^2 + 2 Σ c_k (x·delta_k) + Σ_k,l c_k c_l (delta_k · delta_l)
# Minimize: Σ_l (delta_k · delta_l) c_l = -(x · delta_k) for each k
# This is a 3x3 linear system: G c = -b where G_kl = delta_k · delta_l, b_k = x · delta_k

G = mpm(3, 3)
bvec = mpm(3, 1)
for k in range(3):
    for l in range(3):
        s = mpf(0)
        for i in range(3):
            for j in range(3):
                s += Deltas[k][i,j] * Deltas[l][i,j]
        G[k,l] = s
    s = mpf(0)
    for i in range(3):
        for j in range(3):
            s += R0[i,j] * Deltas[k][i,j]
    bvec[k,0] = -s

c_opt = lu_solve(G, bvec)
print(f"\n=== Optimal (c0, c1, c2) at n={n_target} ===")
for k in range(3):
    print(f"  c{k} = {nstr(c_opt[k,0], 30)}")

# Compute R at n_target with optimal c
R_opt = mpm(3, 3)
for i in range(3):
    for j in range(3):
        R_opt[i,j] = R0[i,j]
        for k in range(3):
            R_opt[i,j] += c_opt[k,0] * Deltas[k][i,j]

print(f"\n  R(n={n_target}) with optimal c:")
for i in range(3):
    for j in range(3):
        print(f"    R[{i},{j}] = {nstr(R_opt[i,j], 20)}")
maxR = max(fabs(R_opt[i,j]) for i in range(3) for j in range(3))
print(f"  max|R| = {nstr(maxR, 10)}")

# Now compute the optimal third solution
s_opt_init = [c_opt[0,0], c_opt[1,0], c_opt[2,0]]
s_opt = p27_fwd(s_opt_init, N_MAX+5)
s_opt_hat = [mpf(64)**n * s_opt[n] for n in range(N_MAX+5)]

# Check R(n) growth for the optimal choice across all n
print(f"\n=== R(n) growth with optimal third solution ===")
for n in range(0, min(N_MAX, 25)):
    try:
        R_n = R_matrix(n, s_opt_hat)
        maxR = max(fabs(R_n[i,j]) for i in range(3) for j in range(3))
        logR = log10(maxR) if maxR > 0 else mpf(-999)
        print(f"  n={n:3d}: max|R| = {nstr(maxR, 8):>20s}  log10 = {nstr(logR, 5)}")
    except Exception as e:
        print(f"  n={n:3d}: ERROR: {e}")

# === Multi-point optimization: minimize max growth across n=5,...,25 ===
print("\n\n=== Multi-point optimization: minimize sum of log|R(n)| for n=5,...,20 ===")

# For each n, compute R0(n) and Delta_k(n)
# Then form a big least-squares system

n_range = list(range(5, 21))
# Accumulate A^T A and A^T b for the full system
GG = mpm(3, 3)
bb = mpm(3, 1)

for n in n_range:
    M_n = Dh(n) * Phi_Z(n)
    Minv_n = M_n**(-1)

    # R0(n)
    R0_n = mpm(3, 3)
    for i in range(3):
        qh_i = [qhat[n+2], qhat[n+1], qhat[n]][i]
        ph_i = [phat[n+2], phat[n+1], phat[n]][i]
        for j in range(3):
            R0_n[i,j] = qh_i * Minv_n[0,j] + ph_i * Minv_n[1,j]

    # Delta_k(n)
    Dk_n = []
    for k in range(3):
        D = mpm(3, 3)
        sk = basis_seqs[k]
        for i in range(3):
            s_i = [sk[n+2], sk[n+1], sk[n]][i]
            for j in range(3):
                D[i,j] = s_i * Minv_n[2,j]
        Dk_n.append(D)

    # Accumulate
    for kk in range(3):
        for ll in range(3):
            s = mpf(0)
            for i in range(3):
                for j in range(3):
                    s += Dk_n[kk][i,j] * Dk_n[ll][i,j]
            GG[kk,ll] += s
        s = mpf(0)
        for i in range(3):
            for j in range(3):
                s += R0_n[i,j] * Dk_n[kk][i,j]
        bb[kk,0] -= s

c_multi = lu_solve(GG, bb)
print(f"Optimal c (multi-point):")
for k in range(3):
    print(f"  c{k} = {nstr(c_multi[k,0], 30)}")

# Compute optimal s sequence
s_multi_init = [c_multi[0,0], c_multi[1,0], c_multi[2,0]]
s_multi = p27_fwd(s_multi_init, N_MAX+5)
s_multi_hat = [mpf(64)**n * s_multi[n] for n in range(N_MAX+5)]

print(f"\n  R(n) growth with multi-point optimal:")
for n in range(0, min(N_MAX, 25)):
    try:
        R_n = R_matrix(n, s_multi_hat)
        maxR = max(fabs(R_n[i,j]) for i in range(3) for j in range(3))
        logR = log10(maxR) if maxR > 0 else mpf(-999)
        print(f"  n={n:3d}: max|R| = {nstr(maxR, 8):>20s}  log10 = {nstr(logR, 5)}")
    except Exception as e:
        print(f"  n={n:3d}: ERROR: {e}")

# Try PSLQ on the optimal c values — are they simple rationals?
print("\n=== PSLQ on optimal c values ===")
from mpmath import pslq as mp_pslq
for k in range(3):
    val = c_multi[k,0]
    # Check if it's rational with small denominator
    for d in range(1, 500):
        frac = val * d
        rounded = round(float(frac))
        if fabs(frac - rounded) < mpf(10)**(-50):
            print(f"  c{k} ≈ {rounded}/{d}")
            break
    else:
        print(f"  c{k} = {nstr(val, 40)} (not a simple rational)")
