#!/usr/bin/env python3
"""Extract intermediate-mode ratio r_int(n) via dominant-mode cancellation.

1. Compute two forward solutions u1, u2 in high precision
2. Ratio u1(n)/u2(n) → A1/A2 exponentially fast (rate ~(16/543)^n)
3. Form w = u1 - (A1/A2)*u2 → no dominant mode
4. w(n+1)/w(n) → r_int(n) exponentially fast (rate ~(0.47/16)^n)
5. Reconstruct rational function from integer evaluations
"""
from mpmath import mp, mpf, nstr, matrix as mp_matrix
mp.dps = 800

# Recurrence coefficients from exact QQ (re-derive via polynomial evaluation at mpf)
# Using the EXACT polynomial coefficients from the null-space computation
# I'll use the 3x3 matrix recurrence directly for better numerical stability

def M_mat(n):
    """3x3 matrix M(n) over mpf."""
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

# The scalar sequence q_N is defined as A * T_N where T_N = M(0)*...*M(N-1)
# A = [33750, -36000, 9000] (second row of the initial matrix)
# q_N = A * T_N * e1 = A * prod(M(k), k=0..N-1) * [1,0,0]^T

# Compute using three independent initial conditions for the scalar recurrence
# via the matrix product. The scalar recurrence is:
# c3(n)u(n+3) + c2(n)u(n+2) + c1(n)u(n+1) + c0(n)u(n) = 0

# But I'll use the matrix directly. The scalar sequence u(N) = v * T_N * e1
# where v is any row vector. Different v's give different solutions.

# Use v1 = [1,0,0], v2 = [0,1,0], v3 = [0,0,1]
N_max = 55

# Compute T_N column by column
# T_N e1 is the first column of T_N = M(0)*...*M(N-1)
# We compute T_N e_j for j=1,2,3

print("Computing three forward solutions via matrix product...")

# We track col = [a, b, c] representing T_N * e_j
# u_j(N) = v_i * col = col[i] (for v_i = e_i)

# Actually, let's track the FULL 3x3 matrix T_N
# T(N+1) = T(N) * M(N)

T = [[mpf(1) if i==j else mpf(0) for j in range(3)] for i in range(3)]
# T is identity at N=0

# Store the scalar sequences: u_ij(N) = e_i^T * T_N * e_j
# We need u for all i, and we'll use the first column (j=0)
# u_i(N) = T[i][0] (the (i,0) entry of T_N)

u = [[] for _ in range(3)]  # u[i][n] = T_n[i][0]

for N in range(N_max + 3):
    for i in range(3):
        u[i].append(T[i][0])
    
    if N < N_max + 2:
        M = M_mat(N)
        T_new = [[mpf(0)]*3 for _ in range(3)]
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    T_new[i][j] += T[i][k] * M[k][j]
        T = T_new
    
    if N % 10 == 0:
        print(f"  N={N} done, |u0|={nstr(abs(u[0][-1]),5)}")

print(f"Computed {len(u[0])} terms of each solution")

# Now compute scalar recurrence coefficients from the matrix
# The scalar recurrence for each u_i uses the SAME recurrence (different initial conditions)
# Let me compute the recurrence from the matrix product directly

# Actually, the key property: all three u_i satisfy the SAME scalar recurrence
# (they are the three components of the matrix product, which all satisfy the
#  minimal polynomial of the transfer matrix)

# Step 1: Compute dominant ratio u0(n)/u1(n) for large n
print("\nStep 1: Dominant cancellation")
print("u0(n)/u1(n) for large n:")
for n in range(35, 51):
    if abs(u[1][n]) > 0:
        ratio = u[0][n] / u[1][n]
        print(f"  n={n}: {nstr(ratio, 40)}")

# The ratio should converge. Let R = u0(N)/u1(N) at large N
R = u[0][50] / u[1][50]
print(f"\nUsing R = u0(50)/u1(50) = {nstr(R, 50)}")

# Step 2: Form w = u0 - R*u1 (dominant mode cancelled)
w = [u[0][n] - R * u[1][n] for n in range(len(u[0]))]

# Check: w(n+1)/w(n) should approach r_int for large enough n
print("\nStep 2: Checking w(n+1)/w(n) [should → r_int(n)]:")
for n in range(5, 45):
    if abs(w[n]) > 0:
        ratio = w[n+1] / w[n]
        # Compare with -16*n^7
        r_over_n7 = ratio / (mpf(-16) * mpf(n)**7)
        print(f"  n={n}: w(n+1)/w(n) = {nstr(ratio, 20)}  ratio/(-16*n^7) = {nstr(r_over_n7, 15)}")

# The w(n+1)/w(n) should be r_int(n) once the recessive component is negligible
# Recessive decays like (c_rec/c_int)^n ≈ (0.47/16)^n ≈ 0.029^n
# At n=30: 0.029^30 ≈ 10^{-46}, so convergence to ~46 digits

# Step 3: Further cancellation to remove recessive mode
# Use ANOTHER ratio: form w2 = u0 - R2*u2 and combine w, w2
print("\n\nStep 3: Double cancellation")
R2 = u[0][50] / u[2][50]  # dominant ratio with u2
w2 = [u[0][n] - R2 * u[2][n] for n in range(len(u[0]))]

# Now w and w2 both have no dominant mode (approximately)
# w = B_w * h_int + C_w * h_rec
# w2 = B_w2 * h_int + C_w2 * h_rec
# Cancel recessive: ww = w - (C_w/C_w2)*w2 = pure h_int

# The recessive coefficient ratio C_w/C_w2 can be found from the ratio
# w(n)/w2(n) at the SMALLEST n (where recessive is largest relative to intermediate)
# Actually at large n, w/w2 → B_w/B_w2. At small n, it's mixed.

# Better: use w(n+1)/w(n) at two different n-ranges to separate
# At n and n+1: w(n+1)/w(n) = (B h_int(n+1) + C h_rec(n+1))/(B h_int(n) + C h_rec(n))
# = r_int(n) * (1 + (C/B)(h_rec/h_int)(n+1)/(1 + (C/B)(h_rec/h_int)(n)))

# Since h_rec/h_int ~ (c_rec/c_int)^n → 0:
# w(n+1)/w(n) ≈ r_int(n) * (1 - (C/B)(c_rec/c_int)^n * δ(n))
# where δ(n) captures the polynomial corrections

# At n=30: correction ~ 0.029^30 ~ 10^{-46}
# At n=40: correction ~ 0.029^40 ~ 10^{-62}

# So r_int(n) = w(n+1)/w(n) to about min(46, 800-digits-lost) digits at n=30

# Let me print r_int values with error estimate
print("\nr_int(n) = w(n+1)/w(n) for n=0..40:")
for n in range(41):
    if abs(w[n]) > mpf('1e-500'):
        r_int = w[n+1] / w[n]
        print(f"  r_int({n}) = {nstr(r_int, 60)}")

