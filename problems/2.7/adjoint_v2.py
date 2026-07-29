"""
Adjoint slow solution and bilinear concomitant for P2.7, v2.
Uses MONIC rational operator form (a_3 = 1).

P2.7 monic form (after shifting n -> n+2):
u_{n+3} + a_2(n)*u_{n+2} + a_1(n)*u_{n+1} + a_0(n)*u_n = 0

where a_0(n) = -D(n)/A(n), a_1(n) = C(n+1)/A(n+1), a_2(n) = -B(n+2)/A(n+2)

Adjoint: w_{n-3} + a_2(n-2)*w_{n-2} + a_1(n-1)*w_{n-1} + a_0(n)*w_n = 0
=> w_{n-3} = -a_0(n)*w_n - a_1(n-1)*w_{n-1} - a_2(n-2)*w_{n-2}
"""
import mpmath
mpmath.mp.dps = 200

from mpmath import mpf, mp, zeta
from fractions import Fraction as Q

def A_c(n):
    return int(1024*(2*n+5)**4*(2*n+7)**3*(2*n+9)**3*(946*n*n+6407*n+10860))

def B_c(n):
    return int(128*(2*n+7)**3*(2*n+9)**3*(104060*n**6+1745370*n**5+12145238*n**4+44886481*n**3+92943995*n**2+102256019*n+46709052))

def C_c(n):
    return int(16*(n+3)**4*(2*n+9)**3*(3784*n**5+57792*n**4+351019*n**3+1059230*n**2+1587211*n+944620))

def D_c(n):
    return int((n+3)**4*(n+4)**6*(946*n*n+4515*n+5399))

# Compute q_n, p_n using exact rationals
N = 80
q = [Q(-215040420000), Q(-167282265043404, 905), Q(-964185327658080, 6071)]
p = [Q(-612218384750), Q(-9525021973931919, 18100), Q(-29561828382772029, 65380)]

for i in range(3, N):
    n = i - 1
    Ac, Bc, Cc_prev, Dc_prev2 = A_c(n), B_c(n), C_c(n-1), D_c(n-2)
    new_q = Q(Bc, Ac) * q[-1] - Q(Cc_prev, A_c(n-1)) * q[-2] + Q(Dc_prev2, A_c(n-2)) * q[-3]
    new_p = Q(Bc, Ac) * p[-1] - Q(Cc_prev, A_c(n-1)) * p[-2] + Q(Dc_prev2, A_c(n-2)) * p[-3]
    q.append(new_q)
    p.append(new_p)

# Verify monic standard form
print("=== Verifying monic standard form ===")
for n in range(5):
    a0_val = Q(-D_c(n), A_c(n))
    a1_val = Q(C_c(n+1), A_c(n+1))
    a2_val = Q(-B_c(n+2), A_c(n+2))
    res = q[n+3] + a2_val*q[n+2] + a1_val*q[n+1] + a0_val*q[n]
    print(f"  n={n}: {float(res):.3e}")

# Adjoint backward iteration using mpmath
# w_{n-3} = D(n)/A(n)*w_n - C(n)/A(n)*w_{n-1} + B(n)/A(n)*w_{n-2}
N_back = 300
w = [mpf(0)] * (N_back + 1)
w[N_back] = mpf(1)
w[N_back-1] = mpf(0)
w[N_back-2] = mpf(0)

print(f"\n=== Backward iteration for adjoint slow solution (N_back={N_back}) ===")
for n in range(N_back, 2, -1):
    Dn = mpf(D_c(n))
    Cn = mpf(C_c(n))
    Bn = mpf(B_c(n))
    An = mpf(A_c(n))
    w[n-3] = (Dn*w[n] - Cn*w[n-1] + Bn*w[n-2]) / An

# Normalize w[0] = 1
if w[0] != 0:
    norm = w[0]
    for i in range(len(w)):
        w[i] /= norm

print(f"  w_0 = {w[0]}")
print(f"  w_1 = {w[1]}")
print(f"  w_2 = {w[2]}")
print(f"  w_3 = {w[3]}")

# Verify growth rate
print("\nGrowth ratios w_{n+1}/w_n:")
for n in [5, 10, 20, 30, 50, 70, 100, 150]:
    if n+1 < len(w) and w[n] != 0:
        r = w[n+1]/w[n]
        print(f"  n={n}: {r}")

# Verify adjoint recurrence
print("\n=== Verifying adjoint recurrence ===")
for n in range(3, 10):
    a0_val = mpf(-D_c(n)) / mpf(A_c(n))
    a1_val = mpf(C_c(n)) / mpf(A_c(n))
    a2_val = mpf(-B_c(n)) / mpf(A_c(n))
    res = w[n-3] + a2_val*w[n-2] + a1_val*w[n-1] + a0_val*w[n]
    print(f"  n={n}: {res}")

# Compute e_n = p_n - L*q_n
L_val = zeta(2) + zeta(3)
print(f"\nL = zeta(2)+zeta(3) = {L_val}")

q_mp = [mpf(int(q[i].numerator))/mpf(int(q[i].denominator)) for i in range(N)]
p_mp = [mpf(int(p[i].numerator))/mpf(int(p[i].denominator)) for i in range(N)]
e_mp = [p_mp[i] - L_val*q_mp[i] for i in range(N)]

# Check e_n decay
print("\n=== e_n decay ===")
for n in [0,1,2,5,10,15,20,25]:
    if n < N:
        print(f"  e_{n} = {e_mp[n]}")

# Now compute the bilinear concomitant J(w, u)
# For monic order-3 operator L = E^3 + a_2 E^2 + a_1 E + a_0:
# The concomitant satisfying Δ_n J = w_n(Lu)_n - u_n(L*w)_n is:
#
# J(w,u)_n = w_{n-1}*u_{n+1} + (a_2(n-1)*w_{n-1} + w_{n-2})*u_n
#
# Wait, I need to derive this properly. Let me compute by telescoping.
#
# S(n) = w_n(Lu)_n - u_n(L*w)_n
# = w_n[u_{n+3} + a_2(n)u_{n+2} + a_1(n)u_{n+1} + a_0(n)u_n]
# - u_n[w_{n-3} + a_2(n-2)w_{n-2} + a_1(n-1)w_{n-1} + a_0(n)w_n]
#
# = w_n u_{n+3} - u_n w_{n-3}
# + a_2(n) w_n u_{n+2} - a_2(n-2) u_n w_{n-2}
# + a_1(n) w_n u_{n+1} - a_1(n-1) u_n w_{n-1}
#
# We need S(n) = J(n+1) - J(n).
#
# Define P_n = w_n u_{n+2}, Q_n = w_n u_{n+1}. Then:
# w_n u_{n+3} = P_{n+1} - ... nah, let me try the matrix approach.
#
# Actually, I'll just verify constancy numerically by computing partial sums.

# Direct constancy test: for L-solution u and L*-solution w,
# define S_N = sum_{n=0}^{N-1} [w_n(Lu)_n - u_n(L*w)_n]
# Since Lu=0 and L*w=0, S_N = 0 for all N.
# So J(w,u)_N = J(w,u)_0 for all N.

# Instead of finding the explicit J, let me use a DIFFERENT approach:
# The Casoratian bracket for the monic operator.
# For the monic form with companion matrix:
# C(n) = [[0, 1, 0], [0, 0, 1], [-a_0(n), -a_1(n), -a_2(n)]]
# X_n = (u_n, u_{n+1}, u_{n+2})^T
# X_{n+1} = C(n) X_n
# For adjoint: Y_n = (w_n, w_{n-1}, w_{n-2})^T
# Y_n^T C(n) X_n = Y_{n+1}^T X_{n+1}? No...
# The adjoint companion: the transpose of the INVERSE of C(n) relates Y.

# Actually, for monic order-3: det C(n) = -a_0(n) = D(n)/A(n).
# The bilinear concomitant is: J(w,u)_n = Y_n^T · M · X_n
# for some matrix M that is constant (related to the symplectic form).

# For order 3, the correct concomitant is the Lagrangian form:
# J = sum of terms involving w and u at consecutive indices.
# Let me try the CASORATIAN approach directly.

# For two solutions u, v of L, the Casoratian is:
# W(u,v)_n = u_n v_{n+1} - u_{n+1} v_n
# For third-order, the full Casoratian involves a 3x3 determinant.

# The CORRECT bilinear form for the adjoint theory:
# Consider the pairing <w, u> = w_n u_{n+2} - (sum involving a_2) + ...
# Let me just try several candidate formulas and test constancy.

print("\n=== Testing candidate concomitant formulas ===")

# Candidate 1: J1 = w_{n-1}*u_{n+1} - w_n*u_n (order-2 style, likely wrong for order 3)
print("Candidate 1: w_{n-1}*u_{n+1} - w_n*u_n")
for n in range(3, 8):
    J1 = w[n-1]*q_mp[n+1] - w[n]*q_mp[n]
    print(f"  n={n}: {J1}")

# Candidate 2: J2 = w[n]*u[n+2] - w[n+1]*u[n+1] + w[n+2]*u[n]
print("\nCandidate 2: w_n*u_{n+2} - w_{n+1}*u_{n+1} + w_{n+2}*u_n")
for n in range(1, 8):
    J2 = w[n]*q_mp[n+2] - w[n+1]*q_mp[n+1] + w[n+2]*q_mp[n]
    print(f"  n={n}: {J2}")

# Candidate 3: a_0(n)*w_n*u_{n+2}*... more complex
# For the order-3 monic operator L = E^3 + a_2 E^2 + a_1 E + a_0:
# J(w,u)_n = w_{n-2}*u_{n+1} - w_{n-1}*(u_n + a_2(n-1)*u_{n+1}) + ... no

# Let me derive step by step. We need Δ_n J = w_n(Lu)_n - u_n(L*w)_n = 0.
# So J_n must be constant. It's a bilinear form in
# (w_{n-3}, w_{n-2}, w_{n-1}) and (u_n, u_{n+1}, u_{n+2}).
# [These are the "state variables" at step n for L* and L respectively.]

# Actually, the standard approach: define the state vectors:
# X_n = (u_n, u_{n+1}, u_{n+2})  for the forward operator
# The companion matrix C(n): X_{n+1} = C(n) X_n, where
# C(n) = [[0, 1, 0], [0, 0, 1], [-a_0(n), -a_1(n), -a_2(n)]]

# For the adjoint, define Y_n = some vector involving w. The bilinear
# invariant is J_n = Y_n^T S X_n for some fixed matrix S such that
# J_{n+1} = J_n when L u = 0 and L* w = 0.

# For this to work: Y_{n+1}^T S X_{n+1} = Y_n^T S X_n
# => Y_{n+1}^T S C(n) X_n = Y_n^T S X_n for all X_n
# => Y_{n+1}^T S C(n) = Y_n^T S
# => Y_{n+1}^T = Y_n^T S C(n)^{-1} S^{-1}

# The adjoint state vector Y relates to w. For the standard pairing,
# S is the "flip" matrix (anti-diagonal identity):
# S = [[0, 0, 1], [0, -1, 0], [1, 0, 0]]
# Then J_n = Y_n^T S X_n = w_{n-3}*u_{n+2} - w_{n-2}*u_{n+1} + w_{n-1}*u_n
# Let me test this!

print("\nCandidate 3 (anti-diagonal): w_{n-3}*u_{n+2} - w_{n-2}*u_{n+1} + w_{n-1}*u_n")
J3_vals = []
for n in range(3, 15):
    J3 = w[n-3]*q_mp[n+2] - w[n-2]*q_mp[n+1] + w[n-1]*q_mp[n]
    J3_vals.append((n, J3))
    if n <= 10:
        print(f"  n={n}: {J3}")
print(f"  Variation: {max(abs(J3_vals[i][1] - J3_vals[0][1]) for i in range(min(8, len(J3_vals))))}")

# That probably won't work because it ignores the a_2 coefficient.
# The correct formula needs to account for the companion matrix structure.

# Let me try the MATRIX approach directly.
# C(n) = [[0, 1, 0], [0, 0, 1], [-a_0(n), -a_1(n), -a_2(n)]]
# det C(n) = a_0(n) = -D(n)/A(n)
# C(n)^{-1} = (1/det) * adj(C) = ...

# For 3x3: C^{-1} = (1/a_0) * [[-a_2, -1, 0], [a_0, 0, 0], [0, a_0, 0]]... hmm no.
# C = [[0,1,0],[0,0,1],[-a0,-a1,-a2]]
# det C = 0*(0*(-a2)-1*(-a1)) - 1*(0*(-a2)-1*(-a0)) + 0 = -(0+a0) = -a0 WAIT
# det C = 0*(0*(-a2) - 1*(-a1)) - 1*(0*(-a2) - 1*(-a0)) + 0*(0*(-a1) - 0*(-a0))
#       = 0 - 1*(0 - (-a0)) + 0 = -1*a0 = -a0
# Wait: a0 = -D/A, so det C = -(-D/A) = D/A.

# Forget the matrix approach. Let me just numerically find the right bilinear form.
# I know J_n = α * w_{n-3}u_{n+2} + β * w_{n-2}u_{n+1} + γ * w_{n-1}u_n
#            + δ * w_{n-3}u_{n+1} + ε * w_{n-2}u_n + ... etc.
# This has up to 9 terms (3 w-indices × 3 u-indices).
# With the constraint that J is constant, I can solve for the coefficients.

# Set up: for each n, J_n = sum_{i,j} c_{ij} w_{n-3+i} u_{n+j} for 0≤i≤2, 0≤j≤2
# Then J_{n+1} - J_n should be 0 when Lu=0 and L*w=0.
# This gives a system for the c_{ij} (which may depend on n!).

# Actually, for a non-self-adjoint operator, J_n DOES involve n-dependent
# coefficients. The bilinear form has coefficients that are rational functions of n.

# The simplest approach: compute J via the CUMULATIVE SUM directly.
# J_n = J_0 + sum_{m=0}^{n-1} S(m) where S(m) = w_m(Lu)_m - u_m(L*w)_m = 0.
# So J_n = J_0 = const. But I need to know J_0 = sum involving w and u at small indices.

# PRAGMATIC approach: just compute c_0(e) directly from the asymptotic behavior.
# c_0(e) = lim_{n->inf} e_n / (rho_0^n)
# If c_0 = 0, this limit is 0.

# Compute rho_0 (dominant Poincaré root of P2.7)
# Poincaré polynomial: 1048576*rho^3 - 901120*rho^2 + 512*rho - 1 = 0
# Or equivalently: 4*nu^3 - 220*nu^2 + 8*nu - 1 = 0 where rho = nu/64

import numpy as np
p_poly = [4, -220, 8, -1]  # 4x^3 - 220x^2 + 8x - 1
nu_roots = np.roots(p_poly)
print(f"\nPoincaré roots nu: {nu_roots}")
rho_roots = nu_roots / 64
print(f"Poincaré roots rho = nu/64: {rho_roots}")

nu0 = max(nu_roots.real)  # dominant root
rho0 = mpf(nu0) / 64
print(f"\nrho_0 = {rho0}")

# Compute e_n / rho_0^n
print("\n=== c_0 test: e_n / rho_0^n ===")
for n in range(min(70, N)):
    r0n = rho0 ** n
    ratio = e_mp[n] / r0n
    if n <= 10 or n % 10 == 0:
        print(f"  n={n}: e_n/rho_0^n = {ratio}")

# Also compute log|e_n| / n to see the effective decay rate
print("\n=== Effective decay rate: log|e_n|/n ===")
for n in range(1, min(60, N)):
    if e_mp[n] != 0:
        rate = mpmath.log(abs(e_mp[n])) / n
        if n <= 10 or n % 5 == 0:
            print(f"  n={n}: log|e_n|/n = {rate} (rho_1 gives {mpmath.log(abs(rho_roots[1]))})")

# Compute |rho_1| for comparison
rho1_abs = abs(complex(rho_roots[1]))
print(f"\n|rho_1| = {rho1_abs}")
print(f"log|rho_1| = {np.log(rho1_abs)}")
print(f"rho_0 = {float(rho0)}")
print(f"log(rho_0) = {float(mpmath.log(rho0))}")

# If e_n/rho_0^n -> 0 exponentially, then c_0(e) = 0.
# The ratio e_n/rho_0^n should decay at rate (rho_1/rho_0)^n ~ (0.001/0.859)^n ~ 0.00122^n
