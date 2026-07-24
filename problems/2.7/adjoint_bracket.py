"""
Compute the adjoint slow solution w^{(0)} of P2.7 via Miller backward iteration,
then verify c_0(e)=0 through the bilinear concomitant J(w^{(0)}, e).

The P2.7 recurrence (standard form, shifted):
a_0(n)*u_n + a_1(n)*u_{n+1} + a_2(n)*u_{n+2} + a_3(n)*u_{n+3} = 0

where a_0(n) = -D(n), a_1(n) = C(n+1), a_2(n) = -B(n+2), a_3(n) = A(n+2).

The adjoint:
L*[w]_n = a_0(n)*w_n + a_1(n-1)*w_{n-1} + a_2(n-2)*w_{n-2} + a_3(n-3)*w_{n-3} = 0

Bilinear concomitant J(w,u)_n = const when Lu=0 and L*w=0.
Explicit formula for order 3 (derived from Green-Lagrange identity):
J(w,u)_n = a_3(n) * [w_{n+2}*u_n - w_{n+1}*u_{n+1} + w_n*u_{n+2}]
          + a_2(n) * [w_{n+1}*u_n - w_n*u_{n+1}]
          + a_1(n) * w_n * u_n
Wait, this needs to be derived carefully. Let me use a different approach.
"""
import mpmath
mpmath.mp.dps = 100

from mpmath import mpf, mp, pi, zeta
from fractions import Fraction as Q
from math import comb

def A_c(n):
    n = Q(n)
    return (1024*(2*n+5)**4*(2*n+7)**3*(2*n+9)**3*(946*n*n+6407*n+10860))

def B_c(n):
    n = Q(n)
    return (128*(2*n+7)**3*(2*n+9)**3*(104060*n**6+1745370*n**5+12145238*n**4+44886481*n**3+92943995*n**2+102256019*n+46709052))

def C_c(n):
    n = Q(n)
    return (16*(n+3)**4*(2*n+9)**3*(3784*n**5+57792*n**4+351019*n**3+1059230*n**2+1587211*n+944620))

def D_c(n):
    n = Q(n)
    return ((n+3)**4*(n+4)**6*(946*n*n+4515*n+5399))

# Standard form coefficients:
# a_0(n)*u_n + a_1(n)*u_{n+1} + a_2(n)*u_{n+2} + a_3(n)*u_{n+3} = 0
# From: A_c(n)*u_{n+1} - B_c(n)*u_n + C_c(n-1)*u_{n-1} - D_c(n-2)*u_{n-2} = 0
# Shift n -> n+2: A_c(n+2)*u_{n+3} - B_c(n+2)*u_{n+2} + C_c(n+1)*u_{n+1} - D_c(n)*u_n = 0
def a0(n): return -D_c(n)
def a1(n): return C_c(n+1)
def a2(n): return -B_c(n+2)
def a3(n): return A_c(n+2)

# Verify on q_n
N = 60
q = [Q(-215040420000), Q(-167282265043404, 905), Q(-964185327658080, 6071)]
for i in range(3, N):
    n = i - 1
    new_q = Q(B_c(n), A_c(n)) * q[-1] - Q(C_c(n-1), A_c(n-1)) * q[-2] + Q(D_c(n-2), A_c(n-2)) * q[-3]
    q.append(new_q)

p = [Q(-612218384750), Q(-9525021973931919, 18100), Q(-29561828382772029, 65380)]
for i in range(3, N):
    n = i - 1
    new_p = Q(B_c(n), A_c(n)) * p[-1] - Q(C_c(n-1), A_c(n-1)) * p[-2] + Q(D_c(n-2), A_c(n-2)) * p[-3]
    p.append(new_p)

# Verify standard form on q
print("=== Verifying standard form recurrence on q_n ===")
for n in range(min(5, N-3)):
    res = a0(n)*q[n] + a1(n)*q[n+1] + a2(n)*q[n+2] + a3(n)*q[n+3]
    print(f"  n={n}: {res}")

# Adjoint recurrence:
# L*[w]_n = a_0(n)*w_n + a_1(n-1)*w_{n-1} + a_2(n-2)*w_{n-2} + a_3(n-3)*w_{n-3} = 0
# Forward: w_n = -[a_1(n-1)*w_{n-1} + a_2(n-2)*w_{n-2} + a_3(n-3)*w_{n-3}] / a_0(n)
# Backward: w_{n-3} = -[a_0(n)*w_n + a_1(n-1)*w_{n-1} + a_2(n-2)*w_{n-2}] / a_3(n-3)

# The adjoint Poincaré roots are 1/rho_j:
# mu_0 = 1/rho_0 ≈ 1.164 (slow adjoint)
# |mu_1| = 1/|rho_1| ≈ 952 (fast adjoint)
# Forward iteration → dominates by fast modes (mu_1^n ~ 952^n)
# Backward iteration → converges to slow mode (mu_0^n ~ 1.164^n)

print("\n=== Computing adjoint slow solution w^{(0)} via backward iteration ===")
# Use mpmath for precision
N_back = 200  # Start backward from here

# Seed: w_N = 1, w_{N-1} = 0, w_{N-2} = 0 (arbitrary)
w_back = [mpf(0)] * (N_back + 1)
w_back[N_back] = mpf(1)
w_back[N_back - 1] = mpf(0)
w_back[N_back - 2] = mpf(0)

# Backward iteration: w_{n-3} = -[a_0(n)*w_n + a_1(n-1)*w_{n-1} + a_2(n-2)*w_{n-2}] / a_3(n-3)
for n in range(N_back, 2, -1):
    a3_val = mpf(int(a3(n-3)))
    if a3_val == 0:
        print(f"  WARNING: a_3(n-3) = 0 at n={n}")
        continue
    w_nm3 = -(mpf(int(a0(n))) * w_back[n]
              + mpf(int(a1(n-1))) * w_back[n-1]
              + mpf(int(a2(n-2))) * w_back[n-2]) / a3_val
    w_back[n-3] = w_nm3

# Normalize so w_0 = 1
if w_back[0] != 0:
    norm = w_back[0]
    for i in range(len(w_back)):
        w_back[i] /= norm

print(f"  w^(0)_0 = {w_back[0]}")
print(f"  w^(0)_1 = {w_back[1]}")
print(f"  w^(0)_2 = {w_back[2]}")
print(f"  w^(0)_3 = {w_back[3]}")

# Check growth rate: w^{(0)}_{n+1} / w^{(0)}_n should approach mu_0 ≈ 1.164
print("\nGrowth ratios w^{(0)}_{n+1} / w^{(0)}_n:")
for n in [5, 10, 20, 30, 40, 50]:
    if n+1 < len(w_back) and w_back[n] != 0:
        ratio = w_back[n+1] / w_back[n]
        print(f"  n={n}: {ratio}")

# Verify adjoint recurrence
print("\n=== Verifying adjoint recurrence on w^{(0)} ===")
for n in range(3, min(10, len(w_back))):
    res = (mpf(int(a0(n))) * w_back[n]
           + mpf(int(a1(n-1))) * w_back[n-1]
           + mpf(int(a2(n-2))) * w_back[n-2]
           + mpf(int(a3(n-3))) * w_back[n-3])
    print(f"  n={n}: residual = {res}")

# Now compute the bilinear concomitant J(w, u)
# For order-3 operator L = sum a_j E^j, the concomitant is:
# We derive it from the Green-Lagrange identity:
# sum_{m=M}^{N} [w_m*(Lu)_m - u_m*(L*w)_m] = J_{N+1} - J_M
#
# The standard formula for third-order L:
# J(w,u)_n = a_3(n-1)*det2(w,u,n-1,n+1) + a_3(n-2)*det2(w,u,n-2,n)
#          + a_2(n-1)*det2(w,u,n-1,n)
#
# where det2(w,u,i,j) = w_i*u_j - w_j*u_i
#
# Let me derive it explicitly. For L = a_0 E^0 + a_1 E^1 + a_2 E^2 + a_3 E^3:
#
# w_n*(Lu)_n = w_n*[a_0(n)*u_n + a_1(n)*u_{n+1} + a_2(n)*u_{n+2} + a_3(n)*u_{n+3}]
# u_n*(L*w)_n = u_n*[a_0(n)*w_n + a_1(n-1)*w_{n-1} + a_2(n-2)*w_{n-2} + a_3(n-3)*w_{n-3}]
#
# Difference: w_n*a_0*u_n - u_n*a_0*w_n = 0 (cancels)
# + a_1(n)*w_n*u_{n+1} - a_1(n-1)*u_n*w_{n-1}
# + a_2(n)*w_n*u_{n+2} - a_2(n-2)*u_n*w_{n-2}
# + a_3(n)*w_n*u_{n+3} - a_3(n-3)*u_n*w_{n-3}
#
# This should telescope. Let me define:
# T_n = a_1(n)*w_n*u_{n+1} contribution: sum telescopes as
#   sum_n a_1(n)*w_n*u_{n+1} - a_1(n-1)*u_n*w_{n-1}
# = sum_n [a_1(n)*w_n*u_{n+1} - a_1(n)*u_{n+1}*w_n] DOESN'T telescope simply.
#
# Actually, let's just verify constancy numerically.
# The standard result for order r:
# J(w,u)_n = sum_{j=1}^{r} sum_{k=0}^{j-1} a_j(n-r+j)*w_{n-r+k}*u_{n-r+2j-k-1}*(-1)^k
# Hmm, this formula is not standard. Let me just compute J by brute force.

# BRUTE FORCE: compute sum_{m=0}^{n-1} [w_m*(Lu)_m - u_m*(L*w)_m] = J_n - J_0
# Since Lu = 0 and L*w = 0, this sum = 0, so J is constant.
# I'll compute J_n directly as a specific bilinear form.

# For a SECOND order operator L = a_0 + a_1 E + a_2 E^2:
# J(w,u)_n = a_2(n-1)*(w_{n-1}*u_n - w_n*u_{n-1}) + ... nah let me just
# verify constancy empirically.

# Let me define candidate concomitants and test which is constant.
# For order 3, J should involve w_n, w_{n+1} and u_n, u_{n+1}, u_{n+2} or similar.

# Standard formula (from Elaydi's textbook, adapted):
# For L = sum_{j=0}^3 a_j(n) E^j, the bilinear concomitant is:
# J(w,u)_n = sum_{j=1}^3 sum_{i=0}^{j-1} (-1)^i * a_j(n+i-j) * w_{n+i-j} * u_{n+j-1-i}

# Let me expand this:
# j=1: i=0: (-1)^0 * a_1(n-1) * w_{n-1} * u_n = a_1(n-1) * w_{n-1} * u_n
# j=2: i=0: a_2(n-2) * w_{n-2} * u_{n+1}
#       i=1: -a_2(n-1) * w_{n-1} * u_n
# j=3: i=0: a_3(n-3) * w_{n-3} * u_{n+2}
#       i=1: -a_3(n-2) * w_{n-2} * u_{n+1}
#       i=2: a_3(n-1) * w_{n-1} * u_n

# So: J(w,u)_n = a_1(n-1)*w_{n-1}*u_n
#              + a_2(n-2)*w_{n-2}*u_{n+1} - a_2(n-1)*w_{n-1}*u_n
#              + a_3(n-3)*w_{n-3}*u_{n+2} - a_3(n-2)*w_{n-2}*u_{n+1} + a_3(n-1)*w_{n-1}*u_n

# Simplify:
# J = w_{n-1}*u_n * [a_1(n-1) - a_2(n-1) + a_3(n-1)]
#   + w_{n-2}*u_{n+1} * [a_2(n-2) - a_3(n-2)]
#   + w_{n-3}*u_{n+2} * a_3(n-3)

print("\n=== Testing bilinear concomitant constancy ===")

L = float(zeta(2) + zeta(3))

q_mp = [mpf(int(q[i].numerator)) / mpf(int(q[i].denominator)) for i in range(N)]
p_mp = [mpf(int(p[i].numerator)) / mpf(int(p[i].denominator)) for i in range(N)]
e_mp = [p_mp[i] - mpf(L) * q_mp[i] for i in range(N)]

def J_bracket(w, u, n):
    """Compute J(w,u)_n using the candidate formula."""
    t1 = w[n-1]*u[n] * (mpf(int(a1(n-1))) - mpf(int(a2(n-1))) + mpf(int(a3(n-1))))
    t2 = w[n-2]*u[n+1] * (mpf(int(a2(n-2))) - mpf(int(a3(n-2))))
    t3 = w[n-3]*u[n+2] * mpf(int(a3(n-3)))
    return t1 + t2 + t3

# Test constancy on q
print("J(w^{(0)}, q)_n for several n:")
J_q_vals = []
for n in range(3, min(40, N-2)):
    Jval = J_bracket(w_back, q_mp, n)
    J_q_vals.append(Jval)
    if n <= 10 or n % 10 == 0:
        print(f"  n={n}: {Jval}")

# Test constancy on p
print("\nJ(w^{(0)}, p)_n for several n:")
J_p_vals = []
for n in range(3, min(40, N-2)):
    Jval = J_bracket(w_back, p_mp, n)
    J_p_vals.append(Jval)
    if n <= 10 or n % 10 == 0:
        print(f"  n={n}: {Jval}")

# Test on e = p - L*q (should be zero if c_0 = 0)
print("\nJ(w^{(0)}, e)_n for several n:")
J_e_vals = []
for n in range(3, min(40, N-2)):
    Jval = J_bracket(w_back, e_mp, n)
    J_e_vals.append(Jval)
    if n <= 10 or n % 10 == 0:
        print(f"  n={n}: {Jval}")

# Check J(w^{(0)}, p) / J(w^{(0)}, q) = zeta(2)+zeta(3)?
if J_q_vals and J_q_vals[0] != 0:
    ratio = J_p_vals[0] / J_q_vals[0]
    print(f"\nJ(w^(0), p) / J(w^(0), q) = {ratio}")
    print(f"zeta(2) + zeta(3) = {mpf(L)}")
    print(f"Difference: {ratio - mpf(L)}")

# Also: check if the formula is correct by checking constancy
if len(J_q_vals) >= 5:
    print(f"\nConstancy check J(w,q): max variation = {max(abs(J_q_vals[i] - J_q_vals[0]) for i in range(min(10, len(J_q_vals))))}")
    print(f"Constancy check J(w,p): max variation = {max(abs(J_p_vals[i] - J_p_vals[0]) for i in range(min(10, len(J_p_vals))))}")
