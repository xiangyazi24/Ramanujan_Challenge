"""
Bilinear concomitant J(w^{(0)}, u) for the P2.7 monic order-3 operator.

For L = S^3 + a_2(n)S^2 + a_1(n)S + a_0(n), the bilinear concomitant is:

  J(w,u)_n = w_{n-1}·u_{n+2}
           + [w_{n-2} + a_2(n-1)·w_{n-1}]·u_{n+1}
           + [w_{n-3} + a_2(n-2)·w_{n-2} + a_1(n-1)·w_{n-1}]·u_n

For Lu = 0 and L*w = 0, J(w,u)_n is CONSTANT in n.
The value J(w^{(0)}, u) = c_0(u), the coefficient of the dominant mode.
So c_0(e) = 0 iff J(w^{(0)}, e) = 0.
"""

import mpmath
mpmath.mp.dps = 200

from mpmath import mpf, mp, pi, zeta
from fractions import Fraction as Q

def A_c(n):
    return int(1024*(2*n+5)**4*(2*n+7)**3*(2*n+9)**3*(946*n*n+6407*n+10860))

def B_c(n):
    return int(128*(2*n+7)**3*(2*n+9)**3*(104060*n**6+1745370*n**5+12145238*n**4+44886481*n**3+92943995*n**2+102256019*n+46709052))

def C_c(n):
    return int(16*(n+3)**4*(2*n+9)**3*(3784*n**5+57792*n**4+351019*n**3+1059230*n**2+1587211*n+944620))

def D_c(n):
    return int((n+3)**4*(n+4)**6*(946*n*n+4515*n+5399))

# Monic operator coefficients: u_{n+3} + a2(n)u_{n+2} + a1(n)u_{n+1} + a0(n)u_n = 0
# From P2.7: u_{n+1} = (B/A)u_n - (C_{n-1}/A_{n-1})u_{n-1} + (D_{n-2}/A_{n-2})u_{n-2}
# Shift n -> n+2: u_{n+3} = (B(n+2)/A(n+2))u_{n+2} - (C(n+1)/A(n+1))u_{n+1} + (D(n)/A(n))u_n
# So: u_{n+3} - (B(n+2)/A(n+2))u_{n+2} + (C(n+1)/A(n+1))u_{n+1} - (D(n)/A(n))u_n = 0
# Hence: a2(n) = -B(n+2)/A(n+2), a1(n) = C(n+1)/A(n+1), a0(n) = -D(n)/A(n)

def a0(n):
    return mpf(-D_c(n)) / mpf(A_c(n))

def a1(n):
    return mpf(C_c(n+1)) / mpf(A_c(n+1))

def a2(n):
    return mpf(-B_c(n+2)) / mpf(A_c(n+2))

# Forward iteration for P2.7 solutions
N = 120

# q_n
q = [mpf(0)] * (N + 10)
q[0] = mpf(-215040420000)
q[1] = mpf(-167282265043404) / mpf(905)
q[2] = mpf(-964185327658080) / mpf(6071)
for i in range(3, N + 5):
    n = i - 1
    q[i] = mpf(B_c(n))/mpf(A_c(n)) * q[i-1] - mpf(C_c(n-1))/mpf(A_c(n-1)) * q[i-2] + mpf(D_c(n-2))/mpf(A_c(n-2)) * q[i-3]

# p_n
p = [mpf(0)] * (N + 10)
p[0] = mpf(-612218384750)
p[1] = mpf(-9525021973931919) / mpf(18100)
p[2] = mpf(-29561828382772029) / mpf(65380)
for i in range(3, N + 5):
    n = i - 1
    p[i] = mpf(B_c(n))/mpf(A_c(n)) * p[i-1] - mpf(C_c(n-1))/mpf(A_c(n-1)) * p[i-2] + mpf(D_c(n-2))/mpf(A_c(n-2)) * p[i-3]

# e_n = p_n - (zeta(2) + zeta(3)) * q_n
z23 = zeta(2) + zeta(3)
e = [p[n] - z23 * q[n] for n in range(N + 5)]

# Adjoint backward iteration for w^{(0)}
# L* w = 0: w_{n-3} + a2(n-2)w_{n-2} + a1(n-1)w_{n-1} + a0(n)w_n = 0
# => w_{n-3} = -a2(n-2)w_{n-2} - a1(n-1)w_{n-1} - a0(n)w_n
# Backward: given w_N, w_{N-1}, w_{N-2}, compute w_{N-3}, w_{N-4}, ...
# Replace n by n+3: w_n = -a2(n+1)w_{n+1} - a1(n+2)w_{n+2} - a0(n+3)w_{n+3}

M = 350  # start backward from here
w_raw = [mpf(0)] * (M + 5)
w_raw[M] = mpf(1)
w_raw[M-1] = mpf(0)
w_raw[M-2] = mpf(0)

for n in range(M-3, -4, -1):
    idx = n
    if idx < -3:
        break
    w_raw[idx] = (-a2(idx+1)*w_raw[idx+1] - a1(idx+2)*w_raw[idx+2] - a0(idx+3)*w_raw[idx+3])

# Normalize w so that w stays manageable
# w^{(0)} grows at rate ~1/rho_0 ≈ 1.164^n for large n
# Normalize at n=10
norm_val = w_raw[10]
w = [mpf(0)] * (M + 5)
for n in range(M + 1):
    w[n] = w_raw[n] / norm_val if norm_val != 0 else w_raw[n]

print("=== Adjoint solution w^{(0)} (normalized at n=10) ===")
for n in range(15):
    print(f"  w[{n}] = {w[n]}")

# Verify adjoint recurrence: w_{n-3} + a2(n-2)w_{n-2} + a1(n-1)w_{n-1} + a0(n)w_n = 0
print("\n=== Adjoint recurrence residuals ===")
for n in range(6, 20):
    res = w[n-3] + a2(n-2)*w[n-2] + a1(n-1)*w[n-1] + a0(n)*w[n]
    print(f"  n={n}: {res}")

# Bilinear concomitant:
# J(w,u)_n = w_{n-1}·u_{n+2}
#          + [w_{n-2} + a2(n-1)·w_{n-1}]·u_{n+1}
#          + [w_{n-3} + a2(n-2)·w_{n-2} + a1(n-1)·w_{n-1}]·u_n

def J(w, u, n):
    return (w[n-1] * u[n+2]
            + (w[n-2] + a2(n-1)*w[n-1]) * u[n+1]
            + (w[n-3] + a2(n-2)*w[n-2] + a1(n-1)*w[n-1]) * u[n])

# Test constancy on q_n (should be constant = c_0(q))
print("\n=== J(w^{(0)}, q) at various n (should be constant) ===")
J_q_vals = []
for n in range(5, 60):
    Jval = J(w, q, n)
    J_q_vals.append(Jval)
    if n <= 15 or n % 10 == 0:
        print(f"  n={n}: J = {Jval}")

# Test constancy on p_n
print("\n=== J(w^{(0)}, p) at various n (should be constant) ===")
J_p_vals = []
for n in range(5, 60):
    Jval = J(w, p, n)
    J_p_vals.append(Jval)
    if n <= 15 or n % 10 == 0:
        print(f"  n={n}: J = {Jval}")

# THE KEY TEST: J(w^{(0)}, e) at various n
print("\n=== J(w^{(0)}, e) at various n (should be ZERO if c_0(e) = 0) ===")
J_e_vals = []
for n in range(5, 60):
    Jval = J(w, e, n)
    J_e_vals.append(Jval)
    if n <= 15 or n % 10 == 0:
        print(f"  n={n}: J = {Jval}")

# Also check: J(w^{(0)}, e) = J(w^{(0)}, p) - (z2+z3)*J(w^{(0)}, q) directly
print("\n=== Consistency check: J(p) - (z2+z3)*J(q) vs J(e) ===")
for n in range(5, 20):
    direct = J(w, e, n)
    indirect = J(w, p, n) - z23 * J(w, q, n)
    diff = direct - indirect
    print(f"  n={n}: J(e)={direct}, J(p)-z23*J(q)={indirect}, diff={diff}")

# Print ratio J(p)/J(q) to see if it equals z2+z3
print("\n=== J(w^{(0)}, p) / J(w^{(0)}, q) vs zeta(2)+zeta(3) ===")
print(f"  zeta(2)+zeta(3) = {z23}")
for n in [10, 20, 30, 40, 50]:
    if n-5 < len(J_q_vals) and J_q_vals[n-5] != 0:
        ratio = J_p_vals[n-5] / J_q_vals[n-5]
        print(f"  n={n}: J(p)/J(q) = {ratio}")
        print(f"         diff from z23 = {ratio - z23}")
