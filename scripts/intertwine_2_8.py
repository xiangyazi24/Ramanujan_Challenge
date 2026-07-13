#!/usr/bin/env python3
"""Problem 2.8: Find the gauge matrix G(n) intertwining M(n) with Chudnovsky T_k.
Goal: M(n) * G(n+1) = G(n) * T(n+1), where T is the canonical Chudnovsky transfer.

This is the key algebraic bridge making the 2.8 proof complete."""
from sympy import *

n, k = symbols('n k')

C = 640320
A_chud = 545140134
B_chud = 13591409
R = 1 + C**3 // 1728  # = 151931373056001

# Chudnovsky term ratio r_k = h_{k+1}/h_k
# r_k = -24*(2k+1)*(6k+1)*(6k+5) / (C^3 * (k+1)^3)
# With k = n+1: r_{n+1} = -24*(2n+3)*(6n+7)*(6n+11) / (C^3*(n+2)^3)

def r_chud(k_val):
    """Chudnovsky term ratio."""
    return -24*(2*k_val+1)*(6*k_val+1)*(6*k_val+5) / (C**3 * (k_val+1)**3)

# Canonical Chudnovsky transfer matrix T_k (4x4)
# State: y_k = (1, S_k, h_k, k*h_k)
# T_k maps y_k to y_{k+1}
def T_chud(k_val):
    rk = r_chud(k_val)
    return Matrix([
        [1, 0, 0, 0],
        [0, 1, 0, 0],
        [0, B_chud, rk, rk],
        [0, A_chud, 0, rk]
    ])

# Verify T_k is correct by checking: y_{k+1} = y_k * T_k
# y_k = (1, S_k, h_k, k*h_k)
# y_{k+1} should be (1, S_{k+1}, h_{k+1}, (k+1)*h_{k+1})
# S_{k+1} = S_k + (A*k + B)*h_k = S_k + A*(k*h_k) + B*h_k
# h_{k+1} = r_k * h_k
# (k+1)*h_{k+1} = (k+1)*r_k*h_k = r_k*(k*h_k) + r_k*h_k

print("=== Verifying canonical Chudnovsky transfer T_k ===")
h, S, kh = symbols('h S kh')  # h = h_k, S = S_k, kh = k*h_k
y = Matrix([[1, S, h, kh]])
rk_sym = symbols('rk')
T = Matrix([
    [1, 0, 0, 0],
    [0, 1, 0, 0],
    [0, B_chud, rk_sym, rk_sym],
    [0, A_chud, 0, rk_sym]
])
y_next = y * T
print(f"y * T = {y_next}")
print(f"  coord 0: {y_next[0,0]} (should be 1)")
print(f"  coord 1: {y_next[0,1]} (should be S + B*h + A*kh = S_{k+1})")
print(f"  coord 2: {y_next[0,2]} (should be rk*h + rk*kh = rk*(h+kh))")
print(f"  coord 3: {y_next[0,3]} (should be rk*h + rk*kh = rk*(h+kh))")
# Wait, coord 2 should be h_{k+1} = rk*h
# Let me recheck: T row 2 (0-indexed): [0, B, rk, rk]
# y * T column 2: 1*0 + S*0 + h*rk + kh*rk = rk*(h+kh)
# But h_{k+1} = rk * h_k, NOT rk*(h+kh). So the matrix might use right-multiplication
# or a different state convention.
print("\nNote: The state convention may be column vectors with LEFT multiplication.")
print("Let me try: T_k * y_k^T = y_{k+1}^T")

y_col = Matrix([1, S, h, kh])
y_next_col = T * y_col
print(f"\nT * y^T = {y_next_col}")
print(f"  coord 0: {y_next_col[0]} (should be 1)")
print(f"  coord 1: {y_next_col[1]} (should be S)")
print(f"  coord 2: {y_next_col[2]} (should be B*S + rk*h + rk*kh)")
print(f"  coord 3: {y_next_col[3]} (should be A*S + rk*kh)")

# Hmm, neither convention gives the right thing directly.
# Let me reconsider: maybe the state is (1, S_k, h_k, k*h_k) as a ROW,
# and the product is y_{k+1} = y_k * T_k where T_k acts on the right.
# Then:
# (y*T)[1] = S*1 + h*B + kh*A = S + B*h + A*kh = S + (B+Ak)*h = S_{k+1} ✓
# (y*T)[2] = h*rk + kh*rk = (h+kh)*rk = (1+k)*h*rk = (k+1)*h_{k+1}
# But we want h_{k+1} = rk*h, not (k+1)*h_{k+1}!
#
# So the state vector is NOT (1, S, h, k*h). It might be (1, S, h, (k+1)*h).
# Let me try: y = (1, S_k, h_k, (k+1)*h_k)
# Then (y*T)[2] = h*rk + (k+1)*h*rk = (k+2)*h*rk = (k+2)*h_{k+1} ✓ for next k+1 -> k+2? No...

# Actually, I think the issue is simpler. Let me re-derive T_k for the state
# z_k = (1, S_k, h_k, k*h_k) where S_{k+1} = S_k + (Ak+B)*h_k

# z_{k+1} = (1, S_k + (Ak+B)h_k, r_k*h_k, (k+1)*r_k*h_k)
# = (1, S_k + B*h_k + A*(k*h_k), r_k*h_k, r_k*(k*h_k) + r_k*h_k)

# So in matrix form z_{k+1} = z_k * T_k:
# T_k = [[1, 0, 0, 0],
#         [0, 1, 0, 0],
#         [0, B, r_k, r_k],
#         [0, A, 0, r_k]]

# Check: z_k * T_k:
# [0] = 1*1 + 0 + 0 + 0 = 1 ✓
# [1] = 0 + S_k*1 + h_k*B + (k*h_k)*A = S_k + B*h_k + A*k*h_k = S_{k+1} ✓
# [2] = 0 + 0 + h_k*r_k + 0 = r_k*h_k = h_{k+1} ✓
# [3] = 0 + 0 + h_k*r_k + (k*h_k)*r_k = r_k*h_k*(1+k) = (k+1)*h_{k+1} ✓

print("\n=== Verified: T_k is correct with RIGHT multiplication ===")
print("State: z_k = (1, S_k, h_k, k*h_k)")
print("z_{k+1} = z_k * T_k  ✓")

# Now: need M(n) * G(n+1) = G(n) * T_{n+1}
# where M(n) is the challenge matrix (4x4) and G(n) is the unknown gauge
#
# With k = n+1:
# T_{n+1} has r_{n+1} = -24*(2n+3)*(6n+7)*(6n+11) / (C^3*(n+2)^3)
#
# The computation: solve for G(n) entry by entry
# This is a 16-equation polynomial system over Q(n).
# Heavy symbolic algebra — dispatch to uisai2 Sage.

print("\n=== Dispatching gauge computation to uisai2 ===")
print("Need: solve M(n)*G(n+1) = G(n)*T(n+1) for G(n) in GL_4(Q(n,R))")
print("This is 16 polynomial equations in the 16 entries of G(n).")
