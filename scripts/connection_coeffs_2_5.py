#!/usr/bin/env python3
"""
Problem 2.5: Find constant connection coefficients a, b, c such that
    t(N) = a · D_N² + b · D_N · E_N + c · E_N²
where t(N) = q(N) / gauge(N) is the gauged CMF scalar.

Try multiple gauge candidates:
(A) g(N) = (-16)^N · (N!)^7
(B) g(N) = (-16)^N · Γ-product from factored recurrence
(C) Empirically determined from ratios q(N)/(a·D_N²+b·D_N·E_N+c·E_N²)
"""
from fractions import Fraction as F
from math import factorial, gcd
from functools import reduce

def M_int(n):
    m11 = (-2*n-5)*(n+3)**2 * (136*n**4 + 1424*n**3 + 5548*n**2 + 9551*n + 6141)
    m12 = 384*n**6 + 6384*n**5 + 44168*n**4 + 162698*n**3 + 336377*n**2 + 369933*n + 169011
    m13 = -480*n**4 - 4980*n**3 - 19210*n**2 - 32690*n - 20730
    m21 = (n+2)**2*(n+3)**2*(4*n+10)*(48*n**3 + 386*n**2 + 1017*n + 879)
    m22 = (n+2)**2*(-272*n**5 - 3848*n**4 - 21732*n**3 - 61184*n**2 - 85761*n - 47808)
    m23 = (n+2)**2*(320*n**3 + 2540*n**2 + 6610*n + 5640)
    m31 = (-4*n-10)*(n+2)**2*(n+3)**2*(32*n**4 + 302*n**3 + 1037*n**2 + 1530*n + 813)
    m32 = (n+2)**2*(192*n**6 + 2984*n**5 + 19116*n**4 + 64452*n**3 + 120256*n**2 + 117279*n + 46476)
    m33 = (n+2)**2*(-16*n**5 - 408*n**4 - 2912*n**3 - 8884*n**2 - 12254*n - 6240)
    return [[m11, m12, m13], [m21, m22, m23], [m31, m32, m33]]

def mat_mul(A, B):
    n = len(A)
    return [[sum(A[i][k]*B[k][j] for k in range(n)) for j in range(n)] for i in range(n)]

# Compute CMF scalar sequence q(N) = product[0][0]
N_MAX = 40
I3 = [[1,0,0],[0,1,0],[0,0,1]]
prod = [row[:] for row in I3]
q = [1]
for N in range(N_MAX):
    prod = mat_mul(prod, M_int(N))
    q.append(prod[0][0])

# Compute Delannoy sequences (over Q)
D = [F(1), F(3)]
E = [F(0), F(1)]
for n in range(1, N_MAX + 5):
    D.append((F(3*(2*n+1)) * D[n] - F(n) * D[n-1]) / F(n+1))
    E.append((F(3*(2*n+1)) * E[n] - F(n) * E[n-1]) / F(n+1))

# Sym² solutions
v1 = [D[n]**2 for n in range(N_MAX + 5)]
v2 = [D[n]*E[n] for n in range(N_MAX + 5)]
v3 = [E[n]**2 for n in range(N_MAX + 5)]

print("=== Approach C: Solve for gauge empirically ===")
print("If q(N) = gauge(N) * (a·D_N² + b·D_N·E_N + c·E_N²),")
print("then gauge(N) = q(N) / (a·D_N² + b·D_N·E_N + c·E_N²)")
print()

# Use N=0,1,2 to determine a,b,c from gauge_candidate
# But first, try WITHOUT gauge: q(N) = a·D² + b·DE + c·E²
# This gives a=1 from N=0, then solve for b,c from N=1,2

# System:
# N=0: 1 = a*1 + b*0 + c*0  =>  a = 1
# N=1: q(1) = 9*1 + 3*b + 1*c  =>  3b + c = q(1) - 9
# N=2: q(2) = 169*1 + (117/2)*b + (81/4)*c  =>  117b/2 + 81c/4 = q(2) - 169

q_F = [F(x) for x in q]

a = F(1)
# From N=1: 3b + c = q(1) - 9
# From N=2: 117b/2 + 81c/4 = q(2) - 169
A11, A12, B1 = F(3), F(1), q_F[1] - 9
A21, A22, B2 = F(117, 2), F(81, 4), q_F[2] - 169

det_A = A11*A22 - A12*A21
b_no_gauge = (B1*A22 - B2*A12) / det_A
c_no_gauge = (A11*B2 - A21*B1) / det_A

print(f"Without gauge: a=1, b={b_no_gauge}, c={c_no_gauge}")
print(f"  b numerator digits: {len(str(abs(b_no_gauge.numerator)))}")
print(f"  c numerator digits: {len(str(abs(c_no_gauge.numerator)))}")

# Check at N=3
test3 = a*v1[3] + b_no_gauge*v2[3] + c_no_gauge*v3[3]
print(f"\n  Verify N=3: a*D²+b*DE+c*E² = {test3}")
print(f"  q(3) = {q[3]}")
print(f"  Match: {test3 == q_F[3]}")

# Now try with factorial gauge
print("\n\n=== Gauge (A): g(N) = (-16)^N · (N!)^7 ===")
t = [q_F[N] / (F((-16)**N) * F(factorial(N))**7) if N > 0 else F(1) for N in range(min(N_MAX+1, 25))]

print(f"t(0) = {t[0]}")
print(f"t(1) = {t[1]} = {float(t[1]):.6f}")
print(f"t(2) = {t[2]}")

# Solve for a, b, c with gauged sequence
a_g = t[0]  # = 1
B1_g = t[1] - 9
B2_g = t[2] - 169
b_g = (B1_g*A22 - B2_g*A12) / det_A
c_g = (A11*B2_g - A21*B1_g) / det_A

print(f"\na={a_g}, b={b_g}, c={c_g}")
print(f"b = {float(b_g):.10f}")
print(f"c = {float(c_g):.10f}")

# Check: is this constant across N?
print("\nVerification (should match t(N) for all N):")
for N in range(min(len(t), 10)):
    pred = a_g * v1[N] + b_g * v2[N] + c_g * v3[N]
    diff = pred - t[N]
    if diff == 0:
        print(f"  N={N}: EXACT MATCH")
    else:
        print(f"  N={N}: MISMATCH, diff = {float(diff):.6e}")

# Since factorial gauge likely doesn't work perfectly, try the RATIO approach
print("\n\n=== Approach: Solve gauge(N) from connection ===")
print("Assume q(N) = gauge(N) · (D_N² + β·D_N·E_N + γ·E_N²)")
print("Then gauge(N) = q(N) / (D_N² + β·D_N·E_N + γ·E_N²)")
print("The ratio gauge(N+1)/gauge(N) should be a rational function of N.\n")

# First, determine β, γ from the ASYMPTOTIC behavior.
# For large N: D_N ~ C₊·λ₊^N + C₋·λ₋^N with λ₊ = 3+2√2, λ₋ = 3-2√2
# E_N ~ C₊'·λ₊^N + C₋'·λ₋^N (with different constants)
# D_N² ~ C₊²·λ₊^{2N} + 2C₊C₋ + C₋²·λ₋^{2N}
# The DOMINANT term of v₁ is C₊²·λ₊^{2N}

# If q(N) ~ K·λ₊^{2N}·N^σ·gauge(N), the dominant behavior fixes the ratio.

# Actually, let me try a direct approach. Pick THREE values of N (say 5,6,7)
# and solve for a, b, c AND gauge(5), gauge(6), gauge(7).
# This is underdetermined. Better: fix a=1, then for each N compute:
# gauge(N) = q(N) / (D_N² + b·D_N·E_N + c·E_N²)
# and require gauge(N+1)/gauge(N) = R(N) for a rational function R.

# The RATIO of consecutive gauges:
# gauge(N+1)/gauge(N) = [q(N+1)/q(N)] / [(D_{N+1}² + b·D_{N+1}·E_{N+1} + c·E_{N+1}²) / (D_N² + b·D_N·E_N + c·E_N²)]

# For the SPECIFIC case b=0, c=0 (just a=1):
print("Gauge ratios for q(N)/D_N²:")
for N in range(1, 15):
    g = q_F[N] / v1[N]
    if N >= 2:
        g_prev = q_F[N-1] / v1[N-1]
        ratio = g / g_prev
        print(f"  N={N}: gauge = {float(g):.6e}, ratio = {float(ratio):.6f}")

# Those ratios should stabilize to show what R(N) looks like
print("\nFitting ratio = α·N^7 (expected from degree analysis):")
for N in range(3, 15):
    g = q_F[N] / v1[N]
    g_prev = q_F[N-1] / v1[N-1]
    ratio = g / g_prev
    normalized = ratio / F(N)**7
    print(f"  N={N}: ratio/N^7 = {float(normalized):.6f}")

# The asymptotic ratio should be (-16) * product_of_shifts / product_of_shifts
# Let's see what (-16)·N^7 gives us
print("\nRatio / [(-16)·N^7]:")
for N in range(3, 15):
    g = q_F[N] / v1[N]
    g_prev = q_F[N-1] / v1[N-1]
    ratio = g / g_prev
    normalized = ratio / (F(-16) * F(N)**7)
    print(f"  N={N}: {float(normalized):.10f}")

print("\nDone.")
