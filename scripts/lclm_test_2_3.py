#!/usr/bin/env python3
"""Test LCLM hypothesis for Problem 2.3 NUMERICALLY.
If the order-4 operator is LCLM of two order-2 operators,
then the 4-dimensional solution space splits into two 2-dimensional subspaces,
each annihilated by an order-2 operator.

Strategy: compute 4 independent solutions, then search for order-2 sub-recurrences."""
from mpmath import mp, mpf, matrix, lu_solve, pi, e as euler_e, euler
import numpy as np

mp.dps = 50

# Recurrence coefficients
def c0(n): return -n**3 + 2*n**2 + 7*n + 3
def c1(n): return (n+2)*(2*n**4 + n**3 - 26*n**2 - 48*n - 19)
def c2(n): return (n+2)*(n**6 + 9*n**5 + 8*n**4 - 87*n**3 - 249*n**2 - 234*n - 68)
def c3(n): return (n+1)**2*(n+2)*(2*n**5 + 3*n**4 - 13*n**3 - 21*n**2 + 4)
def c4(n): return -n**3*(n+1)**2*(n+2)*(n**3 + n**2 - 8*n - 11)

def compute_solution(init, N):
    """Compute u_0,...,u_N given 4 initial values u_{-3},...,u_0."""
    seq = list(init)
    for nn in range(1, N):
        un = -(c1(nn)*seq[-1] + c2(nn)*seq[-2] + c3(nn)*seq[-3] + c4(nn)*seq[-4]) / c0(nn)
        seq.append(un)
    return seq

# Four linearly independent solutions
# Solution 1: the p-sequence from the challenge
sol_p = compute_solution([mpf(1), mpf(1), mpf(20), mpf(296)], 60)
# Solution 2: the q-sequence from the challenge
sol_q = compute_solution([mpf(1), mpf(0), mpf(4), mpf(48)], 60)
# Solution 3: random initial conditions
sol_r = compute_solution([mpf(0), mpf(1), mpf(0), mpf(0)], 60)
# Solution 4: another random
sol_s = compute_solution([mpf(0), mpf(0), mpf(1), mpf(0)], 60)

print("=== LCLM Test for Problem 2.3 ===")
print(f"4 solutions computed to N=60, {mp.dps} digits")

# If LCLM(L1, L2) with L1 order-2 and L2 order-2:
# There exist 2 solutions annihilated by L1 and 2 by L2.
# An order-2 recurrence: alpha(n)*u_n + beta(n)*u_{n-1} + gamma(n)*u_{n-2} = 0
# For a given solution u, at each n: u_n/u_{n-2} and u_{n-1}/u_{n-2} satisfy a linear relation.

# Test: for each solution, check if it satisfies an order-2 recurrence
# by fitting alpha, beta, gamma polynomially.

print("\n=== Testing if individual solutions satisfy order-2 recurrences ===")

def test_order2(sol, name, N_test=40):
    """Test if sol satisfies an order-2 recurrence.
    At each n, we need: a(n)*sol[n] + b(n)*sol[n-1] + c(n)*sol[n-2] = 0
    where a,b,c are polynomials in n.

    If this holds, then sol[n]/sol[n-2] = -(b(n)*sol[n-1]/sol[n-2] + c(n))/a(n)
    which means sol[n]/sol[n-2] and sol[n-1]/sol[n-2] are linearly related.

    We test by checking the 3x3 Casorati determinant:
    det | u_n   u_{n-1}   u_{n-2} |
        | u_{n+1} u_n     u_{n-1} |
        | u_{n+2} u_{n+1} u_n     |
    If this vanishes for all n, the solution satisfies an order-2 recurrence."""

    print(f"\n  {name}:")
    for n in [10, 20, 30]:
        i = n + 3  # offset for indexing (initial values at -3,-2,-1,0)
        det = (sol[i] * (sol[i+1]*sol[i+4] - sol[i+2]*sol[i+3])
             - sol[i+1] * (sol[i]*sol[i+4] - sol[i+2]**2)  # Wait, this isn't right
             + sol[i+2] * (sol[i]*sol[i+3] - sol[i+1]**2))
        # Actually the Casorati (Wronskian) determinant for order-2 test:
        # det | u(n)   u(n+1)  u(n+2) |
        #     | v(n)   v(n+1)  v(n+2) |  for TWO solutions u,v
        # If a single solution u satisfies order-2, then
        # u(n)*u(n+2) - u(n+1)^2 = ... (this is a Casorati condition)
        # Not quite right. Let me think differently.
        pass

# Better approach: test pairwise if two solutions span an order-2 subspace
print("\n=== Testing PAIRS of solutions for order-2 subspaces ===")

def casorati_2(u, v, n):
    """2x2 Casorati determinant: u(n)*v(n+1) - u(n+1)*v(n)"""
    return u[n]*v[n+1] - u[n+1]*v[n]

def casorati_3(u, v, w, n):
    """3x3 Casorati determinant"""
    return (u[n] * (v[n+1]*w[n+2] - v[n+2]*w[n+1])
          - u[n+1] * (v[n]*w[n+2] - v[n+2]*w[n])
          + u[n+2] * (v[n]*w[n+1] - v[n+1]*w[n]))

# If two solutions u, v span an order-2 invariant subspace,
# then for any third solution w NOT in that subspace,
# the 3x3 Casorati det(u,v,w) should NOT vanish.
# But the 3x3 Casorati det(u,v,u) = 0 trivially.
#
# Better test: an order-2 subspace means there exists a(n), b(n), c(n)
# such that a*u_n + b*u_{n-1} + c*u_{n-2} = 0 for BOTH u and v.
# Equivalently, the 3x2 matrix [u_n u_{n-1} u_{n-2}; v_n v_{n-1} v_{n-2}]
# has a kernel.

# Actually, the correct test is:
# 3 solutions of order-2 are dependent → 3x3 Casorati = 0
# 2 solutions of order-2 → 2x2 Casorati non-zero (they're independent within order-2)
# but a THIRD solution NOT in the subspace would give non-zero 3x3 Casorati

# The REAL test:
# For order-2 subspace spanned by {u, v}, the Casorati det for the
# ORDER-2 recurrence should be a specific rational function of n.
# We can just try all (4 choose 2) = 6 pairs and check if any pair
# generates an order-2 subspace by testing if their 3x3 Casorati with
# any other solution is "linearly dependent" (i.e., the third solution
# is a linear combination of the first two plus an order-2 independent solution).

# Simplest approach: for each pair (u,v), compute the rational function
# r(n) = (u_n*v_{n+2} - u_{n+2}*v_n) / (u_n*v_{n+1} - u_{n+1}*v_n)
# If u,v satisfy a common order-2 recurrence alpha*S^2 + beta*S + gamma = 0,
# then r(n) = -beta(n)/alpha(n) should be a rational function of n.

sols = {"p": sol_p, "q": sol_q, "r": sol_r, "s": sol_s}
pairs = [("p","q"), ("p","r"), ("p","s"), ("q","r"), ("q","s"), ("r","s")]

for name1, name2 in pairs:
    u = sols[name1]
    v = sols[name2]
    print(f"\n  Pair ({name1}, {name2}):")
    ratios = []
    for n in range(5, 45):
        i = n + 3
        cas2 = u[i]*v[i+1] - u[i+1]*v[i]
        if cas2 != 0:
            cas2_shifted = u[i+1]*v[i+2] - u[i+2]*v[i+1]
            r = cas2_shifted / cas2
            ratios.append((n, r))

    # Check if the ratio r(n) is a simple rational function of n
    # by looking at r(n)*n, r(n)/n, r(n)*n^2, etc.
    if ratios:
        print(f"    r(n) = C_{n+1}/C_n at selected n:")
        for n, r in ratios[:3]:
            print(f"      n={n}: r = {r}")
        for n, r in ratios[-3:]:
            print(f"      n={n}: r = {r}")

        # Check if r(n)/n^2 stabilizes (Poincaré root test)
        print(f"    r(n)/n^2 at large n:")
        for n, r in ratios[-5:]:
            print(f"      n={n}: r/n^2 = {r/n**2}")

# Also check the actual limit p_n/q_n
print(f"\n=== Limit check ===")
target = pi + euler_e
print(f"pi + e = {target}")
idx = len(sol_p) - 1
ratio = sol_p[idx] / sol_q[idx]
print(f"p_n/q_n (n={idx-3}) = {ratio}")
print(f"diff = {ratio - target}")
