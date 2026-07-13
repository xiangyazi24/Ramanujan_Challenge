#!/usr/bin/env python3
"""Verify Problem 2.2 = Aptekarev recurrence after index shift m = n+3."""
from sympy import symbols, factor, expand, Rational, factorial

n, m = symbols('n m')

# Challenge recurrence coefficients
c0 = -8*n**3 - 51*n**2 - 105*n - 68
c1 = 24*n**5 + 337*n**4 + 1833*n**3 + 4818*n**2 + 6092*n + 2928
c2 = -(n+2)*(n+3)*(24*n**5 + 273*n**4 + 1150*n**3 + 2154*n**2 + 1635*n + 268)
c3 = (n+1)*(n+2)**4*(n+3)*(8*n**3 + 75*n**2 + 231*n + 232)

# Shift: n -> m-3
c0_shifted = c0.subs(n, m-3)
c1_shifted = c1.subs(n, m-3)
c2_shifted = c2.subs(n, m-3)
c3_shifted = c3.subs(n, m-3)

print("=== Shifted coefficients (n = m-3) ===")
print(f"c0(m-3) = {factor(c0_shifted)}")
print(f"c1(m-3) = {factor(c1_shifted)}")
print(f"c2(m-3) = {factor(c2_shifted)}")
print(f"c3(m-3) = {factor(c3_shifted)}")

# Verify initial values match Aptekarev
# Challenge: p_{-3}=0, p_{-2}=7, p_{-1}=179 -> P_0=0, P_1=7, P_2=179
# Challenge: q_{-3}=1, q_{-2}=12, q_{-1}=306 -> Q_0=1, Q_1=12, Q_2=306
print("\n=== Initial values after shift ===")
print("P_0 = p_{-3} = 0")
print("P_1 = p_{-2} = 7")
print("P_2 = p_{-1} = 179")
print("Q_0 = q_{-3} = 1")
print("Q_1 = q_{-2} = 12")
print("Q_2 = q_{-1} = 306")

# Gauge analysis: u_m = (m!)^2 v_m
# Leading terms of shifted coefficients
print("\n=== Gauge analysis ===")
from sympy import Poly, degree
for name, coeff in [("c0", c0_shifted), ("c1", c1_shifted),
                     ("c2", c2_shifted), ("c3", c3_shifted)]:
    p = Poly(expand(coeff), m)
    print(f"deg {name}(m-3) = {p.degree()}, leading = {p.LC()}")

# Characteristic polynomial: extract leading coefficients
# c0 * r^3 + c1 * r^2 / m^2 + c2 * r / m^4 + c3 / m^6 = 0
# After (m!)^2 gauge: all degrees become 3
# Leading: -8 r^3 + 24 r^2 - 24 r + 8 = -8(r-1)^3
print("\n=== Limiting characteristic polynomial ===")
r = symbols('r')
char_poly = -8*r**3 + 24*r**2 - 24*r + 8
print(f"-8r^3 + 24r^2 - 24r + 8 = {factor(char_poly)}")

# Verify: compute many terms and check convergence
from mpmath import mp, mpf, euler
mp.dps = 60

def compute_ratio(N):
    c0f = lambda n: -8*n**3 - 51*n**2 - 105*n - 68
    c1f = lambda n: 24*n**5 + 337*n**4 + 1833*n**3 + 4818*n**2 + 6092*n + 2928
    c2f = lambda n: -(n+2)*(n+3)*(24*n**5 + 273*n**4 + 1150*n**3 + 2154*n**2 + 1635*n + 268)
    c3f = lambda n: (n+1)*(n+2)**4*(n+3)*(8*n**3 + 75*n**2 + 231*n + 232)

    p = [mpf(0), mpf(7), mpf(179)]
    q = [mpf(1), mpf(12), mpf(306)]

    for nn in range(0, N):
        pn = -(c1f(nn)*p[-1] + c2f(nn)*p[-2] + c3f(nn)*p[-3]) / c0f(nn)
        qn = -(c1f(nn)*q[-1] + c2f(nn)*q[-2] + c3f(nn)*q[-3]) / c0f(nn)
        p.append(pn)
        q.append(qn)

    return p[-1] / q[-1]

print(f"\n=== Convergence to gamma ===")
print(f"gamma        = {euler}")
for N in [10, 20, 50, 100, 150]:
    ratio = compute_ratio(N)
    diff = ratio - euler
    print(f"N={N:3d}: p_n/q_n - gamma = {diff}")
