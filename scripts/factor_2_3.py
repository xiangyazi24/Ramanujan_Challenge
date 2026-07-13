#!/usr/bin/env python3
"""Problem 2.3: Factor the order-4 recurrence operator.
Test Fable's prediction: LCLM of two order-2 operators (π-part + e-part)."""
from sympy import *

n = symbols('n')

# Recurrence coefficients
c0 = -n**3 + 2*n**2 + 7*n + 3
c1 = (n+2)*(2*n**4 + n**3 - 26*n**2 - 48*n - 19)
c2 = (n+2)*(n**6 + 9*n**5 + 8*n**4 - 87*n**3 - 249*n**2 - 234*n - 68)
c3 = (n+1)**2*(n+2)*(2*n**5 + 3*n**4 - 13*n**3 - 21*n**2 + 4)
c4 = -n**3*(n+1)**2*(n+2)*(n**3 + n**2 - 8*n - 11)

print("=== Factor all coefficient polynomials ===")
print(f"c0 = {factor(c0)}")
print(f"c1 = {factor(c1)}")
print(f"c2 = {factor(c2)}")
print(f"c3 = {factor(c3)}")
print(f"c4 = {factor(c4)}")

# Check leading polynomial factorizations
print("\n=== Roots of irreducible factors ===")
# c0 = -n^3 + 2n^2 + 7n + 3
print(f"Roots of c0: {solve(c0, n)}")
# cubic in c4: n^3 + n^2 - 8n - 11
cubic_c4 = n**3 + n**2 - 8*n - 11
print(f"Roots of cubic in c4: {solve(cubic_c4, n)}")

# Gauge analysis: u_n = (n!)^k v_n
print("\n=== Degree analysis for gauge ===")
for name, coeff in [("c0", c0), ("c1", c1), ("c2", c2), ("c3", c3), ("c4", c4)]:
    d = degree(Poly(expand(coeff), n))
    print(f"deg {name} = {d}")

# Degrees: 3, 5, 7, 7, 9 -> gauge k: 3, 5-k, 7-2k, 7-3k, 9-4k
# Equal when 3 = 5-k -> k=2. Check: 3, 3, 3, 1, 1. Not balanced!
# Try: after factoring out common n+2 from c1,c2,c3,c4
# Actually let me compute the leading coefficients for Poincaré analysis
print("\n=== Poincaré analysis ===")
# Leading term of each coefficient as n -> infinity
# c0 ~ -n^3
# c1 ~ (n)(2n^4) = 2n^5
# c2 ~ (n)(n^6) = n^7
# c3 ~ n^2 * n * 2n^5 = 2n^8
# c4 ~ -n^3 * n^2 * n * n^3 = -n^9
# After gauge u_n = (n!)^k:
# effective degrees: 3, 5-k, 7-2k, 8-3k, 9-4k
# Can't make all equal! But for order-4 Poincaré with factorial gauge:
# The recurrence is sum_{j=0}^4 c_j(n) u_{n-j} = 0
# With gauge u_n = (n!)^k v_n, u_{n-j} = (n!)^k / (n(n-1)...(n-j+1))^k * v_{n-j}
# So effective: c_j(n) / (n^{jk} * leading) for large n

# Let me just compute the Poincaré characteristic polynomial directly
# Recurrence: c0 u_n + c1 u_{n-1} + c2 u_{n-2} + c3 u_{n-3} + c4 u_{n-4} = 0
# Leading coefficient of c_j is the coefficient of n^{deg(c_j)} in c_j
# c0 leading: -1 (n^3)
# c1: expand (n+2)(2n^4+...) -> leading 2n^5
# c2: expand (n+2)(n^6+...) -> leading n^7
# c3: expand (n+1)^2(n+2)(2n^5+...) -> leading 2n^8 (check: n^2 * n * 2n^5 = 2n^8)

# Actually for Poincaré we need the gauge. Let me try k=2 first.
# With u_n = (n!)^2 v_n: divisions by n^{2j} for j-step back
# effective leading: c0 ~ -n^3, c1/n^2 ~ 2n^3, c2/n^4 ~ n^3, c3/n^6 ~ 2n^2, c4/n^8 ~ -n
# Not balanced. Try k=1:
# c0 ~ -n^3, c1/n ~ 2n^4, c2/n^2 ~ n^5, c3/n^3 ~ 2n^5, c4/n^4 ~ -n^5
# Still not balanced. The issue is the degrees aren't arithmetic.

# Let me re-examine exact degrees
print("\nExact degree of each c_j:")
for j, (name, coeff) in enumerate(zip(["c0","c1","c2","c3","c4"],
                                        [c0, c1, c2, c3, c4])):
    p = Poly(expand(coeff), n)
    print(f"  c_{j}: degree {p.degree()}, leading coeff {p.LC()}")

# c0: 3, c1: 5, c2: 7, c3: 8, c4: 9
# For Poincaré with gauge (n!)^k:
# effective leading: n^{3}, n^{5-k}, n^{7-2k}, n^{8-3k}, n^{9-4k}
# To have at most two dominant: we need the TOP TWO to match
# With k=2: 3, 3, 3, 2, 1 -> top 3 match at degree 3
# Characteristic poly from the top-3 terms: -r^4 + 2r^3 + r^2 = 0
# Wait, let me be more careful.

# Standard Poincaré: c0(n) u_n + c1(n) u_{n-1} + ... + c4(n) u_{n-4} = 0
# With gauge u_n = Gamma(n+1)^k * v_n:
# c_j(n) * Gamma(n-j+1)^k / Gamma(n+1)^k * v_{n-j}
# = c_j(n) / (n(n-1)...(n-j+1))^k * v_{n-j}
# ~ c_j(n) / n^{jk} * v_{n-j}  for large n

# With k=2:
# j=0: c_0 / n^0 ~ -n^3
# j=1: c_1 / n^2 ~ 2n^3
# j=2: c_2 / n^4 ~ n^3
# j=3: c_3 / n^6 ~ 2n^2
# j=4: c_4 / n^8 ~ -n^1

# The top 3 terms dominate: -r^0 + 2r^1 + r^2 = 0 ... no wait
# Standard: sum alpha_j * r^{4-j} = 0 where alpha_j = leading of c_j/n^{jk}
# With k=2: alpha_0=-1, alpha_1=2, alpha_2=1, alpha_3~0, alpha_4~0
# Char poly: -r^4 + 2r^3 + r^2 = r^2(-r^2 + 2r + 1) = 0

# Hmm, that gives r^2 * (-(r-1)^2 + 2) = 0 -> r^2(2-(r-1)^2)=0
# r=0 (double), r = 1±√2

print("\n=== Poincaré characteristic polynomial (gauge k=2) ===")
r = symbols('r')
char = -r**4 + 2*r**3 + r**2
print(f"Char poly: {char} = {factor(char)}")
print(f"Roots: {solve(char, r)}")

# With k=1:
# alpha_0=-n^3, alpha_1=2n^4, alpha_2=n^5, alpha_3=2n^5, alpha_4=-n^5
# Not same degree. The max degree is 5 (from c2,c3,c4)
# Char poly from degree-5 terms: 0 + 0 + r^2 + 2r - 1 = 0 ???
# This is getting messy. Let me just compute numerically.

print("\n=== Numerical Poincaré root estimation ===")
from mpmath import mp, mpf
mp.dps = 30

# Compute many terms and look at ratios u_n / u_{n-1}
def compute_2_3_terms(N, use_p=True):
    def c0f(n): return -n**3 + 2*n**2 + 7*n + 3
    def c1f(n): return (n+2)*(2*n**4 + n**3 - 26*n**2 - 48*n - 19)
    def c2f(n): return (n+2)*(n**6 + 9*n**5 + 8*n**4 - 87*n**3 - 249*n**2 - 234*n - 68)
    def c3f(n): return (n+1)**2*(n+2)*(2*n**5 + 3*n**4 - 13*n**3 - 21*n**2 + 4)
    def c4f(n): return -n**3*(n+1)**2*(n+2)*(n**3 + n**2 - 8*n - 11)

    if use_p:
        seq = [mpf(1), mpf(1), mpf(20), mpf(296)]
    else:
        seq = [mpf(1), mpf(0), mpf(4), mpf(48)]

    for nn in range(1, N):
        un = -(c1f(nn)*seq[-1] + c2f(nn)*seq[-2] + c3f(nn)*seq[-3] + c4f(nn)*seq[-4]) / c0f(nn)
        seq.append(un)

    return seq

p_seq = compute_2_3_terms(30, True)
q_seq = compute_2_3_terms(30, False)

print("Ratios p_n / p_{n-1}:")
for i in range(5, 15):
    if p_seq[i-1] != 0:
        print(f"  n={i}: {p_seq[i]/p_seq[i-1]}")

print("\nRatios q_n / q_{n-1}:")
for i in range(5, 15):
    if q_seq[i-1] != 0:
        print(f"  n={i}: {q_seq[i]/q_seq[i-1]}")

# Check if numerator and denominator terms suggest the same Poincaré root
print("\nRatios p_n/p_{n-1} at large n:")
for i in [20, 25, 29]:
    if p_seq[i-1] != 0:
        ratio = p_seq[i] / p_seq[i-1]
        # Divide by n^2 to see the base Poincaré root
        print(f"  n={i}: raw={ratio}, /n^2={ratio/i**2}")
