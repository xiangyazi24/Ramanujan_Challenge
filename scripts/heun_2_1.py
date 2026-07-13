#!/usr/bin/env python3
"""Problem 2.1: Derive and analyze the order-3 Heun carrier equation.
The gamma gauge g_n splits b_n into two cubic factors.
Goal: explicit ODE, Frobenius bases, connection coefficients → π."""
from sympy import *

n, x, theta = symbols('n x theta')

# The recurrence: u_n = a(n)*u_{n-1} + b(n)*u_{n-2}
# a(n) = -220n^3 - 484n^2 - 301n - 42
# b(n) = 4n^2(2n+1)^2(5n-4)(5n+6) = 100*n^2*(n+1/2)^2*(n-4/5)*(n+6/5)

# The cubic gamma gauge (from ChatGPT Q4663):
# g_n = n * (n+1/2) * (n-4/5) (the "source" cubic from b_n)
# OR g_n = n * (n+1/2) * (n+6/5) (the "target" cubic)
# The natural split: b(n) = g_n * h_n where g_n and h_n are degree 3 each

# From Q4663: the natural split uses g_n = n(2n+1)(5n-4)/10 roughly
# Let me derive the ODE directly.

# With the gauge rho_n = product_{j=1}^{n} something, define F(x) = sum rho_n^{-1} Q_n x^n

# Actually, let me just compute the order-3 ODE operator L3 symbolically.
# From the recurrence Q_n = a(n)Q_{n-1} + b(n)Q_{n-2}, with n >= 1:
# n^3(n-1)^3 c_n - a(n)(n-1)^3 c_{n-1} - b(n) c_{n-2} = 0 (after (n!)^3 gauge)
# This is the ORDER-6 equation.

# For the order-3 carrier, use the cubic gamma gauge:
# Set rho_n = (2n+1)(5n-4)n or similar. Then Q_n = rho_0*rho_1*...*rho_n * d_n
# and d_n satisfies a recurrence where the coefficient degrees are all 3.

# From Q4663: the natural gauge is
# rho_n = n * (2n+1) * (5n-4) / (some normalization)
# giving d_{n+1} = A(n)*d_n + B(n)*d_{n-1} with deg(A)=deg(B)=3

# Let me compute this gauge explicitly
a_n = -220*n**3 - 484*n**2 - 301*n - 42
b_n = 4*n**2*(2*n+1)**2*(5*n-4)*(5*n+6)

# b_n = [n*(2n+1)*(5n-4)] * [4*n*(2n+1)*(5n+6)]
# = [n*(2n+1)*(5n-4)] * [4*n*(2n+1)*(5n+6)]
# Hmm, this doesn't split into two clean factors of degree 3 each.
# b_n has factors: n^2, (2n+1)^2, (5n-4), (5n+6) = degree 2+2+1+1 = 6

# The split: g_n = n*(2n+1)*(5n+6), h_n = 4*n*(2n+1)*(5n-4)
# Then b_n = g_n * h_n / (something)... no, g_n*h_n = 4*n^2*(2n+1)^2*(5n+6)*(5n-4) = b_n ✓

# So split: source = n*(2n+1)*(5n+6), target = 4*n*(2n+1)*(5n-4)
# Or: source = n*(2n+1)*(5n-4), target = 4*n*(2n+1)*(5n+6)

# Let's use the Q4663 split: source is the denominator, target is the numerator
# of the gauged recurrence.

# Define rho_n = product_{j=1}^{n} [j*(2j+1)*(5j+6)]
# Then Q_n = rho_n * d_n
# Q_{n-1} = rho_{n-1} * d_{n-1}, and rho_n = n*(2n+1)*(5n+6) * rho_{n-1}
# So Q_n = n*(2n+1)*(5n+6) * rho_{n-1} * d_n

# Substituting into Q_n = a(n)*Q_{n-1} + b(n)*Q_{n-2}:
# n*(2n+1)*(5n+6)*rho_{n-1}*d_n = a(n)*rho_{n-1}*d_{n-1} + b(n)*rho_{n-2}*d_{n-2}

# rho_{n-1}/rho_{n-2} = (n-1)*(2n-1)*(5n+1)

# n*(2n+1)*(5n+6)*d_n = a(n)*d_{n-1} + b(n)/((n-1)*(2n-1)*(5n+1)) * d_{n-2}

# b(n) / [(n-1)*(2n-1)*(5n+1)] = 4*n^2*(2n+1)^2*(5n-4)*(5n+6) / [(n-1)*(2n-1)*(5n+1)]

# This doesn't simplify nicely. Let me try the OTHER split.
# source = n*(2n+1)*(5n-4), target = 4*n*(2n+1)*(5n+6)

# rho_n = product_{j=1}^{n} [j*(2j+1)*(5j-4)]
# rho_n/rho_{n-1} = n*(2n+1)*(5n-4)

# n*(2n+1)*(5n-4)*d_n = a(n)*d_{n-1} + [4*n*(2n+1)*(5n+6)] * [(n-1)*(2n-1)*(5n-9)]^{-1} * d_{n-2}
# Still messy.

# Let me try the Q4663 approach directly: use g_n with the "natural cubic split"
# rho_0 = 1, rho_n = n * (2n+1) * (5n+6) for n >= 1
# Then rho_1 = 1*3*11 = 33, rho_2 = 2*5*16 = 160, etc.

# Q_n / (product_{j=1}^{n} rho_j) = d_n
# Recurrence for d_n: c3(n)*d_n = c1(n)*d_{n-1} + c0(n)*d_{n-2}
# where c3(n) = rho_n = n*(2n+1)*(5n+6)
# c1(n) = a(n) (the original a coefficient)
# c0(n) = b(n) / product of two consecutive rho's... this is getting circular

# Let me just compute the coefficients NUMERICALLY for specific n and verify
print("Computing gauged recurrence numerically...")
from mpmath import mp, mpf, fac
mp.dps = 30

# Compute Q_n
Q = [mpf(0), mpf(1)]  # Q_{-1}=0, Q_0=1
for nn in range(1, 20):
    a_val = -220*nn**3 - 484*nn**2 - 301*nn - 42
    b_val = 4*nn**2*(2*nn+1)**2*(5*nn-4)*(5*nn+6)
    Q.append(a_val*Q[-1] + b_val*Q[-2])

# Compute d_n = Q_n / product_{j=1}^{n} [j*(2j+1)*(5j+6)]
def rho(j): return j*(2*j+1)*(5*j+6)

d = [mpf(0), mpf(1)]  # d_{-1}=0, d_0=1 (Q_0=1, product is empty)
for nn in range(1, 18):
    prod = mpf(1)
    for j in range(1, nn+1):
        prod *= rho(j)
    d.append(Q[nn+1] / prod)

print("d_n values:")
for nn in range(10):
    print(f"  d_{nn} = {mp.nstr(d[nn+1], 15)}")

# Check if d_n satisfies a degree-3 recurrence
print("\nGauged recurrence check:")
for nn in range(2, 8):
    # d_n should satisfy: c3(n)*d_n + c1(n)*d_{n-1} + c0(n)*d_{n-2} = 0
    # with all c_j degree 3
    # Let's compute the ratio d_n / d_{n-1} and see the structure
    if d[nn] != 0:
        ratio = d[nn+1] / d[nn]
        print(f"  d_{nn}/d_{nn-1} = {mp.nstr(ratio, 12)}")
