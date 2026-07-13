#!/usr/bin/env sage
"""Problem 2.3: Factor the order-4 recurrence operator using ore_algebra.
This is the CRITICAL test: if it factors as LCLM of two order-2 operators,
then pi+e splits into pi-part + e-part and the problem collapses."""

from ore_algebra import *

# Set up the Ore algebra for shift operators
R.<n> = QQ['n']
A.<Sn> = OreAlgebra(R, 'Sn')

# The recurrence: c0(n) u_n + c1(n) u_{n-1} + c2(n) u_{n-2} + c3(n) u_{n-3} + c4(n) u_{n-4} = 0
# Rewrite as operator acting on u: sum c_j(n) S^{-j} u_n = 0
# Or equivalently, shift everything: sum c_j(n+4) S^{4-j} u_{n} = 0 for forward shift

# Define coefficients (as in the challenge)
c0 = -n^3 + 2*n^2 + 7*n + 3
c1_coeff = (n+2)*(2*n^4 + n^3 - 26*n^2 - 48*n - 19)
c2_coeff = (n+2)*(n^6 + 9*n^5 + 8*n^4 - 87*n^3 - 249*n^2 - 234*n - 68)
c3_coeff = (n+1)^2*(n+2)*(2*n^5 + 3*n^4 - 13*n^3 - 21*n^2 + 4)
c4_coeff = -n^3*(n+1)^2*(n+2)*(n^3 + n^2 - 8*n - 11)

# The recurrence is: c0(n)*u(n) + c1(n)*u(n-1) + c2(n)*u(n-2) + c3(n)*u(n-3) + c4(n)*u(n-4) = 0
# In forward-shift form with S = shift-by-1:
# Replace n with n+4 and multiply through:
# c0(n+4)*S^4 + c1(n+4)*S^3 + c2(n+4)*S^2 + c3(n+4)*S + c4(n+4) = 0

L = c0(n=n+4)*Sn^4 + c1_coeff(n=n+4)*Sn^3 + c2_coeff(n=n+4)*Sn^2 + c3_coeff(n=n+4)*Sn + c4_coeff(n=n+4)

print("=== The order-4 operator ===")
print(f"L = {L}")
print(f"Order: {L.order()}")

print("\n=== Attempting right factorization ===")
try:
    factors = L.factor()
    print(f"Factored: {factors}")
except Exception as e:
    print(f"Direct factorization failed: {e}")

print("\n=== Attempting to find right factors of order 1 ===")
try:
    rf1 = L.right_factors(order=1)
    print(f"Right factors of order 1: {rf1}")
except Exception as e:
    print(f"No order-1 right factors: {e}")

print("\n=== Attempting to find right factors of order 2 ===")
try:
    rf2 = L.right_factors(order=2)
    print(f"Right factors of order 2: {rf2}")
except Exception as e:
    print(f"No order-2 right factors: {e}")

print("\n=== Symmetric product / LCLM decomposition ===")
# If L = LCLM(L1, L2) where L1, L2 are order 2,
# then any solution of L1 or L2 is a solution of L.
# Equivalently, L divides L1*L2 from the left.
# We can try: compute the symmetric square and look for order-2 right divisors

# Alternative: try directly with the recurrence in backward-shift form
print("\n=== Alternative: backward shift form ===")
# u_n = -(c1/c0)*u_{n-1} - (c2/c0)*u_{n-2} - (c3/c0)*u_{n-3} - (c4/c0)*u_{n-4}
# Let S^{-1} be the backward shift. The operator is:
# c0*I + c1*S^{-1} + c2*S^{-2} + c3*S^{-3} + c4*S^{-4}

# In ore_algebra, we can also try the differential form
# or just compute generalized series solutions

print("\n=== Generalized series solutions (for identification) ===")
try:
    sols = L.generalized_series_solutions(n=10)
    print("Generalized series solutions:")
    for i, s in enumerate(sols):
        print(f"  sol[{i}]: {s}")
except Exception as e:
    print(f"Series solutions: {e}")

print("\n=== Hypergeometric solutions ===")
try:
    hyp = L.hypergeometric_solutions()
    print(f"Hypergeometric solutions: {hyp}")
except Exception as e:
    print(f"Hypergeometric: {e}")
