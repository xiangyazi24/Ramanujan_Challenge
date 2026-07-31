#!/usr/bin/env python3
"""Problem 2.1: Polynomial Continued Fraction for π.

a_n = -220n³ - 484n² - 301n - 42
b_n = 4n²(2n+1)²(5n-4)(5n+6)

Prove: a_0 + K(b_n/a_n) = 6/(3-π)

Compute convergents and analyze the recurrence.
"""
from mpmath import mp, mpf, nstr, pi
mp.dps = 100

def a(n):
    return -220*n**3 - 484*n**2 - 301*n - 42

def b(n):
    return 4*n**2 * (2*n+1)**2 * (5*n-4) * (5*n+6)

# Convergents via the 3-term recurrence
# p_n = a_n * p_{n-1} + b_n * p_{n-2}
# q_n = a_n * q_{n-1} + b_n * q_{n-2}
target = mpf(6) / (3 - pi)
print(f"Target: 6/(3-π) = {nstr(target, 50)}")
print(f"a_0 = {a(0)}")
print()

p_prev, p_curr = mpf(1), mpf(a(0))  # p_{-1}=1, p_0=a_0
q_prev, q_curr = mpf(0), mpf(1)      # q_{-1}=0, q_0=1

print(f"{'n':>4} {'p_n/q_n':>60} {'error':>20}")
for n in range(1, 61):
    an = mpf(a(n))
    bn = mpf(b(n))
    p_new = an * p_curr + bn * p_prev
    q_new = an * q_curr + bn * q_prev
    p_prev, p_curr = p_curr, p_new
    q_prev, q_curr = q_curr, q_new
    
    conv = p_curr / q_curr
    err = abs(conv - target)
    if n <= 10 or n % 10 == 0:
        log_err = float(mp.log10(err)) if err > 0 else -999
        print(f"{n:>4} {nstr(conv, 50):>60} {log_err:>10.2f}")

# Factor structure analysis
print("\n=== Factor analysis ===")
print("\na_n factored:")
for n in range(6):
    an = a(n)
    print(f"  a({n}) = {an}")

print("\nb_n factored:")
for n in range(1, 6):
    bn = b(n)
    print(f"  b({n}) = {bn} = 4*{n}^2*{2*n+1}^2*{5*n-4}*{5*n+6}")

# Check: is a_n related to a factored form?
print("\n=== Checking a_n roots ===")
# a_n = -220n³ - 484n² - 301n - 42 = 0
# Dividing by -1: 220n³ + 484n² + 301n + 42 = 0
# Try rational roots: ±{1,2,3,6,7,14,21,42}/{1,2,4,5,10,11,20,22,44,55,110,220}
from fractions import Fraction
print("Rational roots of 220n³ + 484n² + 301n + 42:")
for num in [1, 2, 3, 6, 7, 14, 21, 42]:
    for den in [1, 2, 4, 5, 10, 11, 20, 22, 44, 55, 110, 220]:
        for sign in [1, -1]:
            r = Fraction(sign * num, den)
            val = 220*r**3 + 484*r**2 + 301*r + 42
            if val == 0:
                print(f"  n = {r} is a root!")

# Also try: -1/n factors  
# 220 = 4*5*11, 42 = 2*3*7
# Try n = -1/4, -1/5, -3/4, -7/20, -2/5, -3/10, -7/10, -7/4, -7/20, -21/20
for r in [Fraction(-1,4), Fraction(-1,5), Fraction(-3,4), Fraction(-7,20), 
          Fraction(-2,5), Fraction(-3,10), Fraction(-7,10), Fraction(-7,4),
          Fraction(-7,20), Fraction(-21,20), Fraction(-1,11), Fraction(-2,11),
          Fraction(-3,11), Fraction(-6,11), Fraction(-7,11), Fraction(-14,11),
          Fraction(-21,11), Fraction(-42,11), Fraction(-1,20), Fraction(-3,20),
          Fraction(-7,44), Fraction(-7,55), Fraction(-6,55), Fraction(-21,44),
          Fraction(-21,110), Fraction(-42,55), Fraction(-6,220), Fraction(-42,220)]:
    val = 220*r**3 + 484*r**2 + 301*r + 42
    if val == 0:
        print(f"  n = {r} is a root!")

