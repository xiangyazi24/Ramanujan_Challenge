#!/usr/bin/env python3
"""High-precision verification of Problem 2.1: PCF for π.

Verify: a_0 + b_1/(a_1 + b_2/(a_2 + ...)) = 6/(3 - π)

where a_n = -220n³ - 484n² - 301n - 42, b_n = 4n²(2n+1)²(5n-4)(5n+6).

Also: investigate the underlying recurrence and attempt hypergeometric identification.
"""
from mpmath import mp, mpf, pi, sqrt, log, gamma, hyper, nstr, fac, rf

mp.dps = 250

def a(n):
    return -220*n**3 - 484*n**2 - 301*n - 42

def b(n):
    return 4*n**2 * (2*n+1)**2 * (5*n-4) * (5*n+6)

# Compute CF value by backward recurrence (Lentz/Steed or direct tail evaluation)
def cf_value(N=500):
    """Compute a_0 + b_1/(a_1 + b_2/(...)) by backward recurrence from n=N."""
    val = mpf(a(N))
    for n in range(N, 0, -1):
        val = a(n-1) + b(n) / val
    return val

target = mpf(6) / (3 - pi)
print(f"Target 6/(3-π) = {nstr(target, 50)}")

cf = cf_value(500)
print(f"CF(500)        = {nstr(cf, 50)}")
print(f"Difference     = {nstr(cf - target, 15)}")
print(f"|diff|         = {float(abs(cf - target)):.3e}")

# Also verify with N=200 and N=1000 to check convergence rate
cf200 = cf_value(200)
cf1000 = cf_value(1000)
print(f"\nCF(200)  diff  = {float(abs(cf200 - target)):.3e}")
print(f"CF(500)  diff  = {float(abs(cf - target)):.3e}")
print(f"CF(1000) diff  = {float(abs(cf1000 - target)):.3e}")

# Now investigate the 3-term recurrence:
# P_n = a_n * P_{n-1} + b_n * P_{n-2}
# Q_n = a_n * Q_{n-1} + b_n * Q_{n-2}
# with P_{-1}=1, P_0=a_0 and Q_{-1}=0, Q_0=1
print("\n=== Convergent analysis ===")
P_prev, P_curr = mpf(1), mpf(a(0))
Q_prev, Q_curr = mpf(0), mpf(1)
for n in range(1, 30):
    P_new = a(n) * P_curr + b(n) * P_prev
    Q_new = a(n) * Q_curr + b(n) * Q_prev
    P_prev, P_curr = P_curr, P_new
    Q_prev, Q_curr = Q_curr, Q_new
    ratio = P_curr / Q_curr
    err = float(abs(ratio - target))
    if n <= 15 or n % 5 == 0:
        print(f"  n={n:3d}: P/Q - target = {err:.3e}")

# Poincaré analysis: characteristic equation of the recurrence
# y_{n+1} = a_n * y_n + b_n * y_{n-1}
# For large n: a_n ~ -220n³, b_n ~ 4·4·25·n⁶ = 400n⁶
# So characteristic: c = -220 + 400/c => c² + 220c - 400 = 0
# c = (-220 ± √(220² + 4·400))/2 = (-220 ± √(48400 + 1600))/2 = (-220 ± √50000)/2
# = (-220 ± 100√5)/2 = -110 ± 50√5
print("\n=== Poincaré roots ===")
r1 = -110 + 50*sqrt(5)
r2 = -110 - 50*sqrt(5)
print(f"r1 = -110 + 50√5 = {nstr(r1, 30)}")
print(f"r2 = -110 - 50√5 = {nstr(r2, 30)}")
print(f"r1*r2 = {nstr(r1*r2, 15)} (should be -400)")
print(f"r1+r2 = {nstr(r1+r2, 15)} (should be -220)")

# Golden ratio connection
phi = (1 + sqrt(5)) / 2
print(f"\nφ = {nstr(phi, 20)}")
print(f"φ⁵ = {nstr(phi**5, 20)}")
print(f"r1/20 = {nstr(r1/20, 20)}")
print(f"-φ⁻⁵ = {nstr(-phi**(-5), 20)}")
# r1 = -110 + 50√5 = 20(-5.5 + 2.5√5) = 20·(-11/2 + 5√5/2)
# φ⁵ = (1+√5)⁵/32 = (11+5√5)/2 ≈ 11.09
# So r1 = 20·(5√5-11)/2 = 10(5√5-11) and φ⁵ = (11+5√5)/2
# r1 = 20·(φ⁵ - 11) = 20φ⁵ - 220... no.
# Let me check: -r2 = 110 + 50√5 = 10(11+5√5) = 20·(11+5√5)/2 = 20φ⁵
print(f"-r2 = {nstr(-r2, 20)}")
print(f"20φ⁵ = {nstr(20*phi**5, 20)}")
print(f"r1 = {nstr(r1, 20)}")
print(f"-20φ⁻⁵ = {nstr(-20*phi**(-5), 20)}")
# So Poincaré roots are -20φ⁻⁵ and -20φ⁵

# Now let's look for hypergeometric identification.
# The gauge: if we set y_n = Γ-product · ỹ_n, the recurrence simplifies.
# Factor b_n = 4n²(2n+1)²(5n-4)(5n+6)
# = 4 · n · n · (2n+1) · (2n+1) · (5n-4) · (5n+6)
# Pochhammer form: rising factorials at half-integers and fifths
# (5n-4) = 5(n - 4/5) and (5n+6) = 5(n + 6/5)
# So b_n = 100 · n² · (n+1/2)² · (n - 4/5) · (n + 6/5)

# The ₃F₂ connection:
# For a ₃F₂(a,b,c; d,e; z) recurrence, the Poincaré roots at z determine the CF.
# With roots ±20φ^{±5}, z = 1/(product of roots) or similar.

# Let's check: is there a ₃F₂ with these parameters?
# The recurrence y_{n+1} - a_n y_n - b_n y_{n-1} = 0 has degree pattern:
# a_n = O(n³), b_n = O(n⁶) - this is characteristic of a balanced ₃F₂
# after appropriate gauge transform (multiply by (n!)^k etc.)

# Try gauge: Q_n = n!³ · (something) · q̃_n
# With b_n/((n+1)³ · n³) = 4(2n+1)²(5n-4)(5n+6) / ((n+1)³ · n)
# Hmm, let me try to compute the ratio q_{n+1}/q_n for the recessive solution.

print("\n=== Recessive solution ratio ===")
# Compute forward recurrence with high precision
# y_{n+1} = a_n * y_n + b_n * y_{n-1}
# Start with dominant + recessive mixture, then the ratio stabilizes to dominant.
# For recessive, use backward recurrence.
# Let's compute backward: y_{N} = 1, y_{N-1} = 0, and go back.
N = 300
y = [mpf(0)] * (N+1)
y[N] = mpf(1)
y[N-1] = mpf(0)
for n in range(N-1, 0, -1):
    # y_{n+1} = a_n * y_n + b_n * y_{n-1}
    # => y_{n-1} = (y_{n+1} - a_n * y_n) / b_n
    y[n-1] = (y[n+1] - a(n) * y[n]) / b(n)

# Normalize so y[0] = 1
if y[0] != 0:
    norm = y[0]
    for i in range(N+1):
        y[i] /= norm

print("Recessive solution (normalized y[0]=1):")
for n in range(min(10, N)):
    print(f"  y[{n}] = {nstr(y[n], 30)}")

# Check ratio y[n+1]/y[n] for large n
print("\nRatio y[n+1]/y[n] for recessive:")
for n in [5, 10, 20, 50, 100, 150]:
    if n < N and y[n] != 0:
        ratio = y[n+1]/y[n]
        print(f"  n={n:3d}: ratio = {nstr(ratio, 20)}, ratio/n³ ~ {nstr(ratio/n**3, 10)}")

print("\n=== Searching for ₃F₂ parameters ===")
# The CF a_0 + K(b_n/a_n) with the specific degrees suggests
# a ₃F₂ or ₄F₃ structure. Let me check the factorization of a_n and b_n.
# a_n = -220n³ - 484n² - 301n - 42
# Factor: a_n at n=0: -42. n=1: -220-484-301-42 = -1047.
# a_n = -(220n³ + 484n² + 301n + 42)
# Try to factor over Q: 220 = 4·5·11, 42 = 2·3·7
# Rational roots: ±1, ±2, ±3, ±6, ±7, ±14, ±21, ±42 and fractions
# a(n) = 0: 220n³ + 484n² + 301n + 42 = 0
# n = -1/5: 220(-1/125) + 484(1/25) + 301(-1/5) + 42
#         = -1.76 + 19.36 - 60.2 + 42 = -0.6 ≠ 0
# n = -2/5: 220(-8/125) + 484(4/25) + 301(-2/5) + 42
#         = -14.08 + 77.44 - 120.4 + 42 = -15.04 ≠ 0
# n = -7/10: 220(-343/1000) + 484(49/100) + 301(-7/10) + 42
#         = -75.46 + 237.16 - 210.7 + 42 = -7.0 ≠ 0
# n = -3/5: 220(-27/125) + 484(9/25) + 301(-3/5) + 42
#         = -47.52 + 174.24 - 180.6 + 42 = -11.88 ≠ 0
# n = -6/5: 220(-216/125) + 484(36/25) + 301(-6/5) + 42
#         = -380.16 + 696.96 - 361.2 + 42 = -2.4 ≠ 0
# Hmm. Let me just find the roots numerically.
from mpmath import polyroots
roots_a = polyroots([-220, -484, -301, -42])
print("Roots of a_n = 0:")
for r in roots_a:
    print(f"  {nstr(r, 20)}")
# So a_n probably doesn't factor nicely over Q.

# Factor b_n more carefully:
# b_n = 4n²(2n+1)²(5n-4)(5n+6)
# In Pochhammer notation:
# n² = (n)_1 · n, (2n+1)² = (2n+1)²
# (5n-4)(5n+6) = 25n² + 30n - 20n - 24 = 25n² + 10n - 24
# = 25(n + 1/5)(n - 24/25)... no
# (5n-4) = 5(n - 4/5), (5n+6) = 5(n + 6/5)
# So b_n = 4 · 5 · 5 · n² · (n+1/2)² · (n-4/5) · (n+6/5)
# = 100 n² (n+1/2)² (n-4/5)(n+6/5)

# For a generalized hypergeometric recurrence of the form
# c_2(n) y_{n+1} + c_1(n) y_n + c_0(n) y_{n-1} = 0
# with c_2, c_0 being products of linear factors and c_1 polynomial,
# the ₃F₂ parameters come from the ratio c_0(n)/c_2(n).

# Here: y_{n+1} = a_n y_n + b_n y_{n-1}, i.e., c_2 = 1, c_1 = -a_n, c_0 = -b_n.
# This is NOT in the standard form yet. We need a gauge.

# Standard ₃F₂ recurrence (contiguous relation):
# (n+d₁)(n+d₂) y_{n+1} - [(2n+...)(stuff)] y_n + (n+a₁)(n+a₂)(n+a₃)/(stuff) y_{n-1} = 0
# This gets complicated. Let me try a different approach:
# compute the first ~20 terms of Q_n, look for pattern.

print("\n=== Forward recurrence (dominant solution Q_n) ===")
# Q_{-1} = 0, Q_0 = 1
Q = [mpf(0)] * 50
Q[0] = mpf(0)  # Q_{-1}
Q[1] = mpf(1)  # Q_0
for n in range(1, 49):
    # Q_{n} = a_{n-1} * Q_{n-1} + b_{n-1} * Q_{n-2}
    # Wait, we need to be careful with indexing.
    # The CF recurrence is: H_n = a_n H_{n-1} + b_n H_{n-2}
    # with H_{-1} = 1, H_0 = a_0 for numerator P
    # and H_{-1} = 0, H_0 = 1 for denominator Q
    pass

# Let me redo with standard indexing
P = {}
Q2 = {}
P[-1] = mpf(1)
P[0] = mpf(a(0))
Q2[-1] = mpf(0)
Q2[0] = mpf(1)
for n in range(1, 40):
    P[n] = a(n) * P[n-1] + b(n) * P[n-2]
    Q2[n] = a(n) * Q2[n-1] + b(n) * Q2[n-2]

print("Convergents P_n/Q_n:")
for n in range(0, 20):
    if Q2[n] != 0:
        conv = P[n] / Q2[n]
        err = float(abs(conv - target))
        print(f"  n={n:2d}: err={err:.3e}  log10(|Q_n|)={float(log(abs(Q2[n]),10)):.1f}")

# Check the growth rate of Q_n
print("\nQ_n growth (log10):")
for n in range(0, 35):
    if Q2[n] != 0:
        print(f"  n={n:2d}: log10|Q_n| = {float(log(abs(Q2[n]),10)):.2f}")
