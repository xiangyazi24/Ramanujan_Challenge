#!/usr/bin/env python3
"""Problem 2.1: Identify the ₃F₂ underlying the PCF for π.
Strategy: compute the minimal solution via backward recurrence,
then try to match individual terms with known hypergeometric sums."""
from mpmath import mp, mpf, sqrt, pi, fac, binomial, hyper, gamma, mpf as F
from fractions import Fraction

mp.dps = 80

def a(n): return mpf(-220)*n**3 - 484*n**2 - 301*n - 42
def b(n): return 4*n**2*(2*n+1)**2*(5*n-4)*(5*n+6)

# Compute minimal solution by backward recurrence
N = 500
tail = [mpf(0)] * (N + 2)
tail[N+1] = mpf(0)
tail[N] = mpf(1)
for n in range(N, 0, -1):
    tail[n-1] = (tail[n+1] - a(n)*tail[n]) / b(n)

# Normalize: tail[0] = 1
norm = tail[0]
for i in range(N+1):
    tail[i] /= norm

# The CF value = a(0) + b(1)*tail[0]/tail_minus1
# tail_minus1 = (tail[1] - a(1)*tail[0]) / b(1)
tail_m1 = (tail[1] - a(1)*tail[0]) / b(1)
cf_value = a(0) + b(1) / (tail[1]/tail[0])  # wrong, let me use the direct backward

# Actually, just recompute the CF directly
cf = a(N)
for n in range(N-1, -1, -1):
    cf = a(n) + b(n+1) / cf
target = 6/(3-pi)
print(f"CF = {cf}")
print(f"6/(3-pi) = {target}")
print(f"diff = {cf - target}")

# Now try to identify the minimal solution terms
# After gauge: m_n / (n!)^3 should be the hypergeometric-like sequence
phi = (1+sqrt(mpf(5)))/2
rp = 20*phi**(-5)  # minimal Poincaré root

print(f"\nMinimal solution (gauged by (n!)^3 * r+^n):")
print(f"r+ = {rp}")
for n in range(1, 20):
    gauged = tail[n] / (fac(n)**3 * rp**n)
    print(f"  n={n:2d}: tail[n]/((n!)^3 * r+^n) = {mp.nstr(gauged, 20)}")

# Try to identify: is the gauged sequence related to binomial sums?
# Apéry-style: sum_k (-1)^k binom(n,k)^a binom(n+k,k)^b ...

# Compute the RATIO m_{n+1}/m_n after removing (n!)^3 * r+^n gauge
print(f"\nRatios of gauged sequence:")
for n in range(1, 15):
    g_n = tail[n] / (fac(n)**3 * rp**n)
    g_n1 = tail[n+1] / (fac(n+1)**3 * rp**(n+1))
    ratio = g_n1 / g_n
    print(f"  g_{n+1}/g_n = {mp.nstr(ratio, 20)}")

# Try specific ₃F₂ candidates
# The b_n factors suggest parameters involving 0, 1/2, -4/5, 6/5
# Try: ₃F₂(-n, n+1/2, ?; ?, ?; z) for various z

# Check: does ₃F₂(1/2, 4/5, -1/5; 1, 1; z) at some z give something pi-related?
print(f"\n=== Testing ₃F₂ candidates ===")

# Candidate 1: ₃F₂(1/2, -4/5, 6/5; 1, 1; z)
# At z = -1/(20*phi^5) (related to Poincaré root)
z_test = -1/(20*phi**5)
val = hyper([F(1)/2, F(-4)/5, F(6)/5], [1, 1], z_test)
print(f"₃F₂(1/2,-4/5,6/5; 1,1; {mp.nstr(z_test,8)}) = {mp.nstr(val,15)}")

# Candidate 2: at z = 1
# ₃F₂ at z=1 needs Saalschütz or balanced conditions
# Balanced: a+b+c+1 = d+e
# 1/2 + (-4/5) + 6/5 + 1 = d + e => 1/2 + 2/5 + 1 = d+e => 19/10 = d+e

# Try: ₃F₂(1/10, 1/2, 9/10; 1, 1; z) — parameters with fifths
for params_top in [(F(1)/10, F(1)/2, F(9)/10),
                    (F(1)/5, F(1)/2, F(4)/5),
                    (F(2)/5, F(1)/2, F(3)/5),
                    (F(1)/2, F(-4)/5, F(6)/5)]:
    for z in [F(1), F(-1), F(1)/2, -F(1)/4]:
        try:
            val = hyper(list(params_top), [1, 1], z)
            # Check if val is related to pi
            for c in [pi, 1/pi, pi**2, 1/(3-pi), 6/(3-pi)]:
                ratio = val / c
                if abs(ratio - round(float(ratio))) < 0.01:
                    print(f"  HIT: ₃F₂({params_top}; 1,1; {z}) = {mp.nstr(val,10)} ≈ {round(float(ratio))}*{c}")
        except:
            pass

# Try Gauss ₂F₁ connection: π = 4*arctan(1) = 4*₂F₁(1/2, 1; 3/2; -1)... no
# Actually π/4 = ₂F₁(1/2, 1; 3/2; 1) = ... well, through arctan(1)

# The value 6/(3-π) = -6/(π-3). Let's see if (π-3)/6 has a nice ₃F₂ form.
print(f"\n(π-3)/6 = {(pi-3)/6}")
print(f"(π-3) = {pi-3}")
print(f"3-π = {3-pi}")
# 3-π = 3 - 4*arctan(1) = integral_0^1 (3x²-1)/(1+x²) dx

# Key: the PCF may be related to a RATIO of contiguous ₃F₂ values
# CF = ₃F₂(a+1,b,c;d,e;z) / ₃F₂(a,b,c;d,e;z) after equivalence transform
# The contiguous relation generates a 3-term recurrence in a parameter

print(f"\n=== Contiguous ₃F₂ ratio test ===")
# If CF = F(a+1)/F(a) for some parameter, then shifting a → a+n gives the CF
# The PCF converges to -m_0/m_{-1} = 6/(3-π)
# So we need: ₃F₂(a_0+1,...)/₃F₂(a_0,...) = 6/(3-π) for some a_0 and parameters

# With b_n having factors (5n-4)(5n+6), the shift is in a parameter by 1/5 steps
# Actually the shift is by 1 in n, which corresponds to shift by 5 in the 5n parameter

# Try: F(n) = ₃F₂(something involving n; ...; z)
# and the recurrence matches our a_n, b_n

# For now, let me check: is a_0 = -42 = -6*7 related to any ₃F₂ at n=0?
print(f"a_0 = {int(a(0))} = -6*7")
print(f"b_1 = {int(b(1))} = 4*1*9*1*11 = 396")
print(f"a_1 = {int(a(1))} = -1047")
print(f"CF = a_0 + b_1/(a_1 + ...) = -42 + 396/(-1047 + ...)")
print(f"Note: -42 = a_0, and 6/(3-π) ≈ -42.375")
print(f"So the CF correction beyond a_0 is ≈ -0.375 = -3/8")
