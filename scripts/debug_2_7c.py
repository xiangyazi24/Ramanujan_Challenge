#!/usr/bin/env python3
"""Debug 2.7: Identify the actual limit and test initial condition variants."""
from mpmath import mp, mpf, zeta, log10, fabs, identify, pi

mp.dps = 150

def A(n):
    n = mpf(n)
    return 1024*(2*n+5)**4*(2*n+7)**3*(2*n+9)**3*(946*n**2+6407*n+10860)

def B(n):
    n = mpf(n)
    P6 = 104060*n**6+1745370*n**5+12145238*n**4+44886481*n**3+92943995*n**2+102256019*n+46709052
    return 128*(2*n+7)**3*(2*n+9)**3*P6

def C(n):
    n = mpf(n)
    P5 = 3784*n**5+57792*n**4+351019*n**3+1059230*n**2+1587211*n+944620
    return 16*(n+3)**4*(2*n+9)**3*P5

def D(n):
    n = mpf(n)
    return (n+3)**4*(n+4)**6*(946*n**2+4515*n+5399)

def gen(p0, p1, p2, N):
    u = [p0, p1, p2]
    for n in range(2, N):
        u_next = B(n)/A(n)*u[n] - C(n-1)/A(n-1)*u[n-1] + D(n-2)/A(n-2)*u[n-2]
        u.append(u_next)
    return u

target = zeta(2) + zeta(3)
N = 100

# Original
p0 = mpf(-612218384750)
p1 = mpf(-9525021973931919) / mpf(18100)
p2 = mpf(-29561828382772029) / mpf(65380)
q0 = mpf(-215040420000)
q1 = mpf(-1672822650043404) / mpf(905)
q2 = mpf(-964185327658080) / mpf(6071)

p = gen(p0, p1, p2, N)
q = gen(q0, q1, q2, N)
L = p[80]/q[80]
print(f"Current limit:  {mp.nstr(L, 50)}")
print(f"ζ(2)+ζ(3):      {mp.nstr(target, 50)}")
print(f"Difference:     {mp.nstr(L - target, 20)}")

# Try to identify the limit
print(f"\nTrying to identify {mp.nstr(L, 30)}:")
result = identify(L, tol=1e-15)
if result:
    print(f"  Identified: {result}")

# Check if q₁ numerator might have a different digit count
# Original: -1672822650043404/905
# Let me try: maybe it's -167282265043404/905 (one less digit — drop one '0' or '4')
print("\n--- Testing q₁ variants ---")
for q1_num in [
    -1672822650043404,  # original
    -167282265043404,   # one fewer digit (remove leading '1'? no, remove internal)
    -1672822650434040,  # different ending
    -16728226500434040, # extra digit
    -1672822265043404,  # swap 650→265
    -1672822650043440,  # last digits different
]:
    q1_test = mpf(q1_num) / mpf(905)
    q_test = gen(q0, q1_test, q2, N)
    p_test = gen(p0, p1, p2, N)
    ratio = p_test[80] / q_test[80]
    err = fabs(ratio - target)
    if err > 0:
        digits = -log10(err)
    else:
        digits = mpf('inf')
    print(f"  q1_num={q1_num}: p/q = {mp.nstr(ratio, 15)}, digits = {mp.nstr(digits, 5)}")

# Also test p₁ variants
print("\n--- Testing p₁ variants ---")
for p1_num in [
    -9525021973931919,  # original
    -9525021973931199,  # swap last digits
    -9525021937931919,  # internal swap
    -952502197393191,   # fewer digits
    -95250219739319190, # extra digit
]:
    p1_test = mpf(p1_num) / mpf(18100)
    p_test = gen(p0, p1_test, p2, N)
    q_test = gen(q0, q1, q2, N)
    ratio = p_test[80] / q_test[80]
    err = fabs(ratio - target)
    if err > 0:
        digits = -log10(err)
    else:
        digits = mpf('inf')
    print(f"  p1_num={p1_num}: p/q = {mp.nstr(ratio, 15)}, digits = {mp.nstr(digits, 5)}")

# The key insight: if both p and q satisfy the SAME recurrence, then p/q → constant.
# The constant is determined by the initial condition vectors (p0,p1,p2) and (q0,q1,q2).
# If ANY digit is wrong, the limit changes.

# Let me try a REVERSE approach: find what p₁ should be to get limit = ζ(2)+ζ(3)
# p_n/q_n → L means p_n ≈ L·q_n for large n.
# p_n = c₁ f₁(n) + c₂ f₂(n) + c₃ f₃(n)
# q_n = d₁ f₁(n) + d₂ f₂(n) + d₃ f₃(n)
# p_n/q_n → c₁/d₁  (dominant mode)
# So we need c₁/d₁ = ζ(2)+ζ(3)

# For the REMAINDER r_n = p_n - (ζ(2)+ζ(3))·q_n, the dominant mode must cancel:
# r_n = (c₁ - L·d₁)f₁(n) + ... = 0 · f₁(n) + ...
# So r_n is dominated by the subdominant modes.

# Compute r_n = p_n - target*q_n
print("\n--- Remainder r_n = p_n - (ζ(2)+ζ(3))·q_n ---")
for n in [5, 10, 20, 30, 40, 50, 60]:
    r = p[n] - target * q[n]
    if n > 0:
        r_prev = p[n-1] - target * q[n-1]
        ratio_r = r / r_prev if r_prev != 0 else mpf(0)
        print(f"  n={n}: r = {mp.nstr(r, 12)}, r/r_{{{n-1}}} = {mp.nstr(ratio_r, 8)}")
    else:
        print(f"  n={n}: r = {mp.nstr(r, 12)}")

# If the initial conditions are correct, r_n/r_{n-1} should approach the Poincaré root 0.8588
# If wrong, it approaches 0.8588 (same dominant mode) but with wrong proportionality constant
# The ratio r_n/r_{n-1} → dominant Poincaré root means r has a DOMINANT component, i.e.,
# p and q DON'T have the correct connection coefficient!
