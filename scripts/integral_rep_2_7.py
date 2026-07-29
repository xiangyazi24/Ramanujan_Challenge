#!/usr/bin/env python3
"""Problem 2.7: Investigate the integral representation.

The recurrence for q_n has Poincaré roots:
  r_0 = 0.8588... (real dominant)
  r_± = complex conjugate, |r_±| = 0.00105...

Since p_n/q_n → ζ(2)+ζ(3), and the recurrence is connected to Cooper's
level-11 system, let's investigate:

1. What ODE does the generating function satisfy?
2. Is there a Beukers-type integral for the error p_n - (ζ(2)+ζ(3))q_n?
3. The Poincaré polynomial 4μ³ - 220μ² + 8μ - 1 = 0 (where μ = 64c)
   - what Picard-Fuchs equation has this as its indicial data?

Cooper's level-11 ODE is:
  θ³ - 2x(5θ²+5θ+2)(2θ+1) + 8x²(7θ²+7θ+2)(2θ+1)²
  - 22x³(2θ+1)³ = 0
where θ = x d/dx.

The spectral factorization 16P₂₇((t-2)²/4) = -C₁₁(t)C₁₁(4-t)
suggests the Problem 2.7 ODE is a "diagonal" of the Cooper ODE.
"""
from mpmath import mp, mpf, zeta, log10, fabs, polyroots, identify

mp.dps = 50

# Poincaré polynomial of Problem 2.7: 4μ³ - 220μ² + 8μ - 1 = 0
# where μ = 64c, so the char poly for c is:
# A_∞ c³ - B_∞ c² + C_∞ c - D_∞ = 0
# with A_∞=991952896, B_∞=852459520, C_∞=484352, D_∞=946

# Normalized: divide by A_∞
a = mpf(991952896)
b = mpf(852459520)
c_coeff = mpf(484352)
d = mpf(946)

roots = polyroots([a, -b, c_coeff, -d])
print("Poincaré roots (c):")
for r in roots:
    print(f"  {mp.nstr(r, 20)}  |r| = {mp.nstr(abs(r), 20)}")

# The μ roots
print("\nPoincaré roots (μ = 64c):")
for r in roots:
    mu = 64 * r
    print(f"  μ = {mp.nstr(mu, 20)}  |μ| = {mp.nstr(abs(mu), 20)}")

# The irrationality measure
r0 = abs(roots[0])
r_sub = abs(roots[1])
if r0 < r_sub:
    r0, r_sub = r_sub, r0
print(f"\nDominant |r| = {mp.nstr(r0, 15)}")
print(f"Subdominant |r| = {mp.nstr(r_sub, 15)}")
print(f"Ratio = {mp.nstr(r_sub/r0, 15)}")
print(f"Digits/step = {mp.nstr(-log10(r_sub/r0), 6)}")
print(f"Irrationality measure bound = {mp.nstr(1 + log10(r0)/log10(r_sub), 8)}")

# Cooper's level-11 Poincaré polynomial
# θ³ at x→0: (n+1)³T_{n+1} ≈ ...
# The recurrence is: (n+1)³ T_{n+1} = 2(2n+1)(5n²+5n+2)T_n - 8n(7n²+1)T_{n-1} + 22n(2n-1)(n-1)T_{n-2}
# Leading degrees: 1, 20, 56, 22 → Poincaré poly: c³ - 20c² + 56c - 22 = 0
print("\n--- Cooper level-11 ---")
cooper_roots = polyroots([1, -20, 56, -22])
print("Cooper Poincaré roots:")
for r in cooper_roots:
    print(f"  {mp.nstr(r, 20)}")

# Check spectral factorization: 16 P_27((t-2)²/4) = -C_11(t) C_11(4-t)
# P_27(c) = c³ - (B/A)c² + (C/A)c - (D/A)
# C_11(t) = t³ - 20t² + 56t - 22
print("\n--- Spectral factorization check ---")
# C_11(t) * C_11(4-t)
# Let's compute numerically at a few points
def C11(t):
    return t**3 - 20*t**2 + 56*t - 22

def P27(c):
    return c**3 - (b/a)*c**2 + (c_coeff/a)*c - d/a

for t_val in [mpf(1), mpf(2), mpf(3), mpf('0.5'), mpf('-1')]:
    c_val = (t_val - 2)**2 / 4
    lhs = 16 * P27(c_val)
    rhs = -C11(t_val) * C11(4 - t_val)
    print(f"  t={mp.nstr(t_val,3)}: 16*P27(({t_val}-2)²/4) = {mp.nstr(lhs,15)}, -C11(t)*C11(4-t) = {mp.nstr(rhs,15)}, diff = {mp.nstr(lhs-rhs,5)}")

# The spectral map t ↦ (t-2)²/4 sends:
# Cooper root t₁ → c₁ = (t₁-2)²/4
print("\n--- Spectral map Cooper → P27 ---")
for r in cooper_roots:
    c_val = (r - 2)**2 / 4
    p_val = P27(c_val)
    print(f"  Cooper root t={mp.nstr(r, 12)} → c=(t-2)²/4={mp.nstr(c_val, 12)}, P27(c)={mp.nstr(p_val, 8)}")

# Check: do the P27 roots come from Cooper roots via the spectral map?
print("\nP27 roots via spectral map:")
for c_root in roots:
    # c = (t-2)²/4 → t = 2 ± 2√c
    t_plus = 2 + 2*mp.sqrt(c_root)
    t_minus = 2 - 2*mp.sqrt(c_root)
    print(f"  c={mp.nstr(c_root,12)} → t+ = {mp.nstr(t_plus,12)}, t- = {mp.nstr(t_minus,12)}")
    print(f"    C11(t+) = {mp.nstr(C11(t_plus),8)}, C11(t-) = {mp.nstr(C11(t_minus),8)}")

# Target identification
target = zeta(2) + zeta(3)
print(f"\nζ(2) + ζ(3) = {mp.nstr(target, 40)}")
print(f"ζ(2) = {mp.nstr(zeta(2), 40)}")
print(f"ζ(3) = {mp.nstr(zeta(3), 40)}")

# Try to identify components
# ζ(2) = π²/6 is a period (weight 2)
# ζ(3) is a period (weight 3, Apéry's constant)
# Their sum mixes weights — unusual for a single CMF
# But this is exactly what the 4-term recurrence does

# Check: the dominant Poincaré root
r0_exact = roots[0]  # should be real
print(f"\nDominant root: {mp.nstr(r0_exact, 30)}")
print(f"1/dominant: {mp.nstr(1/r0_exact, 30)}")

# Try to identify 1/r0
for desc, val in [
    ("1/r0", 1/r0_exact),
    ("64/r0", 64/r0_exact),
    ("1/(64r0)", 1/(64*r0_exact)),
]:
    result = identify(val, tol=1e-10)
    print(f"  identify({desc} = {mp.nstr(val, 15)}): {result}")
