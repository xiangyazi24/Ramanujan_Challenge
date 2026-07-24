#!/usr/bin/env python3
"""
Test the KEY HYPOTHESIS:

P2.7 sequences = Zudilin sequences evaluated at n → n + 17/22

Specifically: the P2.7 error is proportional to the Zudilin Barnes integral
J_{n+17/22}^Z, where:

J_ν = (1/2πi) ∫ R_ν(t) (π/sin πt)² (1 - π cot πt) dt

R_ν(t) = Γ(t)⁴ / [Γ(t-ν)³ · Γ(ν+1)² · Γ(t+ν+1)]

The 17/22 shift comes from the coboundary: Q_209(n+83/22) vs Q_209(n+3),
where 83/22 - 3 = 17/22.

Strategy:
1. Compute the "shifted a_n" := a_{n+17/22} using the Gamma-function
   representation of a_n as a residue.
2. Compare to q_n^{P2.7}
"""
import mpmath
mpmath.mp.dps = 50

s = mpmath.mpf(17) / mpmath.mpf(22)  # The shift

# ============================================================
# Step 1: Compute a_n for integer n (verification)
# ============================================================
def a_int(n):
    """a_n = Σ C(n,k)² C(n+k,k) C(n+2k,n)"""
    s = mpmath.mpf(0)
    for k in range(n+1):
        term = mpmath.binomial(n,k)**2 * mpmath.binomial(n+k,k) * mpmath.binomial(n+2*k,n)
        s += term
    return s

print("=== Integer a_n values ===")
for n in range(6):
    print(f"  a_{n} = {a_int(n)}")

# ============================================================
# Step 2: Compute the "shifted a" via Gamma representation
# ============================================================
# From Zudilin's construction:
# R_ν(t) = Γ(t)⁴ / [Γ(t-ν)³ · Γ(ν+1)² · Γ(t+ν+1)]
#
# The "a" coefficient comes from the partial fraction at t=0:
# a_ν = (-1)^ν · Γ(ν+1)⁴ / [Γ(ν+1)² · Γ(-ν)³ · Γ(ν+1)] ... hmm
#
# Actually, a_n corresponds to the coefficient of Li_1 in the partial
# fraction decomposition. For integer n:
# a_n = Σ_{k=0}^n C(n,k)² C(n+k,k) C(n+2k,n)
#
# This has a Gamma function representation that EXTENDS to non-integer ν:
# a_ν = Σ_{k=0}^∞ [Γ(ν+1)/(Γ(k+1)Γ(ν-k+1))]² [Γ(ν+k+1)/(Γ(k+1)Γ(ν+1))]
#        × [Γ(ν+2k+1)/(Γ(2k+1)Γ(ν+1))]
# But this requires ν-k+1 > 0, i.e., k < ν+1, which for non-integer ν
# means the sum is infinite (all terms are nonzero).
#
# However, the ALTERNATIVE identity:
# a_n = (-1)^n Σ (-1)^k C(n,k) C(n+k,k)³
# has terms C(n,k) = Γ(n+1)/[Γ(k+1)Γ(n-k+1)] which vanish for k > n
# (because Γ(n-k+1) has a pole when n-k is a negative integer).
#
# For non-integer ν: C(ν,k) = Γ(ν+1)/[Γ(k+1)Γ(ν-k+1)]
# This does NOT vanish for k > ν (since ν-k is not a negative integer).
# So the sum is infinite.
#
# But it might converge! Let's check.

def R_nu(nu, t):
    """R_ν(t) = Γ(t)⁴ / [Γ(t-ν)³ · Γ(ν+1)² · Γ(t+ν+1)]"""
    try:
        return (mpmath.gamma(t)**4 /
                (mpmath.gamma(t-nu)**3 * mpmath.gamma(nu+1)**2 * mpmath.gamma(t+nu+1)))
    except:
        return mpmath.mpf(0)

# Check R_n(t) at integer n, t values
print("\n=== R_n(t) check at integer n ===")
for n in range(4):
    for t_val in [0.5, 1.5, 2.5]:
        t = mpmath.mpf(t_val)
        # Direct formula
        num = 1
        for j in range(1, n+1):
            num *= (t - j)**3
        den = mpmath.fac(n)**2 * t
        for j in range(1, n+1):
            den *= (t + j)
        r_direct = num/den if num != 0 else 0
        r_gamma = R_nu(n, t)
        print(f"  R_{n}({t_val}) = {r_direct:.10f} (direct), {r_gamma:.10f} (Gamma)")

# ============================================================
# Step 3: Compute the shifted "a_ν" using the convergent sum
# ============================================================
# a_ν = (-1)^ν Σ_{k=0}^∞ (-1)^k C(ν,k) C(ν+k,k)³
# = (-1)^ν Σ_{k=0}^∞ (-1)^k Γ(ν+1)/[Γ(k+1)Γ(ν-k+1)]
#   × [Γ(ν+k+1)/(Γ(k+1)Γ(ν+1))]³
# = (-1)^ν Σ_k (-1)^k Γ(ν+k+1)³ / [Γ(k+1)⁴ Γ(ν-k+1) Γ(ν+1)²]

# Check convergence: for large k, C(ν,k) ~ k^{-ν-1}/Γ(-ν) (asymptotic)
# and C(ν+k,k)³ ~ k^{3ν}/[Γ(ν+1)³] (growing like k^{3ν})
# So the general term ~ (-1)^k k^{3ν-ν-1} = (-1)^k k^{2ν-1}
# For ν = n + 17/22 > 0, this DIVERGES!

# So the alternating series doesn't converge for non-integer ν > 0.
# We need a different approach.

# ============================================================
# Step 4: Numerical Barnes integral via contour integration
# ============================================================
# J_ν = (1/2πi) ∫_{C-i∞}^{C+i∞} R_ν(t) (π/sin πt)² (1-π cot πt) dt
#
# We close the contour to the RIGHT and pick up residues at positive integers.
# For non-integer ν, R_ν(t) has NO zeros at positive integer t (unlike
# integer n where t=1,...,n are zeros).
#
# So the residue sum is:
# J_ν = -Σ_{m=1}^∞ Res_{t=m} [R_ν(t) (π/sin πt)² (1-π cot πt)]

# Near t = m (integer), sin(πt) ≈ (-1)^m π(t-m), so (π/sin πt)² ≈ 1/(t-m)²
# Also, π cot πt ≈ 1/(t-m) - π²(t-m)/3 + ...
# So (1 - π cot πt) ≈ 1 - 1/(t-m) + O(t-m)
# And R_ν(t) (π/sin πt)² (1-π cot πt) ≈ R_ν(m)/(t-m)² - R_ν(m)/(t-m)³ + R_ν'(m)/(t-m) + ...

# The residue at t=m is:
# Res = R_ν'(m) - R_ν(m) [coefficient of 1/(t-m) in (π/sin πt)² (1 - π cot πt)]
# ... this is getting complicated.

# Let me use a cleaner approach. From Q5048:
# r̃_ν = -Σ R_ν'(m) (contribution to ζ(2))
# r̃̃_ν = ½ Σ R_ν''(m) (contribution to ζ(3))
# J_ν = r̃_ν + r̃̃_ν = -Σ R_ν'(m) + ½ Σ R_ν''(m)

# For non-integer ν, these sums have terms for ALL m ≥ 1 (no cancellation).
# But R_ν(m) decays rapidly for m > ν, so the sums converge.

# Let me compute this numerically.

def R_nu_deriv1(nu, t):
    """dR_ν/dt using numerical differentiation"""
    h = mpmath.mpf(1e-10)
    return (R_nu(nu, t+h) - R_nu(nu, t-h)) / (2*h)

def R_nu_deriv2(nu, t):
    """d²R_ν/dt² using numerical differentiation"""
    h = mpmath.mpf(1e-8)
    return (R_nu(nu, t+h) - 2*R_nu(nu, t) + R_nu(nu, t-h)) / h**2

# Test at integer ν first
print("\n=== Barnes integral at integer ν (should give Zudilin sequences) ===")
for n in range(4):
    nu = mpmath.mpf(n)
    # Sum residues at m=1, 2, ...
    r_tilde = mpmath.mpf(0)  # -Σ R'(m)
    r_tilde2 = mpmath.mpf(0)  # ½ Σ R''(m)
    a_from_res = mpmath.mpf(0)  # Σ R(m)/m² ... hmm

    max_m = 200
    for m in range(1, max_m):
        t_m = mpmath.mpf(m)
        Rm = R_nu(nu, t_m)
        Rm_d1 = R_nu_deriv1(nu, t_m)
        Rm_d2 = R_nu_deriv2(nu, t_m)

        r_tilde -= Rm_d1
        r_tilde2 += Rm_d2 / 2

    J_n = r_tilde + r_tilde2
    a_n = a_int(n)
    zeta2 = mpmath.zeta(2)
    zeta3 = mpmath.zeta(3)
    L = zeta2 + zeta3

    # J_n should be a_n * L - (b̃ + b̃̃)
    # So J_n / a_n should converge to something...
    if n > 0:
        print(f"  ν={n}: J = {J_n:.10f}, a = {a_n:.1f}, J/a = {J_n/a_n:.10f}")
    else:
        print(f"  ν={n}: J = {J_n:.10f}, a = {a_n:.1f}")
    print(f"         a*L = {a_n * L:.10f}, J - a*L = {J_n - a_n*L:.10f}")

# ============================================================
# Step 5: Barnes integral at shifted ν = n + 17/22
# ============================================================
print("\n=== Barnes integral at ν = n + 17/22 ===")
s_shift = mpmath.mpf(17) / mpmath.mpf(22)

# P2.7 initial values
q_p27 = [
    mpmath.mpf(-215040420000),
    mpmath.mpf(-167282265043404) / mpmath.mpf(905),
    mpmath.mpf(-964185327658080) / mpmath.mpf(6071),
]
p_p27 = [
    mpmath.mpf(-612218384750),
    mpmath.mpf(-9525021973931919) / mpmath.mpf(18100),
    mpmath.mpf(-29561828382772029) / mpmath.mpf(65380),
]

L_val = mpmath.zeta(2) + mpmath.zeta(3)
e_p27 = [p_p27[i] - L_val * q_p27[i] for i in range(3)]

for n in range(4):
    nu = mpmath.mpf(n) + s_shift

    r_tilde = mpmath.mpf(0)
    r_tilde2 = mpmath.mpf(0)

    max_m = 300
    for m in range(1, max_m):
        t_m = mpmath.mpf(m)
        Rm_d1 = R_nu_deriv1(nu, t_m)
        Rm_d2 = R_nu_deriv2(nu, t_m)

        r_tilde -= Rm_d1
        r_tilde2 += Rm_d2 / 2

    J_shifted = r_tilde + r_tilde2

    # Also compute "a" at shifted ν (the residue at t=0 analogue)
    # For non-integer ν, a_ν comes from summing R_ν(m) with appropriate weights
    a_shifted = mpmath.mpf(0)
    for m in range(1, max_m):
        a_shifted += R_nu(nu, mpmath.mpf(m))

    print(f"  ν = {n}+17/22 = {float(nu):.6f}:")
    print(f"    r̃ = {r_tilde:.10e}")
    print(f"    r̃̃ = {r_tilde2:.10e}")
    print(f"    J = {J_shifted:.10e}")
    print(f"    a_shifted = {a_shifted:.10e}")

    if n < len(e_p27):
        e = e_p27[n]
        print(f"    e_n^P27 = {e:.10e}")
        if abs(J_shifted) > 1e-30:
            print(f"    e/J = {e/J_shifted:.10e}")
        if abs(a_shifted) > 1e-30:
            print(f"    q_p27/a_shifted = {q_p27[n]/a_shifted:.10e}")
