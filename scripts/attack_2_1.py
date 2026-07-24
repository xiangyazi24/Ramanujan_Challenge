#!/usr/bin/env python3
"""Problem 2.1: Attack the PCF identity via gauge transform + hypergeometric identification.

Key facts:
- a_n = -220n³ - 484n² - 301n - 42
- b_n = 4n²(2n+1)²(5n-4)(5n+6) = 100·n²·(n+1/2)²·(n-4/5)·(n+6/5)
- Poincaré roots: r₁ = -20φ⁻⁵, r₂ = -20φ⁵
- CF = 6/(3-π) ≈ -42.375...

Strategy: Find gauge g_n such that y_n = g_n · z_n transforms the recurrence
y_{n+1} = a_n y_n + b_n y_{n-1} into a ₃F₂ contiguous relation.
"""
from mpmath import (mp, mpf, pi, sqrt, gamma, rf, fac, nstr, log,
                    hyper, polyroots, matrix, binomial, beta, quad, exp)
from mpmath import mp as MP

mp.dps = 100

phi = (1 + sqrt(5)) / 2

# ---- Recurrence ----
def a(n):
    return -220*n**3 - 484*n**2 - 301*n - 42

def b(n):
    return 4*n**2 * (2*n+1)**2 * (5*n-4) * (5*n+6)

# ---- Error term analysis ----
# Compute P_n, Q_n, and the error ε_n = P_n - (6/(3-π)) Q_n
target = mpf(6) / (3 - pi)
print(f"Target = {nstr(target, 40)}")

P, Q = {}, {}
P[-1] = mpf(1); P[0] = mpf(a(0))
Q[-1] = mpf(0); Q[0] = mpf(1)
for n in range(1, 60):
    P[n] = a(n)*P[n-1] + b(n)*P[n-2]
    Q[n] = a(n)*Q[n-1] + b(n)*Q[n-2]

# Error term
eps = {}
for n in range(-1, 60):
    eps[n] = P[n] - target * Q[n]

print("\n=== Error term ε_n = P_n - (6/(3-π))·Q_n ===")
for n in range(0, 15):
    print(f"  ε[{n:2d}] = {nstr(eps[n], 25)}")

# Ratio ε_n / ε_{n-1} should approach recessive Poincaré root
print("\nRatio ε_n/ε_{n-1}:")
for n in range(2, 20):
    if eps[n-1] != 0:
        ratio = eps[n] / eps[n-1]
        print(f"  ε[{n:2d}]/ε[{n-1:2d}] = {nstr(ratio, 20)}")

# The error satisfies the SAME recurrence: ε_{n+1} = a_n ε_n + b_n ε_{n-1}
# but lives in the recessive subspace (geometric decay).
# The ratio should approach -20/φ⁵ = -20φ⁻⁵ ≈ 1.803...
recessive_root = -20 * phi**(-5)
print(f"\nRecessive Poincaré root = -20φ⁻⁵ = {nstr(recessive_root, 15)}")

# ---- Gauge transform ----
# Try gauge: y_n = g_n · z_n where g_n normalizes b_n
# b_n = 100 n² (n+1/2)² (n-4/5)(n+6/5)
# = 100 (n)² (n+1/2)² (n-4/5)(n+6/5)

# After gauge g_n = ∏_{k=1}^n f(k), the transformed recurrence is:
# z_{n+1} = (a_n g_n / g_{n+1}) z_n + (b_n g_{n-1} / g_{n+1}) z_{n-1}
# = ã_n z_n + b̃_n z_{n-1}

# Want b̃_n = b_n g_{n-1}/g_{n+1} = b_n / (f(n) f(n+1))

# Try 1: f(n) = n(2n+1)(5n+1)
# f(n)f(n+1) = n(2n+1)(5n+1)(n+1)(2n+3)(5n+6)
# b_n / [f(n)f(n+1)] = 4n(2n+1)(5n-4) / [(5n+1)(n+1)(2n+3)]

print("\n=== Gauge transform: f(n) = n(2n+1)(5n+1) ===")
def f1(n): return n * (2*n+1) * (5*n+1)
for n in range(1, 10):
    b_tilde = b(n) / (f1(n) * f1(n+1))
    a_tilde = a(n) * f1(n) / f1(n+1)
    print(f"  n={n}: ã = {nstr(a_tilde, 15)}, b̃ = {nstr(b_tilde, 15)}")

# Hmm, ã is still cubic in n (degree 3/3 = degree 0 only if f has same leading term).
# f(n) = n(2n+1)(5n+1) ~ 10n³. a_n ~ -220n³. So ã_n = a_n · f(n)/f(n+1) ~ -220n³ · 10n³/(10(n+1)³) ~ -220n³/1 → still cubic.

# Try different gauge to also reduce a_n.
# For ₃F₂ contiguous relation, we want a 3-term recurrence with:
# c₂(n) = (n+d₁)(n+d₂) [degree 2]
# c₀(n) = -z·(n+a₁-1)(n+a₂-1)(n+a₃-1) [degree 3]
# c₁(n) = polynomial [degree 2 or 3]
# So total is: deg 2 for leading, deg 3 for trailing, and the ratio c₀/c₂ has degree 1.

# Let me first check: for large n, the recurrence behaves as
# z_{n+1} ≈ (-220/coefficient) n^k z_n + (100/coefficient²) n^{2k} z_{n-1}
# where k = degree of f.
# For ₃F₂, we need the ratio z_{n+1}/z_n → constant (Poincaré roots are constants).
# Currently the ratio is ~n³. So we need f to be degree 3.
# After gauge of degree 3: ã ~ constant, b̃ ~ constant.
# This gives constant Poincaré roots = the gauge-transformed φ-powers.

# The gauge should be chosen so that f(n)·f(n+1) | b_n for all n.
# b_n = 4n²(2n+1)²(5n-4)(5n+6) = 100n²(n+1/2)²(n-4/5)(n+6/5)

# Try: f(n) = n(n+1/2)(n+6/5) = n(2n+1)(5n+6)/10
# f(n) = n(2n+1)(5n+6)/10
# f(n+1) = (n+1)(2n+3)(5n+11)/10
# f(n)f(n+1) = n(2n+1)(5n+6)(n+1)(2n+3)(5n+11)/100
# b_n/[f(n)f(n+1)] = 100n²(n+1/2)²(n-4/5)(n+6/5) · 100 / [n(2n+1)(5n+6)(n+1)(2n+3)(5n+11)]
# = 10000 n(n+1/2)(n-4/5)(n+6/5) / [(5n+6)(n+1)(2n+3)(5n+11)]
# Simplify (n+6/5)/(5n+6) = (5n+6)/(5(5n+6)) = 1/5
# n/(n+1) stays, (n+1/2)/(2n+3) = (2n+1)/(2(2n+3)), (n-4/5)/(5n+11) doesn't simplify
# This is messy.

# Let me try a systematic approach: use the Ore algebra viewpoint.
# The 2nd-order recurrence y_{n+1} - a_n y_n - b_n y_{n-1} = 0
# can be rewritten as: y_{n+1} + 220n³ y_n + ... = 0 (with appropriate signs)

# Actually, let me try a DIFFERENT approach entirely:
# Directly compute the error ε_n = P_n - α Q_n and try to express it as a sum/integral.

print("\n=== Trying to identify error with a known integral ===")
# The error satisfies the recurrence and decays geometrically.
# For Apéry-like proofs, ε_n often equals an integral:
# ε_n = ∫₀¹ ∫₀¹ x^n y^n R(x,y) dx dy / (stuff)
# or ε_n = ∫₀¹ P_n(x) f(x) dx for a specific f.

# Let's check if ε_n has a nice closed form.
# First, normalize: set ε̃_n = ε_n / Q_n (ratio of recessive to dominant)
print("Normalized error ε̃_n = ε_n / Q_n:")
for n in range(0, 15):
    if Q[n] != 0:
        eps_tilde = eps[n] / Q[n]
        print(f"  ε̃[{n:2d}] = {nstr(eps_tilde, 20)}")

# The error ratio ε̃_n should decay like (φ⁻⁵/φ⁵)^n = φ^{-10n}
print(f"\nφ⁻¹⁰ = {nstr(phi**(-10), 15)} (expected decay base)")
print("Ratio ε̃_{n+1}/ε̃_n:")
for n in range(2, 15):
    if Q[n] != 0 and Q[n-1] != 0 and eps[n-1] != 0:
        e1 = eps[n] / Q[n]
        e0 = eps[n-1] / Q[n-1]
        ratio = e1 / e0
        print(f"  n={n:2d}: {nstr(ratio, 20)}")

# ---- Alternative: contiguous ₃F₂ check ----
# Check if the CF can be written as a ratio of ₃F₂ values.
# Many known π formulas come from ₃F₂(a,b,c; d,e; 1) evaluations.
# With golden-ratio parameters, check Ebisu (2015):
# "Three-term relations for ₃F₂(1)" classifies all CFs of this type.

# The 1/5 parameters in b_n suggest a = -4/5, b = 6/5 or similar.
# Try: ₃F₂(1/2, -4/5, 6/5; d, e; z) for various d, e, z.

# From the Poincaré analysis: the argument z satisfies
# z · (leading coeff of b_n after gauge) = product of Poincaré roots
# After appropriate gauge: z = 1/φ^{10} = φ^{-10} ≈ 0.00813...

z_arg = phi**(-10)
print(f"\n=== Testing ₃F₂ at z = φ⁻¹⁰ = {nstr(z_arg, 15)} ===")

# Try some candidate ₃F₂ values and see if they relate to 6/(3-π)
candidates = [
    (('1/2', '1/5', '4/5'), ('1', '1')),
    (('1/2', '-4/5', '6/5'), ('1', '1')),
    (('1/2', '1/5', '4/5'), ('3/2', '1')),
    (('1/2', '1/10', '9/10'), ('1', '1')),
    (('1/2', '2/5', '3/5'), ('1', '1')),
]
for (a_params, b_params), _ in zip(candidates, range(len(candidates))):
    a_vals = [mpf(eval(x)) for x in a_params]
    b_vals = [mpf(eval(x)) for x in b_params]
    try:
        val = hyper(a_vals, b_vals, z_arg)
        print(f"  ₃F₂({a_params}; {b_params}; φ⁻¹⁰) = {nstr(val, 20)}")
    except:
        print(f"  ₃F₂({a_params}; {b_params}; φ⁻¹⁰) = FAILED")

# Also try at z=1 with different parameters
print("\n=== Testing ₃F₂ at z = 1 (Saalschutz, Dixon, etc.) ===")
# For balanced ₃F₂(a,b,c; d,e; 1), need a+b+c = d+e-1 (Saalschutz)
# a = 1/2, b = -4/5, c = 6/5: sum = 1/2 - 4/5 + 6/5 = 1/2 + 2/5 = 9/10
# Need d + e = 9/10 + 1 = 19/10

# But also, the PCF value -42.375... = 6/(3-π) doesn't look like a simple ₃F₂(1).

# ---- BREAKTHROUGH ATTEMPT: Recurrence factorization ----
# The 2nd order recurrence: y_{n+1} = a_n y_n + b_n y_{n-1}
# Equivalently: y_{n+1} - a_n y_n - b_n y_{n-1} = 0
# If a_n = p_n + q_n and b_n = p_n · q_{n-1} (where p, q are degree 3 polynomials),
# then the operator factors: (S - p_n)(S - q_n) = S² - (p_n+q_n)S + p_n q_{n-1}
# Wait, that gives S² - (p_n + q_{n-1})S + p_n q_{n-1}, not what we want.
# Actually: (S - p_n)(y) = y_{n+1} - p_n y_n. So (S - p_n)(S - q_n)(y) =
# y_{n+2} - p_{n+1} y_{n+1} - q_n(y_{n+1} - p_n y_n) =
# y_{n+2} - (p_{n+1} + q_n) y_{n+1} + p_n q_n y_n
# This has the WRONG sign for the last term (positive, not matching our negative b_n).

# Our recurrence: y_{n+1} = a_n y_n + b_n y_{n-1}, i.e.,
# y_{n+1} - a_n y_n - b_n y_{n-1} = 0
# As operator: S² - a_{n-1} S - b_{n-1} = 0 (applied at n → n+1)

# Hmm, let me try: can I find polynomials p(n), q(n) of degree 3 such that
# p(n) + q(n) = a_n and p(n) · q(n-1) = -b_n ?
# (because the operator (S-p_n)(S-q_n) = S² - (p_n+q_n)S + p_n·q_{n-1} = 0
#  gives S² - a_n S + p_n q_{n-1}, and we need p_n q_{n-1} = -b_n)

# System: p + q = a_n (as polynomials), p(n) · q(n-1) = -b_n
# With a_n = -220n³ - 484n² - 301n - 42
# and b_n = 4n²(2n+1)²(5n-4)(5n+6)

# Since p and q are degree 3, write:
# p(n) = α₃n³ + α₂n² + α₁n + α₀
# q(n) = (-220-α₃)n³ + (-484-α₂)n² + (-301-α₁)n + (-42-α₀)
# Then p(n)·q(n-1) must equal -b_n = -4n²(2n+1)²(5n-4)(5n+6)

# This is a system of polynomial identities. The product p(n)·q(n-1) is degree 6.
# We need it to equal a specific degree-6 polynomial.
# We have 4 free parameters (α₃, α₂, α₁, α₀) and 7 equations (coefficients of n⁶,...,n⁰).
# Generally overdetermined, but if such a factorization exists, it would be remarkable.

# Let's check numerically: for n=1,2,3,...,7, solve the linear system.
print("\n=== Checking operator factorization ===")
# p(n) + q(n) = a_n for all n
# p(n) * q(n-1) = -b_n for all n
# 4 unknowns α₃,α₂,α₁,α₀. Use 4 equations from n=1,2,3,4.

from mpmath import lu_solve

# p(n) = α₃n³ + α₂n² + α₁n + α₀
# q(n) = a_n - p(n) = (-220-α₃)n³ + (-484-α₂)n² + (-301-α₁)n + (-42-α₀)
# q(n-1) = (-220-α₃)(n-1)³ + (-484-α₂)(n-1)² + (-301-α₁)(n-1) + (-42-α₀)

# Constraint: p(n)·q(n-1) + b_n = 0 for each n.
# This is nonlinear in the α's. Let me try to solve numerically.

from mpmath import findroot

def factorization_equations(a3, a2, a1, a0):
    """Return residuals of p(n)*q(n-1) + b(n) = 0 for n=1,2,3,4."""
    residuals = []
    for n_val in [1, 2, 3, 4]:
        n = mpf(n_val)
        p_val = a3*n**3 + a2*n**2 + a1*n + a0
        q_coeff = [(-220-a3), (-484-a2), (-301-a1), (-42-a0)]
        nm1 = n - 1
        q_val = q_coeff[0]*nm1**3 + q_coeff[1]*nm1**2 + q_coeff[2]*nm1 + q_coeff[3]
        residuals.append(p_val * q_val + b(n_val))
    return residuals

try:
    result = findroot(factorization_equations, (-110, -242, -150, -21), tol=1e-80)
    a3, a2, a1, a0 = result
    print(f"Factorization found!")
    print(f"  p(n) = {nstr(a3,15)}n³ + {nstr(a2,15)}n² + {nstr(a1,15)}n + {nstr(a0,15)}")
    q3, q2, q1, q0 = -220-a3, -484-a2, -301-a1, -42-a0
    print(f"  q(n) = {nstr(q3,15)}n³ + {nstr(q2,15)}n² + {nstr(q1,15)}n + {nstr(q0,15)}")

    # Verify at more points
    print("  Verification:")
    for n_val in range(1, 10):
        n = mpf(n_val)
        p_val = a3*n**3 + a2*n**2 + a1*n + a0
        nm1 = n - 1
        q_val = q3*nm1**3 + q2*nm1**2 + q1*nm1 + q0
        residual = p_val * q_val + b(n_val)
        print(f"    n={n_val}: p*q(-1) + b = {nstr(residual, 10)}")
except Exception as e:
    print(f"No factorization found: {e}")

    # Try different initial guesses
    for guess in [(-100, -200, -100, -10), (-150, -300, -200, -30), (-110, -242, -151, -21)]:
        try:
            result = findroot(factorization_equations, guess, tol=1e-50)
            a3, a2, a1, a0 = result
            print(f"\nAlternative factorization with guess {guess}:")
            print(f"  p = {nstr(a3,15)}n³ + {nstr(a2,15)}n² + {nstr(a1,15)}n + {nstr(a0,15)}")
            q3, q2, q1, q0 = -220-a3, -484-a2, -301-a1, -42-a0
            print(f"  q = {nstr(q3,15)}n³ + {nstr(q2,15)}n² + {nstr(q1,15)}n + {nstr(q0,15)}")
            for n_val in [5, 6, 7, 8]:
                n = mpf(n_val)
                p_val = a3*n**3 + a2*n**2 + a1*n + a0
                nm1 = n - 1
                q_val = q3*nm1**3 + q2*nm1**2 + q1*nm1 + q0
                residual = p_val * q_val + b(n_val)
                print(f"    n={n_val}: residual = {nstr(residual, 10)}")
        except:
            pass

# ---- Alternative: check if a_n² + 4b_n is a perfect square (discriminant) ----
print("\n=== Discriminant a_n² + 4b_n ===")
for n_val in range(0, 10):
    disc = a(n_val)**2 + 4*b(n_val)
    print(f"  n={n_val}: a²+4b = {disc}, √ = {nstr(sqrt(abs(disc)), 15)}, sign={'+'if disc>0 else '-'}")
