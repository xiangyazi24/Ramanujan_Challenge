#!/usr/bin/env python3
"""Problem 2.5: Extract L₁ EXACTLY from the order-3 recurrence.

Step 1: Derive α₀, α₁, α₂, α₃ as sympy polynomials
Step 2: Verify α₀+α₁+α₂+α₃ ≡ 0 (constants-killing condition)
Step 3: Read off L₁ = α₃ S² + (α₂+α₃) S + (-α₀)
Step 4: Find the Poincaré polynomial of L₁
"""
from sympy import symbols, Matrix, Poly, factor, simplify, Rational, expand, gcd
from sympy import cancel, quo, rem

n = symbols('n')

def M_sym(N):
    m11 = (-2*N-5)*(N+3)**2 * (136*N**4 + 1424*N**3 + 5548*N**2 + 9551*N + 6141)
    m12 = 384*N**6 + 6384*N**5 + 44168*N**4 + 162698*N**3 + 336377*N**2 + 369933*N + 169011
    m13 = -480*N**4 - 4980*N**3 - 19210*N**2 - 32690*N - 20730
    m21 = (N+2)**2*(N+3)**2*(4*N+10)*(48*N**3 + 386*N**2 + 1017*N + 879)
    m22 = (N+2)**2*(-272*N**5 - 3848*N**4 - 21732*N**3 - 61184*N**2 - 85761*N - 47808)
    m23 = (N+2)**2*(320*N**3 + 2540*N**2 + 6610*N + 5640)
    m31 = (-4*N-10)*(N+2)**2*(N+3)**2*(32*N**4 + 302*N**3 + 1037*N**2 + 1530*N + 813)
    m32 = (N+2)**2*(192*N**6 + 2984*N**5 + 19116*N**4 + 64452*N**3 + 120256*N**2 + 117279*N + 46476)
    m33 = (N+2)**2*(-16*N**5 - 408*N**4 - 2912*N**3 - 8884*N**2 - 12254*N - 6240)
    return Matrix([[expand(m11), expand(m12), expand(m13)],
                   [expand(m21), expand(m22), expand(m23)],
                   [expand(m31), expand(m32), expand(m33)]])

print("Building symbolic matrices...")
M0 = M_sym(n)
M1 = M_sym(n+1)
M2 = M_sym(n+2)

# === Step 1: Derive the order-3 recurrence ===
print("Deriving order-3 recurrence coefficients...")

# 2x2 system for (p_N, r_N) elimination
A11 = M0[1,0]; A12 = M0[2,0]
A21 = expand(M0[1,1]*M1[1,0] + M0[1,2]*M1[2,0])
A22 = expand(M0[2,1]*M1[1,0] + M0[2,2]*M1[2,0])
det_sys = expand(A11*A22 - A12*A21)

B1_q0 = expand(-M0[0,0]); B1_q1 = Rational(1)
B2_q0 = expand(-(M0[0,1]*M1[1,0] + M0[0,2]*M1[2,0]))
B2_q1 = expand(-M1[0,0])

# (p_N, r_N) coefficients
pN_q0_num = expand(A22*B1_q0 - A12*B2_q0)
pN_q1_num = expand(A22*B1_q1 - A12*B2_q1)
pN_q2_num = expand(-A12)
rN_q0_num = expand(A11*B2_q0 - A21*B1_q0)
rN_q1_num = expand(A11*B2_q1 - A21*B1_q1)
rN_q2_num = expand(A11)

# p_{N+1}, r_{N+1}
pN1_q0_num = expand(M0[0,1]*det_sys + M0[1,1]*pN_q0_num + M0[2,1]*rN_q0_num)
pN1_q1_num = expand(M0[1,1]*pN_q1_num + M0[2,1]*rN_q1_num)
pN1_q2_num = expand(M0[1,1]*pN_q2_num + M0[2,1]*rN_q2_num)
rN1_q0_num = expand(M0[0,2]*det_sys + M0[1,2]*pN_q0_num + M0[2,2]*rN_q0_num)
rN1_q1_num = expand(M0[1,2]*pN_q1_num + M0[2,2]*rN_q1_num)
rN1_q2_num = expand(M0[1,2]*pN_q2_num + M0[2,2]*rN_q2_num)

# p_{N+2}, r_{N+2}
m01p = M1[0,1]; m11p = M1[1,1]; m21p = M1[2,1]
m02p = M1[0,2]; m12p = M1[1,2]; m22p = M1[2,2]
pN2_q0_num = expand(m11p*pN1_q0_num + m21p*rN1_q0_num)
pN2_q1_num = expand(m01p*det_sys + m11p*pN1_q1_num + m21p*rN1_q1_num)
pN2_q2_num = expand(m11p*pN1_q2_num + m21p*rN1_q2_num)
rN2_q0_num = expand(m12p*pN1_q0_num + m22p*rN1_q0_num)
rN2_q1_num = expand(m02p*det_sys + m12p*pN1_q1_num + m22p*rN1_q1_num)
rN2_q2_num = expand(m12p*pN1_q2_num + m22p*rN1_q2_num)

# q_{N+3} coefficients (×det_sys)
alpha0_raw = expand(M2[1,0]*pN2_q0_num + M2[2,0]*rN2_q0_num)
alpha1_raw = expand(M2[1,0]*pN2_q1_num + M2[2,0]*rN2_q1_num)
alpha2_raw = expand(M2[0,0]*det_sys + M2[1,0]*pN2_q2_num + M2[2,0]*rN2_q2_num)
alpha3_raw = det_sys

# GCD
print("Computing GCD...")
g01 = gcd(alpha0_raw, alpha1_raw)
print(f"  gcd(α₀,α₁) computed, deg = {Poly(g01, n).degree()}")
g23 = gcd(alpha2_raw, alpha3_raw)
print(f"  gcd(α₂,α₃) computed, deg = {Poly(g23, n).degree()}")
g_all = gcd(g01, g23)
p_gcd = Poly(g_all, n)
print(f"  gcd all = {factor(g_all)}")
print(f"  deg = {p_gcd.degree()}")

# Exact polynomial division using cancel
print("\nDividing out GCD (using cancel)...")
a0 = cancel(alpha0_raw / g_all)
a1 = cancel(alpha1_raw / g_all)
a2 = cancel(alpha2_raw / g_all)
a3 = cancel(alpha3_raw / g_all)

# Convert to Poly
p0 = Poly(a0, n); p1 = Poly(a1, n); p2 = Poly(a2, n); p3 = Poly(a3, n)
print(f"  Reduced degrees: α₃={p3.degree()}, α₂={p2.degree()}, α₁={p1.degree()}, α₀={p0.degree()}")

# === Step 2: Check α₀+α₁+α₂+α₃ ≡ 0 ===
print("\n=== Checking constants-killing condition ===")
total = expand(a0 + a1 + a2 + a3)
p_total = Poly(total, n) if total != 0 else None
if total == 0 or p_total is None:
    print("  α₀+α₁+α₂+α₃ = 0  ✓  (S-1) is a RIGHT FACTOR!")
    has_s_minus_1 = True
else:
    print(f"  α₀+α₁+α₂+α₃ ≠ 0, deg = {p_total.degree()}")
    print(f"  Sum = {total}")
    has_s_minus_1 = False

# === Step 3: Extract L₁ ===
if has_s_minus_1:
    print("\n=== Extracting L₁ ===")
    # L = α₃ S³ + α₂ S² + α₁ S + α₀
    # L = L₁ · (S - 1)
    # L₁ = β₂ S² + β₁ S + β₀
    # where β₂ = α₃, β₁ = α₂+α₃, β₀ = -α₀
    beta2 = a3
    beta1 = expand(a2 + a3)
    beta0 = expand(-a0)

    pb2 = Poly(beta2, n); pb1 = Poly(beta1, n); pb0 = Poly(beta0, n)
    print(f"  L₁ = β₂(n) S² + β₁(n) S + β₀(n)")
    print(f"  deg(β₂) = {pb2.degree()}")
    print(f"  deg(β₁) = {pb1.degree()}")
    print(f"  deg(β₀) = {pb0.degree()}")

    # Factor each coefficient
    print(f"\n  β₂ = {factor(beta2)}")
    print(f"\n  β₁ = {factor(beta1)}")
    print(f"\n  β₀ = {factor(beta0)}")

    # === Step 4: Poincaré polynomial ===
    print("\n=== Poincaré polynomial of L₁ ===")
    lc2 = pb2.LC()
    lc1 = pb1.LC()
    lc0 = pb0.LC()
    print(f"  Leading coefficients: β₂→{lc2}, β₁→{lc1}, β₀→{lc0}")
    print(f"  Poincaré: {lc2}t² + ({lc1})t + ({lc0}) = 0")
    # Normalize
    r1 = Rational(lc1, lc2)
    r0 = Rational(lc0, lc2)
    print(f"  Normalized: t² + ({r1})t + ({r0}) = 0")
    print(f"  Expected: t² - 34t + 1 = 0")

    # Check if it matches t² - 34t + 1
    if r1 == -34 and r0 == 1:
        print("  ✓ MATCHES t² - 34t + 1 = (t - (17+12√2))(t - (17-12√2))!")
    else:
        print(f"  ✗ Does NOT match: got t² + {r1}t + {r0}")

    # GCD of β₀, β₁, β₂
    print("\n=== Simplifying L₁ ===")
    g_L1 = gcd(gcd(beta0, beta1), beta2)
    if g_L1 != 1:
        p_gL1 = Poly(g_L1, n)
        print(f"  Common factor in L₁: {factor(g_L1)}")
        b2_red = cancel(beta2/g_L1)
        b1_red = cancel(beta1/g_L1)
        b0_red = cancel(beta0/g_L1)
        print(f"  Reduced: β₂ = {factor(b2_red)}")
        print(f"           β₁ = {factor(b1_red)}")
        print(f"           β₀ = {factor(b0_red)}")
    else:
        print("  No common factor")
        b2_red = beta2; b1_red = beta1; b0_red = beta0

    # Verify: L₁ · (S-1) = L
    print("\n=== Verification ===")
    # (β₂ S² + β₁ S + β₀)(S - 1) = β₂ S³ + (β₁-β₂) S² + (β₀-β₁) S - β₀
    check_a3 = expand(beta2 - a3)
    check_a2 = expand((beta1 - beta2) - a2)
    check_a1 = expand((beta0 - beta1) - a1)
    check_a0 = expand(-beta0 - a0)
    print(f"  α₃ match: {check_a3 == 0}")
    print(f"  α₂ match: {check_a2 == 0}")
    print(f"  α₁ match: {check_a1 == 0}")
    print(f"  α₀ match: {check_a0 == 0}")

else:
    # (S-1) is NOT a factor. Check numerically.
    print("\n(S-1) is not an exact factor. Checking numerically...")
    from sympy import N as symN
    for k in range(5, 15):
        val = symN(total.subs(n, k))
        print(f"  n={k}: α₀+α₁+α₂+α₃ = {val}")
