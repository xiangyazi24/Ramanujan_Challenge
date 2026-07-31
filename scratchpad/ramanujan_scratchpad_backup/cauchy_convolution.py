#!/usr/bin/env python3
"""
Implement the Cauchy convolution from Q5090 §6.

C(z) = Π_j (1 - x_j z)^{Δ(x_j)}
where x_j are roots of 4x³ - 220x² + 8x - 1 = 0
and Δ(x) = σ₂₇(x) + 3/2

B(z) = C(z) · A(z) where A(z) = Σ a_n z^n (AESZ #209 OGF)

Check if b_n = [z^n] B(z) matches P2.7 (after possible normalization and correction).
"""
import mpmath
mpmath.mp.dps = 100

# ============================================================
# 1. Find roots of 4x³ - 220x² + 8x - 1 = 0
# ============================================================
poly_coeffs = [mpmath.mpf(4), mpmath.mpf(-220), mpmath.mpf(8), mpmath.mpf(-1)]
roots = mpmath.polyroots(poly_coeffs)
# Sort: real dominant root first, then complex conjugate pair
roots = sorted(roots, key=lambda r: -r.real if abs(r.imag) < 1e-30 else -100)
print("Poincaré roots (4μ³ - 220μ² + 8μ - 1 = 0):")
for j, x in enumerate(roots):
    print(f"  x_{j} = {x}")
    # Verify
    val = 4*x**3 - 220*x**2 + 8*x - 1
    print(f"       P(x_{j}) = {val}")

x0 = roots[0]  # dominant real root ~54.96
print(f"\nDominant root x₀ = {float(x0.real):.15f}")

# ============================================================
# 2. Compute σ₂₇(x) and Δ(x) = σ₂₇(x) + 3/2
# ============================================================
def sigma27(x):
    return 24*(4*x - 1) / (220*x**2 - 16*x + 3)

deltas = []
print("\nFormal powers and shifts:")
for j, x in enumerate(roots):
    s = sigma27(x)
    d = s + mpmath.mpf(3)/2
    deltas.append(d)
    print(f"  x_{j}: σ₂₇ = {s}, Δ = σ₂₇ + 3/2 = {d}")

# ============================================================
# 3. Compute ψ(z) = C'/C
# Two sign choices: C = Π(1-xz)^{+Δ} vs C = Π(1-xz)^{-Δ}
# The CORRECT one is C = Π(1-xz)^{-Δ} (so ψ = +Σ Δ_j x_j / (1-x_j z))
# because: b_n ~ x^n n^{-δ-3/2} and we need -δ-3/2 = σ₂₇ → δ = -Δ
# ============================================================
N_TERMS = 40

print("\n=== Trying BOTH signs ===")
for sign_label, sign in [("C = Π(1-xz)^{+Δ}", -1), ("C = Π(1-xz)^{-Δ}", +1)]:
    print(f"\n{'='*60}")
    print(f"  {sign_label}")
    print(f"{'='*60}")

    psi = []
    for k in range(N_TERMS):
        val = sign * sum(deltas[j] * roots[j]**(k+1) for j in range(3))
    psi.append(val)
    if k < 10:
        print(f"  ψ_{k} = {val}")
        # Check if close to rational
        if abs(val.imag) < 1e-50:
            rv = float(val.real)
            # Try to identify as p/q with small denominator
            from fractions import Fraction
            frac = Fraction(rv).limit_denominator(10000)
            print(f"       ≈ {frac} = {float(frac):.15f}")

# ============================================================
# 4. Compute C(z) coefficients via C' = ψ · C, c_0 = 1
# ============================================================
c = [mpmath.mpf(0)] * N_TERMS
c[0] = mpmath.mpf(1)

for n in range(N_TERMS - 1):
    # (n+1) c_{n+1} = Σ_{k=0}^n ψ_k c_{n-k}
    val = sum(psi[k] * c[n-k] for k in range(n+1))
    c[n+1] = val / (n+1)

print("\nC(z) coefficients (should be rational):")
for n in range(15):
    v = c[n]
    print(f"  c_{n} = {v}")
    if abs(v.imag) < 1e-50:
        from fractions import Fraction
        rv = float(v.real)
        if abs(rv) > 1e-20:
            frac = Fraction(rv).limit_denominator(100000)
            print(f"       ≈ {frac}")

# ============================================================
# 5. AESZ #209 sequence
# ============================================================
def binom(n, k):
    if k < 0 or k > n:
        return 0
    result = 1
    for i in range(k):
        result = result * (n - i) // (i + 1)
    return result

a = [0] * N_TERMS
for n in range(N_TERMS):
    val = 0
    for k in range(n+1):
        val += binom(n,k)**2 * binom(n+k,n) * binom(n+2*k,n)
    a[n] = val

print("\nAESZ #209 first values:")
for n in range(10):
    print(f"  a_{n} = {a[n]}")

# ============================================================
# 6. Cauchy convolution b_n = Σ_{m=0}^n c_m · a_{n-m}
# ============================================================
b = [mpmath.mpf(0)] * N_TERMS
for n in range(N_TERMS):
    b[n] = sum(c[m] * a[n-m] for m in range(n+1))

print("\nConvolution b_n = (c * a)_n:")
for n in range(15):
    print(f"  b_{n} = {b[n]}")

# ============================================================
# 7. P2.7 reference sequence
# ============================================================
def A_p27(n): return 1024 * (2*n+5)**4 * (2*n+7)**3 * (2*n+9)**3 * (946*n**2+6407*n+10860)
def B_p27(n): return 128 * (2*n+7)**3 * (2*n+9)**3 * (104060*n**6 + 1745370*n**5 + 12145238*n**4 + 44886481*n**3 + 92943995*n**2 + 102256019*n + 46709052)
def C_p27_coeff(n): return 16 * (n+3)**4 * (2*n+9)**3 * (3784*n**5 + 57792*n**4 + 351019*n**3 + 1059230*n**2 + 1587211*n + 944620)
def D_p27(n): return (n+3)**4 * (n+4)**6 * (946*n**2 + 4515*n + 5399)

q = [mpmath.mpf(0)] * 20
q[0] = mpmath.mpf('-215040420000')
q[1] = mpmath.mpf('-167282265043404') / mpmath.mpf('905')
q[2] = mpmath.mpf('-964185327658080') / mpmath.mpf('6071')
for n in range(2, 19):
    q[n+1] = (mpmath.mpf(B_p27(n))/A_p27(n) * q[n]
              - mpmath.mpf(C_p27_coeff(n-1))/A_p27(n-1) * q[n-1]
              + mpmath.mpf(D_p27(n-2))/A_p27(n-2) * q[n-2])

# Normalize to 64^n scaling
print("\nP2.7 sequence (divided by 64^n):")
q_scaled = [q[n] / mpmath.mpf(64)**n for n in range(15)]
for n in range(10):
    print(f"  q_{n}/64^n = {q_scaled[n]}")

# ============================================================
# 8. Compare ratios b_{n+1}/b_n vs q_{n+1}/q_n
# ============================================================
print("\n=== Growth comparison ===")
print("  n  | b_{n+1}/b_n       | q_{n+1}/q_n")
print("  ---+--------------------|------------------")
for n in range(min(12, N_TERMS-1)):
    if abs(b[n]) > 1e-50:
        rb = b[n+1] / b[n]
        rq = q_scaled[n+1] / q_scaled[n] if n < 14 else 0
        rb_r = float(rb.real) if hasattr(rb, 'real') else float(rb)
        rq_r = float(rq.real) if hasattr(rq, 'real') else float(rq)
        print(f"  {n:2d} | {rb_r:+.12f} | {rq_r:+.12f}")

# ============================================================
# 9. Check if b_n satisfies P2.7 (after 64^n normalization)
# ============================================================
print("\n=== P2.7 recurrence residual for b_n ===")
# P2.7 in 64^n normalization: the ratio u_{n+1}/u_n → x₀ ≈ 54.96
# Since q_n = 64^n · q_scaled_n, the recurrence for q_scaled uses 64-scaled coefficients
# u_{n+1} = (B/A · 64) u_n - (C_{n-1}/A_{n-1} · 64²) u_{n-1} + (D_{n-2}/A_{n-2} · 64³) u_{n-2}
# Actually no: the recurrence coefficients already encode the growth.
# Let me just check directly.
for n in range(2, min(12, N_TERMS-1)):
    pred = (mpmath.mpf(B_p27(n))/A_p27(n) * b[n]
            - mpmath.mpf(C_p27_coeff(n-1))/A_p27(n-1) * b[n-1]
            + mpmath.mpf(D_p27(n-2))/A_p27(n-2) * b[n-2])
    res = b[n+1] - pred
    rel = abs(res / b[n+1]) if abs(b[n+1]) > 1e-50 else 0
    rel_f = float(abs(rel))
    tag = " *** MATCH ***" if rel_f < 1e-10 else ""
    print(f"  n={n:2d}: rel = {rel_f:.6e}{tag}")

# Wait -- b_n uses the x_j roots directly (not x_j/64), so the recurrence
# must also use the un-scaled form. Let me check with the AESZ recurrence too.
print("\n=== AESZ rec check for a_n ===")
# AESZ #209 recurrence: (n+1)^3 a_{n+1} = ...
# Actually order 3: let me use ore_algebra later. For now check growth.

print("\n=== b_n / (normalization · q_n/64^n) ===")
if abs(b[0]) > 0 and abs(q_scaled[0]) > 0:
    norm = b[0] / q_scaled[0]
    print(f"  Normalization factor: b_0/q_0_scaled = {norm}")
    for n in range(min(12, N_TERMS)):
        if abs(q_scaled[n]) > 1e-50:
            ratio = b[n] / (norm * q_scaled[n])
            ratio_f = float(ratio.real) if hasattr(ratio, 'real') else float(ratio)
            print(f"  n={n:2d}: b_n / (norm · q_scaled_n) = {ratio_f:.15f}")
