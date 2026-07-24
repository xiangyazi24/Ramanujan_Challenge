#!/usr/bin/env python3
"""
Cauchy convolution test — try BOTH sign conventions for C(z).
Sign A: C = Π(1-x_j z)^{+Δ_j}  (Q5090 as stated)
Sign B: C = Π(1-x_j z)^{-Δ_j}  (corrected for transfer theorem)
"""
import mpmath
mpmath.mp.dps = 100

# Roots of 4μ³ - 220μ² + 8μ - 1 = 0
roots = mpmath.polyroots([mpmath.mpf(4), mpmath.mpf(-220), mpmath.mpf(8), mpmath.mpf(-1)])
# Identify dominant real root
real_root = [r for r in roots if abs(r.imag) < 1e-30]
complex_roots = [r for r in roots if abs(r.imag) > 1e-30]
x0 = real_root[0].real
print(f"Dominant root x₀ = {float(x0):.12f}")
print(f"Complex roots |x±| = {float(abs(complex_roots[0])):.12f}")

# σ₂₇ and Δ
def sigma27(x):
    return 24*(4*x - 1) / (220*x**2 - 16*x + 3)

deltas = [sigma27(r) + mpmath.mpf(3)/2 for r in roots]
print(f"\nΔ(x₀) = σ₂₇+3/2 = {deltas[roots.index(real_root[0])]}")

# AESZ #209
N = 30
def binom(n, k):
    if k < 0 or k > n: return 0
    r = 1
    for i in range(k): r = r * (n - i) // (i + 1)
    return r

a = [0] * N
for n in range(N):
    a[n] = sum(binom(n,k)**2 * binom(n+k,n) * binom(n+2*k,n) for k in range(n+1))

# P2.7 reference
def A_p(n): return 1024*(2*n+5)**4*(2*n+7)**3*(2*n+9)**3*(946*n**2+6407*n+10860)
def B_p(n): return 128*(2*n+7)**3*(2*n+9)**3*(104060*n**6+1745370*n**5+12145238*n**4+44886481*n**3+92943995*n**2+102256019*n+46709052)
def C_p(n): return 16*(n+3)**4*(2*n+9)**3*(3784*n**5+57792*n**4+351019*n**3+1059230*n**2+1587211*n+944620)
def D_p(n): return (n+3)**4*(n+4)**6*(946*n**2+4515*n+5399)

q = [mpmath.mpf(0)] * 20
q[0] = mpmath.mpf('-215040420000')
q[1] = mpmath.mpf('-167282265043404') / mpmath.mpf('905')
q[2] = mpmath.mpf('-964185327658080') / mpmath.mpf('6071')
for n in range(2, 19):
    q[n+1] = mpmath.mpf(B_p(n))/A_p(n)*q[n] - mpmath.mpf(C_p(n-1))/A_p(n-1)*q[n-1] + mpmath.mpf(D_p(n-2))/A_p(n-2)*q[n-2]

# Test BOTH signs
for sign_name, sign in [("Π(1-xz)^{+Δ}", -1), ("Π(1-xz)^{-Δ}", +1)]:
    print(f"\n{'='*60}")
    print(f"  C(z) = {sign_name}")
    print(f"{'='*60}")

    # ψ_k = sign * Σ Δ_j x_j^{k+1}
    psi = []
    for k in range(N):
        val = sign * sum(deltas[j] * roots[j]**(k+1) for j in range(3))
        psi.append(val)

    print(f"  ψ₀ = {float(psi[0].real):.6f} (should be rational)")
    print(f"  ψ₁ = {float(psi[1].real):.6f}")

    # C(z) via C' = ψ·C
    c = [mpmath.mpf(0)] * N
    c[0] = mpmath.mpf(1)
    for n in range(N-1):
        c[n+1] = sum(psi[k] * c[n-k] for k in range(n+1)) / (n+1)

    # Convolution b = c * a
    b = [mpmath.mpf(0)] * N
    for n in range(N):
        b[n] = sum(c[m] * a[n-m] for m in range(n+1))

    # Growth ratios
    print(f"\n  Growth ratios b_{'{n+1}'}/b_n (target: x₀ ≈ {float(x0):.4f}):")
    for n in range(min(15, N-1)):
        if abs(b[n]) > 1e-50:
            r = b[n+1] / b[n]
            print(f"    n={n:2d}: {float(r.real):+.8f}")

    # Compare with P2.7
    print(f"\n  Ratio b_n / q_n:")
    for n in range(min(10, N)):
        if abs(q[n]) > 1e-50 and abs(b[n]) > 1e-50:
            ratio = b[n] / q[n]
            print(f"    n={n:2d}: {float(ratio.real):+.10e}")

    # Check recurrence
    print(f"\n  P2.7 recurrence residual:")
    for n in range(2, min(12, N-1)):
        pred = (mpmath.mpf(B_p(n))/A_p(n)*b[n] - mpmath.mpf(C_p(n-1))/A_p(n-1)*b[n-1]
                + mpmath.mpf(D_p(n-2))/A_p(n-2)*b[n-2])
        res = b[n+1] - pred
        rel = float(abs(res/b[n+1])) if abs(b[n+1]) > 1e-50 else 0
        tag = " ***" if rel < 1e-8 else ""
        print(f"    n={n:2d}: {rel:.6e}{tag}")

    # Also check: does b_n satisfy the AESZ recurrence? (it shouldn't)
    # AESZ: for order 3, the recurrence of a_n
    # Actually let me just check the ratio b_n / a_n
    print(f"\n  Ratio b_n / a_n:")
    for n in range(min(10, N)):
        if a[n] != 0:
            ratio = b[n] / a[n]
            print(f"    n={n:2d}: {float(ratio.real):+.10e}")

# Also try: C = Π(1-xz)^{δ} with δ = -(σ₂₇+3/2) instead
print(f"\n{'='*60}")
print(f"  C(z) = Π(1-xz)^{{-Δ}} with CUSTOM normalization")
print(f"  (multiply by 64^n to match P2.7 scaling)")
print(f"{'='*60}")

# In this case, b_n = Σ c_m a_{n-m} with the -Δ sign
# Then b_n ~ x₀^n n^{σ₂₇(x₀)} by construction
# But q_n ~ (x₀/64)^n · 64^n · n^{σ₂₇}
# So q_n / 64^n ~ (x₀/64)^n · n^σ, while b_n ~ x₀^n n^σ
# These differ by 64^n. So b_n · (1/64)^n should match q_n/64^n.
# i.e., b_n should match q_n (unscaled).

# Actually: P2.7 recurrence coefficients encode q_n (not q_n/64^n).
# The recurrence: A(n) q_{n+1} = B(n) q_n - C(n-1)/A(n-1) q_{n-1} + D(n-2)/A(n-2) q_{n-2}
# and q_n ~ λ₀^n where λ₀ = x₀/64 ≈ 0.859.
# But b_n ~ x₀^n = (64λ₀)^n = 64^n λ₀^n.
# So b_n = 64^n · (something with growth λ₀^n).
# The "something" should match q_n / C_norm.

# Let me check: b_n / 64^n and compare with q_n
print(f"\n  b_n / 64^n compared to q_n (sign B):")
# Use sign B (+1)
psi2 = []
for k in range(N):
    psi2.append(sum(deltas[j] * roots[j]**(k+1) for j in range(3)))
c2 = [mpmath.mpf(0)] * N
c2[0] = mpmath.mpf(1)
for n in range(N-1):
    c2[n+1] = sum(psi2[k] * c2[n-k] for k in range(n+1)) / (n+1)
b2 = [mpmath.mpf(0)] * N
for n in range(N):
    b2[n] = sum(c2[m] * a[n-m] for m in range(n+1))

for n in range(min(12, N)):
    b_scaled = b2[n] / mpmath.mpf(64)**n
    if abs(q[n]) > 1e-50:
        ratio = b_scaled / q[n]
        print(f"    n={n:2d}: b_n/64^n = {float(b_scaled.real):+.6e}, q_n = {float(q[n]):+.6e}, ratio = {float(ratio.real):+.10e}")
