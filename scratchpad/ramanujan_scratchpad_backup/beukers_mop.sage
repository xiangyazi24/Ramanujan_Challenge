#!/usr/bin/env sage
"""
Multiple Orthogonal Polynomial (MOP) computation for the mixed Beukers kernel.

K(x,y) = (1 - ½ log(xy)) / (1-xy)

Markov functions:
  f₂(z) = Σ_{m≥1} ψ₁(m) z^{-m}   where ψ₁(m) = Σ_{ℓ≥m} 1/ℓ²  (tail of ζ(2))
  f₃(z) = Σ_{m≥1} ψ₂(m) z^{-m}   where ψ₂(m) = Σ_{ℓ≥m} 1/ℓ³  (tail of ζ(3))

Type II MOP conditions:
  Σ_{k=0}^n c_{n,k} · ψ_j(k+s) = 0   for s = 1,...,n_j   (j=1,2)

with n₁ + n₂ = n (step-line).

This determines Q_n(t) = Σ c_{n,k} t^k up to normalization.

Goal: compute Q_n(1) and check if it matches P2.7 growth rate ~0.859^n.
"""
import mpmath
mpmath.mp.dps = 300  # high precision for rational reconstruction

def psi_tail(j, m, prec=300):
    """ψ_j(m) = Σ_{ℓ≥m} 1/ℓ^j = ζ(j) - H_{m-1}^{(j)}"""
    old_dps = mpmath.mp.dps
    mpmath.mp.dps = prec
    if j == 2:
        val = mpmath.zeta(2) - sum(mpmath.mpf(1)/k**2 for k in range(1, m))
    elif j == 3:
        val = mpmath.zeta(3) - sum(mpmath.mpf(1)/k**3 for k in range(1, m))
    mpmath.mp.dps = old_dps
    return val

# ============================================================
# Build the MOP for several step-lines and degrees
# ============================================================
def compute_mop(n, n1, n2, prec=300):
    """
    Compute the type-II MOP Q_n for step-line (n1, n2).
    Returns the coefficient vector c = (c_0, ..., c_n) and Q_n(1).
    """
    assert n1 + n2 == n
    old_dps = mpmath.mp.dps
    mpmath.mp.dps = prec

    # Build the moment matrix M (n rows, n+1 columns)
    M = mpmath.matrix(int(n), int(n+1))

    # First n1 rows: ψ₁ conditions
    for row in range(int(n1)):
        s = int(row + 1)
        for k in range(int(n+1)):
            M[row, k] = psi_tail(2, int(k + s), prec)

    # Next n2 rows: ψ₂ conditions
    for row in range(int(n2)):
        s = int(row + 1)
        for k in range(int(n+1)):
            M[int(n1) + row, k] = psi_tail(3, int(k + s), prec)

    # Find the kernel (null space) of M
    # SVD to find the smallest singular value
    U, S, V = mpmath.svd(M)

    # The kernel vector is the last row of V (corresponding to smallest singular value)
    c = [V[n, j] for j in range(n+1)]

    # Normalize so c_n = 1 (monic)
    if abs(c[n]) > mpmath.mpf('1e-50'):
        norm = c[n]
        c = [c[j] / norm for j in range(n+1)]

    # Q_n(1)
    Qn1 = sum(c)

    mpmath.mp.dps = old_dps
    return c, Qn1

# ============================================================
# P2.7 recurrence coefficients
# ============================================================
def A_p27(n):
    return 1024 * (2*n+5)**4 * (2*n+7)**3 * (2*n+9)**3 * (946*n**2+6407*n+10860)
def B_p27(n):
    return 128 * (2*n+7)**3 * (2*n+9)**3 * (104060*n**6 + 1745370*n**5 +
        12145238*n**4 + 44886481*n**3 + 92943995*n**2 + 102256019*n + 46709052)
def C_p27(n):
    return 16 * (n+3)**4 * (2*n+9)**3 * (3784*n**5 + 57792*n**4 +
        351019*n**3 + 1059230*n**2 + 1587211*n + 944620)
def D_p27(n):
    return (n+3)**4 * (n+4)**6 * (946*n**2 + 4515*n + 5399)

# Compute P2.7 ratios for reference
print("=== P2.7 reference ratios ===")
mpmath.mp.dps = 50
q = [mpmath.mpf(0)] * 20
q[0] = mpmath.mpf('-215040420000')
q[1] = mpmath.mpf('-167282265043404') / mpmath.mpf('905')
q[2] = mpmath.mpf('-964185327658080') / mpmath.mpf('6071')
for n in range(2, 19):
    q[n+1] = (mpmath.mpf(B_p27(n))/A_p27(n) * q[n]
              - mpmath.mpf(C_p27(n-1))/A_p27(n-1) * q[n-1]
              + mpmath.mpf(D_p27(n-2))/A_p27(n-2) * q[n-2])

print("P2.7 dominant root λ₀ ≈ 0.858808")
for n in range(min(10, len(q)-1)):
    if q[n] != 0:
        r = q[n+1] / q[n]
        print("  q_{n+1}/q_n at n=%d: %.12f" % (n, float(r)))

# ============================================================
# Main computation: try different step-lines
# ============================================================
print("\n" + "="*70)
print("MOP with CORRECT polygamma moments")
print("="*70)

max_n = 12

for step_name, n1_fn, n2_fn in [
    ("balanced (⌈n/2⌉, ⌊n/2⌋)", lambda n: (n+1)//2, lambda n: n//2),
    ("ζ(2)-heavy (n-1, 1)", lambda n: n-1, lambda n: 1),
    ("ζ(3)-heavy (1, n-1)", lambda n: 1, lambda n: n-1),
    ("2:1 ratio (⌈2n/3⌉, ⌊n/3⌋)", lambda n: (2*n+2)//3, lambda n: n - (2*n+2)//3),
]:
    print("\n--- Step-line: %s ---" % step_name)
    Qn1_list = []
    for n in range(2, max_n + 1):
        n1 = n1_fn(n)
        n2 = n - n1
        if n2 < 1:
            continue
        try:
            c, Qn1 = compute_mop(n, n1, n2, prec=200)
            Qn1_list.append((n, Qn1))
            print("  n=%d (n1=%d,n2=%d): Q_n(1) = %.15e" % (n, n1, n2, float(Qn1)))
        except Exception as e:
            print("  n=%d: FAILED: %s" % (n, str(e)[:60]))

    # Compute ratios
    if len(Qn1_list) >= 2:
        print("  Ratios Q_{n+1}(1)/Q_n(1):")
        for i in range(len(Qn1_list) - 1):
            n1_idx, v1 = Qn1_list[i]
            n2_idx, v2 = Qn1_list[i+1]
            if v1 != 0 and n2_idx == n1_idx + 1:
                r = v2 / v1
                print("    n=%d→%d: ratio = %.12f  (P2.7 target: ~0.859)" % (n1_idx, n2_idx, float(r)))

# ============================================================
# Also compute the error: e_n = Q_n(1)·L - p_n
# where p_n = Σ c_{n,k} (H_k^{(2)} + H_k^{(3)})
# ============================================================
print("\n" + "="*70)
print("Error decay for balanced step-line")
print("="*70)

mpmath.mp.dps = 200
L = mpmath.zeta(2) + mpmath.zeta(3)

for n in range(2, min(max_n + 1, 10)):
    n1 = (n+1)//2
    n2 = n - n1
    try:
        c, Qn1 = compute_mop(n, n1, n2, prec=200)

        # Compute p_n = Σ c_k (H_k^{(2)} + H_k^{(3)})
        p_n = mpmath.mpf(0)
        for k in range(n+1):
            Hk2 = sum(mpmath.mpf(1)/j**2 for j in range(1, k+1))  # H_k^{(2)}
            Hk3 = sum(mpmath.mpf(1)/j**3 for j in range(1, k+1))  # H_k^{(3)}
            p_n += c[k] * (Hk2 + Hk3)

        e_n = Qn1 * L - p_n
        print("  n=%d: |e_n| = %.6e, |e_n/Q_n(1)| = %.6e" % (n, float(abs(e_n)), float(abs(e_n/Qn1)) if Qn1 != 0 else 0))
    except Exception as e:
        print("  n=%d: %s" % (n, str(e)[:60]))

# Error decay ratios
print("\nError decay ratios |e_{n+1}|/|e_n|:")
errors = []
for n in range(2, min(max_n + 1, 10)):
    n1 = (n+1)//2
    n2 = n - n1
    try:
        c, Qn1 = compute_mop(n, n1, n2, prec=200)
        p_n = mpmath.mpf(0)
        for k in range(n+1):
            Hk2 = sum(mpmath.mpf(1)/j**2 for j in range(1, k+1))
            Hk3 = sum(mpmath.mpf(1)/j**3 for j in range(1, k+1))
            p_n += c[k] * (Hk2 + Hk3)
        e_n = Qn1 * L - p_n
        errors.append((n, e_n))
    except:
        pass

for i in range(len(errors)-1):
    n1_idx, e1 = errors[i]
    n2_idx, e2 = errors[i+1]
    if abs(e1) > 0 and n2_idx == n1_idx + 1:
        r = abs(e2) / abs(e1)
        print("  n=%d→%d: |e_{n+1}|/|e_n| = %.12f  (P2.7 subdominant: ~0.067)" % (n1_idx, n2_idx, float(r)))
