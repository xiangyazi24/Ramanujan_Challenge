#!/usr/bin/env sage
"""
MOP for mixed Beukers kernel — fixed version using direct linear solve.

Type II MOP conditions:
  Σ_{k=0}^n c_{n,k} · ψ_j(k+s) = 0   for s = 1,...,n_j   (j=1 for ζ(2), j=2 for ζ(3))
with n₁ + n₂ = n.

Set c_n = 1 (monic) → solve n×n system for c_0,...,c_{n-1}.
"""
import mpmath
mpmath.mp.dps = 200

def psi_tail(j, m):
    """ψ_j(m) = ζ(j) - H_{m-1}^{(j)} = Σ_{ℓ≥m} 1/ℓ^j"""
    if j == 2:
        return mpmath.zeta(2) - sum(mpmath.mpf(1)/mpmath.mpf(k)**2 for k in range(1, int(m)))
    elif j == 3:
        return mpmath.zeta(3) - sum(mpmath.mpf(1)/mpmath.mpf(k)**3 for k in range(1, int(m)))

def compute_mop(n, n1, n2):
    """Compute monic MOP Q_n. Returns coefficients and Q_n(1)."""
    n = int(n); n1 = int(n1); n2 = int(n2)
    assert n1 + n2 == n and n1 >= 0 and n2 >= 0

    # Build n×n system: A * [c_0,...,c_{n-1}]^T = -b
    # where b[row] = M[row, n] (the column for c_n = 1)
    A = mpmath.matrix(n, n)
    b = mpmath.matrix(n, 1)

    # First n1 rows: ψ₁ conditions (ζ(2) moments)
    for row in range(n1):
        s = row + 1
        for k in range(n):
            A[row, k] = psi_tail(2, k + s)
        b[row, 0] = -psi_tail(2, n + s)

    # Next n2 rows: ψ₂ conditions (ζ(3) moments)
    for row in range(n2):
        s = row + 1
        for k in range(n):
            A[n1 + row, k] = psi_tail(3, k + s)
        b[n1 + row, 0] = -psi_tail(3, n + s)

    # Solve
    c_vec = mpmath.lu_solve(A, b)
    c = [c_vec[k, 0] for k in range(n)] + [mpmath.mpf(1)]

    Qn1 = sum(c)
    return c, Qn1

# P2.7 reference
print("=== P2.7 reference ===")
print("Dominant root λ₀ ≈ 0.858808")
print("Subdominant |λ±| ≈ 0.0674")
print("Subdominant/dominant = 0.0674/0.8588 ≈ 0.0785")

# ============================================================
# Main computation
# ============================================================
max_n = 15

for step_name, n1_fn in [
    ("balanced (⌈n/2⌉, ⌊n/2⌋)", lambda n: (n+1)//2),
    ("ζ(2)-heavy (n-1, 1)", lambda n: n-1),
    ("ζ(3)-heavy (1, n-1)", lambda n: 1),
    ("2:1 ratio (⌈2n/3⌉, ⌊n/3⌋)", lambda n: (2*n+2)//3),
]:
    print("\n--- Step-line: %s ---" % step_name)
    results = []
    for n in range(2, max_n + 1):
        n1 = int(n1_fn(n))
        n2 = int(n) - n1
        if n2 < 1 or n1 < 1:
            continue
        try:
            c, Qn1 = compute_mop(n, n1, n2)
            results.append((n, Qn1, c))
            print("  n=%2d (n1=%d,n2=%d): Q_n(1) = %+.12e" % (n, n1, n2, float(Qn1)))
        except Exception as e:
            print("  n=%2d: FAILED: %s" % (n, str(e)[:80]))

    if len(results) >= 2:
        print("  Growth ratios Q_{n+1}(1)/Q_n(1):")
        for i in range(len(results) - 1):
            na, va, _ = results[i]
            nb, vb, _ = results[i+1]
            if abs(va) > 1e-50 and nb == na + 1:
                r = float(vb / va)
                print("    n=%d→%d: ratio = %+.10f  (target: ~0.859)" % (na, nb, r))

# ============================================================
# Error computation for the best step-line
# ============================================================
print("\n" + "="*70)
print("Error decay analysis")
print("="*70)

L = mpmath.zeta(2) + mpmath.zeta(3)
print("L = ζ(2)+ζ(3) = %.30f" % float(L))

for step_name, n1_fn in [
    ("balanced", lambda n: (n+1)//2),
    ("ζ(2)-heavy", lambda n: n-1),
]:
    print("\n--- %s ---" % step_name)
    errors = []
    for n in range(2, min(max_n + 1, 12)):
        n1 = int(n1_fn(n))
        n2 = int(n) - n1
        if n2 < 1 or n1 < 1:
            continue
        try:
            c, Qn1 = compute_mop(n, n1, n2)

            # p_n = Σ c_k (H_k^{(2)} + H_k^{(3)})
            p_n = mpmath.mpf(0)
            for k in range(int(n)+1):
                Hk2 = sum(mpmath.mpf(1)/mpmath.mpf(j)**2 for j in range(1, int(k)+1))
                Hk3 = sum(mpmath.mpf(1)/mpmath.mpf(j)**3 for j in range(1, int(k)+1))
                p_n += c[k] * (Hk2 + Hk3)

            e_n = Qn1 * L - p_n
            errors.append((n, e_n, Qn1))
            rel = float(abs(e_n / Qn1)) if abs(Qn1) > 1e-100 else 0
            print("  n=%2d: |e_n| = %.6e, |e_n/Q_n(1)| = %.6e" % (n, float(abs(e_n)), rel))
        except Exception as e:
            print("  n=%2d: %s" % (n, str(e)[:60]))

    if len(errors) >= 2:
        print("  Error decay ratios:")
        for i in range(len(errors) - 1):
            na, ea, _ = errors[i]
            nb, eb, _ = errors[i+1]
            if abs(ea) > 1e-200 and nb == na + 1:
                r = float(abs(eb) / abs(ea))
                print("    n=%d→%d: |e_{n+1}|/|e_n| = %.10f  (subdominant target: ~0.067)" % (na, nb, r))

# ============================================================
# Check: does Q_n(1) satisfy the P2.7 recurrence?
# ============================================================
print("\n" + "="*70)
print("P2.7 recurrence check for MOP Q_n(1)")
print("="*70)

def A_p27(n): return 1024 * (2*n+5)**4 * (2*n+7)**3 * (2*n+9)**3 * (946*n**2+6407*n+10860)
def B_p27(n): return 128 * (2*n+7)**3 * (2*n+9)**3 * (104060*n**6 + 1745370*n**5 + 12145238*n**4 + 44886481*n**3 + 92943995*n**2 + 102256019*n + 46709052)
def C_p27(n): return 16 * (n+3)**4 * (2*n+9)**3 * (3784*n**5 + 57792*n**4 + 351019*n**3 + 1059230*n**2 + 1587211*n + 944620)
def D_p27(n): return (n+3)**4 * (n+4)**6 * (946*n**2 + 4515*n + 5399)

for step_name, n1_fn in [("balanced", lambda n: (n+1)//2)]:
    print("\n--- %s ---" % step_name)
    vals = []
    for n in range(2, min(max_n + 1, 12)):
        n1 = int(n1_fn(n))
        n2 = int(n) - n1
        if n2 < 1 or n1 < 1:
            continue
        try:
            c, Qn1 = compute_mop(n, n1, n2)
            vals.append((n, Qn1))
        except:
            pass

    for i in range(3, len(vals)):
        n_curr = vals[i][0]
        if vals[i-1][0] != n_curr-1 or vals[i-2][0] != n_curr-2 or vals[i-3][0] != n_curr-3:
            continue
        u3 = vals[i][1]
        u2 = vals[i-1][1]
        u1 = vals[i-2][1]
        u0 = vals[i-3][1]
        nn = n_curr - 3  # shift so u0 = Q_{nn}, ..., u3 = Q_{nn+3}
        pred = (mpmath.mpf(B_p27(nn+2))/A_p27(nn+2) * u2
                - mpmath.mpf(C_p27(nn+1))/A_p27(nn+1) * u1
                + mpmath.mpf(D_p27(nn))/A_p27(nn) * u0)
        res = u3 - pred
        rel = float(abs(res / u3)) if abs(u3) > 1e-100 else float(abs(res))
        match = "MATCH!" if rel < 1e-10 else ""
        print("  n=%d: residual rel = %.6e %s" % (n_curr, rel, match))
