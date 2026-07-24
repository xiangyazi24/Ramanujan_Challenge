#!/usr/bin/env sage
"""
HP Type II for Li₂(1/z) and Li₃(1/z), evaluated at z=1.

Markov functions:
  f₂(z) = Li₂(1/z) = Σ_{k≥1} z^{-k}/k²
  f₃(z) = Li₃(1/z) = Σ_{k≥1} z^{-k}/k³

HP Type II: find polynomial A_n(z) of degree n such that:
  A_n(z) f₂(z) - B_n^{(2)}(z) = R_n^{(2)}(z) = O(z^{-(n₂+1)})
  A_n(z) f₃(z) - B_n^{(3)}(z) = R_n^{(3)}(z) = O(z^{-(n₃+1)})

Orthogonality conditions:
  Σ_{j=0}^n a_j / (j+m)^s = 0  for m=1,...,n_s; s=2,3

At z=1:
  q_n = A_n(1), p_n = B_n^{(2)}(1) + B_n^{(3)}(1)
  q_n(ζ(2)+ζ(3)) - p_n = R_n^{(2)}(1) + R_n^{(3)}(1) → 0
"""
from ore_algebra import *
from ore_algebra import guess

R_poly.<x> = PolynomialRing(QQ)
R_ore.<n> = PolynomialRing(QQ)
A_ore.<Sn> = OreAlgebra(R_ore)

N_max = 15

def compute_HP_z1(nn, n2, n3):
    """Compute A_n, evaluate at z=1, compute B and R at z=1."""
    assert n2 + n3 == nn

    # Build orthogonality matrix
    M = matrix(QQ, nn, nn+1)
    for row_idx in range(n2):
        m = row_idx + 1
        for j in range(nn+1):
            M[row_idx, j] = QQ(1) / QQ(j+m)^2
    for row_idx in range(n3):
        m = row_idx + 1
        for j in range(nn+1):
            M[n2 + row_idx, j] = QQ(1) / QQ(j+m)^3

    ker = M.right_kernel()
    if ker.dimension() != 1:
        return None

    v = ker.basis()[0]
    # Clear denominators to get integer coefficients
    denoms = [v[j].denominator() for j in range(nn+1)]
    L = lcm(denoms)
    a = [ZZ(v[j] * L) for j in range(nn+1)]

    # A_n(z) = Σ a_j z^j
    # A_n(1) = Σ a_j
    q_val = sum(a)

    # B_n^{(s)}(z) = polynomial part of A_n(z) * f_s(z)
    # At z=1: B_n^{(s)}(1) = polynomial part of A_n(z) * f_s(z) evaluated at z=1
    # f_s(z) = Σ_{k≥1} z^{-k}/k^s
    # A_n(z) f_s(z) = Σ_j a_j z^j · Σ_k z^{-k}/k^s
    # Coefficient of z^m (m = n-1, n-2, ..., 0) is part of B_n^{(s)}
    # Coefficient of z^{-m} (m ≥ 1) is part of R_n^{(s)}

    # B_n^{(s)}(z) = Σ_{m=0}^{n-1} [Σ_{j=max(0,m)}^{n} a_j/(j-m+1)^s ...] z^m
    # Actually: the product A·f has terms a_j z^j · z^{-k}/k^s = a_j z^{j-k}/k^s
    # So the coefficient of z^m in A·f is Σ_{j-k=m, k≥1} a_j/k^s = Σ_{k≥1, j=m+k≤n} a_{m+k}/k^s

    # For m ≥ 0 (polynomial part):
    # [z^m] A·f_s = Σ_{k=1}^{n-m} a_{m+k}/k^s
    # B_n^{(s)}(1) = Σ_{m=0}^{n-1} Σ_{k=1}^{n-m} a_{m+k}/k^s

    b2_at_1 = QQ(0)
    for m in range(nn):  # m=0,...,n-1
        for k in range(1, nn-m+1):
            if m+k <= nn:
                b2_at_1 += QQ(a[m+k]) / QQ(k)^2

    b3_at_1 = QQ(0)
    for m in range(nn):
        for k in range(1, nn-m+1):
            if m+k <= nn:
                b3_at_1 += QQ(a[m+k]) / QQ(k)^3

    p_val = b2_at_1 + b3_at_1

    # R_n^{(s)}(1) = Σ_{m≥1} [z^{-m}] A·f_s
    # [z^{-m}] A·f_s = Σ_{k≥m, j=k-m≥0, j≤n} a_{k-m}/k^s = Σ_{j=0}^{n} a_j/(j+m)^s
    # R_n^{(s)}(1) = Σ_{m≥1} Σ_{j=0}^n a_j/(j+m)^s
    # = Σ_j a_j [Σ_{m≥1} 1/(j+m)^s]
    # = Σ_j a_j [ζ(s) - H_j^{(s)}]  where H_j^{(s)} = Σ_{k=1}^j 1/k^s

    # But we want truncated version for computation. Let's compute directly:
    # R_n^{(s)}(1) = Σ_j a_j [ζ(s) - Σ_{m=1}^{n_s} 1/(j+m)^s - ...]
    # Actually, by the orthogonality conditions:
    # Σ_j a_j/(j+m)^s = 0 for m=1,...,n_s
    # So R_n^{(s)}(1) = Σ_{m≥n_s+1} Σ_j a_j/(j+m)^s

    # For exact computation: R^{(s)}(1) = q_val · ζ(s) - p_val_s (just the single-measure version)
    # where p_val_s = B_n^{(s)}(1) + Σ_{m=1}^{n_s} Σ_j a_j/(j+m)^s
    # But the second term is 0 by orthogonality!
    # So R^{(s)}(1) = q_val · ζ(s) - B_n^{(s)}(1)

    return {
        'a': a,
        'q': q_val,  # A_n(1)
        'b2': b2_at_1,
        'b3': b3_at_1,
        'p': p_val,  # B_n^{(2)}(1) + B_n^{(3)}(1)
    }

print("="*70)
print("HP Type II approximation to ζ(2)+ζ(3)")
print("="*70)

# Try different step-line conventions
for label, get_n23 in [
    ('balanced (n2>=n3)', lambda n: ((n+1)//2, n-(n+1)//2)),
    ('balanced (n3>=n2)', lambda n: (n-(n+1)//2, (n+1)//2)),
    ('equal', lambda n: (n//2, n-n//2)),
]:
    print(f"\n--- {label} ---")
    q_vals = []
    p_vals = []
    for nn in range(1, N_max+1):
        n2, n3 = get_n23(nn)
        res = compute_HP_z1(nn, n2, n3)
        if res is None:
            print(f"  n={nn}: degenerate")
            q_vals.append(None)
            p_vals.append(None)
            continue
        q_vals.append(res['q'])
        p_vals.append(res['p'])
        approx = float(res['p']) / float(res['q']) if res['q'] != 0 else float('inf')
        target = float(pi^2/6 + zeta(3))
        err = approx - target if abs(approx) < 100 else float('inf')
        print(f"  n={nn}: n2={n2}, n3={n3}, q={res['q']}, p/q={approx:.15f}, err={err:.6e}")

    # Guess recurrence for q_n = A_n(1)
    clean = [v for v in q_vals if v is not None]
    if len(clean) >= 8:
        print(f"\n  Guessing recurrence for q_n = A_n(1):")
        try:
            L = guess(clean, A_ore)
            print(f"    Order: {L.order()}, Degree: {max(c.degree() for c in L.list())}")
        except Exception as e:
            print(f"    Could not guess: {e}")

# Also try: the COMBINED approach where we require
# A_n(z)(f₂(z) + f₃(z)) - B_n(z) = O(z^{-(n+1)})
# This is standard single-function Padé for f₂+f₃ = Li₂(1/z) + Li₃(1/z)
print(f"\n{'='*70}")
print("Single-function Padé for f₂+f₃ = Li₂(1/z)+Li₃(1/z)")
print("="*70)

q_combined = []
for nn in range(1, N_max+1):
    M = matrix(QQ, nn, nn+1)
    for m in range(nn):
        for j in range(nn+1):
            M[m, j] = QQ(1)/QQ(j+m+1)^2 + QQ(1)/QQ(j+m+1)^3

    ker = M.right_kernel()
    if ker.dimension() != 1:
        q_combined.append(None)
        continue

    v = ker.basis()[0]
    denoms = [v[j].denominator() for j in range(nn+1)]
    L = lcm(denoms)
    a = [ZZ(v[j] * L) for j in range(nn+1)]
    q_val = sum(a)
    q_combined.append(q_val)

    # Compute p
    b_at_1 = QQ(0)
    for m in range(nn):
        for k in range(1, nn-m+1):
            if m+k <= nn:
                b_at_1 += QQ(a[m+k]) * (QQ(1)/QQ(k)^2 + QQ(1)/QQ(k)^3)
    p_val = b_at_1

    approx = float(p_val)/float(q_val) if q_val != 0 else float('inf')
    target = float(pi^2/6 + zeta(3))
    err = approx - target
    print(f"  n={nn}: q={q_val}, p/q={approx:.15f}, err={err:.6e}")

# Guess recurrence for the combined Padé
clean_comb = [v for v in q_combined if v is not None]
if len(clean_comb) >= 8:
    print(f"\n  Guessing recurrence for combined q_n:")
    try:
        L = guess(clean_comb, A_ore)
        print(f"    Order: {L.order()}, Degree: {max(c.degree() for c in L.list())}")
    except Exception as e:
        print(f"    Could not guess: {e}")
