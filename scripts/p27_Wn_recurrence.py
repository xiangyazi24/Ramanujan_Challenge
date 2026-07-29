#!/usr/bin/env python3
"""P2.7: Find the recurrence for W_n = b_{2n}/256^n (level-11 binomial transform).
Then compare with the P2.7 recurrence to identify the gauge.
"""
import mpmath as mp
mp.mp.dps = 150

from fractions import Fraction
import math

# Cooper's level-11 sequence
def compute_T(N):
    T = [Fraction(1), Fraction(4), Fraction(28)]
    for k in range(2, N):
        num = (2*(2*k+1)*(5*k**2+5*k+2)*T[k]
               - 8*k*(7*k**2+1)*T[k-1]
               + 22*k*(2*k-1)*(k-1)*T[k-2])
        T.append(num / Fraction((k+1)**3))
    return T

T = compute_T(250)

# Binomial transform: b_m = Σ C(m,j)(-2)^{m-j} T_j
def compute_bm(m, T):
    val = Fraction(0)
    for j in range(min(m+1, len(T))):
        val += math.comb(m, j) * Fraction(-2)**(m-j) * T[j]
    return val

# W_n = b_{2n} / 256^n
print("Computing W_n...", flush=True)
W_frac = []
W_mp = []
for n in range(105):
    b2n = compute_bm(2*n, T)
    Wn = b2n / Fraction(256)**n
    W_frac.append(Wn)
    W_mp.append(mp.mpf(Wn.numerator) / mp.mpf(Wn.denominator))
    if n <= 5 or n % 20 == 0:
        print(f"  W[{n}] = {mp.nstr(W_mp[-1], 15)}", flush=True)

# Search for 4-term recurrence: P₃(n)W_{n+3} + P₂(n)W_{n+2} + P₁(n)W_{n+1} + P₀(n)W_n = 0
# where P_j(n) are polynomials of degree d
print(f"\n=== Searching for 4-term recurrence of W_n ===", flush=True)

nvals = len(W_mp)
for deg in range(1, 20):
    nparams = 4 * (deg + 1) - 1
    neq = nvals - 3
    if neq < nparams + 5:
        print(f"  degree {deg}: not enough equations ({neq} < {nparams+5})")
        continue

    # Build system
    A_rows = []
    b_vec = []
    for n in range(neq):
        row = []
        for j in range(4):
            for k in range(deg + 1):
                if j == 3 and k == deg:
                    continue
                row.append(mp.mpf(n)**k * W_mp[n + j])
        A_rows.append(row)
        b_vec.append(-mp.mpf(n)**deg * W_mp[n + 3])

    # Solve first nparams equations
    A_mat = mp.matrix([r[:nparams] for r in A_rows[:nparams]])
    b_mat = mp.matrix([b_vec[i] for i in range(nparams)])

    try:
        sol = mp.lu_solve(A_mat, b_mat)
    except:
        print(f"  degree {deg}: singular", flush=True)
        continue

    # Verify on held-out
    max_res = mp.mpf(0)
    for n in range(nparams, min(neq, nparams + 8)):
        res = b_vec[n]
        for i in range(nparams):
            res -= A_rows[n][i] * sol[i]
        max_res = max(max_res, abs(res))

    print(f"  degree {deg}: holdout residual = {mp.nstr(max_res, 6)}", flush=True)

    if max_res < mp.mpf(10)**(-80):
        print(f"  *** FOUND: degree {deg} ***")

        # Extract Poincaré polynomial (leading coefficients of each P_j)
        idx = 0
        coeffs_by_j = {j: [] for j in range(4)}
        for j in range(4):
            for k in range(deg + 1):
                if j == 3 and k == deg:
                    coeffs_by_j[j].append(mp.mpf(1))
                else:
                    coeffs_by_j[j].append(sol[idx])
                    idx += 1

        lc = [coeffs_by_j[j][deg] for j in range(4)]
        print(f"  Leading coefficients: {[mp.nstr(c, 15) for c in lc]}")

        char_poly = [lc[j]/lc[3] for j in range(4)]
        print(f"  Poincaré: λ³ + ({mp.nstr(char_poly[2],15)})λ² + ({mp.nstr(char_poly[1],15)})λ + ({mp.nstr(char_poly[0],15)}) = 0")

        try:
            roots = mp.polyroots([char_poly[3], char_poly[2], char_poly[1], char_poly[0]])
            print(f"  Roots: {[mp.nstr(r, 12) for r in roots]}")
            print(f"  |Roots|: {[mp.nstr(abs(r), 12) for r in roots]}")
        except:
            print("  (root-finding failed)")

        # Also extract ALL polynomial coefficients for comparison with P2.7
        print(f"\n  Full recurrence coefficients (ascending n-powers):")
        for j in range(4):
            cs = coeffs_by_j[j]
            # Try to recognize as integers/rationals
            print(f"    P_{j}(n): {[mp.nstr(c, 12) for c in cs]}")

        # Compare with P2.7 recurrence Poincaré: 4μ³-220μ²+8μ-1
        # Normalized: μ³-55μ²+2μ-1/4
        print(f"\n  P2.7 Poincaré for comparison: μ³-55μ²+2μ-1/4")
        print(f"  W_n Poincaré (normalized): λ³+({mp.nstr(char_poly[2],15)})λ²+({mp.nstr(char_poly[1],15)})λ+({mp.nstr(char_poly[0],15)})")

        break

# Also compute W_n ratios
print(f"\n=== W_n ratios ===")
for n in [10, 20, 50, 80, 100]:
    if n < len(W_mp) and n > 0:
        r = W_mp[n] / W_mp[n-1]
        print(f"  W[{n}]/W[{n-1}] = {mp.nstr(r, 15)}")

print(f"\n  Expected dominant ratio: ((t₀-2)²/256) where t₀ ≈ 16.8 (root of H₁₁)")
# Compute t₀ numerically
coeffs_H11 = [mp.mpf(-44), mp.mpf(56), mp.mpf(-20), mp.mpf(1)]
t_roots = mp.polyroots(coeffs_H11[::-1])
print(f"  H₁₁ roots: {[mp.nstr(r, 12) for r in t_roots]}")
for r in t_roots:
    if mp.im(r) == 0:
        mu = (r - 2)**2 / 4
        print(f"  Real root t₀ = {mp.nstr(r, 12)} → μ₀ = (t₀-2)²/4 = {mp.nstr(mu, 12)}")
        c0 = mu / 64
        print(f"  c₀ = μ₀/64 = {mp.nstr(c0, 12)}")
        # Dominant W_n ratio should be c₀ = (t₀-2)²/256
        wrat = (r - 2)**2 / 256
        print(f"  Expected W ratio → {mp.nstr(wrat, 12)}")

print("\nDone.")
