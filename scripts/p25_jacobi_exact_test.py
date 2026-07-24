#!/usr/bin/env python3
"""P2.5: Exact residual test — does Q_n^J satisfy the CMF recurrence?

Q4877 corrected the Jacobi convention:
  J_n(X; ε) = P_n^{(ε-1/2, 0)}(1-2X)
At ε=0: B_n = J_n(-1; 0) = P_n^{(-1/2, 0)}(3)
  Q_n^J = (4n+1)/2 · B_n²

Also test the original (wrong?) convention:
  Q_n^{old} = (4n+1)/2 · [P_n(3)]² = (4n+1)/2 · D_n²

And test plain D_n² for comparison.
"""
from fractions import Fraction
import sys

def eval_c(coeffs, n):
    """Evaluate polynomial with given coefficients at integer n."""
    val = Fraction(0)
    npow = Fraction(1)
    for c in coeffs:
        val += Fraction(c) * npow
        npow *= Fraction(n)
    return val

c0_coeffs = [-170972650800, -826494925500, -1792449886332, -2317972607944, -2000297648936, -1219354055500, -541255279788, -177419351856, -43002662976, -7620091136, -960400960, -81589760, -4190208, -98304]
c1_coeffs = [8781630505200, 38850314624124, 78557994908508, 96136040496551, 79442239242197, 46814452218572, 20241514501104, 6502490145168, 1552168938336, 271943188864, 33995217088, 2871763456, 146952192, 3440640]
c2_coeffs = [-21132458248680, -87529225645944, -165451256319618, -189073879129764, -145809619841418, -80164318460172, -32338316008004, -9694892892592, -2160716677664, -353683596544, -41340724928, -3268370944, -156684288, -3440640]
c3_coeffs = [587448626688, 2442715444224, 4635428285664, 5317694979920, 4116150568664, 2270943978716, 919036676572, 276298241680, 61721801728, 10120470656, 1184128064, 93632000, 4485120, 98304]

# --- Compute Jacobi P_n^{(-1/2, 0)}(x) via three-term recurrence ---
def jacobi_neg_half_zero(N, x):
    """Compute P_n^{(-1/2, 0)}(x) for n = 0, 1, ..., N using exact arithmetic.
    Uses the standard Jacobi three-term recurrence with α=-1/2, β=0.
    """
    x = Fraction(x)
    vals = [Fraction(1)]  # P_0 = 1
    if N == 0:
        return vals

    # P_1^{(-1/2, 0)}(x) = (α+β+2)/2 · x + (α-β)/2 = (3/4)x - 1/4
    vals.append(Fraction(3, 4) * x - Fraction(1, 4))

    alpha = Fraction(-1, 2)
    beta = Fraction(0)

    for n in range(1, N):
        # Standard Jacobi recurrence:
        # 2(n+1)(n+α+β+1)(2n+α+β) P_{n+1} =
        #   (2n+α+β+1)[(2n+α+β)(2n+α+β+2)x + α²-β²] P_n
        #   - 2(n+α)(n+β)(2n+α+β+2) P_{n-1}
        ab = alpha + beta
        A1 = 2 * (n + 1) * (n + ab + 1) * (2*n + ab)
        B1 = (2*n + ab + 1) * ((2*n + ab) * (2*n + ab + 2) * x + alpha**2 - beta**2)
        C1 = 2 * (n + alpha) * (n + beta) * (2*n + ab + 2)

        if A1 == 0:
            print(f"WARNING: A1=0 at n={n}")
            break

        p_next = (B1 * vals[n] - C1 * vals[n-1]) / A1
        vals.append(p_next)

    return vals

# --- Compute Legendre P_n(x) (= Delannoy when x=3) ---
def legendre_at(N, x):
    """P_n(x) for n = 0, ..., N via recurrence."""
    x = Fraction(x)
    vals = [Fraction(1)]
    if N == 0:
        return vals
    vals.append(x)
    for n in range(1, N):
        # (n+1)P_{n+1} = (2n+1)x P_n - n P_{n-1}
        p_next = ((2*n + 1) * x * vals[n] - n * vals[n-1]) / (n + 1)
        vals.append(p_next)
    return vals

NMAX = 25

print("=== Computing sequences ===", flush=True)

# Jacobi P_n^{(-1/2, 0)}(3)
jac_vals = jacobi_neg_half_zero(NMAX + 3, 3)
print(f"  P_n^{{(-1/2,0)}}(3) for n=0..5: {[float(v) for v in jac_vals[:6]]}")

# Legendre P_n(3) = D_n (Delannoy)
leg_vals = legendre_at(NMAX + 3, 3)
print(f"  P_n(3) (Delannoy) for n=0..5: {[int(v) for v in leg_vals[:6]]}")

# Build test sequences
Q_jac = []  # (4n+1)/2 * [P_n^{(-1/2,0)}(3)]^2
Q_del = []  # D_n^2
Q_kdel = []  # (4n+1)/2 * D_n^2

for n in range(NMAX + 4):
    kappa = Fraction(4*n + 1, 2)
    Q_jac.append(kappa * jac_vals[n]**2)
    Q_del.append(leg_vals[n]**2)
    Q_kdel.append(kappa * leg_vals[n]**2)

print(f"\n  Q_jac[0..5] = {[float(v) for v in Q_jac[:6]]}")
print(f"  Q_del[0..5] = {[int(v) for v in Q_del[:6]]}")
print(f"  Q_kdel[0..5] = {[float(v) for v in Q_kdel[:6]]}")

# --- Test each sequence against the CMF recurrence ---
def test_recurrence(name, seq, max_n=None):
    if max_n is None:
        max_n = len(seq) - 4
    print(f"\n=== Testing {name} against CMF recurrence ===", flush=True)
    all_zero = True
    for n in range(min(max_n, len(seq) - 3)):
        res = (eval_c(c0_coeffs, n) * seq[n] +
               eval_c(c1_coeffs, n) * seq[n+1] +
               eval_c(c2_coeffs, n) * seq[n+2] +
               eval_c(c3_coeffs, n) * seq[n+3])
        if res != 0:
            all_zero = False
            if n < 5 or n == max_n - 1:
                print(f"  n={n}: residual = {float(res):.6e} (NONZERO)")
        else:
            if n < 5:
                print(f"  n={n}: residual = 0 (EXACT)")

    if all_zero:
        print(f"  *** {name} SATISFIES the CMF recurrence exactly! ***")
    else:
        print(f"  {name} does NOT satisfy the CMF recurrence.")
    return all_zero

test_recurrence("Q_jac = (4n+1)/2 * [P_n^{(-1/2,0)}(3)]^2", Q_jac, 20)
test_recurrence("Q_del = D_n^2 = [P_n(3)]^2", Q_del, 20)
test_recurrence("Q_kdel = (4n+1)/2 * D_n^2", Q_kdel, 20)

# --- Also find what recurrence D_n^2 and Q_jac satisfy ---
print("\n=== Finding the recurrence for D_n^2 (Sym^2 Delannoy) ===", flush=True)
# D_n satisfies (n+1)^2 D_{n+1} = (6n^2+6n+2) D_n + n^2 D_{n-1}
# Equivalently: (n+1)^2 D_{n+1} - (6n^2+6n+2) D_n - n^2 D_{n-1} = 0
# The Sym^2 recurrence for D_n^2 has order 3.
# Ansatz: a(n) D_n^2 + b(n) D_{n+1}^2 + c(n) D_{n+2}^2 + d(n) D_{n+3}^2 = 0

# Use the Delannoy recurrence to express D_{n+3}^2 etc. in terms of D_n, D_{n+1}
# Then find the minimal annihilating recurrence

# More directly: search for order-3 recurrence with polynomial coeffs for D_n^2
# Try degree d = 0, 1, 2, 3, ...

for deg in range(6):
    # 4(deg+1)-1 parameters
    nparams = 4*(deg+1) - 1
    ntest = min(len(Q_del) - 3, 30)
    if ntest < nparams + 3:
        continue

    # Build system
    A_mat = []
    b_vec = []
    for n in range(ntest):
        row = []
        for j in range(4):
            for k in range(deg+1):
                if j == 3 and k == deg:
                    continue
                row.append(Fraction(n)**k * Q_del[n+j])
        A_mat.append(row)
        b_vec.append(-Fraction(n)**deg * Q_del[n+3])

    # Solve first nparams equations
    import numpy as np
    A_np = np.array([[float(x) for x in row[:nparams]] for row in A_mat[:nparams]])
    b_np = np.array([float(x) for x in b_vec[:nparams]])

    try:
        sol = np.linalg.solve(A_np, b_np)
    except:
        continue

    # Check residual on held-out
    max_res = 0.0
    for n in range(nparams, min(ntest, nparams + 5)):
        res = b_vec[n]
        for i in range(nparams):
            res -= Fraction(A_mat[n][i]) * Fraction(sol[i]).limit_denominator(10**15)
        max_res = max(max_res, abs(float(res)))

    print(f"  degree {deg}: holdout residual = {max_res:.3e}", flush=True)

    if max_res < 1e-10:
        print(f"  *** D_n^2 has order-3 recurrence with polynomial degree {deg} ***")
        # Extract Poincaré polynomial
        idx = 0
        lc = []
        for j in range(4):
            for k in range(deg+1):
                if j == 3 and k == deg:
                    lc.append(1.0)
                    break
                if k == deg:
                    lc.append(sol[idx])
                idx += 1 if not (j == 3 and k == deg) else 0
        # Actually let me just get the leading coefficients
        lc = []
        idx = 0
        for j in range(4):
            lc_j = None
            for k in range(deg+1):
                if j == 3 and k == deg:
                    lc_j = 1.0
                else:
                    if k == deg:
                        lc_j = sol[idx]
                    idx += 1
            if lc_j is None:
                lc_j = sol[idx - 1]  # last one for this j
            lc.append(lc_j)
        print(f"  Leading coefficients: {lc}")
        break

print("\nDone.", flush=True)
