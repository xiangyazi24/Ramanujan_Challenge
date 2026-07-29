#!/usr/bin/env python3
"""P2.5: Test Jacobi-Padé family R_n^J against the CMF recurrence.

From Q4886: The correct carrier phase is Φ(u,v) = φ(u)φ(v) where
φ(u) = (1+u)(u+2)/u = 3 + u + 2/u, with critical values (3±2√2)²
matching the CMF Poincaré roots.

The Jacobi polynomial family J_n(X;ε) = P_n^{(ε-1/2,0)}(1-2X)
gives an exact Catalan moment identity:
  ∫₀¹ (-log t)/(1+t²) R_n^J(t²) dt = G·Q_n^J - P_n^J

Test: does R_n^J satisfy the CMF recurrence?
If not, search for a contiguous band combination.
"""
from fractions import Fraction
import math

def rising_factorial_frac(a, k):
    result = Fraction(1)
    for i in range(k):
        result *= (a + i)
    return result

def odd_harmonic(m):
    return sum(Fraction(1, 2*r+1) for r in range(m))

def J_and_dJ(N):
    """Compute J_N(X;0) and dJ_N(X;0) as coefficient lists."""
    pref = rising_factorial_frac(Fraction(1,2), N) / Fraction(math.factorial(N))
    J_coeffs = [Fraction(0)] * (N + 1)
    dJ_coeffs = [Fraction(0)] * (N + 1)

    for k in range(N + 1):
        c = (Fraction(-1)**k * Fraction(math.comb(N, k)) * pref
             * rising_factorial_frac(Fraction(N) + Fraction(1,2), k)
             / rising_factorial_frac(Fraction(1,2), k))

        dlog = (2 * odd_harmonic(N)
                + 2 * sum(Fraction(1, 2*N+2*r+1) for r in range(k))
                - 2 * odd_harmonic(k))

        J_coeffs[k] = c
        dJ_coeffs[k] = c * dlog

    return J_coeffs, dJ_coeffs

def poly_eval(coeffs, x):
    result = Fraction(0)
    xpow = Fraction(1)
    for c in coeffs:
        result += c * xpow
        xpow *= x
    return result

def poly_deriv(coeffs):
    return [Fraction(k+1) * coeffs[k+1] for k in range(len(coeffs)-1)] if len(coeffs) > 1 else [Fraction(0)]

def poly_mul_scalar(coeffs, s):
    return [c * s for c in coeffs]

def poly_add(a, b):
    n = max(len(a), len(b))
    result = [Fraction(0)] * n
    for i in range(len(a)):
        result[i] += a[i]
    for i in range(len(b)):
        result[i] += b[i]
    return result

def poly_sub(a, b):
    return poly_add(a, poly_mul_scalar(b, Fraction(-1)))

def poly_mul(a, b):
    if not a or not b:
        return [Fraction(0)]
    n = len(a) + len(b) - 1
    result = [Fraction(0)] * n
    for i in range(len(a)):
        for j in range(len(b)):
            result[i+j] += a[i] * b[j]
    return result

def poly_divmod(num, den):
    """Polynomial division: num = q*den + r."""
    num = list(num)
    dend = len(den) - 1
    q = [Fraction(0)] * max(0, len(num) - dend)
    for i in range(len(q) - 1, -1, -1):
        q[i] = num[i + dend] / den[dend]
        for j in range(dend + 1):
            num[i + j] -= q[i] * den[j]
    r = num[:dend] if dend > 0 else []
    return q, r

def R_jacobi(N):
    """Compute R_N^J(X) as coefficient list."""
    J, dJ = J_and_dJ(N)
    B = poly_eval(J, Fraction(-1))
    dB = poly_eval(dJ, Fraction(-1))

    # numerator = B*dJ - dB*J
    num = poly_sub(poly_mul_scalar(dJ, B), poly_mul_scalar(J, dB))

    # Divide by (X+1)
    C, rem = poly_divmod(num, [Fraction(1), Fraction(1)])
    assert all(r == 0 for r in rem), f"Remainder nonzero at N={N}: {rem}"

    kappa = Fraction(4*N+1, 2)

    # R = kappa * [B*J - (X+1)/2 * (2X*C' + C)]
    BJ = poly_mul_scalar(J, B)

    Cprime = poly_deriv(C)
    # 2X*C' = shift Cprime by 1 and multiply by 2
    twoXCp = [Fraction(0)] + poly_mul_scalar(Cprime, Fraction(2))
    inner = poly_add(twoXCp, C)  # 2X*C' + C

    # (X+1)/2 * inner
    half_inner = poly_mul_scalar(inner, Fraction(1,2))
    xp1_half = poly_mul(half_inner, [Fraction(1), Fraction(1)])  # multiply by (X+1)

    R = poly_mul_scalar(poly_sub(BJ, xp1_half), kappa)
    return R

def catalan_monomial_pair(k):
    """∫₀¹ (-log t)/(1+t²) t^{2k} dt = q·G - p"""
    q = Fraction((-1)**k)
    partial = sum(Fraction((-1)**j, (2*j+1)**2) for j in range(k))
    p = q * partial
    return q, p

def moment_pair(R_coeffs):
    """Compute (q, p) such that ∫₀¹ (-log t)/(1+t²) R(t²) dt = q·G - p."""
    q = Fraction(0)
    p = Fraction(0)
    for k, a in enumerate(R_coeffs):
        if a == 0:
            continue
        qk, pk = catalan_monomial_pair(k)
        q += a * qk
        p += a * pk
    return q, p

# CMF matrix and normalization
def M_entries(n):
    n = Fraction(n)
    m11 = (-2*n-5)*(n+3)**2*(136*n**4+1424*n**3+5548*n**2+9551*n+6141)
    m12 = 384*n**6+6384*n**5+44168*n**4+162698*n**3+336377*n**2+369933*n+169011
    m13 = -(480*n**4+4980*n**3+19210*n**2+32690*n+20730)
    m21 = (n+2)**2*(n+3)**2*(4*n+10)*(48*n**3+386*n**2+1017*n+879)
    m22 = (n+2)**2*(-272*n**5-3848*n**4-21732*n**3-61184*n**2-85761*n-47808)
    m23 = (n+2)**2*(320*n**3+2540*n**2+6610*n+5640)
    m31 = (-4*n-10)*(n+2)**2*(n+3)**2*(32*n**4+302*n**3+1037*n**2+1530*n+813)
    m32 = (n+2)**2*(192*n**6+2984*n**5+19116*n**4+64452*n**3+120256*n**2+117279*n+46476)
    m33 = (n+2)**2*(-16*n**5-408*n**4-2912*n**3-8884*n**2-12254*n-6240)
    return [[m11,m12,m13],[m21,m22,m23],[m31,m32,m33]]

def delta_H(n):
    n = Fraction(n)
    return Fraction(-2)*(n+2)**2*(n+3)**2*(2*n+5)*(2*n+7)**2

def MH_at(n):
    M = M_entries(n)
    d = delta_H(n)
    return [[M[i][j] / d for j in range(3)] for i in range(3)]

def mat_vec(M, v):
    return [sum(M[i][j]*v[j] for j in range(3)) for i in range(3)]

# Compute CMF normalized initial pairs
print("Computing CMF pairs (Q̂_n, P̂_n)...", flush=True)
p_row = [Fraction(30921), Fraction(-32972), Fraction(8240)]
q_row = [Fraction(33750), Fraction(-36000), Fraction(9000)]

NMAX = 40

cmf_pairs = []  # (Q̂_n[col0], P̂_n[col0])
for N in range(NMAX + 4):
    cmf_pairs.append((q_row[0], p_row[0]))
    if N < NMAX + 3:
        MH = MH_at(N)
        p_row = mat_vec([[MH[j][i] for j in range(3)] for i in range(3)], p_row)
        q_row = mat_vec([[MH[j][i] for j in range(3)] for i in range(3)], q_row)

# Wait — the CMF iteration is: row · M_H(n), i.e., left multiplication
# p_{n+1} = p_n · M_H(n), so p_{n+1}[j] = Σ_i p_n[i] · M_H(n)[i][j]
# Let me redo this correctly

p_row = [Fraction(30921), Fraction(-32972), Fraction(8240)]
q_row = [Fraction(33750), Fraction(-36000), Fraction(9000)]

cmf_pairs = []
for N in range(NMAX + 4):
    cmf_pairs.append((q_row[0], p_row[0]))
    if N < NMAX + 3:
        MH = MH_at(N)
        p_new = [sum(p_row[i] * MH[i][j] for i in range(3)) for j in range(3)]
        q_new = [sum(q_row[i] * MH[i][j] for i in range(3)) for j in range(3)]
        p_row = p_new
        q_row = q_new

print(f"  CMF: Q̂_0 = {cmf_pairs[0][0]}, P̂_0 = {cmf_pairs[0][1]}", flush=True)
print(f"  CMF: Q̂_1 = {cmf_pairs[1][0]}, P̂_1 = {cmf_pairs[1][1]}", flush=True)
print(f"  CMF: Q̂_2 = {cmf_pairs[2][0]}, P̂_2 = {cmf_pairs[2][1]}", flush=True)

# Derive scalar recurrence coefficients from matrix iteration
# For column 0, the scalar recurrence Σ_{j=0}^3 ell_j(n) u_{n+j} = 0
# We compute these by finding the nullspace of [e_0, M_H(n)e_0, M_H(n)M_H(n+1)e_0, ...]
# For now, just use the recurrence by checking the CMF pairs directly

# Test the Jacobi family
print("\n=== Testing R_n^J against CMF recurrence ===", flush=True)

print("Computing R_n^J and moment pairs...", flush=True)
jacobi_pairs = []
for N in range(NMAX + 4):
    R = R_jacobi(N)
    q, p = moment_pair(R)
    jacobi_pairs.append((q, p))
    if N < 5:
        print(f"  R_{N}^J: deg={len(R)-1}, Q^J_{N}={q}, eval(-1)={poly_eval(R, Fraction(-1))}", flush=True)
    elif N % 10 == 0:
        print(f"  R_{N}^J computed, Q^J={float(q):.6e}", flush=True)

# Q4886 says Q_n^J = κ_n · B_n² where κ_n = (4n+1)/2
# Verify
for N in range(5):
    J, _ = J_and_dJ(N)
    B = poly_eval(J, Fraction(-1))
    kappa = Fraction(4*N+1, 2)
    R = R_jacobi(N)
    Rm1 = poly_eval(R, Fraction(-1))
    expected = kappa * B**2
    print(f"  N={N}: R(-1)={Rm1}, κB²={expected}, match={Rm1==expected}")

# To test recurrence: need the scalar recurrence coefficients ell_j(n)
# The CMF has a characteristic polynomial (degree 13) recurrence.
# Instead of deriving it symbolically, test numerically:
# If R_n^J satisfies the same recurrence as the CMF column,
# then Σ_j ell_j(n) · (q_j^J · G - p_j^J) = 0 for all n
# ⟺ Σ_j ell_j(n) · q_{n+j}^J = 0 AND Σ_j ell_j(n) · p_{n+j}^J = 0

# We can find the recurrence from the CMF pairs directly
# The column-0 sequence satisfies: Σ_{j=0}^3 ell_j(n) · Q̂_{n+j} = 0
# Find ell_j by solving from 4 consecutive Q̂ values

# Better: test the ratio q_n^J / Q̂_n (should be constant if same recurrence)
print("\n=== Ratio test: q_n^J / Q̂_n ===", flush=True)
for N in range(min(15, len(jacobi_pairs))):
    q_cmf = cmf_pairs[N][0]
    q_jac = jacobi_pairs[N][0]
    if q_cmf != 0 and q_jac != 0:
        ratio = q_jac / q_cmf
        print(f"  N={N}: q^J/Q̂ = {float(ratio):.10e} ({ratio})")
    else:
        print(f"  N={N}: q_cmf={q_cmf}, q_jac={q_jac}")

# Check if the recurrence for the CMF column annihilates the Jacobi pairs
# First find the recurrence from CMF data
print("\n=== Finding CMF column-0 recurrence ===", flush=True)

# Order-3 recurrence: ell_0(n) u_n + ell_1(n) u_{n+1} + ell_2(n) u_{n+2} + ell_3(n) u_{n+3} = 0
# ell_j(n) are polynomials in n of degree d
# Search for d

q_vals = [cmf_pairs[N][0] for N in range(NMAX + 4)]

for d in range(20):
    n_unknowns = 4 * (d + 1)
    n_train = n_unknowns + 5
    n_holdout = 5

    if n_train + n_holdout + 3 >= len(q_vals):
        break

    A_rows = []
    for n_idx in range(n_train + n_holdout):
        row = []
        for j in range(4):
            for k in range(d + 1):
                row.append(Fraction(n_idx)**k * q_vals[n_idx + j])
        A_rows.append(row)

    # Find null space
    m = n_train
    aug = [list(A_rows[i]) for i in range(m)]

    pivot_cols = []
    row_idx = 0
    for col in range(n_unknowns):
        found = -1
        for rr in range(row_idx, m):
            if aug[rr][col] != 0:
                found = rr
                break
        if found == -1:
            continue
        aug[row_idx], aug[found] = aug[found], aug[row_idx]
        piv = aug[row_idx][col]
        for j2 in range(n_unknowns):
            aug[row_idx][j2] /= piv
        for rr in range(m):
            if rr == row_idx:
                continue
            if aug[rr][col] == 0:
                continue
            factor = aug[rr][col]
            for j2 in range(n_unknowns):
                aug[rr][j2] -= factor * aug[row_idx][j2]
        pivot_cols.append(col)
        row_idx += 1

    rank = len(pivot_cols)
    nullity = n_unknowns - rank

    if nullity > 0:
        free_cols = [c for c in range(n_unknowns) if c not in pivot_cols]
        x = [Fraction(0)] * n_unknowns
        x[free_cols[0]] = Fraction(1)
        for pi in range(rank - 1, -1, -1):
            pc = pivot_cols[pi]
            val = Fraction(0)
            for j2 in range(n_unknowns):
                if j2 != pc:
                    val += aug[pi][j2] * x[j2]
            x[pc] = -val

        # Verify on holdout
        all_ok = True
        for n_idx in range(n_train, n_train + n_holdout):
            check = Fraction(0)
            for j in range(4):
                for k in range(d + 1):
                    check += x[j*(d+1)+k] * Fraction(n_idx)**k * q_vals[n_idx + j]
            if check != 0:
                all_ok = False
                break

        if all_ok:
            print(f"  CMF recurrence found: order 3, coeff degree {d}", flush=True)

            # Now test Jacobi pairs against this recurrence
            print(f"\n=== Testing R_n^J against CMF recurrence (d={d}) ===", flush=True)
            q_jac = [jacobi_pairs[N][0] for N in range(len(jacobi_pairs))]
            p_jac = [jacobi_pairs[N][1] for N in range(len(jacobi_pairs))]

            max_n = min(len(jacobi_pairs) - 3, NMAX)
            fail_q = None
            fail_p = None

            for n_idx in range(max_n):
                rq = Fraction(0)
                rp = Fraction(0)
                for j in range(4):
                    for k in range(d + 1):
                        coeff = x[j*(d+1)+k] * Fraction(n_idx)**k
                        rq += coeff * q_jac[n_idx + j]
                        rp += coeff * p_jac[n_idx + j]
                if rq != 0 and fail_q is None:
                    fail_q = (n_idx, rq)
                if rp != 0 and fail_p is None:
                    fail_p = (n_idx, rp)

            if fail_q is None and fail_p is None:
                print(f"  *** R_n^J SATISFIES the CMF recurrence! ***")
            else:
                if fail_q:
                    print(f"  R_n^J q-component FAILS at n={fail_q[0]}")
                if fail_p:
                    print(f"  R_n^J p-component FAILS at n={fail_p[0]}")
                print(f"  R_n^J does NOT satisfy the CMF recurrence directly")

            # Store recurrence for band search
            ell_coeffs = x
            ell_degree = d
            break

    if d % 5 == 0:
        print(f"  d={d}: rank={rank}, nullity={nullity}", flush=True)

print("\nDone.")
