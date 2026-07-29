#!/usr/bin/env python3
"""P2.5: Contiguous band search — Q̂_n = Σ P_r(n) Q^J_{n+r}.

Tests whether the CMF Q̂_n can be expressed as a polynomial-coefficient
linear combination of Jacobi moment pairs Q^J_{n+r}.
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
    J, dJ = J_and_dJ(N)
    B = poly_eval(J, Fraction(-1))
    dB = poly_eval(dJ, Fraction(-1))
    num = poly_sub(poly_mul_scalar(dJ, B), poly_mul_scalar(J, dB))
    C, rem = poly_divmod(num, [Fraction(1), Fraction(1)])
    kappa = Fraction(4*N+1, 2)
    BJ = poly_mul_scalar(J, B)
    Cprime = poly_deriv(C)
    twoXCp = [Fraction(0)] + poly_mul_scalar(Cprime, Fraction(2))
    inner = poly_add(twoXCp, C)
    half_inner = poly_mul_scalar(inner, Fraction(1,2))
    xp1_half = poly_mul(half_inner, [Fraction(1), Fraction(1)])
    R = poly_mul_scalar(poly_sub(BJ, xp1_half), kappa)
    return R

def catalan_monomial_pair(k):
    q = Fraction((-1)**k)
    partial = sum(Fraction((-1)**j, (2*j+1)**2) for j in range(k))
    p = q * partial
    return q, p

def moment_pair(R_coeffs):
    q = Fraction(0)
    p = Fraction(0)
    for k, a in enumerate(R_coeffs):
        if a == 0:
            continue
        qk, pk = catalan_monomial_pair(k)
        q += a * qk
        p += a * pk
    return q, p

# CMF pairs
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

NMAX = 50

print("Computing CMF pairs...", flush=True)
p_row = [Fraction(30921), Fraction(-32972), Fraction(8240)]
q_row = [Fraction(33750), Fraction(-36000), Fraction(9000)]
cmf_q = []
cmf_p = []
for N in range(NMAX):
    cmf_q.append(q_row[0])
    cmf_p.append(p_row[0])
    MH = MH_at(N)
    p_new = [sum(p_row[i] * MH[i][j] for i in range(3)) for j in range(3)]
    q_new = [sum(q_row[i] * MH[i][j] for i in range(3)) for j in range(3)]
    p_row = p_new
    q_row = q_new
print(f"  {len(cmf_q)} CMF pairs", flush=True)

print("Computing Jacobi pairs...", flush=True)
jac_q = []
jac_p = []
for N in range(NMAX):
    R = R_jacobi(N)
    q, p = moment_pair(R)
    jac_q.append(q)
    jac_p.append(p)
    if N % 10 == 0:
        print(f"  N={N} done", flush=True)
print(f"  {len(jac_q)} Jacobi pairs", flush=True)

# Band search: cmf_q[n] = Σ_{r=0}^m Σ_{k=0}^d c_{r,k} n^k jac_q[n+r]
# AND similarly for p component

print("\n=== Band search: Q̂_n = Σ P_r(n) Q^J_{n+r} ===", flush=True)

for m in range(6):  # band width
    for d in range(15):  # polynomial degree
        n_unknowns = (m + 1) * (d + 1)
        n_train = n_unknowns + 5
        n_holdout = 5

        if n_train + n_holdout + m >= len(cmf_q):
            break

        # Build system using BOTH q and p components
        A_rows = []
        b_vec = []
        for n in range(n_train + n_holdout):
            # q equation
            row_q = []
            for r in range(m + 1):
                for k in range(d + 1):
                    row_q.append(Fraction(n)**k * jac_q[n + r])
            A_rows.append(row_q)
            b_vec.append(cmf_q[n])

            # p equation
            row_p = []
            for r in range(m + 1):
                for k in range(d + 1):
                    row_p.append(Fraction(n)**k * jac_p[n + r])
            A_rows.append(row_p)
            b_vec.append(cmf_p[n])

        # Solve train system (using first 2*n_train equations)
        n_eq = 2 * n_train
        aug = [list(A_rows[i]) + [b_vec[i]] for i in range(n_eq)]
        n_cols = n_unknowns

        pivot_cols = []
        row_idx = 0
        for col in range(n_cols):
            found = -1
            for rr in range(row_idx, n_eq):
                if aug[rr][col] != 0:
                    found = rr
                    break
            if found == -1:
                continue
            aug[row_idx], aug[found] = aug[found], aug[row_idx]
            piv = aug[row_idx][col]
            for j2 in range(n_cols + 1):
                aug[row_idx][j2] /= piv
            for rr in range(n_eq):
                if rr == row_idx:
                    continue
                if aug[rr][col] == 0:
                    continue
                factor = aug[rr][col]
                for j2 in range(n_cols + 1):
                    aug[rr][j2] -= factor * aug[row_idx][j2]
            pivot_cols.append(col)
            row_idx += 1

        rank = len(pivot_cols)
        if rank < n_unknowns:
            if d <= 2:
                print(f"  m={m}, d={d}: rank={rank}<{n_unknowns}, underdetermined", flush=True)
            continue

        # Check consistency (remaining rows should be zero)
        consistent = True
        for i in range(rank, n_eq):
            if aug[i][n_cols] != 0:
                consistent = False
                break

        if not consistent:
            if d <= 2 or d % 5 == 0:
                print(f"  m={m}, d={d}: INCONSISTENT (q and p can't both be satisfied)", flush=True)
            continue

        # Extract solution
        x = [Fraction(0)] * n_unknowns
        for pi, pc in enumerate(pivot_cols):
            x[pc] = aug[pi][n_cols]

        # Verify on holdout
        all_ok = True
        for n in range(n_train, n_train + n_holdout):
            pred_q = Fraction(0)
            pred_p = Fraction(0)
            for r in range(m + 1):
                for k in range(d + 1):
                    coeff = x[r * (d + 1) + k] * Fraction(n)**k
                    pred_q += coeff * jac_q[n + r]
                    pred_p += coeff * jac_p[n + r]
            if pred_q != cmf_q[n] or pred_p != cmf_p[n]:
                all_ok = False
                break

        if all_ok:
            print(f"\n  *** EXACT MATCH: m={m}, d={d} ***", flush=True)
            for r in range(m + 1):
                coeffs = [x[r*(d+1)+k] for k in range(d+1)]
                nonzero = [(k,c) for k,c in enumerate(coeffs) if c != 0]
                if nonzero:
                    terms = [f"({c})·n^{k}" if k > 0 else f"{c}" for k, c in nonzero]
                    print(f"  P_{r}(n) = {' + '.join(terms)}")

            # Full verification
            ok = True
            for n in range(n_train + n_holdout, len(cmf_q) - m):
                pred_q = Fraction(0)
                for r in range(m + 1):
                    for k in range(d + 1):
                        pred_q += x[r*(d+1)+k] * Fraction(n)**k * jac_q[n+r]
                if pred_q != cmf_q[n]:
                    ok = False
                    print(f"  FAIL at n={n}")
                    break
            if ok:
                print(f"  Verified for ALL n=0..{len(cmf_q)-m-1}")
            break
        else:
            if d <= 2 or d % 5 == 0:
                print(f"  m={m}, d={d}: holdout fails", flush=True)

    else:
        continue
    break

print("\nDone.")
