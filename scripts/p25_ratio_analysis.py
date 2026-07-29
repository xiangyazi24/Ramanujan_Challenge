#!/usr/bin/env python3
"""P2.5: Deep analysis of Q̂_n / D_n².

If Q̂_n = Σ_k F_D(n,k) · R(n,k), what is R(n,k)?
We know Q̂_n / D_n² ~ 159000·n. Let's find the exact structure.
"""
from fractions import Fraction
import math

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
    M = M_entries(n); d = delta_H(n)
    return [[M[i][j]/d for j in range(3)] for i in range(3)]

NMAX = 35
print("Computing CMF Q̂_n and P̂_n values...", flush=True)
q_row = [Fraction(33750), Fraction(-36000), Fraction(9000)]
p_row = [Fraction(30921), Fraction(-32972), Fraction(8240)]
cmf_q = []; cmf_p = []
for N in range(NMAX):
    cmf_q.append(q_row[0]); cmf_p.append(p_row[0])
    MH = MH_at(N)
    q_row = [sum(q_row[i]*MH[i][j] for i in range(3)) for j in range(3)]
    p_row = [sum(p_row[i]*MH[i][j] for i in range(3)) for j in range(3)]

def F_D(n, k):
    return Fraction(2)**k * Fraction(math.comb(2*k, k)) * Fraction(math.comb(n, k)) * Fraction(math.comb(n+k, k))

def D_n_sq(n):
    return sum(F_D(n, k) for k in range(n+1))

def H(m):
    return sum(Fraction(1,j) for j in range(1, m+1))

def H2(m):
    return sum(Fraction(1,j**2) for j in range(1, m+1))

print("=== Ratio Q̂_n / D_n² ===", flush=True)
ratios = []
for n in range(NMAX):
    dn2 = D_n_sq(n)
    r = cmf_q[n] / dn2 if dn2 != 0 else None
    ratios.append(r)
    if n < 12:
        print(f"  n={n}: ratio = {r} = {float(r):.10f}")

# Check if ratio = polynomial + harmonic correction
print("\n=== Is ratio linear in n? ===", flush=True)
for n in range(1, 10):
    dr = ratios[n] - ratios[n-1]
    print(f"  Δ(n={n}) = {dr} = {float(dr):.10f}")

# Second differences
print("\n=== Second differences ===", flush=True)
for n in range(1, 9):
    d2 = ratios[n+1] - 2*ratios[n] + ratios[n-1]
    print(f"  Δ²(n={n}) = {float(d2):.15e}")

# Check if Q̂_n / D_n² - c·n - d involves H_n
print("\n=== Subtracting linear part ===", flush=True)
# Fit c, d from large n
c = ratios[NMAX-1] - ratios[NMAX-2]  # slope at end
d = ratios[0]  # intercept
print(f"  Approx slope c = {float(c):.10f}")
print(f"  d = ratios[0] = {float(d):.10f}")

# Better: use two points
c_exact = (ratios[20] - ratios[10]) / 10
d_exact = ratios[10] - c_exact * 10
print(f"  Better slope: {float(c_exact):.12f}")
print(f"  Better intercept: {float(d_exact):.12f}")

residuals = []
for n in range(NMAX):
    res = ratios[n] - c_exact * n - d_exact
    residuals.append(res)
    if n < 12:
        print(f"  n={n}: residual = {float(res):.15e}")

# Check if residuals ~ 1/n
print("\n=== Residuals * n ===", flush=True)
for n in range(1, 12):
    print(f"  n={n}: res*n = {float(residuals[n]*n):.12f}")

# Now try: ratio = a + b*n + c*H_n + d*H_n^(2)
print("\n=== Fitting ratio = a + b*n + c*H_n + d*H2_n ===", flush=True)
# Build system
A = []
b_vec = []
for n in range(4, 14):
    row = [Fraction(1), Fraction(n), H(n), H2(n)]
    A.append(row)
    b_vec.append(ratios[n])

# Solve 4x4
from copy import deepcopy
n_eq = len(A)
n_unk = 4
aug = [list(A[i]) + [b_vec[i]] for i in range(n_eq)]
pivot = []
ri = 0
for col in range(n_unk):
    found = -1
    for rr in range(ri, n_eq):
        if aug[rr][col] != 0:
            found = rr; break
    if found == -1: continue
    aug[ri], aug[found] = aug[found], aug[ri]
    piv = aug[ri][col]
    for j in range(n_unk+1): aug[ri][j] /= piv
    for rr in range(n_eq):
        if rr == ri: continue
        if aug[rr][col] == 0: continue
        f = aug[rr][col]
        for j in range(n_unk+1): aug[rr][j] -= f*aug[ri][j]
    pivot.append(col); ri += 1

if len(pivot) >= 4:
    x = [aug[i][n_unk] for i in range(4)]
    print(f"  a = {x[0]} = {float(x[0]):.12f}")
    print(f"  b = {x[1]} = {float(x[1]):.12f}")
    print(f"  c = {x[2]} = {float(x[2]):.12f}")
    print(f"  d = {x[3]} = {float(x[3]):.12f}")

    # Verify
    ok = True
    for n in range(15, NMAX):
        pred = x[0] + x[1]*n + x[2]*H(n) + x[3]*H2(n)
        if pred != ratios[n]:
            rel_err = abs(float(pred - ratios[n]) / float(ratios[n]))
            print(f"  n={n}: pred={float(pred):.10f}, actual={float(ratios[n]):.10f}, rel_err={rel_err:.2e}")
            ok = False
            break
    if ok:
        print("  ALL VERIFIED!")

# Try more terms: a + b*n + c*n² + d*H_n + e*n*H_n + f*H2_n
print("\n=== Fitting ratio = a + b*n + c*n² + d*H_n + e*n*H_n + f*H2_n ===", flush=True)
n_unk2 = 6
A2 = []
b_vec2 = []
for n in range(2, 2 + n_unk2 + 8):
    row = [Fraction(1), Fraction(n), Fraction(n)**2, H(n), Fraction(n)*H(n), H2(n)]
    A2.append(row)
    b_vec2.append(ratios[n])

aug2 = [list(A2[i]) + [b_vec2[i]] for i in range(len(A2))]
pivot2 = []; ri2 = 0
for col in range(n_unk2):
    found = -1
    for rr in range(ri2, len(aug2)):
        if aug2[rr][col] != 0:
            found = rr; break
    if found == -1: continue
    aug2[ri2], aug2[found] = aug2[found], aug2[ri2]
    piv = aug2[ri2][col]
    for j in range(n_unk2+1): aug2[ri2][j] /= piv
    for rr in range(len(aug2)):
        if rr == ri2: continue
        if aug2[rr][col] == 0: continue
        f = aug2[rr][col]
        for j in range(n_unk2+1): aug2[rr][j] -= f*aug2[ri2][j]
    pivot2.append(col); ri2 += 1

if len(pivot2) >= n_unk2:
    x2 = [aug2[i][n_unk2] for i in range(n_unk2)]
    for i, name in enumerate(['a','b','c','d(H_n)','e(n*H_n)','f(H2_n)']):
        print(f"  {name} = {float(x2[i]):.15e}")

    ok2 = True
    for n in range(2 + n_unk2 + 8, NMAX):
        pred = x2[0] + x2[1]*n + x2[2]*n**2 + x2[3]*H(n) + x2[4]*n*H(n) + x2[5]*H2(n)
        if pred != ratios[n]:
            rel_err = float(abs(pred - ratios[n]) / abs(ratios[n]))
            if rel_err > 1e-20:
                print(f"  FAIL n={n}: rel_err={rel_err:.2e}")
                ok2 = False
                break
    if ok2:
        print("  ALL VERIFIED! Q̂_n = D_n² · (a + b·n + c·n² + d·H_n + e·n·H_n + f·H2_n)")

# Also check: does P̂_n (the G-approximation numerator) have similar structure?
print("\n=== P̂_n / D_n² ===", flush=True)
for n in range(8):
    dn2 = D_n_sq(n)
    r = cmf_p[n] / dn2 if dn2 != 0 else None
    print(f"  n={n}: P̂_n/D_n² = {float(r):.10f}" if r else f"  n={n}: D²=0")

# Try: Q̂_n = D_n² · (A(n) + B(n)·H_n) where A,B are rational in n
# This is: look at residual after subtracting the polynomial part, see if it's H_n times something
print("\n=== Check Q̂_n against column-0 of A·M(0)·...·M(N-1) with alternative initial row ===", flush=True)
# The initial row A = [[30921,-32972,8240],[33750,-36000,9000]]
# P = A[0], Q = A[1]
# Note: Q̂_0 = 33750, D_0² = 1, so ratio_0 = 33750
# Q̂_1 = -36000 / delta_H(0) * ... actually computed via matrix iteration

# Let me look at the CMF recurrence coefficients
print("\n=== CMF scalar recurrence check ===", flush=True)
# From p25_scalar_rec.py, order 3, degree 13
# ell_0 leading = -98304·n^13, ell_3 leading = 98304·n^13
# This means Σ_{j=0}^3 ell_j(n) q_{n+j} = 0

# Verify numerically that our Q̂_n values satisfy this
print("  (Verification would require the full ell_j — skipping, was already verified)")

# Instead, let me check the D_n recurrence
print("\n=== D_n values and recurrence check ===", flush=True)
D_vals = [1, 3, 13, 63, 321, 1683, 8989, 48639, 265729, 1462563]
for n in range(8):
    dn2 = D_n_sq(n)
    dn = int(dn2.numerator**0.5) if dn2.denominator == 1 else None
    print(f"  D_{n}² = {dn2}, D_{n} = {dn}" if dn else f"  D_{n}² = {dn2}")

# Key insight: check if Q̂_n / D_n² is in Q(n) + Q(n)·H_n + Q(n)·H2_n
# by computing enough terms and fitting rational functions
print("\n=== Rational function fit for Q̂_n / D_n² - linear_part ===", flush=True)
# Remove the leading linear growth
# Compute the "fractional part" more carefully

# Use Padé-like approach: find p(n)/q(n) matching the ratio
# For degree (d,d) Padé:
for d_num in range(1, 6):
    for d_den in range(d_num, d_num+2):
        n_unk_pade = d_num + 1 + d_den  # numerator coeffs + denominator coeffs (leading=1)
        n_pts = n_unk_pade + 3

        if n_pts + 2 > len(ratios): break

        # ratio[n] ≈ p(n)/q(n) where p has degree d_num, q has degree d_den with leading coeff 1
        # So ratio[n] · q(n) = p(n)
        # ratio[n] · (n^d_den + c_{d_den-1}·n^{d_den-1} + ... + c_0) = a_{d_num}·n^{d_num} + ... + a_0
        # This gives a linear system in the unknowns a_0,...,a_{d_num}, c_0,...,c_{d_den-1}

        rows = []; bv = []
        for n in range(2, 2 + n_pts):
            # LHS of ratio[n]·q(n) = p(n)
            # ratio[n]·n^d_den - a_{d_num}·n^{d_num} - ... = -ratio[n]·(c_{d_den-1}·n^{d_den-1} + ...)
            row = []
            # p coefficients: -n^0, -n^1, ..., -n^{d_num}
            for k in range(d_num + 1):
                row.append(-Fraction(n)**k)
            # q coefficients (c_0 to c_{d_den-1}): ratio[n]·n^0, ..., ratio[n]·n^{d_den-1}
            for k in range(d_den):
                row.append(ratios[n] * Fraction(n)**k)
            rows.append(row)
            bv.append(-ratios[n] * Fraction(n)**d_den)

        # Solve
        m = len(rows)
        nc = len(rows[0])
        aug3 = [list(rows[i]) + [bv[i]] for i in range(m)]
        p3 = []; r3 = 0
        for col in range(nc):
            found = -1
            for rr in range(r3, m):
                if aug3[rr][col] != 0: found = rr; break
            if found == -1: continue
            aug3[r3], aug3[found] = aug3[found], aug3[r3]
            pv = aug3[r3][col]
            for j in range(nc+1): aug3[r3][j] /= pv
            for rr in range(m):
                if rr == r3: continue
                if aug3[rr][col] == 0: continue
                f = aug3[rr][col]
                for j in range(nc+1): aug3[rr][j] -= f*aug3[r3][j]
            p3.append(col); r3 += 1

        if len(p3) >= nc:
            sol = [aug3[i][nc] for i in range(nc)]
            # Check holdout
            ok3 = True
            for n in range(2 + n_pts, min(2 + n_pts + 5, NMAX)):
                p_val = sum(sol[k] * Fraction(n)**k for k in range(d_num+1))
                q_val = Fraction(n)**d_den + sum(sol[d_num+1+k] * Fraction(n)**k for k in range(d_den))
                pred = p_val / q_val
                if pred != ratios[n]:
                    ok3 = False; break

            if ok3:
                print(f"  Padé ({d_num},{d_den}) MATCHES!")
                print(f"  Numerator coeffs: {[float(sol[k]) for k in range(d_num+1)]}")
                print(f"  Denominator coeffs (below leading): {[float(sol[d_num+1+k]) for k in range(d_den)]}")
                break
        else:
            pass  # underdetermined, skip
    else:
        continue
    break

print("\nDone.")
