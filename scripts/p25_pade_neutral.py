#!/usr/bin/env python3
"""P2.5: Padé reconstruction of the neutral multiplier r₀(n) and eigenvector v(n).

Compute r₀(n) and v(n) to high precision for n = 0..80, then try to identify
r₀(n) as a rational function of n using Padé/rational reconstruction.
"""
from mpmath import mp, mpf, matrix, polyroots, nstr
from fractions import Fraction

mp.dps = 80

def M_exact(n):
    m11 = (-2*n-5)*(n+3)**2*(136*n**4+1424*n**3+5548*n**2+9551*n+6141)
    m12 = 384*n**6+6384*n**5+44168*n**4+162698*n**3+336377*n**2+369933*n+169011
    m13 = -480*n**4-4980*n**3-19210*n**2-32690*n-20730
    m21 = (n+2)**2*(n+3)**2*(4*n+10)*(48*n**3+386*n**2+1017*n+879)
    m22 = (n+2)**2*(-272*n**5-3848*n**4-21732*n**3-61184*n**2-85761*n-47808)
    m23 = (n+2)**2*(320*n**3+2540*n**2+6610*n+5640)
    m31 = (-4*n-10)*(n+2)**2*(n+3)**2*(32*n**4+302*n**3+1037*n**2+1530*n+813)
    m32 = (n+2)**2*(192*n**6+2984*n**5+19116*n**4+64452*n**3+120256*n**2+117279*n+46476)
    m33 = (n+2)**2*(-16*n**5-408*n**4-2912*n**3-8884*n**2-12254*n-6240)
    return matrix([[m11,m12,m13],[m21,m22,m23],[m31,m32,m33]])

def delta(n):
    return mpf(-2) * (n+2)**2 * (n+3)**2 * (2*n+5) * (2*n+7)**2

def balanced_matrix(n):
    M = M_exact(n)
    d = delta(n)
    n1 = mpf(n+1)
    n2 = mpf(n+2)
    B = matrix(3, 3)
    for i in range(3):
        for j in range(3):
            scale_inv = [1, 1/n1, 1/n1**2][i]
            scale_next = [1, n2, n2**2][j]
            B[i,j] = M[i,j] * scale_inv * scale_next / d
    return B

def neutral_eigen(n):
    B = balanced_matrix(n)
    tr = B[0,0] + B[1,1] + B[2,2]
    cofsum = (B[0,0]*B[1,1] - B[0,1]*B[1,0]
            + B[0,0]*B[2,2] - B[0,2]*B[2,0]
            + B[1,1]*B[2,2] - B[1,2]*B[2,1])
    det = (B[0,0]*(B[1,1]*B[2,2]-B[1,2]*B[2,1])
         - B[0,1]*(B[1,0]*B[2,2]-B[1,2]*B[2,0])
         + B[0,2]*(B[1,0]*B[2,1]-B[1,1]*B[2,0]))
    roots = polyroots([1, -tr, cofsum, -det])
    neutral_idx = min(range(3), key=lambda i: abs(roots[i] - 1))
    lam = roots[neutral_idx]

    C = matrix(3, 3)
    for i in range(3):
        for j in range(3):
            C[i,j] = B[i,j] - (lam if i==j else 0)
    row0 = [C[0,0], C[0,1], C[0,2]]
    row1 = [C[1,0], C[1,1], C[1,2]]
    v1 = row0[1]*row1[2] - row0[2]*row1[1]
    v2 = row0[2]*row1[0] - row0[0]*row1[2]
    v3 = row0[0]*row1[1] - row0[1]*row1[0]
    if abs(v3) > 1e-50:
        return lam, v1/v3, v2/v3, mpf(1)
    return lam, v1, v2, v3

# Compute r₀(n) for n = 0..60
N_MAX = 60
r0_vals = []
v1_vals = []
v2_vals = []

print("Computing neutral eigendata for n = 0..%d..." % N_MAX)
for n in range(N_MAX + 1):
    lam, v1, v2, v3 = neutral_eigen(n)
    r0_vals.append(lam)
    v1_vals.append(v1)
    v2_vals.append(v2)

print("Done.\n")

# Try rational reconstruction of r₀(n)
# If r₀(n) = P(n)/Q(n) with deg P = p, deg Q = q, then
# we need p + q + 2 values to determine P and Q.
# Try various (p, q) pairs.

def rational_recon(vals, p, q, start=5):
    """Try to fit vals[start], vals[start+1], ..., vals[start+p+q+1] to P(n)/Q(n)
    with deg P = p, deg Q = q. Return (P_coeffs, Q_coeffs, residual)."""
    npts = p + q + 2
    if start + npts > len(vals):
        return None
    # Set up linear system: Q(n)*r₀(n) = P(n)
    # Q(n) = n^q + a_{q-1} n^{q-1} + ... + a_0 (monic)
    # P(n) = b_p n^p + ... + b_0
    # So for each sample point n_i:
    # b_p n_i^p + ... + b_0 - r₀(n_i) * (n_i^q + a_{q-1} n_i^{q-1} + ... + a_0) = 0
    # Unknowns: b_0, ..., b_p, a_0, ..., a_{q-1} (total p+1+q)
    # But we have p+q+2 equations and p+q+1 unknowns, so overdetermined.
    # Use p+q+1 equations, then check the last one.

    num_unknowns = p + 1 + q
    A = matrix(npts, num_unknowns)
    b = matrix(npts, 1)
    for idx in range(npts):
        n = start + idx
        rn = vals[n]
        # P coefficients: b_k contributes n^k
        for k in range(p + 1):
            A[idx, k] = mpf(n)**k
        # Q coefficients: a_k contributes -r₀(n) * n^k
        for k in range(q):
            A[idx, p + 1 + k] = -rn * mpf(n)**k
        # RHS: r₀(n) * n^q (from the monic leading term of Q)
        b[idx, 0] = rn * mpf(n)**q

    # Solve using first num_unknowns equations
    A_sq = A[:num_unknowns, :]
    b_sq = b[:num_unknowns, :]
    try:
        x = mp.lu_solve(A_sq, b_sq)
    except:
        return None

    # Check last equation
    residual = sum(A[npts-1, j] * x[j] for j in range(num_unknowns)) - b[npts-1, 0]

    P_coeffs = [x[k] for k in range(p+1)]
    Q_coeffs = [x[p+1+k] for k in range(q)] + [mpf(1)]

    return P_coeffs, Q_coeffs, residual

print("="*80)
print("Padé reconstruction of r₀(n)")
print("="*80)

for p in range(1, 12):
    for q in range(1, 12):
        result = rational_recon(r0_vals, p, q, start=10)
        if result is None:
            continue
        P, Q, res = result
        if abs(res) < mpf(10)**(-40):
            # Found a candidate! Verify against more points
            max_err = mpf(0)
            for n in range(N_MAX + 1):
                pval = sum(P[k] * mpf(n)**k for k in range(len(P)))
                qval = sum(Q[k] * mpf(n)**k for k in range(len(Q)))
                if abs(qval) > 1e-50:
                    predicted = pval / qval
                    err = abs(predicted - r0_vals[n])
                    max_err = max(max_err, err)
            if max_err < mpf(10)**(-30):
                print(f"  (p,q)=({p},{q}): residual={nstr(res,5)}, max_err={nstr(max_err,5)}")
                # Try to identify rational coefficients
                print(f"    P coeffs: {[nstr(c, 20) for c in P]}")
                print(f"    Q coeffs: {[nstr(c, 20) for c in Q]}")
                # Round to rational
                from fractions import Fraction
                P_rat = [Fraction(float(c)).limit_denominator(10000) for c in P]
                Q_rat = [Fraction(float(c)).limit_denominator(10000) for c in Q]
                rat_err = mpf(0)
                for n in range(N_MAX + 1):
                    pval = sum(mpf(P_rat[k].numerator)/mpf(P_rat[k].denominator) * mpf(n)**k for k in range(len(P_rat)))
                    qval = sum(mpf(Q_rat[k].numerator)/mpf(Q_rat[k].denominator) * mpf(n)**k for k in range(len(Q_rat)))
                    if abs(qval) > 1e-50:
                        predicted = pval / qval
                        rat_err = max(rat_err, abs(predicted - r0_vals[n]))
                print(f"    Rational approx error: {nstr(rat_err, 5)}")
                if rat_err < mpf(10)**(-30):
                    print(f"    P_rat: {[str(c) for c in P_rat]}")
                    print(f"    Q_rat: {[str(c) for c in Q_rat]}")
                print()

# Also try direct: check if r₀(n) = product of known linear factors
print("="*80)
print("Test specific candidate multipliers")
print("="*80)

# Candidate 1: r₀(n) = (n+a)^2(2n+b) / ((n+c)^2(2n+d)) for various a,b,c,d
# From Q4855 §3.1: first candidate is r₀(n) = (2n+3)(n+1)² / ((2n+5)(n+2)²)
# But that didn't match. Let me try the determinant-adapted version.

candidates = [
    ("(n+1)²(2n+3)/((n+2)²(2n+5))", lambda n: (n+1)**2*(2*n+3)/((n+2)**2*(2*n+5))),
    ("(n+1)²(2n+3)(2n+1)/((n+2)²(2n+5)(2n+7))", lambda n: (n+1)**2*(2*n+3)*(2*n+1)/((n+2)**2*(2*n+5)*(2*n+7))),
    ("(2n+1)(2n+3)/((2n+5)(2n+7))", lambda n: (2*n+1)*(2*n+3)/((2*n+5)*(2*n+7))),
    ("(n+1)³/((n+2)³)", lambda n: (n+1)**3/(n+2)**3),
    ("(n+1)²(n+3/2)/((n+2)²(n+5/2))", lambda n: (n+1)**2*(n+1.5)/((n+2)**2*(n+2.5))),
]

for name, f in candidates:
    errs = [abs(f(n) - float(r0_vals[n])) for n in range(1, 30)]
    max_err = max(errs)
    print(f"  {name}: max err (n=1..29) = {max_err:.6e}")

# Check the asymptotic expansion more carefully
print()
print("="*80)
print("Asymptotic expansion of r₀(n)")
print("="*80)
print("n   r₀(n)              1-r₀     (1-r₀)*n     (1-r₀)*n²/(n·c)")

# Compute (1-r₀)·n for large n
for n in [20, 30, 40, 50, 60]:
    r = r0_vals[n]
    x = 1 - r
    print(f"  n={n:2d}: r₀={nstr(r,20)}, 1-r₀={nstr(x,15)}, (1-r₀)·n={nstr(x*n,15)}, (1-r₀)·n²={nstr(x*n*n,15)}")
