#!/usr/bin/env python3
"""
P2.5: Type-II Hermite-Padé comparison using mpmath for accurate moments.

Weight: dσ₀ = dx/(1+x) on [0,1]
Functions: f₀(z) = ∫ dσ₀/(z-x), f₁(z) = ∫ (-log x) dσ₀/(z-x)

Moments (closed form):
  m_k^{(0)} = ∫₀¹ x^k/(1+x) dx = (-1)^k [ln2 - H_k^{alt}]
    where H_k^{alt} = Σ_{j=1}^k (-1)^{j-1}/j

  m_k^{(1)} = ∫₀¹ x^k(-log x)/(1+x) dx (involves G and rationals)
    = (-1)^k [G_k - ln2 · H_k^{alt,1}]  ... actually compute numerically

Type-II HP poly P_{n,n}(x): degree 2n, monic, orthogonal to
both x^j dσ₀ and x^j dσ₁ for j=0,...,n-1.

Compare P_{n,n}(-1) growth rate with CMF Q̂_N growth rate.
"""
from mpmath import mp, mpf, matrix, log, catalan, pi, quad, power, fac

mp.dps = 80

def compute_moments(kmax):
    m0 = []
    m1 = []
    for k in range(kmax + 1):
        # Exact numerical integration with mpmath
        v0 = quad(lambda x: power(x, k) / (1 + x), [0, 1])
        v1 = quad(lambda x: power(x, k) * (-log(x)) / (1 + x), [0, 1])
        m0.append(v0)
        m1.append(v1)
        if k % 20 == 0:
            print(f"  moments k={k} done", flush=True)
    return m0, m1

def type2_poly_eval_at_minus1(n, m0, m1):
    """
    Compute P_{n,n}(-1) using mpmath linear algebra.
    P(x) = x^{2n} + a_{2n-1}x^{2n-1} + ... + a_0 (monic, deg 2n)
    Orthogonality: 2n conditions from two measures.
    """
    deg = 2 * n
    # Build 2n x 2n matrix
    mat = matrix(deg, deg)
    rhs = matrix(deg, 1)

    for j in range(n):
        for i in range(deg):
            mat[j, i] = m0[i + j]
        rhs[j] = -m0[deg + j]

    for j in range(n):
        for i in range(deg):
            mat[n + j, i] = m1[i + j]
        rhs[n + j] = -m1[deg + j]

    # Solve
    try:
        coeffs_vec = mp.lu_solve(mat, rhs)
    except Exception as e:
        print(f"  n={n}: solve failed: {e}")
        return None

    # Evaluate at x = -1
    # P(-1) = (-1)^{2n} + Σ a_i (-1)^i = 1 + Σ a_i (-1)^i
    val = mpf(1)  # a_{2n} = 1, (-1)^{2n} = 1
    for i in range(deg):
        val += coeffs_vec[i] * power(-1, i)
    return val

NMAX = 25
print("Computing moments up to k =", 2*NMAX + NMAX + 5, "...", flush=True)
kmax_needed = 2 * NMAX + NMAX + 5  # deg = 2n, need m_{deg+n-1}
m0, m1 = compute_moments(kmax_needed)

print(f"\nm_0^(0) = {m0[0]}  (should be ln2 = {log(2)})")
print(f"m_0^(1) = {m1[0]}  (should be pi^2/12 = {pi**2/12})")
print(f"  (Note: ∫₀¹ (-log x)/(1+x) dx = pi²/12, NOT G)")
print(f"  G = {catalan}")

pvals = []
print(f"\nComputing type-II HP P_{{n,n}}(-1)...", flush=True)
for n in range(NMAX + 1):
    if n == 0:
        pvals.append(mpf(1))
        continue
    val = type2_poly_eval_at_minus1(n, m0, m1)
    if val is None:
        break
    pvals.append(val)
    if n <= 5 or n % 5 == 0:
        print(f"  P_{{{n},{n}}}(-1) = {mp.nstr(val, 15)}", flush=True)

# CMF computation
print(f"\nComputing CMF Q̂_N...", flush=True)

def M_entries(n):
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
    return -2*(n+2)**2*(n+3)**2*(2*n+5)*(2*n+7)**2

from fractions import Fraction as Fr
rows = [[Fr(1), Fr(0), Fr(0)],
        [Fr(0), Fr(1), Fr(0)],
        [Fr(0), Fr(0), Fr(1)]]
u = {j: [rows[j][0]] for j in range(3)}

for N in range(NMAX):
    M = M_entries(N)
    d = Fr(delta_H(N))
    MH = [[Fr(M[i][j]) / d for j in range(3)] for i in range(3)]
    for j in range(3):
        r = rows[j]
        new_r = [sum(r[i]*MH[i][k] for i in range(3)) for k in range(3)]
        rows[j] = new_r
        u[j].append(new_r[0])

q = [Fr(33750), Fr(-36000), Fr(9000)]
Qhat = [sum(q[j] * u[j][N] for j in range(3)) for N in range(NMAX + 1)]

# Consecutive ratios
print(f"\nConsecutive ratios:")
print(f"  {'N':>3s}  {'Q̂ ratio':>22s}  {'P ratio':>22s}  {'P/16^n ratio':>22s}")
rho = 17 + 12 * 2**0.5
for N in range(1, min(NMAX, len(pvals) - 1)):
    if Qhat[N] == 0 or pvals[N] == 0:
        continue
    rq = float(Qhat[N+1]) / float(Qhat[N])
    rp = float(pvals[N+1]) / float(pvals[N])
    # P_{n,n}(-1) / 16^n growth rate
    rp_norm = rp / 16 if N > 0 else 0
    print(f"  {N:3d}  {rq:22.10f}  {rp:22.10f}  {rp_norm:22.10f}")

print(f"\n  17+12√2 = {17+12*2**0.5:.10f}")
print(f"  16*(17+12√2) = {16*(17+12*2**0.5):.10f}")

# Direct comparison
print(f"\nQ̂_N vs P_{{N,N}}(-1) ratio:")
for N in range(1, min(NMAX + 1, len(pvals))):
    if pvals[N] == 0:
        continue
    ratio = mpf(Qhat[N].numerator) / mpf(Qhat[N].denominator) / pvals[N]
    print(f"  N={N}: Q̂/P = {mp.nstr(ratio, 20)}")

# Also try: does Q̂_N / (16^N * P_{N,N}(-1)) have a limit?
print(f"\nQ̂_N / (16^N * P_{{N,N}}(-1)):")
for N in range(1, min(NMAX + 1, len(pvals))):
    if pvals[N] == 0:
        continue
    ratio = mpf(Qhat[N].numerator) / mpf(Qhat[N].denominator) / (mpf(16)**N * pvals[N])
    print(f"  N={N}: ratio = {mp.nstr(ratio, 20)}")

print("\nDone.")
