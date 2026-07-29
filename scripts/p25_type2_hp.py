#!/usr/bin/env python3
"""
P2.5: Compute confluent type-II Hermite-Padé polynomials
and compare with CMF denominator.

From Q4914: The confluent JP system uses measures:
  μ₀: dσ = dx/(1+x²) on [0,1]  (wait, need to check)

Actually, G = ∫₀¹ (-log t)/(1+t²) dt. With x = t²:
  G = (1/2) ∫₀¹ (-log x)/(x(1+x)) · √x dx ...

Let me use the simpler weight: σ(x) = 1/(1+x) on [0,1].
Then F(z) = ∫₀¹ dσ/(z-x) = ∫₀¹ dx/((1+x)(z-x))

Moments: m_k = ∫₀¹ x^k/(1+x) dx
= ∫₀¹ x^k Σ_{j≥0} (-x)^j dx = Σ_j (-1)^j/(k+j+1)

For the confluent type-II:
Two functions: f₀(z) = ∫₀¹ dσ₀/(z-x), f₁(z) = ∫₀¹ (-log x)dσ₀/(z-x)
where dσ₀ = dx/(1+x)

Moments:
m_k^{(0)} = ∫₀¹ x^k/(1+x) dx = Σ_j (-1)^j/(k+j+1)
m_k^{(1)} = ∫₀¹ x^k (-log x)/(1+x) dx = Σ_j (-1)^j/(k+j+1)²

Note: m_0^{(1)} = ∫₀¹ (-log x)/(1+x) dx = G (Catalan!)
More generally: m_k^{(1)} = Σ_j (-1)^j/(k+j+1)² = S_k

For the type-II HP with multi-index (n,n):
Find P_{n,n}(x) of degree 2n such that:
∫₀¹ x^j P_{n,n}(x)/(1+x) dx = 0  for j = 0, ..., n-1
∫₀¹ x^j P_{n,n}(x)(-log x)/(1+x) dx = 0  for j = 0, ..., n-1

Total: 2n conditions for 2n+1 coefficients → monic polynomial.
"""
from fractions import Fraction as Fr
from math import comb

def compute_moments(KMAX):
    """Compute m_k^{(0)} and m_k^{(1)} for k=0,...,KMAX."""
    # m_k^{(0)} = ∫₀¹ x^k/(1+x) dx = Σ_j (-1)^j/(k+j+1)
    # m_k^{(1)} = ∫₀¹ x^k(-log x)/(1+x) dx = Σ_j (-1)^j/(k+j+1)²

    JMAX = 200  # enough terms for convergence
    m0 = []
    m1 = []
    for k in range(KMAX + 1):
        s0 = Fr(0)
        s1 = Fr(0)
        for j in range(JMAX):
            term0 = Fr((-1)**j, k + j + 1)
            term1 = Fr((-1)**j, (k + j + 1)**2)
            s0 += term0
            s1 += term1
        m0.append(s0)
        m1.append(s1)
    return m0, m1

def type2_poly(n, m0, m1):
    """
    Compute monic type-II HP polynomial P_{n,n}(x) of degree 2n.
    Orthogonality conditions:
      ∫ x^j P(x) dσ₀ = 0 for j = 0,...,n-1
      ∫ x^j P(x) dσ₁ = 0 for j = 0,...,n-1

    P(x) = x^{2n} + a_{2n-1} x^{2n-1} + ... + a_0

    ∫ x^j P(x) dσ₀ = Σ_{i=0}^{2n} a_i · m0[i+j] = 0  (a_{2n} = 1)
    → Σ_{i=0}^{2n-1} a_i · m0[i+j] = -m0[2n+j]

    Similarly for dσ₁.
    """
    deg = 2 * n
    # 2n unknowns: a_0, ..., a_{2n-1}
    # 2n equations

    # Build matrix
    mat = []
    rhs = []
    for j in range(n):
        row = []
        for i in range(deg):
            row.append(m0[i + j])
        mat.append(row)
        rhs.append(-m0[deg + j])
    for j in range(n):
        row = []
        for i in range(deg):
            row.append(m1[i + j])
        mat.append(row)
        rhs.append(-m1[deg + j])

    # Solve mat · a = rhs (exact fractions)
    N = deg
    assert len(mat) == N and len(mat[0]) == N

    # Gaussian elimination
    aug = [mat[i][:] + [rhs[i]] for i in range(N)]
    for col in range(N):
        # Find pivot
        pivot_row = None
        for row in range(col, N):
            if aug[row][col] != 0:
                pivot_row = row
                break
        if pivot_row is None:
            print(f"  Singular at col {col}!")
            return None
        aug[col], aug[pivot_row] = aug[pivot_row], aug[col]
        for row in range(N):
            if row != col and aug[row][col] != 0:
                factor = aug[row][col] / aug[col][col]
                for c in range(N + 1):
                    aug[row][c] -= factor * aug[col][c]

    coeffs = [aug[i][N] / aug[i][i] for i in range(N)]
    coeffs.append(Fr(1))  # monic: a_{2n} = 1
    return coeffs  # P(x) = Σ coeffs[i] x^i

def eval_poly(coeffs, x):
    """Evaluate polynomial at x (Horner)."""
    result = Fr(0)
    for i in range(len(coeffs) - 1, -1, -1):
        result = result * x + coeffs[i]
    return result

# Main computation
NMAX = 30
print("Computing moments...", flush=True)
m0, m1 = compute_moments(4 * NMAX + 10)

print(f"m_0^(0) = {float(m0[0]):.15f} (should be log 2 = 0.6931...)")
print(f"m_0^(1) = {float(m1[0]):.15f} (should be G = 0.9160...)")

# Compute type-II polynomials and evaluate at x = -1
print(f"\nComputing type-II HP polynomials P_{{n,n}}(-1) for n=0,...,{NMAX}...", flush=True)

pnn_vals = []
for n in range(NMAX + 1):
    if n == 0:
        # P_{0,0}(x) = 1
        pnn_vals.append(Fr(1))
        continue
    coeffs = type2_poly(n, m0, m1)
    if coeffs is None:
        print(f"  n={n}: FAILED")
        break
    val = eval_poly(coeffs, Fr(-1))
    pnn_vals.append(val)
    if n <= 5 or n % 5 == 0:
        print(f"  P_{{{n},{n}}}(-1) = {float(val):.10e}", flush=True)

# Now compute CMF Q̂_N for comparison
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

# Compare ratios Q̂_N / P_{n,n}(-1)
print(f"\nComparison Q̂_N / P_{{N,N}}(-1):")
for N in range(1, min(NMAX + 1, len(pnn_vals))):
    if pnn_vals[N] == 0:
        continue
    ratio = Qhat[N] / pnn_vals[N]
    print(f"  N={N}: Q̂/P = {float(ratio):.15e}")

# Check consecutive ratios of Q̂ and P
print(f"\nConsecutive ratios:")
print(f"  {'N':>3s}  {'Q̂_{N+1}/Q̂_N':>20s}  {'P_{N+1}/P_N':>20s}")
for N in range(1, min(NMAX, len(pnn_vals) - 1)):
    if Qhat[N] == 0 or pnn_vals[N] == 0:
        continue
    rq = float(Qhat[N+1] / Qhat[N])
    rp = float(pnn_vals[N+1] / pnn_vals[N])
    print(f"  {N:3d}  {rq:20.10f}  {rp:20.10f}")

print("\nDone.")
