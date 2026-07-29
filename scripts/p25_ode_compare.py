#!/usr/bin/env python3
"""P2.5: Compare the CMF module with the integrated-K module.

The integrated K ODE is: k(1-k²)Y''' + (1-3k²)Y'' - kY' = 0
Under Brafman substitution k = 4√(2z)/(1-z), this becomes a z-ODE.

We verify the module identification by checking:
1. The z-ODE solutions at z=0 include 1, K(k(z)), Y(k(z))
2. The Taylor coefficients of K(k(z)) are related to D_n² (Delannoy squares) via Brafman
3. The CMF solutions span the same space as {1, K(k(z)), Y(k(z))}

Numerical approach: compute the CMF fundamental matrix and check if the
columns match the Legendre model after appropriate gauge transform.
"""
from mpmath import mp, mpf, matrix, catalan, sqrt, pi, log, rf

mp.dps = 100

# ============================================================
# Part 1: Compute the z-ODE from the k-ODE
# ============================================================
# k-ODE: k(1-k²)Y''' + (1-3k²)Y'' - kY' = 0
# Substitution: k = 4√(2z)/(1-z), so k² = 32z/(1-z)²
#
# Rather than doing the algebra, let's verify NUMERICALLY that
# the CMF solutions match the integrated K module solutions.

def legendre_P(n, x):
    if n == 0: return mpf(1)
    if n == 1: return mpf(x)
    p0, p1 = mpf(1), mpf(x)
    for k in range(1, n):
        p2 = ((2*k+1)*x*p1 - k*p0) / (k+1)
        p0, p1 = p1, p2
    return p1

def legendre_Q(n, x):
    x = mpf(x)
    q0 = log((x+1)/(x-1)) / 2
    if n == 0: return q0
    q1 = x * q0 - 1
    if n == 1: return q1
    for k in range(1, n):
        q2 = ((2*k+1)*x*q1 - k*q0) / (k+1)
        q0, q1 = q1, q2
    return q1

# Central Delannoy numbers D_n = P_n(3)
def delannoy(n):
    return legendre_P(n, mpf(3))

# Delannoy squares
def D2(n):
    d = delannoy(n)
    return d * d

# ============================================================
# Part 2: CMF matrix and scalar recurrence
# ============================================================
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

# ============================================================
# Part 3: Check if D_n² satisfies the CMF scalar recurrence
# ============================================================
# If the CMF module contains Sym²(Delannoy) as a sub-quotient,
# then D_n² should satisfy a related equation.

# First, let's compute the CMF scalar recurrence from the matrix.
# For a 3×3 system v(n+1) = (1/δ(n)) M(n) v(n),
# the scalar recurrence for the first component v₁ is obtained
# by eliminating v₂ and v₃.
#
# Alternative: just check if specific sequences satisfy the matrix recurrence.

# Check: does the sequence D_n² satisfy any column of the matrix system?
print("="*60)
print("Part 3: Does D_n² satisfy the CMF matrix equation?")
print("="*60)

# The CMF matrix equation: A·M(0)·M(1)·...·M(N-1)
# gives 2 rows (P and Q). We want to check if D_n² appears
# as one of these sequences for some initial condition.

# For D_n² to be a solution, the vector (D_{n+2}², D_{n+1}², D_n²)
# must satisfy v(n+1) = (1/δ(n)) M(n) v(n).
# i.e., M(n) (D_{n+3}², D_{n+2}², D_{n+1}²)ᵀ = δ(n) (D_{n+2}², D_{n+1}², D_n²)ᵀ

# Wait, the convention matters. Let me check:
# If v(n) = (v₁(n), v₂(n), v₃(n))ᵀ and v(n+1) = (1/δ(n)) M(n) v(n),
# then v₁(n+1) = (1/δ(n)) [m₁₁(n)v₁(n) + m₁₂(n)v₂(n) + m₁₃(n)v₃(n)]

# But the initial matrix A is 2×3, and A·M_N gives 2 rows.
# The three COLUMNS of M_N are the three linearly independent solutions
# in the 3D solution space.

# Let me check: do the COLUMN vectors of the identity matrix I₃,
# evolved by M(0), M(1), ..., give sequences related to D_n²?

# Column 0 of M_N gives the evolution of (1,0,0)
# Column 1 gives (0,1,0), etc.

# Check for column 0:
v = matrix([[1], [0], [0]])
cols = [[v[0,0], v[1,0], v[2,0]]]
for n in range(15):
    v = M_exact(n) * v
    cols.append([v[0,0], v[1,0], v[2,0]])

print("Column 0 of M_N (first component):")
for n, c in enumerate(cols[:10]):
    d2 = D2(n)
    ratio = c[0] / d2 if d2 != 0 else None
    print(f"  n={n}: v1={mp.nstr(c[0], 20)}, D_n²={mp.nstr(d2, 10)}, ratio={mp.nstr(ratio, 15) if ratio else 'N/A'}")

# ============================================================
# Part 4: Pochhammer-normalized CMF comparison
# ============================================================
print()
print("="*60)
print("Part 4: Normalized comparison with H_n twist")
print("="*60)

# H_n = (-16)^n (2)_n^2 (3)_n^2 (5/2)_n (7/2)_n^2
def H(n):
    return (mpf(-16)**n * rf(2,n)**2 * rf(3,n)**2 * rf(mpf(5)/2,n) * rf(mpf(7)/2,n)**2)

# The normalized sequence: Q̂_n = Q_n / H_n
# Check if Q̂_n / D_n² is a rational function of n.

p0 = [mpf(30921), mpf(-32972), mpf(8240)]
q0 = [mpf(33750), mpf(-36000), mpf(9000)]

q_row = list(q0)

print("n    Q_{n,0}/H_n    D_n²    ratio Q_hat/D²    diff")
for N in range(12):
    Q_vals = list(q_row)
    Hn = H(N) if N > 0 else mpf(1)
    d2 = D2(N)
    q_hat = Q_vals[0] / Hn if Hn != 0 and N > 0 else Q_vals[0]
    ratio = q_hat / d2 if d2 != 0 else None

    if N > 0:
        print(f"  n={N:2d}: Q̂={mp.nstr(q_hat, 15):>20s}, D²={mp.nstr(d2, 10):>12s}, "
              f"Q̂/D²={mp.nstr(ratio, 15) if ratio else 'N/A':>18s}")

    if N < 11:
        M = M_exact(N)
        new_q = [mpf(0)]*3
        for col in range(3):
            for k in range(3):
                new_q[col] += q_row[k] * M[k, col]
        q_row = new_q

# ============================================================
# Part 5: Check the Sym²(Delannoy) recurrence
# ============================================================
print()
print("="*60)
print("Part 5: Sym²(Delannoy) recurrence for D_n²")
print("="*60)

# The WZ-proved recurrence from proof.tex:
# (n+2)²(2n+1) D²_{n+2} - (2n+3)P(n) D²_{n+1} + (2n+1)P(n) D²_n - (2n+3)n² D²_{n-1} = 0
# where P(n) = 35n² + 70n + 26

for n in range(1, 10):
    P_n = 35*n**2 + 70*n + 26
    lhs = ((n+2)**2*(2*n+1)*D2(n+2)
           - (2*n+3)*P_n*D2(n+1)
           + (2*n+1)*P_n*D2(n)
           - (2*n+3)*n**2*D2(n-1))
    print(f"  n={n}: recurrence residual = {mp.nstr(lhs, 5)}")

# ============================================================
# Part 6: Numerical ODE check
# ============================================================
print()
print("="*60)
print("Part 6: Verify the integrated-K ODE")
print("="*60)

# The ODE: k(1-k²)Y''' + (1-3k²)Y'' - kY' = 0
# Solutions at k=0: {1, K(k), (1/2)∫₀ᵏ K(s)ds}
#
# K(k) = (π/2) ₂F₁(1/2, 1/2; 1; k²)
# Taylor expansion in k²:
# K(k) = (π/2) Σ C(2m,m)²/4^(2m) · k^(2m)
# = (π/2) [1 + (1/4)k² + (9/64)k⁴ + ...]

# Under z = k²(1-k²)/32 ≈ ... no wait, the Brafman substitution is:
# k = 4√(2z)/(1-z), k² = 32z/(1-z)²

# Let's verify: K(k(z)) = (π(1-z)/2) Σ D_n² z^n by numerical check
z_test = mpf('0.01')
k_test = 4*sqrt(2*z_test)/(1-z_test)
K_direct = mp.ellipk(k_test**2)

brafman_sum = sum(D2(n) * z_test**n for n in range(50))
K_brafman = (pi/2) * (1-z_test) * brafman_sum

print(f"z = {z_test}")
print(f"k(z) = {mp.nstr(k_test, 15)}")
print(f"K(k²) direct = {mp.nstr(K_direct, 30)}")
print(f"K via Brafman = {mp.nstr(K_brafman, 30)}")
print(f"Agreement: {mp.nstr(abs(K_direct - K_brafman), 5)}")

# ============================================================
# Part 7: Key test - does the CMF GF match F(z) = Σ D_n² z^n?
# ============================================================
print()
print("="*60)
print("Part 7: CMF GF vs Brafman GF")
print("="*60)

# The CMF has GF: Σ Q_{N,j} z^N
# We want to check if Σ Q_{N,j}/H_N · z^N is related to F(z) = Σ D_n² z^n
# or to the integrated K module

# The normalized CMF GF should be related to the z-ODE solutions.
# The question is WHAT the relationship is.

# Check: do the normalized CMF sequences Q̂_N/D_N² approach a limiting ratio?
print("\nQ̂_{N,0}/D_N² for various N:")
q_row = list(q0)
cum_delta = mpf(1)
for N in range(25):
    Q_vals = list(q_row)
    Hn = H(N) if N > 0 else mpf(1)
    d2 = D2(N)
    if N > 0 and Hn != 0 and d2 != 0:
        q_hat = Q_vals[0] / Hn
        ratio = q_hat / d2
        if N <= 5 or N % 5 == 0:
            print(f"  N={N:2d}: Q̂/D² = {mp.nstr(ratio, 20)}")

    if N < 24:
        M = M_exact(N)
        new_q = [mpf(0)]*3
        for col in range(3):
            for k in range(3):
                new_q[col] += q_row[k] * M[k, col]
        q_row = new_q
