#!/usr/bin/env sage
"""
Two-measure step-line Hermite-Padé for ζ(2)+ζ(3).

Measures:
  dμ₂ = -log(t) dt on (0,1), moments ψ₂(k) = 1/(k+1)²
  dμ₃ = ½log²(t) dt on (0,1), moments ψ₃(k) = 1/(k+1)³

Step-line type-II HP: find Q_n(t) of degree n such that:
  ∫₀¹ t^j Q_n(t) dμ₂ = 0 for j = 0,...,n₁-1
  ∫₀¹ t^j Q_n(t) dμ₃ = 0 for j = 0,...,n₂-1

where (n₁,n₂) follows a step-line path with n₁+n₂ = n.

The step-line path: (0,0) → (1,0) → (1,1) → (2,1) → (2,2) → ...
"""

N_MAX = 15

# ============================================================
# Moment matrices
# ============================================================
def psi2(m, k):
    """Moment ψ₂(m+k) = 1/(m+k+1)²"""
    return QQ(1) / (m + k + 1)^2

def psi3(m, k):
    """Moment ψ₃(m+k) = 1/(m+k+1)³"""
    return QQ(1) / (m + k + 1)^3

# ============================================================
# Step-line path
# ============================================================
def step_line(n):
    """Return (n1, n2) for the n-th step.
    Path: (0,0),(1,0),(1,1),(2,1),(2,2),(3,2),(3,3),...
    """
    n1 = (n + 1) // 2
    n2 = n // 2
    return (n1, n2)

# ============================================================
# Compute HP polynomial Q_n
# ============================================================
def hp_polynomial(n):
    """Compute the type-II HP polynomial Q_n of degree n.

    Orthogonality conditions:
      Σ_{m=0}^n c_m · ψ₂(m+k) = 0 for k = 0,...,n₁-1
      Σ_{m=0}^n c_m · ψ₃(m+k) = 0 for k = 0,...,n₂-1

    Returns coefficients [c_0, ..., c_n] with c_n = 1 (monic).
    """
    n1, n2 = step_line(n)
    assert n1 + n2 == n

    if n == 0:
        return [QQ(1)]

    # Build the linear system: n equations for n unknowns (c_0,...,c_{n-1})
    # with c_n = 1 fixed
    rows = []
    rhs_vals = []

    # From μ₂: n₁ conditions
    for k in range(n1):
        row = [psi2(m, k) for m in range(n)]
        rows.append(row)
        rhs_vals.append(-psi2(n, k))

    # From μ₃: n₂ conditions
    for k in range(n2):
        row = [psi3(m, k) for m in range(n)]
        rows.append(row)
        rhs_vals.append(-psi3(n, k))

    M = matrix(QQ, rows)
    b = vector(QQ, rhs_vals)

    try:
        x = M.solve_right(b)
        coeffs = list(x) + [QQ(1)]
        return coeffs
    except:
        return None

# ============================================================
# Compute Q_n for n = 0,...,N_MAX
# ============================================================
print("=== Step-line HP polynomials ===")
Q_vals_at_1 = []
Q_polys = []

Rt = PolynomialRing(QQ, 't')
t_var = Rt.gen()

for n in range(N_MAX + 1):
    coeffs = hp_polynomial(n)
    if coeffs is None:
        print("n=%d: SINGULAR" % n)
        Q_vals_at_1.append(None)
        continue

    Q_at_1 = sum(coeffs)
    Q_vals_at_1.append(Q_at_1)
    Q_poly = sum(coeffs[m] * t_var^m for m in range(len(coeffs)))
    Q_polys.append(Q_poly)

    n1, n2 = step_line(n)
    print("n=%2d: (n1=%d,n2=%d), Q_n(1) = %s" % (n, n1, n2, Q_at_1))
    if n <= 5:
        print("       Q_n(t) = %s" % Q_poly)

# ============================================================
# Check: does Q_n satisfy a four-term recurrence?
# ============================================================
print("\n=== Four-term recurrence for Q_n(1) ===")
# Q_{n+1}(1) = alpha_n Q_n(1) + beta_n Q_{n-1}(1) + gamma_n Q_{n-2}(1)

print("Q_n(1) values:")
for n in range(min(10, len(Q_vals_at_1))):
    if Q_vals_at_1[n] is not None:
        print("  Q_%d(1) = %s = %.15f" % (n, Q_vals_at_1[n], float(Q_vals_at_1[n])))

# Try to fit: look at ratios
print("\nRatios Q_{n+1}(1)/Q_n(1):")
for n in range(min(10, len(Q_vals_at_1)-1)):
    if Q_vals_at_1[n] is not None and Q_vals_at_1[n] != 0 and Q_vals_at_1[n+1] is not None:
        print("  n=%d: %.15f" % (n, float(Q_vals_at_1[n+1] / Q_vals_at_1[n])))

# ============================================================
# Compute companion p_n from the HP system
# ============================================================
print("\n=== HP companion for ζ(2) + ζ(3) ===")
# For each n, the HP gives:
# Q_n(1)·f₂(1) - P₂(1) = r₂  (remainder for ζ(2))
# Q_n(1)·f₃(1) - P₃(1) = r₃  (remainder for ζ(3))
#
# But f₂(1) = Σ 1/(k+1)² = ζ(2), f₃(1) = Σ 1/(k+1)³ = ζ(3)
#
# The HP gives SEPARATE approximations to ζ(2) and ζ(3):
# P₂_n(1)/Q_n(1) → ζ(2) and P₃_n(1)/Q_n(1) → ζ(3)
#
# For ζ(2)+ζ(3), we take p_n = P₂_n(1) + P₃_n(1).
#
# Compute P₂ and P₃ from the HP system:
# [t^k](Q_n · f₂) = P₂[k] for k = 0,...,n₁-1
# [t^k](Q_n · f₃) = P₃[k] for k = 0,...,n₂-1

import mpmath
mpmath.mp.dps = 100
L = mpmath.zeta(2) + mpmath.zeta(3)
zeta2 = mpmath.zeta(2)
zeta3 = mpmath.zeta(3)

hp_q = []
hp_p = []
hp_p2 = []
hp_p3 = []
hp_e = []

for n in range(min(12, len(Q_polys))):
    coeffs = list(Q_polys[n])
    while len(coeffs) <= n:
        coeffs.append(QQ(0))

    n1, n2 = step_line(n)

    # P₂(t) coefficients: [t^k](Q_n · f₂) for k = 0,...,n₁-1
    p2_coeffs = []
    for k in range(n1):
        val = sum(coeffs[m] * psi2(m, k) for m in range(n+1))
        p2_coeffs.append(val)
    P2_at_1 = sum(p2_coeffs)

    # P₃(t) coefficients: [t^k](Q_n · f₃) for k = 0,...,n₂-1
    p3_coeffs = []
    for k in range(n2):
        val = sum(coeffs[m] * psi3(m, k) for m in range(n+1))
        p3_coeffs.append(val)
    P3_at_1 = sum(p3_coeffs)

    q_n = sum(coeffs)  # Q_n(1)
    p_n = P2_at_1 + P3_at_1  # combined companion

    hp_q.append(q_n)
    hp_p.append(p_n)
    hp_p2.append(P2_at_1)
    hp_p3.append(P3_at_1)

    # Error: e_n = p_n - (ζ(2)+ζ(3))·q_n
    e_n = mpmath.mpf(p_n) - L * mpmath.mpf(q_n)

    # Also separate errors
    e2 = mpmath.mpf(P2_at_1) - zeta2 * mpmath.mpf(q_n)
    e3 = mpmath.mpf(P3_at_1) - zeta3 * mpmath.mpf(q_n)

    hp_e.append(e_n)

    print("n=%2d: q=%s, p2+p3=%s" % (n, q_n, p_n))
    print("       e = p-(ζ(2)+ζ(3))q = %s" % mpmath.nstr(e_n, 20))
    print("       e₂ = P₂-ζ(2)q = %s" % mpmath.nstr(e2, 15))
    print("       e₃ = P₃-ζ(3)q = %s" % mpmath.nstr(e3, 15))

# ============================================================
# Check error decay and recurrence
# ============================================================
print("\n=== Error decay ===")
for n in range(1, len(hp_e)):
    if abs(hp_e[n-1]) > 0:
        ratio = hp_e[n] / hp_e[n-1]
        print("n=%d: e_n/e_{n-1} = %s" % (n, mpmath.nstr(ratio, 15)))

# Guess recurrence for q_n
print("\n=== Recurrence guess for HP q_n ===")
from ore_algebra import OreAlgebra, guess
Rn = PolynomialRing(QQ, 'nn')
OS = OreAlgebra(Rn, 'Snn')

q_list = hp_q[:12]
for order in [3, 4, 5]:
    try:
        rec = guess(q_list, OS, order=order)
        print("Found order-%d recurrence!" % order)
        for j in range(order+1):
            print("  P_%d: degree %d" % (j, rec[j].degree()))
        print("  P_%d factored: %s" % (order, factor(rec[order])))
        print("  P_0 factored: %s" % factor(rec[0]))
        break
    except:
        print("  order=%d: no relation found" % order)

# Compare P2.7 q_n with HP q_n
print("\n=== Comparison with P2.7 ===")
def A_p27(n):
    return QQ(1024) * (2*n+5)^4 * (2*n+7)^3 * (2*n+9)^3 * (946*n^2+6407*n+10860)
def B_p27(n):
    return QQ(128) * (2*n+7)^3 * (2*n+9)^3 * (104060*n^6 + 1745370*n^5 +
        12145238*n^4 + 44886481*n^3 + 92943995*n^2 + 102256019*n + 46709052)
def C_p27(n):
    return QQ(16) * (n+3)^4 * (2*n+9)^3 * (3784*n^5 + 57792*n^4 +
        351019*n^3 + 1059230*n^2 + 1587211*n + 944620)
def D_p27(n):
    return QQ(1) * (n+3)^4 * (n+4)^6 * (946*n^2 + 4515*n + 5399)

p27_q = [QQ(0)] * 12
p27_q[0] = QQ(-215040420000)
p27_q[1] = QQ(-167282265043404) / QQ(905)
p27_q[2] = QQ(-964185327658080) / QQ(6071)

for n in range(2, 11):
    p27_q[n+1] = B_p27(n)/A_p27(n)*p27_q[n] - C_p27(n-1)/A_p27(n-1)*p27_q[n-1] + D_p27(n-2)/A_p27(n-2)*p27_q[n-2]

print("P2.7 q_n / HP q_n ratio:")
for n in range(min(10, len(hp_q))):
    if hp_q[n] != 0 and p27_q[n] != 0:
        ratio = float(p27_q[n] / hp_q[n])
        print("  n=%d: P2.7/HP = %e" % (n, ratio))
