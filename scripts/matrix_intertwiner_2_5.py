#!/usr/bin/env python3
"""
Problem 2.5: Compute the matrix-level intertwiner V(n) such that
M(n) · V(n+1) = V(n) · S̃(n)

where M(n) is the CMF companion matrix and S̃(n) is the (-16)-twisted
Sym²(Delannoy) companion matrix.

If V(n) is a RATIONAL matrix function of n, then the entries of V
give the Ore intertwiner coefficients.
"""
from fractions import Fraction as F

# Delannoy
D = [F(1), F(3)]
E = [F(0), F(1)]
for n in range(1, 40):
    D.append((F(3*(2*n+1)) * D[n] - F(n) * D[n-1]) / F(n+1))
    E.append((F(3*(2*n+1)) * E[n] - F(n) * E[n-1]) / F(n+1))

# Sym² Delannoy recurrence coefficients (from gauge_from_factorization.py)
# a₀(n) = -(2n+5)(n+1)²
# a₁(n) = (2n+5)(35n²+140n+131)  [verified: 70n³+455n²+980n+655]
# a₂(n) = -(2n+3)(35n²+140n+131) [= -(70n³+385n²+655n+393)]
# a₃(n) = (2n+3)(n+3)²

def a0(n): return -(2*n+5)*(n+1)**2
def a1(n): return 70*n**3 + 455*n**2 + 980*n + 655
def a2(n): return -(70*n**3 + 385*n**2 + 655*n + 393)
def a3(n): return (2*n+3)*(n+3)**2

# Sym² companion matrix: U_{n+1} = [U_n, U_{n+1}, U_{n+2}]^T
# From a₃·U_{n+3} = -a₂·U_{n+2} - a₁·U_{n+1} - a₀·U_n
# Companion: v(n+1) = S(n)·v(n) where v = [U_n, U_{n+1}, U_{n+2}]^T
# S(n) = [[0, 1, 0], [0, 0, 1], [-a₀/a₃, -a₁/a₃, -a₂/a₃]]

def S_sym2(n):
    a3n = F(a3(n))
    return [[F(0), F(1), F(0)],
            [F(0), F(0), F(1)],
            [F(-a0(n), a3(n)), F(-a1(n), a3(n)), F(-a2(n), a3(n))]]

# CMF matrix M(n) — integer entries
def M_int(n):
    m11 = (-2*n-5)*(n+3)**2 * (136*n**4 + 1424*n**3 + 5548*n**2 + 9551*n + 6141)
    m12 = 384*n**6 + 6384*n**5 + 44168*n**4 + 162698*n**3 + 336377*n**2 + 369933*n + 169011
    m13 = -480*n**4 - 4980*n**3 - 19210*n**2 - 32690*n - 20730
    m21 = (n+2)**2*(n+3)**2*(4*n+10)*(48*n**3 + 386*n**2 + 1017*n + 879)
    m22 = (n+2)**2*(-272*n**5 - 3848*n**4 - 21732*n**3 - 61184*n**2 - 85761*n - 47808)
    m23 = (n+2)**2*(320*n**3 + 2540*n**2 + 6610*n + 5640)
    m31 = (-4*n-10)*(n+2)**2*(n+3)**2*(32*n**4 + 302*n**3 + 1037*n**2 + 1530*n + 813)
    m32 = (n+2)**2*(192*n**6 + 2984*n**5 + 19116*n**4 + 64452*n**3 + 120256*n**2 + 117279*n + 46476)
    m33 = (n+2)**2*(-16*n**5 - 408*n**4 - 2912*n**3 - 8884*n**2 - 12254*n - 6240)
    return [[F(m11), F(m12), F(m13)],
            [F(m21), F(m22), F(m23)],
            [F(m31), F(m32), F(m33)]]

def mat_mul(A, B):
    n = len(A)
    m = len(B[0])
    k = len(B)
    return [[sum(A[i][l]*B[l][j] for l in range(k)) for j in range(m)] for i in range(n)]

def mat_inv(M):
    """Invert a 3x3 matrix over Q."""
    a = M[0][0]; b = M[0][1]; c = M[0][2]
    d = M[1][0]; e = M[1][1]; f = M[1][2]
    g = M[2][0]; h = M[2][1]; i = M[2][2]
    det = a*(e*i-f*h) - b*(d*i-f*g) + c*(d*h-e*g)
    if det == 0:
        return None
    inv_det = F(1) / det
    return [[(e*i-f*h)*inv_det, (c*h-b*i)*inv_det, (b*f-c*e)*inv_det],
            [(f*g-d*i)*inv_det, (a*i-c*g)*inv_det, (c*d-a*f)*inv_det],
            [(d*h-e*g)*inv_det, (b*g-a*h)*inv_det, (a*e-b*d)*inv_det]]

# Strategy: compute V(n) = P(n)⁻¹ · V₀ · Q(n) for various V₀
# and find V₀ such that V(n) is rational (polynomial).
#
# P(n) = M(0)·M(1)·...·M(n-1) (CMF fundamental matrix)
# Q(n) = S(0)·S(1)·...·S(n-1) (Sym² fundamental matrix)
#
# The correct V₀ is the one where the exponential growth cancels.

# Compute P(n) and Q(n) for n=0,...,15
print("Computing fundamental matrices...")
P = [[[F(1),F(0),F(0)],[F(0),F(1),F(0)],[F(0),F(0),F(1)]]]
Q = [[[F(1),F(0),F(0)],[F(0),F(1),F(0)],[F(0),F(0),F(1)]]]

for n in range(15):
    P.append(mat_mul(P[-1], M_int(n)))
    Q.append(mat_mul(Q[-1], S_sym2(n)))
    if n < 3:
        print(f"  P({n+1})[0][0] = {P[-1][0][0]}")
        print(f"  Q({n+1})[0][0] = {Q[-1][0][0]}")

# At each n: V(n) = P(n)⁻¹ · V₀ · Q(n)
# For V₀ = I: V(n) = P(n)⁻¹ · Q(n)
# This should show the structure.

print("\n=== V(n) = P(n)⁻¹ · Q(n) (with V₀=I) ===")
print("(Looking for rational matrix entries)")

for n in range(8):
    P_inv = mat_inv(P[n])
    if P_inv is None:
        print(f"n={n}: P singular")
        continue
    V = mat_mul(P_inv, Q[n])
    # Print the matrix entries (as floats for readability)
    print(f"\nV({n}):")
    for i in range(3):
        row = [float(V[i][j]) for j in range(3)]
        print(f"  [{row[0]:>15.6e}, {row[1]:>15.6e}, {row[2]:>15.6e}]")

# The entries of V(n) should grow/shrink exponentially if V₀ is wrong.
# Let's check the growth pattern.
print("\n=== Growth of V(n)[0][0] ===")
v00_vals = []
for n in range(12):
    P_inv = mat_inv(P[n])
    if P_inv is None:
        continue
    V = mat_mul(P_inv, Q[n])
    v00_vals.append(float(V[0][0]))
    print(f"  n={n}: V[0][0] = {v00_vals[-1]:.10e}")

# Check ratio of consecutive V[0][0]
print("\n=== Ratio V(n+1)[0][0] / V(n)[0][0] ===")
for n in range(len(v00_vals)-1):
    if v00_vals[n] != 0:
        ratio = v00_vals[n+1] / v00_vals[n]
        print(f"  n={n}: ratio = {ratio:.10f}")

# The ratio should converge to a Poincaré root ratio:
# λ_CMF / λ_Sym² = -16 (for the dominant modes)
# If ratio → -16, then V(n) ~ (-16)^n, and V₀ = I is wrong.
# We need to find V₀ such that the exponential cancels.

# Alternative approach: work with the GAUGED companion.
# Define M̃(n) = M(n) / (-16·n^7) (normalized CMF companion)
# and S̃(n) = S(n) (Sym² companion)
# Then V(n) should satisfy M̃(n)·V(n+1) = V(n)·S̃(n)·(-16·n^7)^{-1}

# Actually, the right gauge is: if we define the twisted Sym²
# matrix as T(n) = (-16)·(n+1)^7 / [a₃(n)] · S(n), then the
# Poincaré roots of T match those of M.

# Let me try: compute V(n) = P(n)⁻¹ · V₀ · Q_twisted(n)
# where Q_twisted(n) accounts for the -16·n^7 twist.

# But the twist is on the EIGENVALUES, not a uniform scalar.
# The Poincaré roots of M are -16, -16·(17±12√2).
# The Poincaré roots of S are 1, 17±12√2.
# The ratio is -16 for ALL three roots.

# So: M(n) ≈ -16·n^7 · P · diag(1, 17+12√2, 17-12√2) · P⁻¹
#     S(n) ≈ P' · diag(1, 17+12√2, 17-12√2) · P'⁻¹

# And the twist factor is a COMMON scalar: λ_M = -16·n^7 · λ_S

# This means M(n) ≈ (-16·n^7) · S̃(n) where S̃ is basis-changed S.
# So V(n+1) = M(n)⁻¹ · V(n) · S(n) ≈ 1/(-16·n^7) · S̃(n)⁻¹ · V(n) · S(n)

# For large n with the right basis, V converges to the identity (up to formal exponents).

# Let me compute the GAUGED intertwiner:
# V_g(n) = Λ(n)⁻¹ · P(n)⁻¹ · V₀ · Q(n)
# where Λ(n) = ∏_{k=0}^{n-1} (-16·(k+1)^7·...) is the gauge factor.
# Wait, I need to be more precise. The CMF eigenvalues at step n are
# approximately -16n^7 times the Sym² eigenvalues.

# Let me instead compute the per-step gauge:
# d(n) = det(M(n)) / det(S(n))^{1/1} ... no, this doesn't work for 3x3.

# Let me use the Wronskian approach:
# det(P(n)) = ∏_{k=0}^{n-1} det(M(k))
# det(Q(n)) = ∏_{k=0}^{n-1} det(S(k))
# The ratio det(P)/det(Q) is the "scalar gauge" factor.

print("\n\n=== Determinant-based gauge ===")
det_M_prod = F(1)
det_S_prod = F(1)
for n in range(12):
    M = M_int(n)
    S = S_sym2(n)
    det_M = M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1]) - M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0]) + M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0])
    det_S = S[0][0]*(S[1][1]*S[2][2]-S[1][2]*S[2][1]) - S[0][1]*(S[1][0]*S[2][2]-S[1][2]*S[2][0]) + S[0][2]*(S[1][0]*S[2][1]-S[1][1]*S[2][0])
    det_M_prod *= det_M
    det_S_prod *= det_S
    if n < 8:
        ratio = det_M / det_S
        print(f"n={n}: det(M)/det(S) = {float(ratio):.6e}")
        print(f"  det(M) = {float(det_M):.6e}")
        print(f"  det(S) = {float(det_S):.6e}")

# The per-step ratio det(M(n))/det(S(n)) should be related to (-16)^3 · n^21 · ...
# (since M ≈ -16·n^7 · S̃, and det scales as (-16)^3 · n^21)
print("\n=== det(M(n)) / [(-16)³ · det(S(n))] ===")
for n in range(12):
    M = M_int(n)
    S = S_sym2(n)
    det_M = M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1]) - M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0]) + M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0])
    det_S = S[0][0]*(S[1][1]*S[2][2]-S[1][2]*S[2][1]) - S[0][1]*(S[1][0]*S[2][2]-S[1][2]*S[2][0]) + S[0][2]*(S[1][0]*S[2][1]-S[1][1]*S[2][0])
    ratio = det_M / (F(-4096) * det_S)
    print(f"n={n}: {float(ratio):.10e}  (expect ~ n^21 · ...)")

# Now: the gauged intertwiner.
# Define Λ(n) = ∏_{k=0}^{n-1} (-4096)^(1/3) · [det_ratio(k)]^(1/3) · Identity
# This doesn't quite work because the gauge varies per mode.

# ALTERNATIVE: Directly compute V(n) at each step using the recurrence
# V(n+1) = M(n)⁻¹ · V(n) · S(n)
# and normalize V(n) to have unit determinant or fixed first column.

print("\n\n=== Computing V(n) via recurrence ===")
print("V(n+1) = M(n)⁻¹ · V(n) · S(n)")
print("Starting from V(0) = I")

V = [[F(1),F(0),F(0)],[F(0),F(1),F(0)],[F(0),F(0),F(1)]]
for n in range(10):
    M_inv = mat_inv(M_int(n))
    S = S_sym2(n)
    V_new = mat_mul(mat_mul(M_inv, V), S)
    V = V_new
    # Print first row
    print(f"\nV({n+1}) first row: [{float(V[0][0]):.6e}, {float(V[0][1]):.6e}, {float(V[0][2]):.6e}]")
    print(f"  ndigits: [{len(str(abs(V[0][0].numerator)))}, {len(str(abs(V[0][1].numerator)))}, {len(str(abs(V[0][2].numerator)))}]")

# The entries grow because of the gauge factor.
# Let me normalize by the cumulative gauge: Λ(n) = ∏ det(M(k))/det(S(k))
# V_norm(n) = Λ(n)^{-1/3} · V(n)

# Actually, let's normalize V(n) at each step by the det(M)^{1/3} factor.
# det(V(n+1)) = det(M(n)⁻¹) · det(V(n)) · det(S(n))
# = det(V(n)) · det(S(n)) / det(M(n))

print("\n\n=== Normalized V(n) ===")
print("V_norm(n) = V(n) · ∏_{k=0}^{n-1} [det(M(k))/det(S(k))]")

V = [[F(1),F(0),F(0)],[F(0),F(1),F(0)],[F(0),F(0),F(1)]]
gauge = F(1)
for n in range(10):
    M_inv = mat_inv(M_int(n))
    S = S_sym2(n)
    V_new = mat_mul(mat_mul(M_inv, V), S)
    V = V_new

    M = M_int(n)
    det_M = M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1]) - M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0]) + M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0])
    Sn = S_sym2(n)
    det_S = Sn[0][0]*(Sn[1][1]*Sn[2][2]-Sn[1][2]*Sn[2][1]) - Sn[0][1]*(Sn[1][0]*Sn[2][2]-Sn[1][2]*Sn[2][0]) + Sn[0][2]*(Sn[1][0]*Sn[2][1]-Sn[1][1]*Sn[2][0])
    gauge *= det_M / det_S

    # V_norm = V * gauge
    V_norm = [[V[i][j] * gauge for j in range(3)] for i in range(3)]
    print(f"\nV_norm({n+1}) first row: [{float(V_norm[0][0]):.10e}, {float(V_norm[0][1]):.10e}, {float(V_norm[0][2]):.10e}]")

print("\nDone.")
