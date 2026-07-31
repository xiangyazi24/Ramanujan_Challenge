#!/usr/bin/env python3
"""
Compute the P2.7 adjoint recurrence, the adjoint recessive solution w⁰,
and the Lagrange bilinear form J(w⁰, e; 0).

Goal: verify J(w⁰, p; 0) / J(w⁰, q; 0) = ζ(2) + ζ(3) to high precision.
"""
from mpmath import mp, mpf, pi, zeta, log, matrix, power, fac

mp.dps = 300  # 300 digits

# ===== P2.7 recurrence coefficients =====
def A(n):
    n = mpf(n)
    return mpf(1024) * (2*n+5)**4 * (2*n+7)**3 * (2*n+9)**3 * (946*n**2 + 6407*n + 10860)

def B(n):
    n = mpf(n)
    return mpf(128) * (2*n+7)**3 * (2*n+9)**3 * (104060*n**6 + 1745370*n**5 + 12145238*n**4 + 44886481*n**3 + 92943995*n**2 + 102256019*n + 46709052)

def C(n):
    n = mpf(n)
    return mpf(16) * (n+3)**4 * (2*n+9)**3 * (3784*n**5 + 57792*n**4 + 351019*n**3 + 1059230*n**2 + 1587211*n + 944620)

def D(n):
    n = mpf(n)
    return (n+3)**4 * (n+4)**6 * (946*n**2 + 4515*n + 5399)

# Standard form: A(n) u_{n+1} = B(n) u_n - C(n-1) u_{n-1} + D(n-2) u_{n-2}
# Rewrite as order-3 operator (shift by 2):
# A(n+2) u_{n+3} - B(n+2) u_{n+2} + C(n+1) u_{n+1} - D(n) u_n = 0
# So: p3(n) = A(n+2), p2(n) = -B(n+2), p1(n) = C(n+1), p0(n) = -D(n)

def p3(n): return A(n+2)
def p2(n): return -B(n+2)
def p1(n): return C(n+1)
def p0(n): return -D(n)

# ===== Forward solutions =====
from fractions import Fraction as F

# Use exact arithmetic for initial values
q0 = mpf(-215040420000)
q1 = mpf(-167282265043404) / mpf(905)
q2 = mpf(-964185327658080) / mpf(6071)

p0_val = mpf(-612218384750)
p1_val = mpf(-9525021973931919) / mpf(18100)
p2_val = mpf(-29561828382772029) / mpf(65380)

N = 500  # Compute to n=500 for high precision Miller's algorithm

# Forward q_n
q = [mpf(0)] * (N + 3)
q[0], q[1], q[2] = q0, q1, q2
for n in range(2, N):
    q[n+1] = (B(n) * q[n] - C(n-1) * q[n-1] + D(n-2) * q[n-2]) / A(n)

# Forward p_n
p = [mpf(0)] * (N + 3)
p[0], p[1], p[2] = p0_val, p1_val, p2_val
for n in range(2, N):
    p[n+1] = (B(n) * p[n] - C(n-1) * p[n-1] + D(n-2) * p[n-2]) / A(n)

print("Forward solutions computed.")
print(f"q[0] = {q[0]}")
print(f"q[1] = {q[1]}")
print(f"|q[50]| ~ 10^{float(log(abs(q[50]))/log(10)):.1f}")

L_val = zeta(2) + zeta(3)
print(f"\nL = ζ(2)+ζ(3) = {L_val}")

# Check error decay
for n in [5, 10, 20, 50]:
    e_n = p[n] - L_val * q[n]
    if q[n] != 0:
        ratio = abs(e_n / q[n])
        print(f"  |e_{n}/q_{n}| = {float(ratio):.6e}")

# ===== Adjoint recurrence =====
# Forward: p3(n) u_{n+3} + p2(n) u_{n+2} + p1(n) u_{n+1} + p0(n) u_n = 0
# Adjoint: p0(n) w_{n+3} + p1(n-1) w_{n+2} + p2(n-2) w_{n+1} + p3(n-3) w_n = 0
# i.e., w_{n+3} = [-p1(n-1) w_{n+2} - p2(n-2) w_{n+1} - p3(n-3) w_n] / p0(n)

# The adjoint Poincaré roots are 1/μ_j.
# 1/μ₀ ≈ 0.0182 (most recessive)
# 1/|μ±| ≈ 14.84 (dominant for adjoint)

# Miller's algorithm: to find the adjoint recessive solution w⁰,
# start at large M with w_M = 1, w_{M-1} = w_{M-2} = 0 and run BACKWARD.

# The BACKWARD adjoint recurrence:
# p3(n-3) w_n = -p0(n) w_{n+3} - p1(n-1) w_{n+2} - p2(n-2) w_{n+1}
# So: w_n = [-p0(n) w_{n+3} - p1(n-1) w_{n+2} - p2(n-2) w_{n+1}] / p3(n-3)

M = N  # Start at n=M
w = [mpf(0)] * (M + 4)
w[M] = mpf(1)
w[M-1] = mpf(0)
w[M-2] = mpf(0)

print("\nRunning Miller's algorithm for adjoint recessive solution...")
for n in range(M-3, -1, -1):
    # w_n = [-p0(n) w_{n+3} - p1(n-1) w_{n+2} - p2(n-2) w_{n+1}] / p3(n-3)
    w[n] = (-p0(n) * w[n+3] - p1(n-1) * w[n+2] - p2(n-2) * w[n+1]) / p3(n-3)

# Normalize: we'll normalize later based on the bilinear form
print(f"w[0] = {w[0]}")
print(f"w[1] = {w[1]}")
print(f"w[2] = {w[2]}")
print(f"w[3] = {w[3]}")

# Check growth rate
for n in [10, 20, 50, 100]:
    if w[n] != 0 and w[n-1] != 0:
        ratio = abs(w[n] / w[n-1])
        print(f"  |w[{n}]/w[{n-1}]| = {float(ratio):.10f}")

# ===== Lagrange bilinear form =====
# For the recurrence L·u = 0 with L = p3(n)S³ + p2(n)S² + p1(n)S + p0(n),
# the Lagrange identity gives a constant bilinear form J(w,u;n).
#
# The general formula for order-3:
# J(w,u;n) = p3(n)[w_{n+2} u_{n+2} - w_{n+2} u_{n+2}]  -- no, need to derive properly
#
# Actually, for the standard Green-Lagrange identity:
# Σ_{n=a}^{b-1} [w_n (Lu)_n] = boundary terms
#
# For a 3rd order operator L = Σ_{i=0}^3 p_i(n) S^i:
# The adjoint is L† = Σ_{i=0}^3 S^{-i} p_i(n) = Σ_{i=0}^3 p_i(n+i) S^{-i}
# (after multiplying by S³: Σ_{i=0}^3 p_i(n+i) S^{3-i})
# Wait, let me be more careful.
#
# L = p_0(n) + p_1(n)S + p_2(n)S² + p_3(n)S³
# L† = Σ S^{-i} · p_i(n) = p_0(n) + S^{-1}p_1(n) + S^{-2}p_2(n) + S^{-3}p_3(n)
#    = p_0(n) + p_1(n+1)S^{-1} + p_2(n+2)S^{-2} + p_3(n+3)S^{-3}
# (using S^{-k}f(n) = f(n+k)S^{-k})
# Wait no: S^{-1} p(n) f(n) = p(n-1) f(n-1) = S^{-1}[p(n)f(n)]
# Hmm, the adjoint depends on whether we use the formal adjoint or the
# Lagrange adjoint.
#
# The Lagrange adjoint: for a bilinear form ⟨w, Lu⟩ = ⟨L†w, u⟩ + boundary,
# where ⟨f,g⟩ = Σ_n f_n g_n.
#
# For S: ⟨w, Su⟩ = Σ w_n u_{n+1} = Σ w_{n-1} u_n = ⟨S^{-1}w, u⟩
# Actually S^{-1}w means (S^{-1}w)_n = w_{n-1}
# So: ⟨w, p(n)S^k u⟩ = Σ w_n p(n) u_{n+k} = Σ w_{n-k} p(n-k) u_n = ⟨p(n-k) S^{-k} w, u⟩
# Hence L† = Σ_i p_i(n-i) S^{-i}
# Multiply by S³: S³L† = Σ_i p_i(n-i) S^{3-i}
#                       = p_0(n) S³ + p_1(n-1) S² + p_2(n-2) S + p_3(n-3)
#
# So the adjoint forward recurrence is:
# p_0(n) w_{n+3} + p_1(n-1) w_{n+2} + p_2(n-2) w_{n+1} + p_3(n-3) w_n = 0
# This matches what I had above. ✓
#
# Now the Lagrange identity. For the discrete Green formula:
# w_n (Lu)_n = w_n [p_3(n)u_{n+3} + p_2(n)u_{n+2} + p_1(n)u_{n+1} + p_0(n)u_n]
#
# We want to show: Σ_{n=a}^{b} w_n (Lu)_n = [J(w,u;b+1) - J(w,u;a)] + Σ_n u_n (L†w)_n
# where J is the boundary form.
#
# By Abel summation (summation by parts):
# For each term p_i(n) w_n u_{n+i}, we shift:
# Σ_n p_i(n) w_n u_{n+i} = Σ_m p_i(m-i) w_{m-i} u_m + boundary
#
# The boundary form J(w,u;n) for order 3 is:
# J(w,u;n) = p_3(n-1) w_{n-1} u_{n+2} + [p_3(n-2) w_{n-2} + p_2(n-1) w_{n-1}] u_{n+1}
#           + ... (need careful derivation)
#
# Actually, let me use a simpler approach. Since J is constant when Lu=0 and L†w=0,
# I can just VERIFY numerically that a specific bilinear expression is constant.
#
# For order 3, J has the form:
# J(w,u;n) = Σ_{i,j with i+j≤2} A_{ij}(n) w_{n+i} u_{n+j}
# or equivalently it involves w_n, w_{n+1}, w_{n+2} and u_n, u_{n+1}, u_{n+2}.
#
# The simplest approach: compute J(w,u;n) = Σ bilinear terms such that ΔJ = 0
# when Lu=0 and L†w=0.

# Let me use the WRAPPING APPROACH:
# For a system y_{n+1} = M_n y_n where y = (u_n, u_{n+1}, u_{n+2})
# and the adjoint z_{n+1} = (M_n^T)^{-1} z_n,
# the pairing J = z^T y is constant.
#
# But for our recurrence, the companion matrix at step n is:
# u_{n+3} = -(p_2/p_3)(n) u_{n+2} - (p_1/p_3)(n) u_{n+1} - (p_0/p_3)(n) u_n
# So M_n = [[0, 1, 0], [0, 0, 1], [-p0/p3, -p1/p3, -p2/p3]](n)
#
# And the invariant is: z_n^T M_{n-1} M_{n-2} ... M_0 y_0 = z^T y is constant
# more simply: if y_{n+1} = M_n y_n and z_{n+1} = (M_n^T)^{-1} z_n, then z_n^T y_n = const.
#
# But that's just saying (M_n^{-T} z_n)^T (M_n y_n) = z_n^T M_n^{-1} M_n y_n = z_n^T y_n.
# Wait, z_{n+1}^T y_{n+1} = z_n^T M_n^{-1} · M_n y_n = z_n^T y_n. ✓

# The companion matrix: y_n = (u_n, u_{n+1}, u_{n+2})
# y_{n+1} = (u_{n+1}, u_{n+2}, u_{n+3}) = (u_{n+1}, u_{n+2}, -p0/p3·u_n - p1/p3·u_{n+1} - p2/p3·u_{n+2})

# For the adjoint: the "companion dual" z must satisfy z^T y = const.
# If z = (z0, z1, z2) at step n, then z^T y = z0 u_n + z1 u_{n+1} + z2 u_{n+2} = const
# This IS the Lagrange bilinear form:
# J(w,u;n) = z0_n u_n + z1_n u_{n+1} + z2_n u_{n+2}
# where (z0, z1, z2) is the adjoint state vector.

# But z_n is NOT the same as (w_n, w_{n+1}, w_{n+2}). The adjoint companion
# system has a different structure.

# Let me work it out. We need z_{n+1}^T y_{n+1} = z_n^T y_n.
# y_{n+1} = M_n y_n implies z_n^T y_n = z_{n+1}^T M_n y_n for all y_n.
# So z_{n+1}^T M_n = z_n^T, i.e., z_{n+1} = M_n^{-T} z_n.

# M_n = [[0, 1, 0], [0, 0, 1], [-p0(n)/p3(n), -p1(n)/p3(n), -p2(n)/p3(n)]]
# det(M_n) = p0(n)/p3(n) (cofactor expansion along first row... let me compute)
# Actually: det = 0·(0·(-p2/p3) - 1·(-p1/p3)) - 1·(0·(-p2/p3) - 1·(-p0/p3)) + 0·...
# = -1·(p0/p3) = -p0(n)/p3(n)

# So det(M_n) = -p0(n)/p3(n) = D(n)/A(n+2) (since p0 = -D, p3 = A(n+2))
# (Using p0(n) = -D(n), p3(n) = A(n+2))

# M_n^{-1} = (1/det) adj(M_n)
# M_n^{-T} = (1/det) adj(M_n)^T

# This is getting messy. Let me just use the NUMERICAL approach:
# Compute J(w,u;n) = w_n α_n u_n + w_n β_n u_{n+1} + ... numerically.
#
# Actually, the simplest thing: verify that w_n · u_n is NOT constant,
# but some specific combination IS.
#
# For a 3rd order operator, the conserved quantity is the "discrete Wronskian":
# J = p3(n) det [[w_n, w_{n+1}], [u_n, u_{n+1}]] + something
#
# No wait. Let me just compute it directly.
#
# The bilinear form for a 3rd-order operator can be derived as follows.
# Write the recurrence as: p3 u_{n+3} + p2 u_{n+2} + p1 u_{n+1} + p0 u_n = 0
# Adjoint: p0(n) w_{n+3} + p1(n-1) w_{n+2} + p2(n-2) w_{n+1} + p3(n-3) w_n = 0
#
# Define:
# J(w,u;n) = p3(n-1) [w_n u_{n+2} - w_{n+2} u_n]
#           + p3(n-2) [w_n u_{n+1} - w_{n+1} u_n] · ???
#
# Actually, for order 2 the form is just p2(n-1)(w_n u_{n+1} - w_{n+1} u_n).
# For order 3, it's more complex.
#
# Let me use a KNOWN formula. For the operator L = Σ_{k=0}^r a_k(n) S^k,
# the Lagrange identity gives:
# w·(Lu) - (L†w)·u = Δ J(w,u;n)
# where Δf(n) = f(n+1) - f(n).
#
# For order 3: J(w,u;n) = Σ_{j=1}^{3} Σ_{i=0}^{j-1} (-1)^i a_j(n-j+i) w_{n-j+1+i} u_n ???
#
# This is getting confusing. Let me just compute J NUMERICALLY by finding
# a linear combination of products w_{n+i} u_{n+j} that is constant.

# We know J has the form:
# J(w,u;n) = Σ c_{ij}(n) w_{n+i} u_{n+j}  for 0 ≤ i,j and i+j ≤ 2 (or similar)
# Actually, J should involve at most w_{n}, w_{n+1} and u_{n}, u_{n+1}, u_{n+2}
# (or some subset of these).

# For order r=3: J involves products with total shift ≤ r-1 = 2.
# The terms are: w_n u_n, w_n u_{n+1}, w_n u_{n+2}, w_{n+1} u_n, w_{n+1} u_{n+1}.
# (Five terms with polynomial coefficients.)

# Let me compute J by requiring J(w,q;n) = const for n=0,1,...,10.
# This gives linear equations on the coefficient functions c_{ij}(n).

# SIMPLEST APPROACH: For our SPECIFIC w and u=q, just compute
# different linear combinations and find one that's constant.

# Actually, the textbook formula is available.
# For L = a₃ S³ + a₂ S² + a₁ S + a₀ and L† = a₀(n) S³ + a₁(n-1) S² + a₂(n-2) S + a₃(n-3):
#
# J(w,u;n) = a₃(n) w_{n+1} u_{n+2} + [a₃(n-1) + a₂(n)] w_n u_{n+1}
#           - a₃(n) w_{n+2} u_{n+1} - ...
# hmm this is wrong.
#
# Let me try the approach from Elaydi "An Introduction to Difference Equations":
# For the operator L = Σ pₖ(n) Eⁿ where E is the shift:
# The bilinear concomitant is:
# B[w,u](n) = Σⱼ₌₁ʳ Σₖ₌₀ʲ⁻¹ (-1)ᵏ Δᵏ[pⱼ(n) wₙ] · Eʲ⁻¹⁻ᵏ uₙ
#
# For order 3 (r=3, using p₀,p₁,p₂,p₃):
# j=1: k=0: p₁(n) w_n · u_n
# j=2: k=0: p₂(n) w_n · u_{n+1}
#       k=1: -Δ[p₂(n) w_n] · u_n = -(p₂(n+1)w_{n+1} - p₂(n)w_n) · u_n
# j=3: k=0: p₃(n) w_n · u_{n+2}
#       k=1: -Δ[p₃(n)w_n] · u_{n+1} = -(p₃(n+1)w_{n+1} - p₃(n)w_n) · u_{n+1}
#       k=2: Δ²[p₃(n)w_n] · u_n = [p₃(n+2)w_{n+2} - 2p₃(n+1)w_{n+1} + p₃(n)w_n] · u_n
#
# Collecting:
# B[w,u](n) = p₁(n) w_n u_n
#            + p₂(n) w_n u_{n+1} - (p₂(n+1) w_{n+1} - p₂(n) w_n) u_n
#            + p₃(n) w_n u_{n+2} - (p₃(n+1) w_{n+1} - p₃(n) w_n) u_{n+1}
#            + (p₃(n+2) w_{n+2} - 2p₃(n+1) w_{n+1} + p₃(n) w_n) u_n
#
# = w_n u_n [p₁(n) + p₂(n) + p₃(n) + p₃(n)]  -- collecting w_n u_n terms
# Hmm wait, let me be more careful:
#
# Term w_n u_n: p₁(n) - (-p₂(n)) + p₃(n) from j=3,k=2 part
# Wait I need to be more systematic.
#
# From j=1,k=0: +p₁ w u
# From j=2,k=0: +p₂ w u₊₁
# From j=2,k=1: -p₂(n+1)w₊₁ u + p₂(n) w u
# From j=3,k=0: +p₃ w u₊₂
# From j=3,k=1: -p₃(n+1)w₊₁ u₊₁ + p₃(n) w u₊₁
# From j=3,k=2: +p₃(n+2)w₊₂ u - 2p₃(n+1)w₊₁ u + p₃(n) w u
#
# Collecting by (w_{n+i}, u_{n+j}):
# w u:     p₁(n) + p₂(n) + p₃(n)                                    [coeff of w_n u_n]
# w u₊₁:  p₂(n) + p₃(n)                                             [coeff of w_n u_{n+1}]
# w u₊₂:  p₃(n)                                                       [coeff of w_n u_{n+2}]
# w₊₁ u:  -p₂(n+1) - 2p₃(n+1)                                       [coeff of w_{n+1} u_n]
# w₊₁ u₊₁: -p₃(n+1)                                                  [coeff of w_{n+1} u_{n+1}]
# w₊₂ u:  +p₃(n+2)                                                    [coeff of w_{n+2} u_n]
#
# So: B[w,u](n) = [p₁(n)+p₂(n)+p₃(n)] w_n u_n
#               + [p₂(n)+p₃(n)] w_n u_{n+1}
#               + p₃(n) w_n u_{n+2}
#               + [-p₂(n+1)-2p₃(n+1)] w_{n+1} u_n
#               + [-p₃(n+1)] w_{n+1} u_{n+1}
#               + p₃(n+2) w_{n+2} u_n

def J_form(w, u, n):
    """Lagrange bilinear concomitant for our order-3 operator."""
    p3n = p3(n); p2n = p2(n); p1n = p1(n)
    p3n1 = p3(n+1); p2n1 = p2(n+1)
    p3n2 = p3(n+2)

    wn = w[n]; wn1 = w[n+1]; wn2 = w[n+2]
    un = u[n]; un1 = u[n+1]; un2 = u[n+2]

    return ((p1n + p2n + p3n) * wn * un
          + (p2n + p3n) * wn * un1
          + p3n * wn * un2
          + (-p2n1 - 2*p3n1) * wn1 * un
          + (-p3n1) * wn1 * un1
          + p3n2 * wn2 * un)

# Verify constancy of J(w, q; n)
print("\n===== Lagrange bilinear form J(w, q; n) =====")
J_vals = []
for n in range(0, 20):
    Jn = J_form(w, q, n)
    J_vals.append(Jn)
    if n < 5 or n >= 18:
        print(f"  J(w,q;{n}) = {Jn}")

# Check constancy
if len(J_vals) >= 2:
    max_var = max(abs(J_vals[i] - J_vals[0]) for i in range(1, len(J_vals)))
    print(f"  Max variation from J(0): {float(max_var):.6e}")
    rel_var = float(max_var / abs(J_vals[0])) if J_vals[0] != 0 else float('inf')
    print(f"  Relative variation: {rel_var:.6e}")

# Similarly for J(w, p; n)
print("\n===== Lagrange bilinear form J(w, p; n) =====")
Jp_vals = []
for n in range(0, 20):
    Jn = J_form(w, p, n)
    Jp_vals.append(Jn)
    if n < 5 or n >= 18:
        print(f"  J(w,p;{n}) = {Jn}")

if len(Jp_vals) >= 2:
    max_var = max(abs(Jp_vals[i] - Jp_vals[0]) for i in range(1, len(Jp_vals)))
    print(f"  Max variation from J(0): {float(max_var):.6e}")

# ===== The key test: J(w,p;0) / J(w,q;0) = ζ(2)+ζ(3) ? =====
print("\n" + "="*60)
print("KEY TEST: J(w,p;0) / J(w,q;0) = ζ(2)+ζ(3) ?")
print("="*60)

Jwp = J_form(w, p, 0)
Jwq = J_form(w, q, 0)

ratio = Jwp / Jwq
target = L_val

print(f"J(w,p;0) = {Jwp}")
print(f"J(w,q;0) = {Jwq}")
print(f"Ratio = {ratio}")
print(f"ζ(2)+ζ(3) = {target}")
print(f"Difference = {float(ratio - target):.6e}")
print(f"Relative error = {float(abs(ratio - target) / target):.6e}")

# Also compute J(w,e;0) directly
Jwe = Jwp - L_val * Jwq
print(f"\nJ(w,e;0) = J(w,p;0) - L·J(w,q;0) = {float(Jwe):.6e}")
print(f"|J(w,e;0)| / |J(w,q;0)| = {float(abs(Jwe/Jwq)):.6e}")

# ===== Alternative: compute at n=0 using only w_0, w_1, w_2 =====
# The linear functional α₀ u₀ + α₁ u₁ + α₂ u₂ that extracts c₀
# From the form: J(w,u;0) = c₀₀ w_0 u_0 + c₀₁ w_0 u_1 + c₀₂ w_0 u_2
#                           + c₁₀ w_1 u_0 + c₁₁ w_1 u_1 + c₂₀ w_2 u_0
print("\n===== Coefficient structure at n=0 =====")
print(f"p3(0) = A(2) = {A(2)}")
print(f"p2(0) = -B(2) = {-B(2)}")
print(f"p1(0) = C(1) = {C(1)}")
print(f"p0(0) = -D(0) = {-D(0)}")

# The α_i coefficients:
# α_0 = (p1(0)+p2(0)+p3(0))·w_0 + (-p2(1)-2p3(1))·w_1 + p3(2)·w_2
# α_1 = (p2(0)+p3(0))·w_0 + (-p3(1))·w_1
# α_2 = p3(0)·w_0

alpha_2 = p3(0) * w[0]
alpha_1 = (p2(0) + p3(0)) * w[0] + (-p3(1)) * w[1]
alpha_0 = (p1(0) + p2(0) + p3(0)) * w[0] + (-p2(1) - 2*p3(1)) * w[1] + p3(2) * w[2]

J_check = alpha_0 * q[0] + alpha_1 * q[1] + alpha_2 * q[2]
print(f"\nJ(w,q;0) via alphas = {J_check}")
print(f"Direct J(w,q;0)     = {Jwq}")
print(f"Match: {float(abs(J_check - Jwq) / abs(Jwq)):.6e}")

# ===== Growth analysis of w =====
print("\n===== Adjoint solution growth analysis =====")
mu0 = mpf('54.96')  # approximate
for n in [10, 20, 50, 100, 200]:
    if w[n] != 0:
        growth = abs(w[n]) ** (mpf(1)/n)
        print(f"  |w[{n}]|^(1/{n}) = {float(growth):.10f}  (expect 1/μ₀ ≈ {float(1/mu0):.10f})")
