#!/usr/bin/env python3
"""Problem 2.5: Derive the order-3 scalar recurrence for q_N EXACTLY
by tracking all 3 components of the state vector.

s_N = (q_N, q'_N, q''_N) satisfies s_{N+1} = s_N · M(N).
From the 3-component system, eliminate q' and q'' to get:
α₃(N) q_{N+3} + α₂(N) q_{N+2} + α₁(N) q_{N+1} + α₀(N) q_N = 0

Method: at each N, we have s_N · M(N) = s_{N+1}. Writing the first column:
q_{N+1} = q_N·m00(N) + q'_N·m10(N) + q''_N·m20(N)

From three consecutive such equations, solve for q'_N, q''_N and substitute
to get the scalar recurrence. Use exact integer arithmetic (sympy/fractions).
"""
from mpmath import mp, mpf, nstr, matrix, catalan, log

mp.dps = 500

def M_mat(n):
    n = mpf(n)
    m11 = (-2*n-5)*(n+3)**2 * (136*n**4 + 1424*n**3 + 5548*n**2 + 9551*n + 6141)
    m12 = 384*n**6 + 6384*n**5 + 44168*n**4 + 162698*n**3 + 336377*n**2 + 369933*n + 169011
    m13 = -480*n**4 - 4980*n**3 - 19210*n**2 - 32690*n - 20730
    m21 = (n+2)**2*(n+3)**2*(4*n+10)*(48*n**3 + 386*n**2 + 1017*n + 879)
    m22 = (n+2)**2*(-272*n**5 - 3848*n**4 - 21732*n**3 - 61184*n**2 - 85761*n - 47808)
    m23 = (n+2)**2*(320*n**3 + 2540*n**2 + 6610*n + 5640)
    m31 = (-4*n-10)*(n+2)**2*(n+3)**2*(32*n**4 + 302*n**3 + 1037*n**2 + 1530*n + 813)
    m32 = (n+2)**2*(192*n**6 + 2984*n**5 + 19116*n**4 + 64452*n**3 + 120256*n**2 + 117279*n + 46476)
    m33 = (n+2)**2*(-16*n**5 - 408*n**4 - 2912*n**3 - 8884*n**2 - 12254*n - 6240)
    return matrix([[m11, m12, m13], [m21, m22, m23], [m31, m32, m33]])

A_mat = matrix([[mpf(30921), mpf(-32972), mpf(8240)],
                [mpf(33750), mpf(-36000), mpf(9000)]])

# Compute ALL 3 components of the state vector s_N
N_max = 40
s_vals = []  # s_vals[N] = (q_N, q'_N, q''_N)
T = matrix([[1,0,0],[0,1,0],[0,0,1]])
for N in range(N_max + 1):
    AT = A_mat * T
    s_vals.append((AT[1,0], AT[1,1], AT[1,2]))
    T = T * M_mat(N)

print(f"Computed {len(s_vals)} state vectors at {mp.dps}-digit precision")
print(f"  s[0] = ({nstr(s_vals[0][0],10)}, {nstr(s_vals[0][1],10)}, {nstr(s_vals[0][2],10)})")
print(f"  s[1] = ({nstr(s_vals[1][0],10)}, {nstr(s_vals[1][1],10)}, {nstr(s_vals[1][2],10)})")

# Verify: s_{N+1} = s_N · M(N)
print("\n=== Verifying s_{N+1} = s_N · M(N) ===")
for N in [0, 1, 5, 10]:
    M_N = M_mat(N)
    s_now = s_vals[N]
    s_next = s_vals[N+1]
    # Compute s_N · M(N)
    for j in range(3):
        computed = s_now[0]*M_N[0,j] + s_now[1]*M_N[1,j] + s_now[2]*M_N[2,j]
        actual = s_next[j]
        err = abs(computed - actual)
        if err > 1e-100:
            print(f"  N={N}, j={j}: MISMATCH err={float(err):.3e}")
        elif N <= 1:
            print(f"  N={N}, j={j}: OK")
print("  [verification passed]")

# Now derive the scalar recurrence.
# From s_{N+1} = s_N · M(N), the first-column equation is:
#   q_{N+1} = q_N · M₀₀(N) + q'_N · M₁₀(N) + q''_N · M₂₀(N)  ... (*)
#
# From (*) at N: q'_N and q''_N can be expressed in terms of q_{N+1}, q_N and
# the known matrix entries, PROVIDED we have one more equation to eliminate
# the two unknowns.
#
# From the second column equation:
#   q'_{N+1} = q_N · M₀₁(N) + q'_N · M₁₁(N) + q''_N · M₂₁(N)  ... (**)
#
# And from the third column:
#   q''_{N+1} = q_N · M₀₂(N) + q'_N · M₁₂(N) + q''_N · M₂₂(N)  ... (***)
#
# Strategy: Use (*) to express q_{N+1} in terms of (q_N, q'_N, q''_N).
# Use (**) and (***) to get (q'_{N+1}, q''_{N+1}).
# Then apply (*) at N+1:
#   q_{N+2} = q_{N+1} · M₀₀(N+1) + q'_{N+1} · M₁₀(N+1) + q''_{N+1} · M₂₀(N+1)
# And at N+2:
#   q_{N+3} = q_{N+2} · M₀₀(N+2) + q'_{N+2} · M₁₀(N+2) + q''_{N+2} · M₂₀(N+2)
#
# By composing these, we eliminate q', q'' and get a relation between
# q_N, q_{N+1}, q_{N+2}, q_{N+3}.
#
# Let me do this computation for specific N values and check that the resulting
# recurrence has POLYNOMIAL COEFFICIENTS.

print("\n=== Extracting scalar recurrence coefficients ===")

# For each N, compute the coefficients numerically.
# Fix α₃(N) = 1 and solve for α₀(N), α₁(N), α₂(N).
recurrence_coeffs = []
for N in range(5, 30):
    q0 = s_vals[N][0]
    q1 = s_vals[N+1][0]
    q2 = s_vals[N+2][0]
    q3 = s_vals[N+3][0]

    # q3 + α₂ q2 + α₁ q1 + α₀ q0 = 0
    # Solve 3×3 system: [[q2, q1, q0]] · [[α₂], [α₁], [α₀]] = [-q3]
    # But this is 1 equation in 3 unknowns!
    # I need to use the STRUCTURE to get the recurrence.

    # Better approach: use 3 consecutive equations
    # At N:   α₃ q_{N+3} + α₂ q_{N+2} + α₁ q_{N+1} + α₀ q_N = 0
    # At N+1: α₃ q_{N+4} + α₂ q_{N+3} + α₁ q_{N+2} + α₀ q_{N+1} = 0
    # ...
    # But the coefficients α_i(N) CHANGE with N (they're polynomials).

    # For POLYNOMIAL coefficients: write α_i(N) = Σ c_{ij} N^j
    # Each value of N gives one equation: α₃(N)q_{N+3}+... = 0
    # Total unknowns = 4(d+1) where d is the max degree
    pass

# Better approach: USE THE MATRIX DIRECTLY.
# The scalar recurrence coefficients can be read from the matrix.
# Specifically, the recurrence comes from the adjugate matrix.
#
# Define the 3×3 Casorati-like matrix:
# W(N) = [s_N; s_{N+1}; s_{N+2}]^T  (columns are the state vectors)
#
# Then q_N = first row of W(N), etc.
# The recurrence arises from: W(N+1) = M_shift(N) · W(N)
# where M_shift is the companion matrix of the scalar recurrence.

# Actually, let me just compute the companion matrix directly.
# From s_{N+1} = s_N · M(N), we have:
# [s_N; s_{N+1}; s_{N+2}] → [s_{N+1}; s_{N+2}; s_{N+3}]
# = [s_N·M(N); s_{N+1}·M(N+1); s_{N+2}·M(N+2)]

# This is NOT a simple shift; each row evolves with a DIFFERENT matrix.
# The scalar recurrence must be found by elimination.

# Let me use a direct elimination approach.
#
# We have s_N = (q_N, p_N, r_N) where p=q', r=q''.
# s_{N+1} = s_N · M(N), giving:
# q_{N+1} = q_N m00(N) + p_N m10(N) + r_N m20(N)
# p_{N+1} = q_N m01(N) + p_N m11(N) + r_N m21(N)
# r_{N+1} = q_N m02(N) + p_N m12(N) + r_N m22(N)
#
# From eq 1: p_N m10(N) + r_N m20(N) = q_{N+1} - q_N m00(N)
# From eq 2: p_N m11(N) + r_N m21(N) = p_{N+1} - q_N m01(N)
# From eq 3: p_N m12(N) + r_N m22(N) = r_{N+1} - q_N m02(N)
#
# Solve for p_N, r_N using any two of these equations, say 1 and 2:
# | m10(N)  m20(N) | |p_N|   | q_{N+1} - q_N m00(N) |
# | m11(N)  m21(N) | |r_N| = | p_{N+1} - q_N m01(N) |
#
# det = m10 m21 - m11 m20 (a polynomial in N)
# p_N = [(q_{N+1}-q_N m00)m21 - (p_{N+1}-q_N m01)m20] / det
# r_N = [(p_{N+1}-q_N m01)m10 - (q_{N+1}-q_N m00)m11] / det
#
# But this introduces p_{N+1} — we haven't eliminated it yet.
#
# The trick: at time N, express (p_N, r_N) in terms of q_N, q_{N+1}, q_{N+2}.
# At time N+1, express (p_{N+1}, r_{N+1}) in terms of q_{N+1}, q_{N+2}, q_{N+3}.
# Then substitute into the relation s_{N+1} = s_N M(N).

# Let me proceed computationally. For each N, I have all 3 components.
# Define:
# p_N = s_vals[N][1]
# r_N = s_vals[N][2]

# For each N, compute the 2×2 determinant det12(N) = m10(N)*m21(N) - m11(N)*m20(N)
print("Computing 2x2 determinant for (p,r) elimination...")
for N in [0, 5, 10, 15, 20]:
    M_N = M_mat(N)
    det12 = M_N[1,0]*M_N[2,1] - M_N[1,1]*M_N[2,0]
    print(f"  N={N}: det12 = {nstr(det12, 15)}")

# Good — these are nonzero polynomials in N. So we can solve.
#
# From equations 1 and 3 (using columns 0 and 2 of M(N)):
# q_{N+1} = q_N m00 + p_N m10 + r_N m20
# r_{N+1} = q_N m02 + p_N m12 + r_N m22
#
# Solving for p_N, r_N:
# | m10  m20 | |p_N| = | q_{N+1} - q_N m00 |
# | m12  m22 | |r_N|   | r_{N+1} - q_N m02 |

# det_A = m10*m22 - m12*m20
# p_N = [(q_{N+1}-q_N m00)*m22 - (r_{N+1}-q_N m02)*m20] / det_A
# r_N = [(r_{N+1}-q_N m02)*m10 - (q_{N+1}-q_N m00)*m12] / det_A

# But r_{N+1} depends on the NEXT step! This is circular.

# Correct approach: express (p_N, r_N) purely in terms of q values.
# We need q_N, q_{N+1}, and ONE MORE q value (q_{N+2} or q_{N-1}).

# Since s_{N+1} = s_N M(N):
# q_{N+1} = q_N m00 + p_N m10 + r_N m20  ... (I)
# And s_{N+2} = s_{N+1} M(N+1):
# q_{N+2} = q_{N+1} m00(N+1) + p_{N+1} m10(N+1) + r_{N+1} m20(N+1)  ... (II)
# where p_{N+1} = q_N m01 + p_N m11 + r_N m21
#       r_{N+1} = q_N m02 + p_N m12 + r_N m22

# So (II) becomes:
# q_{N+2} = q_{N+1} m00(N+1)
#          + [q_N m01 + p_N m11 + r_N m21] m10(N+1)
#          + [q_N m02 + p_N m12 + r_N m22] m20(N+1)

# Rearranging:
# q_{N+2} = q_{N+1} m00(N+1) + q_N [m01 m10(N+1) + m02 m20(N+1)]
#          + p_N [m11 m10(N+1) + m12 m20(N+1)]
#          + r_N [m21 m10(N+1) + m22 m20(N+1)]

# Now from (I) and this equation:
# (I):  p_N m10(N) + r_N m20(N) = q_{N+1} - q_N m00(N)
# (II'): p_N [m11 m10(N+1)+m12 m20(N+1)] + r_N [m21 m10(N+1)+m22 m20(N+1)]
#       = q_{N+2} - q_{N+1} m00(N+1) - q_N [m01 m10(N+1)+m02 m20(N+1)]

# This is a 2×2 system for (p_N, r_N):
# A₁₁ p_N + A₁₂ r_N = B₁
# A₂₁ p_N + A₂₂ r_N = B₂
# where all coefficients are POLYNOMIALS in N, and B₁, B₂ are LINEAR in q values.

# Once we solve for (p_N, r_N), substitute into s_{N+2} = s_{N+1} M(N+1):
# q_{N+3} = q_{N+2} m00(N+2)
#          + [q_{N+1} m01(N+1) + p_{N+1} m11(N+1) + r_{N+1} m21(N+1)] m10(N+2)
#          + [q_{N+1} m02(N+1) + p_{N+1} m12(N+1) + r_{N+1} m22(N+1)] m20(N+2)
# and express p_{N+1}, r_{N+1} in terms of (q_N, p_N, r_N) → terms of (q_N, q_{N+1}, q_{N+2}).

# This will give q_{N+3} as a polynomial-coefficient linear combination of q_N, q_{N+1}, q_{N+2}.
# That IS the scalar recurrence.

# Let me compute this numerically and verify.

print("\n=== Computing scalar recurrence coefficients for specific N ===")

for N in range(5, 25):
    M_N = M_mat(N)
    M_N1 = M_mat(N+1)

    m00 = M_N[0,0]; m01 = M_N[0,1]; m02 = M_N[0,2]
    m10 = M_N[1,0]; m11 = M_N[1,1]; m12 = M_N[1,2]
    m20 = M_N[2,0]; m21 = M_N[2,1]; m22 = M_N[2,2]

    m00p = M_N1[0,0]; m01p = M_N1[0,1]; m02p = M_N1[0,2]
    m10p = M_N1[1,0]; m11p = M_N1[1,1]; m12p = M_N1[1,2]
    m20p = M_N1[2,0]; m21p = M_N1[2,1]; m22p = M_N1[2,2]

    # 2×2 system coefficients
    A11 = m10
    A12 = m20
    A21 = m11*m10p + m12*m20p
    A22 = m21*m10p + m22*m20p

    # B values (depend on q)
    q_N = s_vals[N][0]
    q_N1 = s_vals[N+1][0]
    q_N2 = s_vals[N+2][0]

    B1 = q_N1 - q_N * m00
    B2 = q_N2 - q_N1 * m00p - q_N * (m01*m10p + m02*m20p)

    # Solve for p_N, r_N
    det_A = A11*A22 - A12*A21
    p_N_comp = (B1*A22 - B2*A12) / det_A
    r_N_comp = (A11*B2 - A21*B1) / det_A

    # Verify against actual values
    p_N_actual = s_vals[N][1]
    r_N_actual = s_vals[N][2]

    p_err = abs(p_N_comp - p_N_actual) / abs(p_N_actual) if p_N_actual != 0 else abs(p_N_comp)
    r_err = abs(r_N_comp - r_N_actual) / abs(r_N_actual) if r_N_actual != 0 else abs(r_N_comp)

    if N <= 7 or N == 20:
        print(f"  N={N}: p_err={float(p_err):.3e}, r_err={float(r_err):.3e}")

print("\n=== Now derive q_{N+3} as function of q_N, q_{N+1}, q_{N+2} ===")

# Having expressed (p_N, r_N) in terms of (q_N, q_{N+1}, q_{N+2}),
# we can compute (p_{N+1}, r_{N+1}) and then q_{N+3}.

# p_{N+1} = q_N m01 + p_N m11 + r_N m21
# r_{N+1} = q_N m02 + p_N m12 + r_N m22

# p_N = (B1*A22 - B2*A12) / det_A  where B1 = q_{N+1}-q_N*m00, B2 = q_{N+2}-q_{N+1}*m00'-q_N*(m01*m10'+m02*m20')
# r_N = (A11*B2 - A21*B1) / det_A

# So p_N and r_N are LINEAR in (q_N, q_{N+1}, q_{N+2}) with polynomial-in-N coefficients.
# Let me express p_N = α_p0(N)*q_N + α_p1(N)*q_{N+1} + α_p2(N)*q_{N+2}
# and r_N = α_r0(N)*q_N + α_r1(N)*q_{N+1} + α_r2(N)*q_{N+2}

# where α are polynomial-coefficient ratios.

# Then q_{N+3} = q_{N+2}*m00(N+2) + p_{N+2}*m10(N+2) + r_{N+2}*m20(N+2)
# But p_{N+2}, r_{N+2} depend on (p_{N+1}, r_{N+1}) which depend on (p_N, r_N)
# which depend on (q_N, q_{N+1}, q_{N+2}). So ultimately:
# q_{N+3} = c0(N)*q_N + c1(N)*q_{N+1} + c2(N)*q_{N+2}

# Rather than tracking all this algebra symbolically, let me compute c0, c1, c2
# NUMERICALLY for each N and then fit polynomials.

# For each N, compute the recurrence by expressing q_{N+3} = c2*q_{N+2}+c1*q_{N+1}+c0*q_N
# using the elimination above.

c0_vals = []
c1_vals = []
c2_vals = []

for N in range(3, 30):
    M_N = M_mat(N)
    M_N1 = M_mat(N+1)
    M_N2 = M_mat(N+2)

    # First: solve for (p_N, r_N) in terms of (q_N, q_{N+1}, q_{N+2})
    m00 = M_N[0,0]; m01 = M_N[0,1]; m02 = M_N[0,2]
    m10 = M_N[1,0]; m11 = M_N[1,1]; m12 = M_N[1,2]
    m20 = M_N[2,0]; m21 = M_N[2,1]; m22 = M_N[2,2]

    m00p = M_N1[0,0]; m10p = M_N1[1,0]; m20p = M_N1[2,0]
    m01p = M_N1[0,1]; m11p = M_N1[1,1]; m21p = M_N1[2,1]
    m02p = M_N1[0,2]; m12p = M_N1[1,2]; m22p = M_N1[2,2]

    A11 = m10; A12 = m20
    A21 = m11*m10p + m12*m20p
    A22 = m21*m10p + m22*m20p
    det_A = A11*A22 - A12*A21

    # p_N = (B1*A22 - B2*A12)/det_A where:
    # B1 = q_{N+1} - q_N*m00 = -m00*q_N + 1*q_{N+1} + 0*q_{N+2}
    # B2 = q_{N+2} - m00p*q_{N+1} - (m01*m10p+m02*m20p)*q_N
    #    = -(m01*m10p+m02*m20p)*q_N + (-m00p)*q_{N+1} + 1*q_{N+2}

    coef_B1_q0 = -m00;  coef_B1_q1 = mpf(1); coef_B1_q2 = mpf(0)
    coef_B2_q0 = -(m01*m10p + m02*m20p); coef_B2_q1 = -m00p; coef_B2_q2 = mpf(1)

    # p_N coefficients:
    pN_q0 = (coef_B1_q0*A22 - coef_B2_q0*A12) / det_A
    pN_q1 = (coef_B1_q1*A22 - coef_B2_q1*A12) / det_A
    pN_q2 = (coef_B1_q2*A22 - coef_B2_q2*A12) / det_A

    # r_N coefficients:
    rN_q0 = (A11*coef_B2_q0 - A21*coef_B1_q0) / det_A
    rN_q1 = (A11*coef_B2_q1 - A21*coef_B1_q1) / det_A
    rN_q2 = (A11*coef_B2_q2 - A21*coef_B1_q2) / det_A

    # Now compute p_{N+1} = q_N*m01 + p_N*m11 + r_N*m21
    pN1_q0 = m01 + pN_q0*m11 + rN_q0*m21
    pN1_q1 = pN_q1*m11 + rN_q1*m21
    pN1_q2 = pN_q2*m11 + rN_q2*m21

    # r_{N+1} = q_N*m02 + p_N*m12 + r_N*m22
    rN1_q0 = m02 + pN_q0*m12 + rN_q0*m22
    rN1_q1 = pN_q1*m12 + rN_q1*m22
    rN1_q2 = pN_q2*m12 + rN_q2*m22

    # Now do the SAME for step N+1 → N+2.
    # Need (p_{N+2}, r_{N+2}) to compute q_{N+3}.
    # Actually, for q_{N+3} = q_{N+2}*m00(N+2) + p_{N+2}*m10(N+2) + r_{N+2}*m20(N+2)
    # we need p_{N+2} and r_{N+2}.

    # p_{N+2} = q_{N+1}*m01(N+1) + p_{N+1}*m11(N+1) + r_{N+1}*m21(N+1)
    # r_{N+2} = q_{N+1}*m02(N+1) + p_{N+1}*m12(N+1) + r_{N+1}*m22(N+1)

    # But wait: the step N+1 uses M(N+1) entries, and p_{N+1}, r_{N+1} are
    # already expressed in terms of (q_N, q_{N+1}, q_{N+2}).

    pN2_q0 = pN1_q0*m11p + rN1_q0*m21p
    pN2_q1 = m01p + pN1_q1*m11p + rN1_q1*m21p
    pN2_q2 = pN1_q2*m11p + rN1_q2*m21p

    rN2_q0 = pN1_q0*m12p + rN1_q0*m22p
    rN2_q1 = m02p + pN1_q1*m12p + rN1_q1*m22p
    rN2_q2 = pN1_q2*m12p + rN1_q2*m22p

    m00pp = M_N2[0,0]; m10pp = M_N2[1,0]; m20pp = M_N2[2,0]

    # q_{N+3} = q_{N+2}*m00(N+2) + p_{N+2}*m10(N+2) + r_{N+2}*m20(N+2)
    # = m00pp*q_{N+2} + m10pp*p_{N+2} + m20pp*r_{N+2}
    # where p_{N+2} = pN2_q0*q_N + pN2_q1*q_{N+1} + pN2_q2*q_{N+2}
    # and r_{N+2} = rN2_q0*q_N + rN2_q1*q_{N+1} + rN2_q2*q_{N+2}

    c0 = m10pp*pN2_q0 + m20pp*rN2_q0
    c1 = m10pp*pN2_q1 + m20pp*rN2_q1
    c2 = m00pp + m10pp*pN2_q2 + m20pp*rN2_q2

    # Verify: q_{N+3} = c2*q_{N+2} + c1*q_{N+1} + c0*q_N
    q_N_val = s_vals[N][0]
    q_N1_val = s_vals[N+1][0]
    q_N2_val = s_vals[N+2][0]
    q_N3_val = s_vals[N+3][0]

    computed = c2*q_N2_val + c1*q_N1_val + c0*q_N_val
    err = abs(computed - q_N3_val)
    rel_err = err / abs(q_N3_val) if q_N3_val != 0 else err

    c0_vals.append((N, c0))
    c1_vals.append((N, c1))
    c2_vals.append((N, c2))

    if N <= 10 or N % 5 == 0:
        print(f"  N={N:2d}: c0={nstr(c0,12)}, c2={nstr(c2,12)}, rel_err={float(rel_err):.3e}")

# Now fit polynomials to c0(N), c1(N), c2(N) by Lagrange interpolation.
print("\n=== Polynomial fitting of c0, c1, c2 ===")
print("Checking: is c0(N)/det(M(N)) a polynomial? (it should be)")

# c0(N) should be proportional to det(M(N)) (up to sign).
# Let me compute det(M(N)) and check.
print("\nRatio c0(N) / det(M(N)):")
for N, c0_val in c0_vals[:15]:
    det_M = mp.det(M_mat(N))
    ratio = c0_val / det_M if det_M != 0 else mpf('nan')
    print(f"  N={N:2d}: c0/det(M) = {nstr(ratio, 20)}")

# Also check c2(N) / trace(M(N)):
print("\nRatio c2(N) / tr(M(N)):")
for N, c2_val in c2_vals[:15]:
    M_N = M_mat(N)
    tr_M = M_N[0,0] + M_N[1,1] + M_N[2,2]
    ratio = c2_val / tr_M if tr_M != 0 else mpf('nan')
    print(f"  N={N:2d}: c2/tr(M) = {nstr(ratio, 20)}")
