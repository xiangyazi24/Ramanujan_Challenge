#!/usr/bin/env python3
"""Problem 2.5: Compare characteristic polynomials of CMF M(n) and D_n² companion.

If the CMF has the same characteristic polynomial as the companion matrix
of the D_n² recurrence (after twist), they're conjugate, which proves
the Ore intertwiner exists.

D_n² recurrence (from Q4802):
  (n+3)²(2n+3) u_{n+3} - (2n+5)(35n²+140n+131) u_{n+2}
  + (2n+3)(35n²+140n+131) u_{n+1} - (2n+5)(n+1)² u_n = 0

Companion matrix C(n) for this recurrence:
  [[0, 1, 0],
   [0, 0, 1],
   [(2n+5)(n+1)²/((n+3)²(2n+3)), -(2n+3)(35n²+140n+131)/((n+3)²(2n+3)), (2n+5)(35n²+140n+131)/((n+3)²(2n+3))]]

But the CMF M(n) is a DIFFERENT matrix (not in companion form).
Their traces, cofactor sums, and determinants should match if they're conjugate.
"""
from fractions import Fraction

def M_entries(n):
    """Return M(n) as a 3×3 matrix of exact integers."""
    m11 = (-2*n-5)*(n+3)**2*(136*n**4+1424*n**3+5548*n**2+9551*n+6141)
    m12 = 384*n**6+6384*n**5+44168*n**4+162698*n**3+336377*n**2+369933*n+169011
    m13 = -480*n**4-4980*n**3-19210*n**2-32690*n-20730
    m21 = (n+2)**2*(n+3)**2*(4*n+10)*(48*n**3+386*n**2+1017*n+879)
    m22 = (n+2)**2*(-272*n**5-3848*n**4-21732*n**3-61184*n**2-85761*n-47808)
    m23 = (n+2)**2*(320*n**3+2540*n**2+6610*n+5640)
    m31 = (-4*n-10)*(n+2)**2*(n+3)**2*(32*n**4+302*n**3+1037*n**2+1530*n+813)
    m32 = (n+2)**2*(192*n**6+2984*n**5+19116*n**4+64452*n**3+120256*n**2+117279*n+46476)
    m33 = (n+2)**2*(-16*n**5-408*n**4-2912*n**3-8884*n**2-12254*n-6240)
    return [[m11,m12,m13],[m21,m22,m23],[m31,m32,m33]]

def trace(M):
    return M[0][0] + M[1][1] + M[2][2]

def cofactor_sum(M):
    """Sum of 2×2 principal minors."""
    return (M[0][0]*M[1][1] - M[0][1]*M[1][0] +
            M[0][0]*M[2][2] - M[0][2]*M[2][0] +
            M[1][1]*M[2][2] - M[1][2]*M[2][1])

def det(M):
    return (M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1])
           -M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0])
           +M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0]))

# The D_n² recurrence: (n+3)²(2n+3) u_{n+3} = (2n+5)(35n²+140n+131) u_{n+2}
#   - (2n+3)(35n²+140n+131) u_{n+1} + (2n+5)(n+1)² u_n
# Companion matrix C(n) has char poly λ³ - tr(C)λ² + cofsum(C)λ - det(C) = 0.
# For the companion matrix:
#   tr(C) = a₂(n)/a₃(n)  where a₃ = (n+3)²(2n+3), a₂ = (2n+5)(35n²+140n+131)
#   cofsum = a₁(n)/a₃(n) where a₁ = -(2n+3)(35n²+140n+131)
#   det(C) = a₀(n)/a₃(n) where a₀ = (2n+5)(n+1)²
# Wait, det of companion matrix [[0,1,0],[0,0,1],[a0,a1,a2]] = -a0
# And cofsum of companion = a1 (for 3×3 companion)
# And trace = a2.

# Actually: companion matrix for u_{n+3} = c2·u_{n+2} + c1·u_{n+1} + c0·u_n:
# C = [[0, 1, 0], [0, 0, 1], [c0, c1, c2]]
# char poly: λ³ - c2·λ² - c1·λ - c0 = 0
# tr = c2, cofsum = -c1 (from (0·0-1·0)+(0·c2-0·0)+(0·c2-1·c0) = -c1... hmm no)

# Let me just compute directly.
# Companion C = [[0,1,0],[0,0,1],[c0,c1,c2]]
# tr(C) = c2
# Principal 2×2 minors:
#   [0,1;0,0] det = 0
#   [0,0;c0,c2] det = 0
#   [0,1;c0,c2] = -c0... no, it's the sum of 2×2 principal minors:
#   M[0,1;0,1] = 0*0-1*0 = 0
#   M[0,2;0,2] = 0*c2-0*c0 = 0
#   M[1,2;1,2] = 0*c2-1*c1 = -c1
# So cofsum = -c1
# det(C) = 0*(0*c2-1*c1) - 1*(0*c2-1*c0) + 0*... = c0
# char poly: λ³ - c2·λ² + (-c1)λ - c0 = 0 ← correct (Newton's identities)

def D2_companion(n):
    """Return the companion matrix coefficients for D_n² recurrence.
    Recurrence: a₃ u_{n+3} = a₂ u_{n+2} - a₁ u_{n+1} + a₀ u_n
    where a₃=(n+3)²(2n+3), a₂=(2n+5)(35n²+140n+131),
    a₁=(2n+3)(35n²+140n+131), a₀=(2n+5)(n+1)²
    """
    a3 = Fraction((n+3)**2 * (2*n+3))
    a2 = Fraction((2*n+5) * (35*n**2 + 140*n + 131))
    a1 = Fraction((2*n+3) * (35*n**2 + 140*n + 131))
    a0 = Fraction((2*n+5) * (n+1)**2)
    # u_{n+3} = (a2/a3)u_{n+2} - (a1/a3)u_{n+1} + (a0/a3)u_n
    c2 = a2 / a3
    c1 = -a1 / a3
    c0 = a0 / a3
    return c0, c1, c2

print("=== Characteristic polynomial comparison: CMF M(n) vs D_n² companion ===\n")

for n in range(8):
    Mn = M_entries(n)
    tr_M = trace(Mn)
    cf_M = cofactor_sum(Mn)
    det_M = det(Mn)

    c0, c1, c2 = D2_companion(n)
    tr_C = c2
    cf_C = -c1
    det_C = c0

    # For the CMF to be conjugate to the D_n² companion (after twist),
    # we need: tr(M) = h(n) · tr(C) or some twist relation.
    # More precisely, if M = P(n)·C(n)·P(n+1)⁻¹ · diag(h), then
    # the characteristic polynomials don't match directly — they differ
    # by the twist factor.

    print(f"n={n}:")
    print(f"  CMF: tr={tr_M}, cofsum={cf_M}, det={det_M}")
    print(f"  D_n²: tr={float(tr_C):.6f}, cofsum={float(cf_C):.6f}, det={float(det_C):.6f}")

    # Check ratios
    if tr_C != 0:
        r_tr = Fraction(tr_M) / tr_C
        r_cf = Fraction(cf_M) / cf_C if cf_C != 0 else None
        r_det = Fraction(det_M) / det_C if det_C != 0 else None
        print(f"  Ratios: tr_M/tr_C = {float(r_tr):.6e}", end="")
        if r_cf is not None:
            print(f", cf_M/cf_C = {float(r_cf):.6e}", end="")
        if r_det is not None:
            print(f", det_M/det_C = {float(r_det):.6e}")
        else:
            print()
    print()

# The characteristic polynomials won't match directly because the CMF
# carries a hypergeometric twist. The TWISTED companion matrix would be:
# M_twisted(n) = diag(h_{n+1}/h_n, h_{n+2}/h_n, h_{n+3}/h_n) · C(n)
# or something similar. Let's check if the RATIOS of the symmetric
# functions have a pattern.

print("\n--- Ratio patterns ---")
for n in range(8):
    Mn = M_entries(n)
    det_M = det(Mn)
    c0, c1, c2 = D2_companion(n)
    det_C = c0
    if det_C != 0:
        r = Fraction(det_M) / det_C
        print(f"n={n}: det(M)/det(C) = {r} = {float(r):.6e}")

# The det ratio should be the product of twist factors:
# det(M(n)) / det(C(n)) = h_{n+1}·h_{n+2}·h_{n+3} / h_n³ or similar
# where h_n is the hypergeometric normalizing factor.
