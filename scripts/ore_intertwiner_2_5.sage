#!/usr/bin/env sage
"""
Problem 2.5: Compute the Ore intertwiner between L_25 and Sym^2(Delannoy).

L_25: order-3 operator with degree pattern (28,21,14,7)
L_D2: Sym^2(Delannoy) = order-3 operator with small coefficients

We seek T = t_0(n) + t_1(n)*S + t_2(n)*S^2 such that L_25 * T = R * L_D2.

Strategy:
1. Reconstruct L_25 from the CMF matrix M(n) using Casorati minors
2. L_D2 is known explicitly (eq 2.3 from ChatGPT)
3. Set up the intertwiner equation as linear algebra over Q(n)
"""

from ore_algebra import *

# Set up the Ore algebra
R.<n> = QQ['n']
A.<Sn> = OreAlgebra(R, 'Sn')

# ---- Sym^2(Delannoy) recurrence ----
# (n+3)^2(2n+3) U_{n+3} - (2n+5)(35n^2+140n+131) U_{n+2}
# + (2n+3)(35n^2+140n+131) U_{n+1} - (2n+5)(n+1)^2 U_n = 0

c0_D2 = -(2*n+5)*(n+1)^2
c1_D2 = (2*n+3)*(35*n^2 + 140*n + 131)
c2_D2 = -(2*n+5)*(35*n^2 + 140*n + 131)
c3_D2 = (n+3)^2*(2*n+3)

L_D2 = c0_D2 + c1_D2*Sn + c2_D2*Sn^2 + c3_D2*Sn^3

print("L_D2 (Sym^2 Delannoy):")
print(L_D2)
print()

# ---- Reconstruct L_25 from the CMF matrix ----
# The CMF matrix M(n) has polynomial entries.
# The scalar recurrence is obtained from Casorati minors.
# We need to compute this from the matrix data.

# M(n) entries:
def M_entry(n_val, i, j):
    """Return M(n)[i,j] as a polynomial in n."""
    n = n_val
    entries = {
        (0,0): (-2*n-5)*(n+3)^2 * (136*n^4 + 1424*n^3 + 5548*n^2 + 9551*n + 6141),
        (0,1): 384*n^6 + 6384*n^5 + 44168*n^4 + 162698*n^3 + 336377*n^2 + 369933*n + 169011,
        (0,2): -480*n^4 - 4980*n^3 - 19210*n^2 - 32690*n - 20730,
        (1,0): (n+2)^2*(n+3)^2*(4*n+10)*(48*n^3 + 386*n^2 + 1017*n + 879),
        (1,1): (n+2)^2*(-272*n^5 - 3848*n^4 - 21732*n^3 - 61184*n^2 - 85761*n - 47808),
        (1,2): (n+2)^2*(320*n^3 + 2540*n^2 + 6610*n + 5640),
        (2,0): (-4*n-10)*(n+2)^2*(n+3)^2*(32*n^4 + 302*n^3 + 1037*n^2 + 1530*n + 813),
        (2,1): (n+2)^2*(192*n^6 + 2984*n^5 + 19116*n^4 + 64452*n^3 + 120256*n^2 + 117279*n + 46476),
        (2,2): (n+2)^2*(-16*n^5 - 408*n^4 - 2912*n^3 - 8884*n^2 - 12254*n - 6240),
    }
    return entries[(i,j)]

# Build M as a matrix over QQ[n]
M_mat = matrix(R, 3, 3, lambda i,j: M_entry(n, i, j))
print("M(n) =")
print(M_mat)
print()
print("det M(n) =", M_mat.det().factor())
print()

# To get the scalar recurrence, we compute the minimal annihilating operator
# for the first column of the matrix product A * M(0) * M(1) * ... * M(N-1).
#
# This is equivalent to finding the operator L such that:
# c_3(N) * q_{N+3} + c_2(N) * q_{N+2} + c_1(N) * q_{N+1} + c_0(N) * q_N = 0
# where q_N = [A * prod(M(k), k=0..N-1)]_{2,1}
#
# The Casorati method: form the vectors v_j(n) = M(n) * M(n+1) * ... * M(n+j-1) * e_1
# and compute the determinantal minors.

# Alternative: use ore_algebra's sequence recognition.
# Compute q_N for many N values, then use guess_rec.

print("Computing scalar sequence for ore_algebra recognition...")
print("(This requires numerical computation with exact arithmetic)")

# We work over QQ for exact results
from sage.all import *

# Compute matrix products exactly over QQ
A_init = matrix(QQ, [[30921, -32972, 8240],
                      [33750, -36000, 9000]])

# Compute q_N = (A * M(0) * ... * M(N-1))_{2,1} for N = 0, 1, ...
# This is the second row, first column of A * product

N_max = 30  # Should be enough for degree (28,21,14,7) recognition
q_vals = []
M_prod = identity_matrix(QQ, 3)

for N in range(N_max + 5):
    AM = A_init * M_prod
    q_val = AM[1, 0]  # second row, first column
    q_vals.append(q_val)
    # Update product: multiply by M(N)
    M_N = matrix(QQ, 3, 3, lambda i,j: M_entry(N, i, j))
    M_prod = M_prod * M_N

print(f"Computed q[0..{len(q_vals)-1}]")
print(f"q[0] = {q_vals[0]}")
print(f"q[1] = {q_vals[1]}")

# Use ore_algebra to guess the recurrence
print("\nGuessing recurrence with ore_algebra...")
from ore_algebra import guess

try:
    L_25 = guess(q_vals, A)
    print("L_25 =", L_25)
    print(f"\nOrder: {L_25.order()}")

    # Check degree pattern
    coeffs_25 = L_25.list()
    print(f"Number of coefficients: {len(coeffs_25)}")
    for i, c in enumerate(coeffs_25):
        if c != 0:
            print(f"  coeff of S^{i}: degree {c.degree()}")

    # Now compute the Ore intertwiner
    print("\n=== Computing Ore intertwiner ===")
    print("Seeking T such that L_25 * T = R * L_D2")

    # Try right division: does L_D2 divide L_25 from the right?
    # In other words, does L_25 = Q * L_D2 for some Q?
    Q, rem = L_25.quo_rem(L_D2)
    print(f"\nL_25 / L_D2: quotient order = {Q.order()}, remainder = {rem}")

    if rem == 0:
        print("EXACT DIVISION! L_25 = Q * L_D2")
        print(f"Q = {Q}")
    else:
        print("L_D2 does not right-divide L_25")
        print(f"Remainder order: {rem.order()}")

        # Try LCLM (least common left multiple)
        print("\nComputing LCLM(L_25, L_D2)...")
        try:
            lclm = L_25.lclm(L_D2)
            print(f"LCLM order: {lclm.order()}")
        except Exception as e:
            print(f"LCLM failed: {e}")

        # Try GCRD (greatest common right divisor)
        print("\nComputing GCRD(L_25, L_D2)...")
        try:
            gcrd = L_25.gcrd(L_D2)
            print(f"GCRD order: {gcrd.order()}")
            print(f"GCRD = {gcrd}")
        except Exception as e:
            print(f"GCRD failed: {e}")

except Exception as e:
    print(f"Guess failed: {e}")
    print("Need more terms or different approach")
