#!/usr/bin/env sage
"""Problem 2.5: Extract order-3 recurrence and factor as L1·(S-1).

Steps:
1. Compute 150 exact QQ terms of scalar sequence from matrix product
2. Use ore_algebra to guess the recurrence
3. Factor it — expecting (S-1) as right factor
4. Extract the gauge r(n) = h(n+1)/h(n) from the hypergeometric solution
"""
from sage.all import *
from sage.matrix.constructor import matrix as sage_matrix

# The 3x3 matrix M(n) from Problem 2.5
R = PolynomialRing(QQ, 'n')
n = R.gen()

def M_mat(nn):
    """Return 3x3 matrix M(nn) over QQ."""
    m11 = (-2*nn-5)*(nn+3)^2 * (136*nn^4 + 1424*nn^3 + 5548*nn^2 + 9551*nn + 6141)
    m12 = 384*nn^6 + 6384*nn^5 + 44168*nn^4 + 162698*nn^3 + 336377*nn^2 + 369933*nn + 169011
    m13 = -480*nn^4 - 4980*nn^3 - 19210*nn^2 - 32690*nn - 20730
    m21 = (nn+2)^2*(nn+3)^2*(4*nn+10)*(48*nn^3 + 386*nn^2 + 1017*nn + 879)
    m22 = (nn+2)^2*(-272*nn^5 - 3848*nn^4 - 21732*nn^3 - 61184*nn^2 - 85761*nn - 47808)
    m23 = (nn+2)^2*(320*nn^3 + 2540*nn^2 + 6610*nn + 5640)
    m31 = (-4*nn-10)*(nn+2)^2*(nn+3)^2*(32*nn^4 + 302*nn^3 + 1037*nn^2 + 1530*nn + 813)
    m32 = (nn+2)^2*(192*nn^6 + 2984*nn^5 + 19116*nn^4 + 64452*nn^3 + 120256*nn^2 + 117279*nn + 46476)
    m33 = (nn+2)^2*(-16*nn^5 - 408*nn^4 - 2912*nn^3 - 8884*nn^2 - 12254*nn - 6240)
    return sage_matrix(QQ, [[m11, m12, m13], [m21, m22, m23], [m31, m32, m33]])

# Initial matrix A (first two rows select p_N and q_N)
A = sage_matrix(QQ, [[30921, -32972, 8240],
                       [33750, -36000, 9000]])

# Compute exact QQ terms of the scalar sequence q_N = (A·T_N)[1,0]
print("Computing 150 exact QQ terms of q_N...")
N_terms = 150
T = identity_matrix(QQ, 3)
q_vals = []

for N in range(N_terms):
    AT = A * T
    q_vals.append(AT[1, 0])  # q_N = second row, first column
    T = T * M_mat(N)
    if N % 30 == 0:
        print(f"  N={N} done, q_N has {len(str(q_vals[-1]))} digits")

print(f"Computed {len(q_vals)} terms")
print(f"First few: {q_vals[:5]}")

# Use ore_algebra to guess the recurrence
print("\nGuessing recurrence with ore_algebra...")
try:
    from ore_algebra import OreAlgebra
    Rn.<n> = QQ[]
    A_ore.<Sn> = OreAlgebra(Rn, 'Sn')

    # Guess the recurrence
    L = guess(q_vals, A_ore)
    print(f"Found operator of order {L.order()}")
    print(f"Degree pattern: {[L[i].degree() for i in range(L.order()+1)]}")

    # Factor it
    print("\nFactoring the operator...")
    factors = L.factor()
    print(f"Factorization: {len(factors)} factors")
    for i, (f, m) in enumerate(factors):
        print(f"  Factor {i}: order {f.order()}, multiplicity {m}")
        print(f"    {f}")

    # Check for (Sn - 1) right factor
    print("\nChecking for (Sn - 1) right factor...")
    # The right factor corresponds to the LAST factor in the factorization
    # since L = L1 * (Sn - 1) means Sn-1 is on the RIGHT

    # Try right division by (Sn - 1)
    R_factor = Sn - 1
    try:
        Q, remainder = L.quo_rem(R_factor)
        print(f"  Quotient order: {Q.order()}")
        print(f"  Remainder: {remainder}")
        if remainder == 0:
            print("  ✓ (Sn - 1) IS a right factor!")
            print(f"\n  Quotient L1 = {Q}")
        else:
            print("  ✗ (Sn - 1) is NOT a right factor of the RAW recurrence")
            print("  Need to find gauge first...")
    except Exception as e:
        print(f"  Division error: {e}")

    # Try to find hypergeometric solutions (Petkovšek)
    print("\nSearching for hypergeometric solutions...")
    try:
        hyp_sols = L.hypergeometric_solutions()
        print(f"Found {len(hyp_sols)} hypergeometric solution(s)")
        for i, sol in enumerate(hyp_sols):
            print(f"  Solution {i}: {sol}")
    except Exception as e:
        print(f"  Hypergeometric search error: {e}")

    # Try right factors of order 1
    print("\nSearching for right factors of order 1...")
    try:
        rfactors = L.right_factors(1)
        print(f"Found {len(rfactors)} right factor(s) of order 1")
        for i, rf in enumerate(rfactors):
            print(f"  Right factor {i}: {rf}")
    except Exception as e:
        print(f"  Right factor search error: {e}")

    # Poincaré polynomial
    print("\nPoincaré polynomial:")
    coeffs = [L[i].leading_coefficient() for i in range(L.order()+1)]
    print(f"  Leading coefficients: {coeffs}")
    c = var('c')
    poincare = sum(coeffs[i] * c^i for i in range(L.order()+1))
    print(f"  Polynomial: {poincare}")
    print(f"  Factored: {factor(poincare)}")

except ImportError:
    print("ore_algebra not installed. Installing...")
    import subprocess
    subprocess.run(["sage", "-pip", "install", "ore_algebra"], check=True)
    print("Installed ore_algebra. Please re-run the script.")
except Exception as e:
    print(f"Error: {e}")
    import traceback
    traceback.print_exc()
