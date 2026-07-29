#!/usr/bin/env python3
"""
Problem 2.5: Check if M(n) = G(n) * Sym^2(A_D(n)) * G(n+1)^{-1}
where A_D(n) is the Delannoy companion matrix.

This directly tests the matrix-level Sym^2 structure.
"""
from fractions import Fraction as F

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
    return [[F(m11), F(m12), F(m13)], [F(m21), F(m22), F(m23)], [F(m31), F(m32), F(m33)]]

def sym2_companion(n):
    """Sym^2 of Delannoy companion [[a, b], [1, 0]] where
    a = 3(2n+1)/(n+1), b = -n/(n+1)."""
    a = F(3*(2*n+1), n+1)
    b = F(-n, n+1)
    return [[a*a,   2*a*b, b*b],
            [a,     a*0+b, b*0],  # ad+bc where d=0, so = bc
            [F(1),  F(0),  F(0)]]
    # Wait: Sym^2([[a,b],[c,d]]) = [[a^2, 2ab, b^2],[ac, ad+bc, bd],[c^2, 2cd, d^2]]
    # With c=1, d=0: [[a^2, 2ab, b^2],[a, b, 0],[1, 0, 0]]

def sym2_companion_correct(n):
    a = F(3*(2*n+1), n+1)
    b = F(-n, n+1)
    return [[a*a,   2*a*b, b*b],
            [a,     b,     F(0)],
            [F(1),  F(0),  F(0)]]

def mat_mul_F(A, B):
    n = len(A)
    return [[sum(A[i][k]*B[k][j] for k in range(n)) for j in range(n)] for i in range(n)]

def mat_inv_F(M):
    """Inverse of 3x3 matrix over Q."""
    a, b, c = M[0]
    d, e, f = M[1]
    g, h, i = M[2]
    det = a*(e*i-f*h) - b*(d*i-f*g) + c*(d*h-e*g)
    if det == 0:
        return None
    inv = [[(e*i-f*h)/det, -(b*i-c*h)/det, (b*f-c*e)/det],
           [-(d*i-f*g)/det, (a*i-c*g)/det, -(a*f-c*d)/det],
           [(d*h-e*g)/det, -(a*h-b*g)/det, (a*e-b*d)/det]]
    return inv

# Compute G(n+1) = M(n)^{-1} * G(n) * Sym^2(A_D(n))
# Start with G(0) = I and see what happens
I3 = [[F(1), F(0), F(0)], [F(0), F(1), F(0)], [F(0), F(0), F(1)]]

print("=== Testing M(n) = G(n) * Sym^2(A_D(n)) * G(n+1)^{-1} ===")
print("Propagating G(n+1) = M(n)^{-1} * G(n) * Sym^2(A_D(n))")
print()

G = [row[:] for row in I3]

for n in range(8):
    Mn = M_int(n)
    Sn = sym2_companion_correct(n)
    Mn_inv = mat_inv_F(Mn)
    if Mn_inv is None:
        print(f"  M({n}) is singular!")
        break
    # G(n+1) = M(n)^{-1} * G(n) * S(n)
    temp = mat_mul_F(G, Sn)
    G = mat_mul_F(Mn_inv, temp)

    print(f"G({n+1}):")
    for i in range(3):
        entries = []
        for j in range(3):
            v = G[i][j]
            if v == 0:
                entries.append("0")
            elif v.denominator == 1:
                entries.append(str(v.numerator))
            else:
                entries.append(f"{v.numerator}/{v.denominator}")
        print(f"  [{', '.join(entries)}]")
    print()

# Check if G(n) stabilizes to a rational function pattern
# If entries grow, the constant-G approach fails
print("\n=== Growth diagnostic ===")
G = [row[:] for row in I3]
for n in range(20):
    Mn = M_int(n)
    Sn = sym2_companion_correct(n)
    Mn_inv = mat_inv_F(Mn)
    if Mn_inv is None:
        break
    temp = mat_mul_F(G, Sn)
    G = mat_mul_F(Mn_inv, temp)

    # Measure size (numerator + denominator digits)
    max_digits = 0
    for i in range(3):
        for j in range(3):
            v = G[i][j]
            if v != 0:
                d = len(str(abs(v.numerator))) + len(str(abs(v.denominator)))
                max_digits = max(max_digits, d)
    print(f"  n={n+1}: max digits = {max_digits}")

print("\n=== Alternative: try G(n) = D_n * H where H is constant ===")
print("If the gauge is G(n) = diag(D_n^2, D_n*E_n, E_n^2) * H...")

# Compute Delannoy
D = [F(1), F(3)]
E = [F(0), F(1)]
for n in range(1, 30):
    D.append((F(3*(2*n+1)) * D[n] - F(n) * D[n-1]) / F(n+1))
    E.append((F(3*(2*n+1)) * E[n] - F(n) * E[n-1]) / F(n+1))

# At each n, the Sym^2 companion maps [D_n^2, D_n*E_n, E_n^2] to [D_{n+1}^2, D_{n+1}*E_{n+1}, E_{n+1}^2]
# So the Sym^2 product is the identity in this basis (it just shifts n)
# Thus: CMF product in the Sym^2 basis should give the gauge matrix.

# Define V(n) = [[D_n^2, D_n*E_n, E_n^2]] (column vectors of Sym^2 solutions)
# Then V(n+1) = Sym^2(A_D(n)) * V(n)
# We want: M(n) * ?? = ?? * Sym^2(A_D(n))

# Let's try: does M(n) * v(n+1) = lambda(n) * v(n) for v = [D^2, DE, E^2]^T?
# If so, M(n) is a rank-1 map in the Sym^2 basis (unlikely for 3x3).

# More general: find C such that M(n) * [D_{n+1}^2; D_{n+1}E_{n+1}; E_{n+1}^2] * C
#   is proportional to [D_n^2; D_n*E_n; E_n^2] * C' ...

# Actually, let me check if the columns of M(n)^{-1} are proportional to
# [D_{n+1}^2, D_{n+1}E_{n+1}, E_{n+1}^2]^T or similar.

print("\nChecking column ratios of M(n)^{-1} vs Sym^2 Delannoy vectors:")
for n in range(5):
    Mn = M_int(n)
    Mn_inv = mat_inv_F(Mn)
    D_vals = [D[n+1]**2, D[n+1]*E[n+1], E[n+1]**2]
    for j in range(3):
        col = [Mn_inv[i][j] for i in range(3)]
        # Check if col is proportional to D_vals
        if D_vals[0] != 0 and col[0] != 0:
            ratio = col[0] / D_vals[0]
            proportional = all(
                (D_vals[i] == 0 and col[i] == 0) or
                (D_vals[i] != 0 and col[i] / D_vals[i] == ratio)
                for i in range(3)
            )
            if proportional:
                print(f"  n={n}, col {j}: PROPORTIONAL to Sym^2 Delannoy! ratio={ratio}")

# Let's check M(n) * [D_{n+1}^2, D_{n+1}E_{n+1}, E_{n+1}^2]^T
print("\nM(n) * [D_{n+1}^2, D_{n+1}E_{n+1}, E_{n+1}^2]^T:")
for n in range(6):
    Mn = M_int(n)
    v = [D[n+1]**2, D[n+1]*E[n+1], E[n+1]**2]
    result = [sum(Mn[i][j] * v[j] for j in range(3)) for i in range(3)]
    # Compare with [D_n^2, D_n*E_n, E_n^2]
    ref = [D[n]**2, D[n]*E[n], E[n]**2]
    print(f"  n={n}:")
    for i in range(3):
        if ref[i] != 0:
            r = result[i] / ref[i]
            print(f"    row {i}: result/ref = {r}")
        else:
            print(f"    row {i}: ref=0, result={result[i]}")
