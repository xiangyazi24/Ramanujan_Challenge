#!/usr/bin/env python3
"""
Problem 2.5: Compute det(M(N)), det(S(N)), and find the right normalization
for the Ore intertwiner.

The key insight: CMF recurrence has degree pattern (28,21,14,7), slope=-7
on the Newton polygon. Solutions grow like (N!)^7 · ρ^N.
Sym² has degree pattern (3,3,3,3), slope=0. Solutions grow like ρ_S^N.

The intertwiner needs to account for this (N!)^7 mismatch.
"""
from fractions import Fraction as F
from sympy import Symbol, factor, Rational, Poly, expand, simplify

N = Symbol('N')

def M_sym(n):
    m11 = (-2*n-5)*(n+3)**2 * (136*n**4 + 1424*n**3 + 5548*n**2 + 9551*n + 6141)
    m12 = 384*n**6 + 6384*n**5 + 44168*n**4 + 162698*n**3 + 336377*n**2 + 369933*n + 169011
    m13 = -480*n**4 - 4980*n**3 - 19210*n**2 - 32690*n - 20730
    m21 = (n+2)**2*(n+3)**2*(4*n+10)*(48*n**3 + 386*n**2 + 1017*n + 879)
    m22 = (n+2)**2*(-272*n**5 - 3848*n**4 - 21732*n**3 - 61184*n**2 - 85761*n - 47808)
    m23 = (n+2)**2*(320*n**3 + 2540*n**2 + 6610*n + 5640)
    m31 = (-4*n-10)*(n+2)**2*(n+3)**2*(32*n**4 + 302*n**3 + 1037*n**2 + 1530*n + 813)
    m32 = (n+2)**2*(192*n**6 + 2984*n**5 + 19116*n**4 + 64452*n**3 + 120256*n**2 + 117279*n + 46476)
    m33 = (n+2)**2*(-16*n**5 - 408*n**4 - 2912*n**3 - 8884*n**2 - 12254*n - 6240)
    return [[m11, m12, m13], [m21, m22, m23], [m31, m32, m33]]

# Compute det(M(N)) numerically at several points to find the polynomial
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
    return [[m11, m12, m13], [m21, m22, m23], [m31, m32, m33]]

def det3(M):
    a,b,c = M[0]; d,e,f = M[1]; g,h,i = M[2]
    return a*(e*i-f*h) - b*(d*i-f*g) + c*(d*h-e*g)

# Compute det(M(n)) for n=0,...,30
det_vals = []
for n in range(30):
    det_vals.append(det3(M_int(n)))

# Try to fit as a polynomial in n
# First check degree by looking at ratios
print("=== det(M(n)) values ===")
for n in range(10):
    print(f"  n={n}: det = {det_vals[n]}")

print(f"\n=== Estimating degree by finite differences ===")
# Apply differences until we get zeros (for polynomial)
diffs = list(det_vals[:25])
for order in range(1, 25):
    diffs = [diffs[i+1] - diffs[i] for i in range(len(diffs)-1)]
    if all(d == 0 for d in diffs[:5]):
        print(f"  Degree = {order - 1} (differences vanish at order {order})")
        break
    if len(diffs) < 5:
        break
    print(f"  Order {order}: diffs[0:3] = {[float(d) for d in diffs[:3]]}")
else:
    print("  Degree > 24 (not polynomial or very high degree)")

# Interpolate the polynomial
from sympy import interpolate as symp_interp
print("\n=== Factoring det(M(N)) ===")
det_points = [(i, det_vals[i]) for i in range(25)]
det_poly = symp_interp(det_points, N)
det_factored = factor(det_poly)
print(f"det(M(N)) = {det_factored}")

# Sym² companion det
a0 = -(2*N+5)*(N+1)**2
a3 = (2*N+3)*(N+3)**2
det_S = -a0/a3  # companion matrix det = (-1)^3 · (-a₀/a₃) = a₀/a₃
# Actually for companion matrix [[0,1,0],[0,0,1],[-a0/a3,-a1/a3,-a2/a3]]:
# det = (-1)^{n+1} · a0/a3 for order 3 companion... let me compute it directly.
# det of [[0,1,0],[0,0,1],[-p,-q,-r]] = 0·(0·(-r)-1·(-q)) - 1·(0·(-r)-1·(-p)) + 0·(...)
# = -1·(p) = -(-a0/a3) = a0/a3... hmm
# Actually: det = expand cofactors along first row:
# 0·det([[0,1],[-q,-r]]) - 1·det([[0,1],[-p,-r]]) + 0·det(...) = -(0·(-r)-1·(-p)) = -p = a0/a3
# Wait: [[0,1,0],[0,0,1],[-a0/a3,-a1/a3,-a2/a3]]
# det = 0·(0·(-a2/a3)-1·(-a1/a3)) - 1·(0·(-a2/a3)-1·(-a0/a3)) + 0·(...)
# = -1·(0 - (-a0/a3)) = -1·(a0/a3) = -a0/a3
# But for standard companion: det(companion) = (-1)^n · a0/a3 for degree n recurrence.
# For n=3: det = -a0/a3 = (2N+5)(N+1)²/[(2N+3)(N+3)²]
det_S_sym = (2*N+5)*(N+1)**2 / ((2*N+3)*(N+3)**2)
print(f"\ndet(S(N)) = {factor(det_S_sym)}")

# Ratio
print("\n=== det(M(N)) / det(S(N)) ===")
det_ratio = det_poly * (2*N+3)*(N+3)**2 / ((2*N+5)*(N+1)**2)
det_ratio_simplified = factor(expand(det_ratio))
print(f"det(M)/det(S) = {det_ratio_simplified}")

# Product ∏ det(M(k))/det(S(k)) for k=0..N-1
# This is the det(G(N)) that grows super-exponentially
print("\n=== ∏ det(M(k))/det(S(k)) for k=0..N-1 ===")
prod = F(1)
for k in range(15):
    dM = det_vals[k]
    # det(S(k)) = (2k+5)(k+1)²/[(2k+3)(k+3)²]
    dS = F((2*k+5)*(k+1)**2, (2*k+3)*(k+3)**2)
    ratio = F(dM) / dS
    prod *= ratio
    print(f"  k={k}: det(M)/det(S) = {float(ratio):.6e}, cumul prod = {float(prod):.6e}")

# Try to identify det(M(N)) / [(2N+5)(N+1)²] as a polynomial times (2N+3)(N+3)²
# i.e., write det(M) = P(N) · (2N+5)(N+1)² / [(2N+3)(N+3)²] where P is polynomial
print("\n\n=== Structure of det(M(N)) ===")
# det(M(N)) evaluated at n as integer
# Check if det(M(N)) is divisible by certain factors
for n in range(5):
    val = det_vals[n]
    factors_to_check = {
        '(n+2)^4': (n+2)**4,
        '(n+3)^2': (n+3)**2,
        '(2n+3)': (2*n+3),
        '(2n+5)': (2*n+5),
        '(n+1)^2': (n+1)**2 if n>0 else 1,
    }
    for name, fac in factors_to_check.items():
        if fac != 0 and val % fac == 0:
            print(f"  n={n}: det(M) divisible by {name}")

print("\nDone.")
