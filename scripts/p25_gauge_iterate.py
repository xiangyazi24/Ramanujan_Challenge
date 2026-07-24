#!/usr/bin/env python3
"""P2.5: Compute V(n) by backwards iteration from V_∞.

From Q4880: M̄(n)·V(n+1) = V(n)·S_D(n)
=> V(n) = M̄(n)·V(n+1)·S_D(n)^{-1}

Start from V(N) = V_∞ for large N, iterate to V(0).
Then look for a rational pattern.
"""
from fractions import Fraction

def M_entries(n):
    n = Fraction(n)
    m11 = (-2*n-5)*(n+3)**2*(136*n**4+1424*n**3+5548*n**2+9551*n+6141)
    m12 = 384*n**6+6384*n**5+44168*n**4+162698*n**3+336377*n**2+369933*n+169011
    m13 = -(480*n**4+4980*n**3+19210*n**2+32690*n+20730)
    m21 = (n+2)**2*(n+3)**2*(4*n+10)*(48*n**3+386*n**2+1017*n+879)
    m22 = (n+2)**2*(-272*n**5-3848*n**4-21732*n**3-61184*n**2-85761*n-47808)
    m23 = (n+2)**2*(320*n**3+2540*n**2+6610*n+5640)
    m31 = (-4*n-10)*(n+2)**2*(n+3)**2*(32*n**4+302*n**3+1037*n**2+1530*n+813)
    m32 = (n+2)**2*(192*n**6+2984*n**5+19116*n**4+64452*n**3+120256*n**2+117279*n+46476)
    m33 = (n+2)**2*(-16*n**5-408*n**4-2912*n**3-8884*n**2-12254*n-6240)
    return [[m11,m12,m13],[m21,m22,m23],[m31,m32,m33]]

def delta_H(n):
    n = Fraction(n)
    return Fraction(-2)*(n+2)**2*(n+3)**2*(2*n+5)*(2*n+7)**2

def M_H(n):
    M = M_entries(n)
    d = delta_H(n)
    return [[M[i][j] / d for j in range(3)] for i in range(3)]

def balanced_M(n):
    MH = M_H(n)
    nn = Fraction(n)
    d_inv = [Fraction(1), Fraction(1, int(nn+1)), Fraction(1, int((nn+1)**2))]
    d_next = [Fraction(1), nn+2, (nn+2)**2]
    return [[d_inv[i] * MH[i][j] * d_next[j] for j in range(3)] for i in range(3)]

def S_D(n):
    """Sym²(B_D(n)) where B_D = [[0,-α],[1,β]], α=(n+1)/(n+2), β=3(2n+3)/(n+2)"""
    nn = Fraction(n)
    alpha = Fraction(nn+1, nn+2)
    beta = Fraction(3*(2*nn+3), nn+2)
    return [
        [Fraction(0), Fraction(0), alpha**2],
        [Fraction(0), -alpha, -2*alpha*beta],
        [Fraction(1), beta, beta**2],
    ]

def mat_mul(A, B):
    n = len(A)
    m = len(B[0])
    k = len(B)
    C = [[Fraction(0)]*m for _ in range(n)]
    for i in range(n):
        for j in range(m):
            for l in range(k):
                C[i][j] += A[i][l] * B[l][j]
    return C

def mat_inv_3(M):
    """Exact 3×3 inverse using Fraction."""
    a, b, c = M[0]
    d, e, f = M[1]
    g, h, k = M[2]
    det = a*(e*k-f*h) - b*(d*k-f*g) + c*(d*h-e*g)
    if det == 0:
        raise ZeroDivisionError("singular matrix")
    di = Fraction(1, 1) / det
    return [
        [di*(e*k-f*h), di*(c*h-b*k), di*(b*f-c*e)],
        [di*(f*g-d*k), di*(a*k-c*g), di*(c*d-a*f)],
        [di*(d*h-e*g), di*(b*g-a*h), di*(a*e-b*d)],
    ]

V_INF = [
    [Fraction(1), Fraction(3), Fraction(17)],
    [Fraction(0), Fraction(-2), Fraction(-12)],
    [Fraction(0), Fraction(0), Fraction(8)],
]

print("=== Backwards iteration of V(n) ===", flush=True)

# V(n) = M̄(n) · V(n+1) · S_D(n)^{-1}
N_START = 30
V = [row[:] for row in V_INF]

V_values = {}

for n in range(N_START, -1, -1):
    Mbar = balanced_M(n)
    SD = S_D(n)
    SD_inv = mat_inv_3(SD)
    V = mat_mul(mat_mul(Mbar, V), SD_inv)
    V_values[n] = [row[:] for row in V]
    if n <= 10 or n == N_START:
        print(f"\nV({n}):", flush=True)
        for row in V:
            print(f"  [{', '.join(str(x) for x in row)}]")

# Verify: M̄(n)·V(n+1) should equal V(n)·S_D(n)
print("\n=== Verification ===", flush=True)
for n in range(10):
    if n in V_values and n+1 in V_values:
        Mbar = balanced_M(n)
        SD = S_D(n)
        LHS = mat_mul(Mbar, V_values[n+1])
        RHS = mat_mul(V_values[n], SD)
        ok = all(LHS[i][j] == RHS[i][j] for i in range(3) for j in range(3))
        print(f"  n={n}: {'EXACT ✓' if ok else 'FAIL'}")

# Analyze the structure of V(n) entries
# Try to find rational expressions
print("\n=== Rational structure analysis ===", flush=True)
print("Looking at V(n)·(n+1)^k for small k...", flush=True)

for k in range(6):
    print(f"\nV(n)·(n+1)^{k}:")
    for n in range(8):
        if n in V_values:
            scaled = [[(Fraction(n+1)**k) * V_values[n][i][j] for j in range(3)] for i in range(3)]
            # Check if entries are polynomial-looking (small denominators)
            max_denom = max(abs(scaled[i][j].denominator) for i in range(3) for j in range(3))
            if max_denom <= 10:
                vals = [[int(scaled[i][j]) if scaled[i][j].denominator == 1 else str(scaled[i][j]) for j in range(3)] for i in range(3)]
                print(f"  n={n}: {vals}  (all integer!)" if max_denom == 1 else f"  n={n}: {vals}")
            else:
                print(f"  n={n}: max_denom = {max_denom}")

# Also look at individual entries
print("\n=== Individual entry analysis ===", flush=True)
for i in range(3):
    for j in range(3):
        print(f"\nV[{i},{j}](n):")
        for n in range(12):
            if n in V_values:
                val = V_values[n][i][j]
                print(f"  n={n}: {val} = {float(val):.8f}")

print("\nDone.")
