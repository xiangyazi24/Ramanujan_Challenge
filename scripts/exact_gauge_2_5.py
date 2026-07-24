#!/usr/bin/env python3
"""
Problem 2.5: Exact rational computation of gauge from Casorati.

Compute r₀(N) exactly and identify r₀(N)/r₀(N-1) as a rational function of N.
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
from sympy import Symbol, factor, Poly, Rational, simplify

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

def mat_mul(A, B):
    n = len(A)
    return [[sum(A[i][k]*B[k][j] for k in range(n)) for j in range(n)] for i in range(n)]

# CMF scalar over Z
N_MAX = 20
prod = [[1,0,0],[0,1,0],[0,0,1]]
q = [1]
for N in range(N_MAX + 3):
    prod = mat_mul(prod, M_int(N))
    q.append(prod[0][0])

# Delannoy over Q
D = [F(1), F(3)]
E = [F(0), F(1)]
for n in range(1, N_MAX + 10):
    D.append((F(3*(2*n+1)) * D[n] - F(n) * D[n-1]) / F(n+1))
    E.append((F(3*(2*n+1)) * E[n] - F(n) * E[n-1]) / F(n+1))

def solve_3x3_exact(M, b):
    det = (M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1])
         - M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0])
         + M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0]))
    if det == 0:
        return None
    x = [F(0)]*3
    for j in range(3):
        Mj = [[M[i][k] if k != j else b[i] for k in range(3)] for i in range(3)]
        det_j = (Mj[0][0]*(Mj[1][1]*Mj[2][2]-Mj[1][2]*Mj[2][1])
               - Mj[0][1]*(Mj[1][0]*Mj[2][2]-Mj[1][2]*Mj[2][0])
               + Mj[0][2]*(Mj[1][0]*Mj[2][1]-Mj[1][1]*Mj[2][0]))
        x[j] = det_j / det
    return x

q_F = [F(x) for x in q]

print("=== Exact Casorati connection coefficients ===\n")
r0_vals = []
gauge_ratios = []

for N in range(N_MAX):
    Phi = [[D[N]**2, D[N]*E[N], E[N]**2],
           [D[N+1]**2, D[N+1]*E[N+1], E[N+1]**2],
           [D[N+2]**2, D[N+2]*E[N+2], E[N+2]**2]]
    b = [q_F[N], q_F[N+1], q_F[N+2]]
    r = solve_3x3_exact(Phi, b)
    if r is None:
        print(f"N={N}: singular")
        r0_vals.append(None)
        continue

    r0 = r[0]
    r0_vals.append(r0)

    if N <= 5:
        print(f"N={N}: r₀ = {r0}")
        print(f"  numerator: {r0.numerator}")
        print(f"  denominator: {r0.denominator}")

    if N >= 1 and r0_vals[N-1] is not None and r0_vals[N-1] != 0:
        ratio = r0 / r0_vals[N-1]
        gauge_ratios.append((N, ratio))
        if N <= 10:
            print(f"N={N}: gauge_ratio = {ratio}")
            print(f"  num: {ratio.numerator}")
            print(f"  den: {ratio.denominator}")

# The gauge ratios should be EXACT rational numbers.
# Let me factor them using sympy.
print("\n\n=== Factoring gauge ratios ===")
Nsym = Symbol('N')

for N, ratio in gauge_ratios[:10]:
    print(f"\nN={N}: ratio = {ratio.numerator}/{ratio.denominator}")

    # Check: is the denominator related to Casorati det?
    # Casorati det(Phi(N)) = Sym² Wronskian
    Phi = [[D[N]**2, D[N]*E[N], E[N]**2],
           [D[N+1]**2, D[N+1]*E[N+1], E[N+1]**2],
           [D[N+2]**2, D[N+2]*E[N+2], E[N+2]**2]]
    det_Phi = (Phi[0][0]*(Phi[1][1]*Phi[2][2]-Phi[1][2]*Phi[2][1])
             - Phi[0][1]*(Phi[1][0]*Phi[2][2]-Phi[1][2]*Phi[2][0])
             + Phi[0][2]*(Phi[1][0]*Phi[2][1]-Phi[1][1]*Phi[2][0]))
    nd = len(str(abs(det_Phi.numerator)))
    dd = len(str(abs(det_Phi.denominator)))
    print(f"  det(Φ({N})): {nd}d num / {dd}d den")

# Now compute Casorati determinant as function of N
print("\n\n=== Casorati determinant values ===")
for N in range(15):
    Phi = [[D[N]**2, D[N]*E[N], E[N]**2],
           [D[N+1]**2, D[N+1]*E[N+1], E[N+1]**2],
           [D[N+2]**2, D[N+2]*E[N+2], E[N+2]**2]]
    det_val = (Phi[0][0]*(Phi[1][1]*Phi[2][2]-Phi[1][2]*Phi[2][1])
             - Phi[0][1]*(Phi[1][0]*Phi[2][2]-Phi[1][2]*Phi[2][0])
             + Phi[0][2]*(Phi[1][0]*Phi[2][1]-Phi[1][1]*Phi[2][0]))
    print(f"  N={N}: det = {det_val}")

# The Casorati det satisfies: det(Phi(N+1))/det(Phi(N)) = det(L_Sym2 companion at N+2)
# = leading coefficient ratio of L_Sym2
print("\n=== Casorati det ratio ===")
det_vals = []
for N in range(15):
    Phi = [[D[N]**2, D[N]*E[N], E[N]**2],
           [D[N+1]**2, D[N+1]*E[N+1], E[N+1]**2],
           [D[N+2]**2, D[N+2]*E[N+2], E[N+2]**2]]
    det_val = (Phi[0][0]*(Phi[1][1]*Phi[2][2]-Phi[1][2]*Phi[2][1])
             - Phi[0][1]*(Phi[1][0]*Phi[2][2]-Phi[1][2]*Phi[2][0])
             + Phi[0][2]*(Phi[1][0]*Phi[2][1]-Phi[1][1]*Phi[2][0]))
    det_vals.append(det_val)
    if N >= 1:
        ratio = det_val / det_vals[N-1]
        print(f"  N={N}: det(N)/det(N-1) = {ratio}")

print("\nDone.")
