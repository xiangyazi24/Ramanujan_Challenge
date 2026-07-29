#!/usr/bin/env python3
"""
Check whether β = r₁/r₀ and γ = r₂/r₀ are EXACTLY constant across N.
Uses exact Fraction arithmetic.
"""
from fractions import Fraction as F

D = [F(1), F(3)]
E = [F(0), F(1)]
for n in range(1, 35):
    D.append((F(3*(2*n+1)) * D[n] - F(n) * D[n-1]) / F(n+1))
    E.append((F(3*(2*n+1)) * E[n] - F(n) * E[n-1]) / F(n+1))

DD = [D[n]**2 for n in range(35)]
DE = [D[n]*E[n] for n in range(35)]
EE = [E[n]**2 for n in range(35)]

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

N_MAX = 15
prod = [[1,0,0],[0,1,0],[0,0,1]]
q = [1]
for N in range(N_MAX + 5):
    prod = mat_mul(prod, M_int(N))
    q.append(prod[0][0])
q_F = [F(x) for x in q]

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

print("=== Exact β(N) = r₁(N)/r₀(N) and γ(N) = r₂(N)/r₀(N) ===\n")
betas = []
gammas = []
for N in range(N_MAX):
    Phi = [[DD[N], DE[N], EE[N]],
           [DD[N+1], DE[N+1], EE[N+1]],
           [DD[N+2], DE[N+2], EE[N+2]]]
    b = [q_F[N], q_F[N+1], q_F[N+2]]
    r = solve_3x3_exact(Phi, b)
    if r is None or r[0] == 0:
        betas.append(None)
        gammas.append(None)
        continue
    beta_N = r[1] / r[0]
    gamma_N = r[2] / r[0]
    betas.append(beta_N)
    gammas.append(gamma_N)
    print(f"N={N}: β = {float(beta_N):.20e}")
    print(f"       γ = {float(gamma_N):.20e}")

# Check if β(N) is EXACTLY constant
print("\n=== Check: β(N) - β(N-1) (should be 0 if constant) ===")
for N in range(1, min(N_MAX, 10)):
    if betas[N] is None or betas[N-1] is None:
        continue
    diff = betas[N] - betas[N-1]
    if diff == 0:
        print(f"N={N}: β(N) - β(N-1) = EXACTLY 0")
    else:
        print(f"N={N}: β(N) - β(N-1) = {float(diff):.6e}  (NOT zero!)")
        print(f"  num digits: {len(str(abs(diff.numerator)))}")

print("\n=== Check: γ(N) - γ(N-1) ===")
for N in range(1, min(N_MAX, 10)):
    if gammas[N] is None or gammas[N-1] is None:
        continue
    diff = gammas[N] - gammas[N-1]
    if diff == 0:
        print(f"N={N}: γ(N) - γ(N-1) = EXACTLY 0")
    else:
        print(f"N={N}: γ(N) - γ(N-1) = {float(diff):.6e}  (NOT zero!)")

# If not constant, check what pattern the differences follow
# diff_β(N) = β(N) - β(N-1), does diff_β(N+1)/diff_β(N) → constant?
if any(betas[N] is not None and betas[N-1] is not None and betas[N] != betas[N-1]
       for N in range(1, min(N_MAX, 10))):
    print("\n=== Ratio of consecutive β-differences (should → convergence rate) ===")
    diffs = []
    for N in range(1, min(N_MAX, 10)):
        if betas[N] is None or betas[N-1] is None:
            diffs.append(None)
            continue
        diffs.append(betas[N] - betas[N-1])
    for N in range(1, len(diffs)):
        if diffs[N] is None or diffs[N-1] is None or diffs[N-1] == 0:
            continue
        ratio = diffs[N] / diffs[N-1]
        print(f"  diff_β({N+1})/diff_β({N}) = {float(ratio):.15f}")

print("\nDone.")
