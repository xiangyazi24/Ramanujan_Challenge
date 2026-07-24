#!/usr/bin/env python3
"""
Problem 2.5: Use Casorati matrix of Sym²(Delannoy) to find connection coefficients.

For each N, solve:
    Φ(N) · r(N) = [q(N), q(N+1), q(N+2)]^T

where Φ(N) is the Casorati matrix of {D², DE, E²}.

If q(N) = gauge(N) · (α·D_N² + β·D_N·E_N + γ·E_N²), then
r₀(N) = gauge(N)·α, r₁(N) = gauge(N+1)·α, r₂(N) = gauge(N+2)·α
(approximately). So the RATIO r₀(N+1)/r₀(N) gives gauge(N+1)/gauge(N).
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
    return [[m11, m12, m13], [m21, m22, m23], [m31, m32, m33]]

def mat_mul(A, B):
    n = len(A)
    return [[sum(A[i][k]*B[k][j] for k in range(n)) for j in range(n)] for i in range(n)]

# Compute CMF scalar q(N)
N_MAX = 30
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

# Sym² solutions
def v1(n): return D[n]**2
def v2(n): return D[n]*E[n]
def v3(n): return E[n]**2

def solve_3x3(M, b):
    """Solve M·x = b over Q (Cramer's rule)."""
    det = (M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1])
         - M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0])
         + M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0]))
    if det == 0:
        return None

    x = [F(0)]*3
    for j in range(3):
        # Replace column j with b
        Mj = [[M[i][k] if k != j else b[i] for k in range(3)] for i in range(3)]
        det_j = (Mj[0][0]*(Mj[1][1]*Mj[2][2]-Mj[1][2]*Mj[2][1])
               - Mj[0][1]*(Mj[1][0]*Mj[2][2]-Mj[1][2]*Mj[2][0])
               + Mj[0][2]*(Mj[1][0]*Mj[2][1]-Mj[1][1]*Mj[2][0]))
        x[j] = det_j / det
    return x

q_F = [F(x) for x in q]

print("=== Casorati Connection Coefficients ===")
print("Solving Φ(N)·r(N) = [q(N), q(N+1), q(N+2)]^T\n")

r_vals = []
for N in range(N_MAX):
    Phi = [[v1(N), v2(N), v3(N)],
           [v1(N+1), v2(N+1), v3(N+1)],
           [v1(N+2), v2(N+2), v3(N+2)]]
    b = [q_F[N], q_F[N+1], q_F[N+2]]
    r = solve_3x3(Phi, b)
    if r is None:
        print(f"N={N}: Casorati singular!")
        r_vals.append(None)
        continue
    r_vals.append(r)

    if N <= 5:
        print(f"N={N}: r = [{float(r[0]):.6e}, {float(r[1]):.6e}, {float(r[2]):.6e}]")
        for k in range(3):
            nd = len(str(abs(r[k].numerator)))
            dd = len(str(abs(r[k].denominator)))
            print(f"  r[{k}]: {nd} num digits, {dd} den digits")

# Key diagnostic: r₀(N+1)/r₀(N) should be gauge(N+1)/gauge(N)
print("\n=== Gauge ratio: r₀(N+1)/r₀(N) ===")
print("If this is a rational function of N with integer coefficients, we win.\n")

gauge_ratios = []
for N in range(1, min(N_MAX, 25)):
    if r_vals[N] is None or r_vals[N-1] is None:
        continue
    if r_vals[N-1][0] == 0:
        continue
    ratio = r_vals[N][0] / r_vals[N-1][0]
    gauge_ratios.append((N, ratio))
    print(f"N={N}: r₀(N)/r₀(N-1) = {float(ratio):.10f}")

# Try to fit ratio as a polynomial/rational function
# First check: ratio / (-16·N^7)
print("\n=== ratio / (-16 · N^7) ===")
for N, ratio in gauge_ratios:
    norm = ratio / (F(-16) * F(N)**7)
    print(f"  N={N}: {float(norm):.10f}")

# Let me also check ratio / [(-16) * (2N+7)!/(2N-1)! * N^k...]
# From the factored recurrence, the gauge ratio should involve
# (N+1)(N+3/2)(N+2)(N+5/2)(N+3)(N+7/2)(N+4) / (...)

# First let me check if ratio is related to c₀/c₃ somehow
print("\n=== ratio vs c₀(N-1)/c₃(N-1) ===")

def c0_coeff(N):
    """Leading coefficient of c₀(N) from factored form."""
    # c₀(N) = (1/768)·(N+1)(N+2)(N+3)^5(N+4)^3(N+5)^2(2N+3)(2N+5)^2(2N+7)^4(2N+9)^3·P₆(N)
    P6 = 3072*N**6 + 74112*N**5 + 738544*N**4 + 3890106*N**3 + 11417947*N**2 + 17696904*N + 11307715
    return F(1, 768) * F(N+1) * F(N+2) * F(N+3)**5 * F(N+4)**3 * F(N+5)**2 * F(2*N+3) * F(2*N+5)**2 * F(2*N+7)**4 * F(2*N+9)**3 * F(P6)

def c3_coeff(N):
    P6p = 3072*N**6 + 55680*N**5 + 414064*N**4 + 1615610*N**3 + 3483853*N**2 + 3929280*N + 1806156
    return F(1, 6144) * F(2*N+7) * F(P6p)

for N, ratio in gauge_ratios:
    c0_val = c0_coeff(N-1)
    c3_val = c3_coeff(N-1)
    if c3_val != 0:
        r = c0_val / c3_val
        print(f"  N={N}: gauge_ratio/c₀_c₃ = {float(ratio / r):.10e}")

# Now check: does r₁/r₀ stay constant?
print("\n=== r₁(N)/r₀(N) — should be constant if β/α ===")
for N in range(min(N_MAX, 15)):
    if r_vals[N] is None or r_vals[N][0] == 0:
        continue
    ratio10 = r_vals[N][1] / r_vals[N][0]
    ratio20 = r_vals[N][2] / r_vals[N][0]
    print(f"  N={N}: r₁/r₀ = {float(ratio10):.10f}, r₂/r₀ = {float(ratio20):.10f}")

print("\n\n=== CRITICAL TEST: Does r₁(N)/r₀(N) converge? ===")
for N in range(5, min(N_MAX, 25)):
    if r_vals[N] is None or r_vals[N][0] == 0:
        continue
    ratio10 = r_vals[N][1] / r_vals[N][0]
    ratio20 = r_vals[N][2] / r_vals[N][0]
    print(f"  N={N}: r₁/r₀ = {float(ratio10):.15f}")

print("\nDone.")
