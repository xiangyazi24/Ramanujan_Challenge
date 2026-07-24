#!/usr/bin/env python3
"""
Factor the exact gauge ratios r₀(N)/r₀(N-1) to identify the rational function.
"""
from fractions import Fraction as F
from sympy import factorint

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
for N in range(N_MAX + 3):
    prod = mat_mul(prod, M_int(N))
    q.append(prod[0][0])

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
r0_vals = []
for N in range(N_MAX):
    Phi = [[D[N]**2, D[N]*E[N], E[N]**2],
           [D[N+1]**2, D[N+1]*E[N+1], E[N+1]**2],
           [D[N+2]**2, D[N+2]*E[N+2], E[N+2]**2]]
    b = [q_F[N], q_F[N+1], q_F[N+2]]
    r = solve_3x3_exact(Phi, b)
    r0_vals.append(r[0] if r else None)

print("=== Factoring gauge ratios ===\n")
for N in range(1, min(N_MAX, 12)):
    if r0_vals[N] is None or r0_vals[N-1] is None or r0_vals[N-1] == 0:
        continue
    ratio = r0_vals[N] / r0_vals[N-1]
    num = abs(ratio.numerator)
    den = ratio.denominator
    sign = -1 if ratio < 0 else 1

    print(f"N={N}: ratio = {'+' if sign > 0 else '-'}{num}/{den}")

    num_factors = factorint(num) if num > 1 else {}
    den_factors = factorint(den) if den > 1 else {}

    # Format nicely
    def format_factors(d):
        parts = []
        for p in sorted(d.keys()):
            if d[p] == 1:
                parts.append(str(p))
            else:
                parts.append(f"{p}^{d[p]}")
        return " · ".join(parts) if parts else "1"

    print(f"  num = {format_factors(num_factors)}")
    print(f"  den = {format_factors(den_factors)}")
    print()

# Also factor r₀(N) directly
print("\n=== Factoring r₀(N) ===\n")
for N in range(min(N_MAX, 8)):
    if r0_vals[N] is None:
        continue
    val = r0_vals[N]
    if val == 0:
        print(f"N={N}: r₀ = 0")
        continue
    v = abs(val.numerator)
    sign = -1 if val < 0 else 1
    factors = factorint(v) if v > 1 else {}

    def format_factors(d):
        parts = []
        for p in sorted(d.keys()):
            if d[p] == 1:
                parts.append(str(p))
            else:
                parts.append(f"{p}^{d[p]}")
        return " · ".join(parts) if parts else "1"

    print(f"N={N}: r₀ = {'+' if sign > 0 else '-'}{format_factors(factors)} (den={val.denominator})")

# Factor the CMF Casorati det (det of 3×3 Casorati of q₀, q₁, q₂)
# First compute q₀, q₁, q₂ (three independent CMF solutions from rows 0,1,2)
print("\n\n=== CMF Casorati determinant ===")
print("(Using q(N) Casorati: det[[q(N),q(N+1),q(N+2)]] — but this is just 1 solution, need 3)")
print()

# Actually, use the L₂₅ Wronskian formula:
# W_{L25}(N+1)/W_{L25}(N) = (-1)^3 · c₃(N)/c₀(N) = -c₃(N)/c₀(N)
# where c₀, c₃ are from the recurrence

def P6(N):
    return 3072*N**6 + 74112*N**5 + 738544*N**4 + 3890106*N**3 + 11417947*N**2 + 17696904*N + 11307715

def P6p(N):
    return 3072*N**6 + 55680*N**5 + 414064*N**4 + 1615610*N**3 + 3483853*N**2 + 3929280*N + 1806156

def c0_exact(N):
    return F(1, 768) * F(N+1) * F(N+2) * F(N+3)**5 * F(N+4)**3 * F(N+5)**2 * F(2*N+3) * F(2*N+5)**2 * F(2*N+7)**4 * F(2*N+9)**3 * F(P6(N))

def c3_exact(N):
    return F(1, 6144) * F(2*N+7) * F(P6p(N))

print("L₂₅ Wronskian ratio = -c₃(N)/c₀(N):")
for N in range(8):
    ratio = -c3_exact(N) / c0_exact(N)
    print(f"  N={N}: -c₃/c₀ = {ratio}")
    # Factor
    num = abs(ratio.numerator)
    den = ratio.denominator
    sign = -1 if ratio < 0 else 1
    num_f = factorint(num) if num > 1 else {}
    den_f = factorint(den) if den > 1 else {}
    def ff(d):
        parts = []
        for p in sorted(d.keys()):
            if d[p] == 1:
                parts.append(str(p))
            else:
                parts.append(f"{p}^{d[p]}")
        return " · ".join(parts) if parts else "1"
    print(f"    = {'+' if sign > 0 else '-'}{ff(num_f)} / {ff(den_f)}")
    print(f"    P₆'(N)/P₆(N) factor = {F(P6p(N), P6(N))}")
    print()

# Compare gauge ratio with Wronskian ratio
print("\n=== Gauge ratio vs Wronskian ratio ===")
print("gauge_ratio(N) = r₀(N)/r₀(N-1)")
print("Wronskian_ratio(N-1) = -c₃(N-1)/c₀(N-1)")
print("Sym² Casorati ratio(N-1) = det_Sym²(N)/det_Sym²(N-1)")
print()
for N in range(2, min(N_MAX, 10)):
    if r0_vals[N] is None or r0_vals[N-1] is None or r0_vals[N-1] == 0:
        continue
    gr = r0_vals[N] / r0_vals[N-1]
    wr = -c3_exact(N-1) / c0_exact(N-1)
    # Sym² det ratio
    det_N = F(3*(2*N+3), (N+1)**2 * (N+2)**2)
    det_Nm1 = F(3*(2*(N-1)+3), N**2 * (N+1)**2)
    sr = det_N / det_Nm1

    print(f"N={N}:")
    print(f"  gauge_ratio = {float(gr):.6e}")
    print(f"  Wronsk_ratio = {float(wr):.6e}")
    print(f"  Sym²_det_ratio = {float(sr):.6f}")
    print(f"  gauge / (Wronsk · Sym²) = {float(gr / (wr * sr)):.10f}")
    print(f"  gauge / Wronsk = {float(gr / wr):.10f}")
    print()

print("Done.")
