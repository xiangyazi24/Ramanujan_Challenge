#!/usr/bin/env python3
"""
Problem 2.5: p-adic analysis of gauge r₀(N).

Compute v_p(r₀(N)) for small primes to detect Pochhammer structure.
Also compute v_p of the gauge RATIO r₀(N)/r₀(N-1) directly.
"""
from fractions import Fraction as F

# Delannoy
D = [F(1), F(3)]
E = [F(0), F(1)]
for n in range(1, 50):
    D.append((F(3*(2*n+1)) * D[n] - F(n) * D[n-1]) / F(n+1))
    E.append((F(3*(2*n+1)) * E[n] - F(n) * E[n-1]) / F(n+1))

DD = [D[n]**2 for n in range(50)]
DE = [D[n]*E[n] for n in range(50)]
EE = [E[n]**2 for n in range(50)]

# CMF scalar
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

N_MAX = 20
prod = [[1,0,0],[0,1,0],[0,0,1]]
q = [1]
for N in range(N_MAX + 5):
    prod = mat_mul(prod, M_int(N))
    q.append(prod[0][0])
q_F = [F(x) for x in q]

# Casorati solve
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

r0_vals = []
for N in range(N_MAX):
    Phi = [[DD[N], DE[N], EE[N]],
           [DD[N+1], DE[N+1], EE[N+1]],
           [DD[N+2], DE[N+2], EE[N+2]]]
    b = [q_F[N], q_F[N+1], q_F[N+2]]
    r = solve_3x3_exact(Phi, b)
    r0_vals.append(r[0] if r else None)

# p-adic valuation
def vp(n, p):
    if n == 0:
        return float('inf')
    n = abs(n)
    v = 0
    while n % p == 0:
        n //= p
        v += 1
    return v

primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]

# Print r₀(N) integer values
print("=== r₀(N) exact integers ===")
for N in range(min(N_MAX, 12)):
    r = r0_vals[N]
    if r is None:
        print(f"N={N}: None")
        continue
    assert r.denominator == 1, f"N={N}: r₀ not integer! den={r.denominator}"
    print(f"N={N}: sign={'+ ' if r >= 0 else '- '} #digits={len(str(abs(r.numerator)))}")

# p-adic valuations of r₀(N)
print("\n=== v_p(r₀(N)) ===")
header = "N   " + "  ".join(f"p={p:>2}" for p in primes)
print(header)
for N in range(min(N_MAX, 12)):
    r = r0_vals[N]
    if r is None or r == 0:
        continue
    val = abs(r.numerator)
    vals = [vp(val, p) for p in primes]
    line = f"{N:<4}" + "  ".join(f"{v:>4}" for v in vals)
    print(line)

# p-adic valuations of gauge RATIO r₀(N)/r₀(N-1)
print("\n=== v_p(r₀(N)/r₀(N-1)) = v_p(r₀(N)) - v_p(r₀(N-1)) ===")
header = "N   " + "  ".join(f"p={p:>2}" for p in primes)
print(header)
for N in range(1, min(N_MAX, 12)):
    r = r0_vals[N]
    rp = r0_vals[N-1]
    if r is None or rp is None or r == 0 or rp == 0:
        continue
    vals_N = [vp(abs(r.numerator), p) for p in primes]
    vals_Nm1 = [vp(abs(rp.numerator), p) for p in primes]
    diffs = [vals_N[i] - vals_Nm1[i] for i in range(len(primes))]
    line = f"{N:<4}" + "  ".join(f"{d:>4}" for d in diffs)
    print(line)

# Now compare with expected Pochhammer patterns
# For (N!)^k: v_p = k * (N - s_p(N))/(p-1)
# For (2N+a)!!: v_2 follows a different pattern
# Let's compute v_p for common Pochhammer products

print("\n=== Comparison: v_p of Pochhammer products ===")
print("\nv_p(N!) for reference:")
from math import factorial
header = "N   " + "  ".join(f"p={p:>2}" for p in primes[:5])
print(header)
for N in range(12):
    vals = [vp(factorial(N), p) for p in primes[:5]]
    line = f"{N:<4}" + "  ".join(f"{v:>4}" for v in vals)
    print(line)

# Compute the Wronskian factor: W_L25(N) = W_L25(0) * prod_{k=0}^{N-1} (-c3(k)/c0(k))
# The Wronskian of L₂₅ is:
# W(N+1)/W(N) = (-1)^{ord} * c_{ord}(N)/c_0(N) = -c₃(N)/c₀(N) (for order 3)
def P6(N):
    return 3072*N**6 + 74112*N**5 + 738544*N**4 + 3890106*N**3 + 11417947*N**2 + 17696904*N + 11307715

def P6p(N):
    return 3072*N**6 + 55680*N**5 + 414064*N**4 + 1615610*N**3 + 3483853*N**2 + 3929280*N + 1806156

def c0_num(N):
    """Numerator of c₀(N) * 768"""
    return (N+1)*(N+2)*(N+3)**5*(N+4)**3*(N+5)**2*(2*N+3)*(2*N+5)**2*(2*N+7)**4*(2*N+9)**3*P6(N)

def c3_num(N):
    """Numerator of c₃(N) * 6144"""
    return (2*N+7)*P6p(N)

# Wronskian ratio = -c₃/c₀ = -(c3_num * 768) / (c0_num * 6144) = -c3_num / (8 * c0_num)
print("\n\n=== v_p of Wronskian ratio numerator/denominator ===")
print("W(N+1)/W(N) = -c₃(N)/c₀(N)")
print("  = -(2N+7)·P₆'(N)·768 / [(N+1)(N+2)(N+3)⁵(N+4)³(N+5)²(2N+3)(2N+5)²(2N+7)⁴(2N+9)³·P₆(N)·6144]")
print("  = -P₆'(N) / [8·(N+1)(N+2)(N+3)⁵(N+4)³(N+5)²(2N+3)(2N+5)²(2N+7)³(2N+9)³·P₆(N)]")
print()

# Since P₆(N) = P₆'(N+1), the product telescopes:
# prod_{k=0}^{N-1} P₆'(k)/P₆(k) = prod P₆'(k)/P₆'(k+1) = P₆'(0)/P₆'(N)
P6p_0 = P6p(0)
print(f"P₆'(0) = {P6p_0}")
print(f"P₆(0) = P₆'(1) = {P6(0)} = {P6p(1)}")
print(f"Telescoping check: P₆'(0) / P₆'(1) = P₆'(0)/P₆(0) = {F(P6p(0), P6(0))}")

# Compute Wronskian product W(N)/W(0) = prod_{k=0}^{N-1} (-c₃(k)/c₀(k))
W_ratio = [F(1)]
for k in range(15):
    c3 = F(1, 6144) * F(2*k+7) * F(P6p(k))
    c0 = F(1, 768) * F(k+1) * F(k+2) * F(k+3)**5 * F(k+4)**3 * F(k+5)**2 * F(2*k+3) * F(2*k+5)**2 * F(2*k+7)**4 * F(2*k+9)**3 * F(P6(k))
    ratio = -c3 / c0
    W_ratio.append(W_ratio[-1] * ratio)

print("\n=== v_p(W(N)/W(0)) = v_p(prod Wronskian ratios) ===")
header = "N   " + "  ".join(f"p={p:>2}" for p in primes[:5])
print(header)
for N in range(12):
    w = W_ratio[N]
    if w == 0:
        continue
    num_vals = [vp(abs(w.numerator), p) for p in primes[:5]]
    den_vals = [vp(w.denominator, p) for p in primes[:5]]
    net = [num_vals[i] - den_vals[i] for i in range(5)]
    line = f"{N:<4}" + "  ".join(f"{v:>4}" for v in net)
    print(line)

# Now compute r₀(N) / W(N) to see what the "extra" factor is
print("\n=== v_p(r₀(N)) - v_p(W(N)/W(0)) ===")
print("(This should reveal the non-Wronskian part of the gauge)")
header = "N   " + "  ".join(f"p={p:>2}" for p in primes[:5])
print(header)
for N in range(1, min(N_MAX, 10)):
    r = r0_vals[N]
    if r is None or r == 0:
        continue
    r_num = abs(r.numerator)
    w = W_ratio[N]
    w_num = abs(w.numerator)
    w_den = w.denominator

    r_vals = [vp(r_num, p) for p in primes[:5]]
    w_net = [vp(w_num, p) - vp(w_den, p) for p in primes[:5]]
    diff = [r_vals[i] - w_net[i] for i in range(5)]
    line = f"{N:<4}" + "  ".join(f"{v:>4}" for v in diff)
    print(line)

# Also compute the Sym² Casorati determinant's contribution
# det_Sym²(N) = 3(2N+3)/[(N+1)²(N+2)²]
# Accumulated: prod_{k=0}^{N-1} det(k+1)/det(k)
print("\n=== Sym² Casorati det product ===")
det_prod = [F(1)]
for k in range(15):
    # det(k) = 3(2k+3) / [(k+1)²(k+2)²]
    det_k = F(3*(2*k+3), (k+1)**2 * (k+2)**2)
    # det(k+1)/det(k)
    det_kp1 = F(3*(2*k+5), (k+2)**2 * (k+3)**2)
    ratio = det_kp1 / det_k
    det_prod.append(det_prod[-1] * ratio)

header = "N   " + "  ".join(f"p={p:>2}" for p in primes[:5])
print(header)
for N in range(12):
    dp = det_prod[N]
    if dp == 0:
        continue
    num_vals = [vp(abs(dp.numerator), p) for p in primes[:5]]
    den_vals = [vp(dp.denominator, p) for p in primes[:5]]
    net = [num_vals[i] - den_vals[i] for i in range(5)]
    line = f"{N:<4}" + "  ".join(f"{v:>4}" for v in net)
    print(line)

# Final attempt: just print the exact gauge ratios in lowest terms
# and try to match them to rational functions by checking specific forms
print("\n\n=== Exact gauge ratios ===")
for N in range(1, min(N_MAX, 8)):
    r = r0_vals[N]
    rp = r0_vals[N-1]
    if r is None or rp is None or rp == 0:
        continue
    ratio = r / rp
    from math import gcd
    g = gcd(abs(ratio.numerator), ratio.denominator)
    num = ratio.numerator // g
    den = ratio.denominator // g
    print(f"\nN={N}: R = {num}")
    print(f"       / {den}")
    print(f"  #digits: num={len(str(abs(num)))}, den={len(str(den))}")

    # Factor out small primes
    for p in [2, 3, 5, 7]:
        vn = vp(abs(num), p)
        vd = vp(den, p)
        if vn > 0 or vd > 0:
            print(f"  v_{p}: num={vn}, den={vd}, net={vn-vd}")

print("\nDone.")
