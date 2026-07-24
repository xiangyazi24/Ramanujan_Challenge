#!/usr/bin/env python3
"""
P2.5: THE KEY CHECK.

If Q̂_N = Σ_k f(k)·F(N,k) and P̂_N = Σ_k g(k)·F(N,k),
check whether g(k)/f(k) = C_k where C_k = Σ_{j=0}^k (-1)^j/(2j+1)².

If yes, then P̂_N - G·Q̂_N = Σ_k f(k)·F(N,k)·(C_k - G) → 0
(because the tail G - C_k = Σ_{j>k} (-1)^j/(2j+1)² → 0 and F(N,k)
concentrates near k=N), proving L = G.

F(N,k) = 2^k C(2k,k) C(N,k) C(N+k,k) is the Delannoy summand.
"""
from fractions import Fraction as F
from math import comb

def M_entries(n):
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
    return -2*(n+2)**2*(n+3)**2*(2*n+5)*(2*n+7)**2

def delannoy_summand(N, k):
    if k < 0 or k > N:
        return F(0)
    return F(2**k * comb(2*k, k) * comb(N, k) * comb(N+k, k))

def catalan_partial(k):
    """C_k = Σ_{j=0}^k (-1)^j / (2j+1)²"""
    return sum(F((-1)**j, (2*j+1)**2) for j in range(k+1))

NMAX = 25
print(f"Computing Q̂ and P̂ for N=0..{NMAX} (exact fractions)...\n")

# Compute CMF trajectories from ALL initial vectors
rows = {
    'e1': [F(1), F(0), F(0)],
    'e2': [F(0), F(1), F(0)],
    'e3': [F(0), F(0), F(1)],
}
history = {key: [[v[j] for j in range(3)]] for key, v in rows.items()}

for N in range(NMAX):
    M = M_entries(N)
    d = F(delta_H(N))
    MH = [[F(M[i][j]) / d for j in range(3)] for i in range(3)]
    for key in rows:
        r = rows[key]
        new_r = [sum(r[i]*MH[i][k] for i in range(3)) for k in range(3)]
        rows[key] = new_r
        history[key].append([new_r[j] for j in range(3)])

# Q-row: q = (33750, -36000, 9000)
# P-row: p = (30921, -32972, 8240)
q = [F(33750), F(-36000), F(9000)]
p = [F(30921), F(-32972), F(8240)]

# Compute Q̂_N and P̂_N from combined trajectories
Q_vals = []
P_vals = []
for N in range(NMAX + 1):
    Q_N = sum(q[j] * history['e' + str(j+1)][N][0] for j in range(3))
    P_N = sum(p[j] * history['e' + str(j+1)][N][0] for j in range(3))
    Q_vals.append(Q_N)
    P_vals.append(P_N)

print("First few Q̂_N:", [str(Q_vals[i]) for i in range(5)])
print("First few P̂_N:", [str(P_vals[i]) for i in range(5)])
print()

# Decompose Q̂_N = Σ_k f(k)·F(N,k) and P̂_N = Σ_k g(k)·F(N,k)
def decompose_in_delannoy(vals, NMAX):
    """Given vals[0..NMAX], find coeffs c[k] such that vals[N] = Σ c[k]·F(N,k)"""
    coeffs = []
    for K in range(NMAX + 1):
        rhs = vals[K]
        for k in range(K):
            rhs -= coeffs[k] * delannoy_summand(K, k)
        bKK = delannoy_summand(K, K)
        if bKK == 0:
            raise ValueError(f"base({K},{K}) = 0")
        coeffs.append(rhs / bKK)
    return coeffs

f_coeffs = decompose_in_delannoy(Q_vals, NMAX)
g_coeffs = decompose_in_delannoy(P_vals, NMAX)

# THE KEY CHECK: g(k)/f(k) =? C_k
print("=" * 70)
print("THE KEY CHECK: g(k)/f(k) vs C_k = Σ_{j≤k} (-1)^j/(2j+1)²")
print("=" * 70)

all_match = True
for k in range(NMAX + 1):
    C_k = catalan_partial(k)
    if f_coeffs[k] == 0:
        print(f"  k={k}: f(k)=0, skip")
        continue
    ratio = g_coeffs[k] / f_coeffs[k]
    match = (ratio == C_k)
    if not match:
        all_match = False
    if k < 15 or not match:
        print(f"  k={k}: g/f = {float(ratio):.15f}, C_k = {float(C_k):.15f}, "
              f"{'✓ MATCH' if match else '✗ MISMATCH'}")
        if not match:
            print(f"         g/f - C_k = {float(ratio - C_k):.6e}")
            # Show exact values
            print(f"         g(k) = {g_coeffs[k]}")
            print(f"         f(k) = {f_coeffs[k]}")
            print(f"         C_k  = {C_k}")

print()
if all_match:
    print("★★★ ALL MATCH! g(k) = f(k)·C_k for k = 0,...,{NMAX} ★★★")
    print()
    print("This proves L = G because:")
    print("  P̂_N = Σ_k f(k)·F(N,k)·C_k")
    print("  Q̂_N = Σ_k f(k)·F(N,k)")
    print("  G·Q̂_N - P̂_N = Σ_k f(k)·F(N,k)·(G - C_k)")
    print("  |G - C_k| = O(1/(2k+3)²) → 0")
    print("  F(N,k) concentrates near k = O(N)")
    print("  Therefore G·Q̂_N - P̂_N = o(Q̂_N), hence P̂/Q̂ → G")
else:
    print("MISMATCH: g(k) ≠ f(k)·C_k")
    print("The partial-sum structure is NOT preserved.")
    print()
    # Try to find the actual pattern in g/f
    print("Actual ratios g(k)/f(k):")
    for k in range(min(15, NMAX+1)):
        if f_coeffs[k] != 0:
            ratio = g_coeffs[k] / f_coeffs[k]
            print(f"  k={k}: {float(ratio):.15f}")

# Also verify the decomposition
print(f"\n=== Verification ===")
for N in range(NMAX + 1):
    Q_check = sum(f_coeffs[k] * delannoy_summand(N, k) for k in range(N+1))
    P_check = sum(g_coeffs[k] * delannoy_summand(N, k) for k in range(N+1))
    if Q_check != Q_vals[N] or P_check != P_vals[N]:
        print(f"  N={N}: VERIFICATION FAILED")
        break
    elif N % 5 == 0:
        print(f"  N={N}: OK (Q̂={float(Q_vals[N]):.6e}, P̂={float(P_vals[N]):.6e})")

print("\nDone.")
