#!/usr/bin/env python3
"""
CORRECTED adjoint bracket computation for P2.7.

Key corrections from v1:
1. Correct adjoint recurrence: p₀(n+3) w_{n+3} + p₁(n+2) w_{n+2} + p₂(n+1) w_{n+1} + p₃(n) w_n = 0
2. Correct Lagrange bilinear form (derived from first principles)
3. Correct Poincaré root: μ₀ ≈ 0.859 (not 54.96 — the 54.96 is the SCALED root ν₀ = 64μ₀)
"""
from mpmath import mp, mpf, pi, zeta, log, fabs

mp.dps = 250

# ===== P2.7 recurrence coefficients =====
def A(n):
    n = mpf(n)
    return mpf(1024) * (2*n+5)**4 * (2*n+7)**3 * (2*n+9)**3 * (946*n**2 + 6407*n + 10860)

def B(n):
    n = mpf(n)
    return mpf(128) * (2*n+7)**3 * (2*n+9)**3 * (104060*n**6 + 1745370*n**5 + 12145238*n**4 + 44886481*n**3 + 92943995*n**2 + 102256019*n + 46709052)

def C(n):
    n = mpf(n)
    return mpf(16) * (n+3)**4 * (2*n+9)**3 * (3784*n**5 + 57792*n**4 + 351019*n**3 + 1059230*n**2 + 1587211*n + 944620)

def D(n):
    n = mpf(n)
    return (n+3)**4 * (n+4)**6 * (946*n**2 + 4515*n + 5399)

# Operator: p₃(n) u_{n+3} + p₂(n) u_{n+2} + p₁(n) u_{n+1} + p₀(n) u_n = 0
# where p₃(n) = A(n+2), p₂(n) = -B(n+2), p₁(n) = C(n+1), p₀(n) = -D(n)
def p3(n): return A(n+2)
def p2(n): return -B(n+2)
def p1(n): return C(n+1)
def p0(n): return -D(n)

# ===== Forward solutions =====
N = 400

q = [mpf(0)] * (N + 10)
q[0] = mpf(-215040420000)
q[1] = mpf(-167282265043404) / mpf(905)
q[2] = mpf(-964185327658080) / mpf(6071)
for n in range(2, N + 5):
    q[n+1] = (B(n) * q[n] - C(n-1) * q[n-1] + D(n-2) * q[n-2]) / A(n)

p = [mpf(0)] * (N + 10)
p[0] = mpf(-612218384750)
p[1] = mpf(-9525021973931919) / mpf(18100)
p[2] = mpf(-29561828382772029) / mpf(65380)
for n in range(2, N + 5):
    p[n+1] = (B(n) * p[n] - C(n-1) * p[n-1] + D(n-2) * p[n-2]) / A(n)

L_val = zeta(2) + zeta(3)
print(f"ζ(2)+ζ(3) = {L_val}")
print(f"q[0] = {q[0]}")
print(f"|q[50]/q[0]|^(1/50) = {float(fabs(q[50]/q[0])**(mpf(1)/50)):.10f}")
print(f"|q[100]/q[0]|^(1/100) = {float(fabs(q[100]/q[0])**(mpf(1)/100)):.10f}")

for n in [5, 10, 20, 50, 100]:
    e_n = p[n] - L_val * q[n]
    ratio = fabs(e_n / q[n])
    print(f"  |e_{n}/q_{n}| = {float(ratio):.6e}")

# ===== Correct adjoint recurrence =====
# L†[w]_n = p₀(n) w_n + p₁(n-1) w_{n-1} + p₂(n-2) w_{n-2} + p₃(n-3) w_{n-3} = 0
# S³-shifted forward form:
# p₀(n+3) w_{n+3} + p₁(n+2) w_{n+2} + p₂(n+1) w_{n+1} + p₃(n) w_n = 0
# i.e. -D(n+3) w_{n+3} + C(n+3) w_{n+2} - B(n+3) w_{n+1} + A(n+2) w_n = 0
# Forward: w_{n+3} = [C(n+3) w_{n+2} - B(n+3) w_{n+1} + A(n+2) w_n] / D(n+3)
# Backward (Miller): w_n = [D(n+3) w_{n+3} - C(n+3) w_{n+2} + B(n+3) w_{n+1}] / A(n+2)

M = N
w = [mpf(0)] * (M + 10)
w[M] = mpf(1)
w[M-1] = mpf(0)
w[M-2] = mpf(0)

print(f"\nMiller's algorithm (M={M})...")
for n in range(M-3, -1, -1):
    w[n] = (D(n+3) * w[n+3] - C(n+3) * w[n+2] + B(n+3) * w[n+1]) / A(n+2)

print(f"w[0] = {w[0]}")
print(f"w[1] = {w[1]}")
print(f"w[2] = {w[2]}")

# Growth analysis
print("\nAdjoint growth analysis (should approach 1/μ₀ ≈ 1.164 for forward, or decay for w[n]/w[n-1]):")
for n in [5, 10, 20, 50, 100, 200]:
    if w[n] != 0 and w[n-1] != 0:
        ratio = fabs(w[n] / w[n-1])
        print(f"  |w[{n}]/w[{n-1}]| = {float(ratio):.10f}")

# Check: |w[n]|^{1/n} should approach growth rate
if w[0] != 0:
    for n in [10, 20, 50, 100, 200]:
        if w[n] != 0:
            growth = fabs(w[n] / w[0]) ** (mpf(1)/n)
            print(f"  |w[{n}]/w[0]|^(1/{n}) = {float(growth):.10f}")

# ===== Verify adjoint recurrence =====
print("\nVerifying adjoint recurrence residuals:")
max_res = mpf(0)
for n in range(0, min(M-3, 50)):
    # p₀(n+3) w_{n+3} + p₁(n+2) w_{n+2} + p₂(n+1) w_{n+1} + p₃(n) w_n = 0
    res = p0(n+3)*w[n+3] + p1(n+2)*w[n+2] + p2(n+1)*w[n+1] + p3(n)*w[n]
    rr = fabs(res)
    if rr > max_res: max_res = rr
    if n < 3:
        print(f"  n={n}: residual = {float(res):.6e}")
print(f"  Max residual (n=0..49): {float(max_res):.6e}")

# ===== Correct Lagrange bilinear form =====
# J(w,u;n) = [p₁(n-1) w_{n-1} + p₂(n-2) w_{n-2} + p₃(n-3) w_{n-3}] u_n
#           + [p₂(n-1) w_{n-1} + p₃(n-2) w_{n-2}] u_{n+1}
#           + p₃(n-1) w_{n-1} u_{n+2}

def J_form(w, u, n):
    """Correct Lagrange bilinear concomitant."""
    t1 = (p1(n-1)*w[n-1] + p2(n-2)*w[n-2] + p3(n-3)*w[n-3]) * u[n]
    t2 = (p2(n-1)*w[n-1] + p3(n-2)*w[n-2]) * u[n+1]
    t3 = p3(n-1)*w[n-1] * u[n+2]
    return t1 + t2 + t3

# Verify constancy of J(w, q; n)
print("\n===== J(w, q; n) constancy check =====")
Jq_vals = []
for n in range(3, 30):
    Jn = J_form(w, q, n)
    Jq_vals.append(Jn)
    if n <= 6 or n >= 27:
        print(f"  J(w,q;{n}) = {Jn}")

if len(Jq_vals) >= 2:
    J0 = Jq_vals[0]
    max_var = max(fabs(J - J0) for J in Jq_vals[1:])
    rel_var = float(max_var / fabs(J0)) if J0 != 0 else float('inf')
    print(f"  Max variation from J(3): {float(max_var):.6e}")
    print(f"  Relative variation: {rel_var:.6e}")

# J(w, p; n) constancy check
print("\n===== J(w, p; n) constancy check =====")
Jp_vals = []
for n in range(3, 30):
    Jn = J_form(w, p, n)
    Jp_vals.append(Jn)
    if n <= 6 or n >= 27:
        print(f"  J(w,p;{n}) = {Jn}")

if len(Jp_vals) >= 2:
    J0p = Jp_vals[0]
    max_var_p = max(fabs(J - J0p) for J in Jp_vals[1:])
    rel_var_p = float(max_var_p / fabs(J0p)) if J0p != 0 else float('inf')
    print(f"  Max variation from J(3): {float(max_var_p):.6e}")
    print(f"  Relative variation: {rel_var_p:.6e}")

# ===== KEY TEST =====
print("\n" + "="*60)
print("KEY TEST: J(w,p;n) / J(w,q;n) = ζ(2)+ζ(3) ?")
print("="*60)

# Use n=3 (first valid index)
Jwp = J_form(w, p, 3)
Jwq = J_form(w, q, 3)

if Jwq != 0:
    ratio = Jwp / Jwq
    diff = ratio - L_val
    rel_err = fabs(diff / L_val)

    print(f"J(w,p;3) = {Jwp}")
    print(f"J(w,q;3) = {Jwq}")
    print(f"Ratio     = {ratio}")
    print(f"ζ(2)+ζ(3) = {L_val}")
    print(f"Difference = {float(diff):.6e}")
    print(f"Relative error = {float(rel_err):.6e}")

    # J(w,e;3) = J(w,p;3) - L * J(w,q;3)
    Jwe = Jwp - L_val * Jwq
    print(f"\nJ(w,e;3) = {float(Jwe):.6e}")
    print(f"|J(w,e;3)/J(w,q;3)| = {float(fabs(Jwe/Jwq)):.6e}")

# Also check at other n values
print("\n--- Ratio at different n values ---")
for n in [4, 5, 10, 20, 50]:
    Jp = J_form(w, p, n)
    Jq = J_form(w, q, n)
    if Jq != 0:
        r = Jp / Jq
        d = r - L_val
        print(f"  n={n}: ratio = {float(r):.15f}, diff = {float(d):.6e}")

# ===== Poincaré polynomial verification =====
print("\n===== Poincaré polynomial of P2.7 =====")
# Original: p₃_∞ μ³ + p₂_∞ μ² + p₁_∞ μ + p₀_∞ = 0
# With p₃ = A(n+2), p₂ = -B(n+2), p₁ = C(n+1), p₀ = -D(n)
# Leading: A_∞ = 1024·1024·946, B_∞ = 128·64·110·946, C_∞ = 128·4·946, D_∞ = 946
# After dividing by 946: 1048576 μ³ - 901120 μ² + 512 μ - 1 = 0
# Substitute ν = 64μ: 4ν³ - 220ν² + 8ν - 1 = 0

# Find exact roots
from mpmath import polyroots
roots = polyroots([4, -220, 8, -1])
print("Roots of 4ν³ - 220ν² + 8ν - 1 = 0:")
for r in roots:
    print(f"  ν = {r}  (|ν| = {float(fabs(r)):.10f})")

mu_roots = [r/64 for r in roots]
print("\nPoincaré roots μ = ν/64:")
for r in mu_roots:
    print(f"  μ = {r}  (|μ| = {float(fabs(r)):.10f})")

# Verify growth rate of q matches μ₀
mu0 = max(mu_roots, key=lambda x: fabs(x))
print(f"\nDominant root μ₀ = {float(mu0.real):.10f}")
print(f"|μ₀| = {float(fabs(mu0)):.10f}")

# Check q growth
for n in [50, 100, 200]:
    if q[n] != 0 and q[n-1] != 0:
        r = q[n] / q[n-1]
        print(f"  q[{n}]/q[{n-1}] = {float(r):.10f}  (should → {float(mu0.real):.10f})")
