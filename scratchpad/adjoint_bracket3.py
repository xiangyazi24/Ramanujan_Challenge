#!/usr/bin/env python3
"""
Adjoint bracket computation for P2.7, using Q5103's exact formulas.

Q5103 eq 2.3 (adjoint): -D(n+2) w_{n+3} + C(n+2) w_{n+2} - B(n+2) w_{n+1} + A(n+1) w_n = 0
Q5103 eq 2.4 (bracket): J_n(w,u) = A(n+1) w_n u_{n+2} + D(n) w_{n+1} u_n
                                    + (D(n+1) w_{n+2} - C(n+1) w_{n+1}) u_{n+1}
Q5103 eq E6 (at n=0):   J_0 = D(0) w_1 u_0 + (D(1) w_2 - C(1) w_1) u_1 + A(1) w_0 u_2

Goal: verify J(w^0, p; 0) / J(w^0, q; 0) = ζ(2) + ζ(3) to 200+ digits.
"""
from mpmath import mp, mpf, pi, zeta, fabs, log, floor

mp.dps = 350

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

# ===== Forward solutions =====
N = 200

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
print(f"zeta(2)+zeta(3) = {L_val}")
print(f"q[0] = {q[0]}")

# ===== Adjoint backward Miller (Q5103 eq 2.3) =====
# -D(n+2) w_{n+3} + C(n+2) w_{n+2} - B(n+2) w_{n+1} + A(n+1) w_n = 0
# w_n = [D(n+2) w_{n+3} - C(n+2) w_{n+2} + B(n+2) w_{n+1}] / A(n+1)

M = N
def miller(seed):
    w = [mpf(0)] * (M + 10)
    w[M] = mpf(seed[0])
    w[M+1] = mpf(seed[1])
    w[M+2] = mpf(seed[2])
    for n in range(M-1, -1, -1):
        w[n] = (D(n+2) * w[n+3] - C(n+2) * w[n+2] + B(n+2) * w[n+1]) / A(n+1)
    return w

# Bracket (Q5103 eq 2.4)
def J(w, u, n):
    return (A(n+1) * w[n] * u[n+2]
            + D(n) * w[n+1] * u[n]
            + (D(n+1) * w[n+2] - C(n+1) * w[n+1]) * u[n+1])

# J at n=0 (E6) — should be identical to J(w,u,0)
def J0(w, u):
    return (D(0) * w[1] * u[0]
            + (D(1) * w[2] - C(1) * w[1]) * u[1]
            + A(1) * w[0] * u[2])

# ===== Verify adjoint recurrence =====
print("\nComputing Miller backward...")
wA = miller([1, 0, 1])
wB = miller([0, 1, 1])

print("Verifying adjoint recurrence residuals (seed A):")
max_res = mpf(0)
for n in range(0, min(M-3, 50)):
    res = -D(n+2)*wA[n+3] + C(n+2)*wA[n+2] - B(n+2)*wA[n+1] + A(n+1)*wA[n]
    rr = fabs(res)
    if rr > max_res:
        max_res = rr
    if n < 3:
        print(f"  n={n}: residual = {float(res):.6e}")
print(f"  Max residual (n=0..49): {float(max_res):.6e}")

# ===== Normalize by J(w, q; 0) = 1 =====
scaleA = J0(wA, q)
scaleB = J0(wB, q)
print(f"\nJ(wA,q;0) before normalization: {float(scaleA):.10e}")
print(f"J(wB,q;0) before normalization: {float(scaleB):.10e}")

wA_norm = [x / scaleA for x in wA]
wB_norm = [x / scaleB for x in wB]

seed_err = max(fabs(wA_norm[i] - wB_norm[i]) for i in range(5))
print(f"Miller seed disagreement (normalized): {float(seed_err):.6e}")

w = wA_norm

# ===== KEY TEST =====
print("\n" + "="*60)
print("KEY TEST: J(w,p;0) / J(w,q;0) = zeta(2)+zeta(3) ?")
print("="*60)

Jq0 = J0(w, q)
Jp0 = J0(w, p)
ratio = Jp0 / Jq0

err = fabs(ratio - L_val)
if err != 0:
    matched = int(floor(-log(err, 10)))
else:
    matched = 'inf'

print(f"J(w,q;0) = {Jq0}")
print(f"J(w,p;0) = {Jp0}")
print(f"ratio     = {ratio}")
print(f"zeta(2)+zeta(3) = {L_val}")
print(f"absolute error  = {float(err):.6e}")
print(f"matching decimal digits = {matched}")

Je0 = Jp0 - L_val * Jq0
print(f"\nJ(w,e;0) = {float(Je0):.6e}")
if Jq0 != 0:
    print(f"|J(w,e;0)/J(w,q;0)| = {float(fabs(Je0/Jq0)):.6e}")

# ===== Bracket constancy =====
print("\nBracket constancy check:")
Jq_ref = J(w, q, 0)
Jp_ref = J(w, p, 0)
print(f"  J(w,q;0) via generic formula = {float(Jq_ref):.15e}")
print(f"  J(w,q;0) via E6 formula      = {float(Jq0):.15e}")
print(f"  Difference: {float(fabs(Jq_ref - Jq0)):.6e}")

for n in [1, 2, 5, 10, 20, 50, 100]:
    Jqn = J(w, q, n)
    Jpn = J(w, p, n)
    dq = fabs(Jqn - Jq0)
    dp = fabs(Jpn - Jp0)
    rn = Jpn / Jqn if Jqn != 0 else mpf(0)
    print(f"  n={n:3d}: |dJ_q| = {float(dq):.3e}, |dJ_p| = {float(dp):.3e}, ratio = {float(rn):.15f}")

# ===== Adjoint growth =====
print("\nAdjoint growth |w[n]|^(1/n):")
for n in [10, 20, 50, 100, 150]:
    if w[n] != 0 and n > 0:
        g = fabs(w[n]) ** (mpf(1)/n)
        print(f"  n={n}: {float(g):.10f}")

# ===== Ratio at multiple n values =====
print("\nRatio J(w,p;n)/J(w,q;n) at various n:")
for n in [0, 1, 2, 5, 10, 20, 50]:
    Jpn = J(w, p, n)
    Jqn = J(w, q, n)
    if Jqn != 0:
        r = Jpn / Jqn
        d = r - L_val
        print(f"  n={n:3d}: diff from L = {float(d):.6e}")
