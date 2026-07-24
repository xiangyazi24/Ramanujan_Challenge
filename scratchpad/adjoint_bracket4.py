#!/usr/bin/env python3
"""
Exact translation of Q5103's Sage code to Python/mpmath.
Uses Q5103's MONIC convention (eq 2.5-2.6, 2.7-2.10).

Key: alpha(n) = B(n)/A(n), beta(n) = -C(n-1)/A(n-1), gamma(n) = D(n-2)/A(n-2)
"""
from mpmath import mp, mpf, pi, zeta, fabs, log, floor, inf

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

# Q5103 monic coefficients (eq 2.5-2.6)
def alpha(n):
    return B(n) / A(n)

def beta(n):
    return -C(n-1) / A(n-1)

def gamma(n):
    return D(n-2) / A(n-2)

# Forward sequences using Q5103's monic convention
NTEST = 200

def forward(init, N):
    u = [mpf(0)] * (N + 3)
    u[0] = mpf(init[0])
    u[1] = mpf(init[1])
    u[2] = mpf(init[2])
    for n in range(2, N + 2):
        u[n+1] = alpha(n) * u[n] + beta(n) * u[n-1] + gamma(n) * u[n-2]
    return u

q_init = [mpf(-215040420000), mpf(-167282265043404)/mpf(905), mpf(-964185327658080)/mpf(6071)]
p_init = [mpf(-612218384750), mpf(-9525021973931919)/mpf(18100), mpf(-29561828382772029)/mpf(65380)]

q = forward(q_init, NTEST)
p = forward(p_init, NTEST)

# Also compute with my "standard monic" (divide everything by A(n))
def forward_std(init, N):
    u = [mpf(0)] * (N + 3)
    u[0] = mpf(init[0])
    u[1] = mpf(init[1])
    u[2] = mpf(init[2])
    for n in range(2, N + 2):
        u[n+1] = (B(n)*u[n] - C(n-1)*u[n-1] + D(n-2)*u[n-2]) / A(n)
    return u

q_std = forward_std(q_init, NTEST)
p_std = forward_std(p_init, NTEST)

# Compare sequences
print("=== Forward sequence comparison ===")
print(f"  q[3] (Q5103 monic):    {q[3]}")
print(f"  q[3] (std monic):      {q_std[3]}")
print(f"  q[3] ratio Q5103/std:  {float(q[3] / q_std[3]):.15f}")
print(f"  q[10] (Q5103 monic):   {float(q[10]):.10e}")
print(f"  q[10] (std monic):     {float(q_std[10]):.10e}")
print(f"  q[10] ratio:           {float(q[10] / q_std[10]):.15f}")

L_val = zeta(2) + zeta(3)
print(f"\nzeta(2)+zeta(3) = {L_val}")

# Miller backward (Q5103 eq 2.7)
# w_n = alpha(n+2)*w_{n+1} + beta(n+3)*w_{n+2} + gamma(n+4)*w_{n+3}
M = 180

def miller(seed):
    w = [mpf(0)] * (M + 3)
    w[M] = mpf(seed[0])
    w[M+1] = mpf(seed[1])
    w[M+2] = mpf(seed[2])
    for n in range(M-1, -1, -1):
        w[n] = alpha(n+2)*w[n+1] + beta(n+3)*w[n+2] + gamma(n+4)*w[n+3]
    return w

# J0 (Q5103 eq 2.10)
def J0(w, u):
    ell0 = gamma(2) * w[1]
    ell1 = beta(2) * w[1] + gamma(3) * w[2]
    ell2 = w[0]
    return ell0 * u[0] + ell1 * u[1] + ell2 * u[2]

# J at m >= 1 (Q5103 eq 2.9)
def J(w, u, m):
    assert m >= 1
    return (gamma(m+2) * w[m+1] * u[m]
            + (w[m-1] - alpha(m+1) * w[m]) * u[m+1]
            + w[m] * u[m+2])

# Normalized Miller
def normalized_miller(seed):
    w = miller(seed)
    scale = J0(w, q)
    assert scale != 0
    return [x / scale for x in w]

print("\nComputing Miller backward (M={})...".format(M))
wA = normalized_miller([1, 0, 1])
wB = normalized_miller([0, 1, 1])

seed_err = max(fabs(wA[i] - wB[i]) for i in range(3))
print(f"Miller seed disagreement: {float(seed_err):.6e}")

w = wA

# KEY TEST
print("\n" + "=" * 60)
print("KEY TEST: J(w,p;0) / J(w,q;0) = zeta(2)+zeta(3) ?")
print("=" * 60)

Jq0 = J0(w, q)
Jp0 = J0(w, p)
ratio = Jp0 / Jq0
err = fabs(ratio - L_val)
matched = int(floor(-log(err, 10))) if err > 0 else 'inf'

ell0 = gamma(2) * w[1]
ell1 = beta(2) * w[1] + gamma(3) * w[2]
ell2 = w[0]

print(f"ell_0 = {ell0}")
print(f"ell_1 = {ell1}")
print(f"ell_2 = {ell2}")
print(f"J(w,q;0) = {Jq0}")
print(f"J(w,p;0) = {Jp0}")
print(f"ratio     = {ratio}")
print(f"zeta(2)+zeta(3) = {L_val}")
print(f"absolute error  = {float(err):.6e}")
print(f"matching decimal digits = {matched}")

assert fabs(Jq0 - 1) < mpf(10)**(-250), f"J(w,q;0) not 1: {Jq0}"
assert err < mpf(10)**(-200), f"Only {matched} matching digits"

# Constancy check
print("\nBracket constancy check:")
for m in [1, 2, 5, 10, 30, 60, 80, 100]:
    if m < M - 5:
        dq = fabs(J(w, q, m) - Jq0)
        dp = fabs(J(w, p, m) - Jp0)
        print(f"  m={m:3d}: |dJ_q| = {float(dq):.3e}, |dJ_p| = {float(dp):.3e}")

# Growth check
print("\nAdjoint growth:")
for m in [20, 40, 80, 120]:
    if m < M:
        print(f"  m={m}: |w_m|^(1/m) = {float(fabs(w[m])**(mpf(1)/m)):.10f}")
