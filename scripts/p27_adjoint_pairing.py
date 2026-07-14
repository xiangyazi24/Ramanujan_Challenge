#!/usr/bin/env python3
"""P2.7: Adjoint certificate for dominant coefficient vanishing.

The Challenge recurrence is MONIC:
  u(n+1) = [B(n)/A(n)] u(n) - [C(n-1)/A(n-1)] u(n-1) + [D(n-2)/A(n-2)] u(n-2)

Adjoint (derived via summation by parts on the monic operator):
  W(n) = [B(n+2) W(n+1) - C(n+2) W(n+2) + D(n+2) W(n+3)] / A(n+2)

Miller backward on L* selects W^(0), the unique recessive adjoint solution
(dual to the dominant forward solution). The Lagrange bracket
  J(m) = [D(m)/A(m)] W(m+1) u(m)
       + [W(m-1) - B(m+1)/A(m+1) W(m)] u(m+1)
       + W(m) u(m+2)
is constant when Lu = 0 and L*W = 0.

If J(W^(0), e) = 0 where e = p - L*q, then c_0(e) = 0, proving e recessive.
"""
from mpmath import mp, mpf, zeta, nstr, fabs, log
import time

mp.dps = 600
L_val = zeta(2) + zeta(3)

def A(n):
    n = mpf(n)
    return 1024*(2*n+5)**4*(2*n+7)**3*(2*n+9)**3*(946*n**2+6407*n+10860)

def B(n):
    n = mpf(n)
    return 128*(2*n+7)**3*(2*n+9)**3*(104060*n**6+1745370*n**5+12145238*n**4+44886481*n**3+92943995*n**2+102256019*n+46709052)

def C(n):
    n = mpf(n)
    return 16*(n+3)**4*(2*n+9)**3*(3784*n**5+57792*n**4+351019*n**3+1059230*n**2+1587211*n+944620)

def D(n):
    n = mpf(n)
    return (n+3)**4*(n+4)**6*(946*n**2+4515*n+5399)

NMAX = 300
SIZE = NMAX + 3

# === Forward recurrence (CORRECT monic form from Challenge) ===
t0 = time.time()
p = [mpf(0)] * SIZE
q = [mpf(0)] * SIZE
p[0] = mpf(-612218384750)
p[1] = mpf(-9525021973931919) / mpf(18100)
p[2] = mpf(-29561828382772029) / mpf(65380)
q[0] = mpf(-215040420000)
q[1] = mpf(-167282265043404) / mpf(905)
q[2] = mpf(-964185327658080) / mpf(6071)

for n in range(2, NMAX + 2):
    An = A(n); Bn = B(n)
    Cn1 = C(n-1); An1 = A(n-1)
    Dn2 = D(n-2); An2 = A(n-2)
    p[n+1] = (Bn/An)*p[n] - (Cn1/An1)*p[n-1] + (Dn2/An2)*p[n-2]
    q[n+1] = (Bn/An)*q[n] - (Cn1/An1)*q[n-1] + (Dn2/An2)*q[n-2]

e = [p[i] - L_val*q[i] for i in range(SIZE)]
t1 = time.time()

print(f"=== Forward recurrence ({t1-t0:.1f}s) ===")
for nn in [50, 100, 150, 200, 250, 300]:
    if nn < SIZE and q[nn] != 0:
        err = fabs(p[nn]/q[nn] - L_val)
        digs = -int(float(log(err, 10))) if err > 0 else 999
        print(f"  p_{nn}/q_{nn} matches L to {digs} digits")

# === Miller backward for adjoint recessive W^(0) ===
t2 = time.time()
W = [mpf(0)] * SIZE
W[NMAX+2] = mpf(1)
W[NMAX+1] = mpf(0)
W[NMAX] = mpf(1)

for n in range(NMAX - 1, -1, -1):
    W[n] = (B(n+2)*W[n+1] - C(n+2)*W[n+2] + D(n+2)*W[n+3]) / A(n+2)

scale = W[0]
if scale != 0:
    for i in range(SIZE):
        W[i] /= scale
t3 = time.time()

print(f"\n=== Adjoint backward ({t3-t2:.1f}s) ===")
print(f"W(0) = 1 (normalized)")
print(f"W(1) = {nstr(W[1], 30)}")
print(f"W(2) = {nstr(W[2], 30)}")

print("\nGrowth rate |W(m)|^(1/m) (expect ~1.164 = 1/rho_0):")
for m in [10, 30, 50, 100, 150, 200, 250]:
    if m < SIZE and W[m] != 0:
        rate = fabs(W[m])**(mpf(1)/m)
        print(f"  m={m:3d}: {nstr(rate, 15)}")

print("\nAdjoint recurrence residual (should be ~0):")
for m in [5, 10, 50, 100, 200]:
    if m+3 < SIZE:
        resid = A(m+2)*W[m] - B(m+2)*W[m+1] + C(m+2)*W[m+2] - D(m+2)*W[m+3]
        print(f"  m={m:3d}: |residual| = {nstr(fabs(resid), 6)}")

# === Lagrange bracket ===
def bracket(W, u, m):
    """J(m) = [D(m)/A(m)] W(m+1) u(m) + [W(m-1) - B(m+1)/A(m+1) W(m)] u(m+1) + W(m) u(m+2)"""
    t1 = (D(m)/A(m)) * W[m+1] * u[m]
    t2 = (W[m-1] - (B(m+1)/A(m+1)) * W[m]) * u[m+1]
    t3 = W[m] * u[m+2]
    return t1 + t2 + t3

test_points = [1, 2, 3, 5, 10, 20, 50, 100, 150, 200, 250]

print("\n=== [W, q]_m (sanity: should be constant) ===")
bq_vals = []
for m in test_points:
    if m+2 < SIZE and m >= 1:
        val = bracket(W, q, m)
        bq_vals.append((m, val))
        print(f"  m={m:3d}: {nstr(val, 40)}")

if len(bq_vals) >= 2:
    r = bq_vals[-1][1] / bq_vals[0][1]
    d = fabs(r - 1)
    digs = -int(float(log(d, 10))) if d > 0 else 999
    print(f"  Constancy: {digs} digits")

print("\n=== [W, p]_m (sanity: should be constant) ===")
bp_vals = []
for m in test_points:
    if m+2 < SIZE and m >= 1:
        val = bracket(W, p, m)
        bp_vals.append((m, val))
        print(f"  m={m:3d}: {nstr(val, 40)}")

if len(bp_vals) >= 2:
    r = bp_vals[-1][1] / bp_vals[0][1]
    d = fabs(r - 1)
    digs = -int(float(log(d, 10))) if d > 0 else 999
    print(f"  Constancy: {digs} digits")

print("\n=== [W, e]_m (KEY TEST: zero iff c_0 = 0) ===")
be_vals = []
for m in test_points:
    if m+2 < SIZE and m >= 1:
        val = bracket(W, e, m)
        be_vals.append((m, val))
        print(f"  m={m:3d}: {nstr(val, 6)}")

# Cross-check: [W,p]/[W,q] should equal L if [W,e] = 0
print("\n=== Cross-check: [W,p]/[W,q] vs zeta(2)+zeta(3) ===")
bq10 = bracket(W, q, 10)
bp10 = bracket(W, p, 10)
if bq10 != 0:
    ratio_pq = bp10 / bq10
    print(f"  [W,p]/[W,q] = {nstr(ratio_pq, 60)}")
    print(f"  zeta(2)+zeta(3) = {nstr(L_val, 60)}")
    diff = ratio_pq - L_val
    if diff != 0:
        d = fabs(diff / L_val)
        digs = -int(float(log(d, 10)))
        print(f"  Relative agreement: {digs} digits")
    else:
        print(f"  EXACT match!")

    be10 = bracket(W, e, 10)
    print(f"\n  [W,e]_10 / [W,q]_10 = {nstr(be10/bq10, 6)}")
    print(f"  (This is c_0(e) up to normalization)")

tf = time.time()
print(f"\nTotal time: {tf-t0:.1f}s")
