#!/usr/bin/env sage
"""
Core computation: compare the P2.7 error to the Zudilin error.

The key hypothesis:
  e_n^{P27} / (64^{-n} · e_n^Z) → polynomial gauge in n

If true, then c₀(e)=0 for P2.7 follows from c₀(e)=0 for Zudilin.

Step 1: Compute rec_a (the AESZ #209 inner recurrence) via guess()
Step 2: Compute p̃_n via rec_a with ICs p̃_0=0, p̃_1=20, p̃_2=7425/16
Step 3: Compute errors e_n^Z = p̃_n - L·a_n and e_n^{P27} = p_n - L·q_n
Step 4: Study the ratio e_n^{P27} / (64^{-n} · e_n^Z)
"""
from ore_algebra import OreAlgebra, guess
import mpmath
mpmath.mp.dps = 200

# ============================================================
# Step 1: Compute a_n and find rec_a
# ============================================================
N = 50
a_list = []
for i in range(N):
    val = QQ(0)
    for k in range(i+1):
        val += binomial(i,k)^2 * binomial(i+k,i) * binomial(i+2*k,i)
    a_list.append(val)

print("a_n:", a_list[:6])

Rn = PolynomialRing(QQ, 'nn')
nn = Rn.gen()
OS = OreAlgebra(Rn, 'Snn')
Snn = OS.gen()

print("\nGuessing rec_a...")
try:
    rec_a_op = guess(a_list[:40], OS, order=3)
    print("Found rec_a!")
    # Extract coefficients: rec_a_op = P0 + P1*Snn + P2*Snn^2 + P3*Snn^3
    P = [rec_a_op[j] for j in range(4)]
    for j in range(4):
        print("  P_%d(nn) = %s  (degree %d)" % (j, P[j], P[j].degree()))
except Exception as e:
    print("guess failed: %s" % e)
    print("Falling back to manual rec_a computation...")
    # Use the AESZ #209 inner recurrence computed elsewhere
    # The recurrence for a_n is the SAME as the Zudilin/rec_a recurrence
    # We know c_3 and c_0 from proof.tex; let me compute the full thing
    # by building the annihilating matrix directly

    # Build kernel matrix: P3(n) a_{n+3} + P2(n) a_{n+2} + P1(n) a_{n+1} + P0(n) a_n = 0
    # We seek P_j as polynomials in n of some degree d
    # For AESZ inner sum, the recurrence has order 3 and degree 6 per coefficient
    max_deg = 8
    nrows = 30
    rows = []
    for n_val in range(nrows):
        row = []
        for j in range(4):  # shifts 0,1,2,3
            for d in range(max_deg+1):  # polynomial degree
                row.append(QQ(n_val^d) * a_list[n_val + j])
        rows.append(row)

    M = matrix(QQ, rows)
    K = M.right_kernel()
    print("Kernel dimension: %d" % K.dimension())

    if K.dimension() >= 1:
        v = K.basis()[0]
        P = []
        for j in range(4):
            poly = sum(v[(max_deg+1)*j + d] * nn^d for d in range(max_deg+1))
            P.append(poly)
            print("  P_%d(nn) = %s" % (j, poly))
    else:
        raise RuntimeError("Cannot find recurrence")

# Verify rec_a
print("\n=== Verify rec_a ===")
for n_val in range(5):
    val = sum(P[j](nn=n_val) * a_list[n_val + j] for j in range(4))
    print("  n=%d: residual = %s" % (n_val, val))

# ============================================================
# Step 2: Compute p̃_n (Zudilin companion)
# ============================================================
# From Q5048: p̃_n = b̃_n + b̃̃_n satisfies rec_a
# Initial conditions: p̃_0 = 0, p̃_1 = 20, p̃_2 = 7425/16
print("\n=== Computing Zudilin companion p̃_n ===")
p_zud = [QQ(0)] * N
p_zud[0] = QQ(0)
p_zud[1] = QQ(20)
p_zud[2] = QQ(7425) / QQ(16)

# rec_a: P0(n)·u_n + P1(n)·u_{n+1} + P2(n)·u_{n+2} + P3(n)·u_{n+3} = 0
# => u_{n+3} = -(P0(n)·u_n + P1(n)·u_{n+1} + P2(n)·u_{n+2}) / P3(n)
for n_val in range(N-3):
    val = sum(P[j](nn=n_val) * p_zud[n_val + j] for j in range(3))
    p_zud[n_val + 3] = -val / P[3](nn=n_val)

# Verify: p̃_n / a_n should converge to ζ(2) + ζ(3)
L_mp = mpmath.zeta(2) + mpmath.zeta(3)
L_float = float(L_mp)

print("\np̃_n / a_n → ζ(2) + ζ(3):")
for n_val in range(min(20, N)):
    if a_list[n_val] != 0:
        r = float(p_zud[n_val] / a_list[n_val])
        print("  n=%2d: p̃/a = %.15f, err = %.3e" % (n_val, r, abs(r - L_float)))

# ============================================================
# Step 3: Compute P2.7 sequences
# ============================================================
def A_p27(n_val):
    return QQ(1024) * (2*n_val+5)^4 * (2*n_val+7)^3 * (2*n_val+9)^3 * (946*n_val^2+6407*n_val+10860)

def B_p27(n_val):
    return QQ(128) * (2*n_val+7)^3 * (2*n_val+9)^3 * (104060*n_val^6 + 1745370*n_val^5 +
        12145238*n_val^4 + 44886481*n_val^3 + 92943995*n_val^2 + 102256019*n_val + 46709052)

def C_p27(n_val):
    return QQ(16) * (n_val+3)^4 * (2*n_val+9)^3 * (3784*n_val^5 + 57792*n_val^4 +
        351019*n_val^3 + 1059230*n_val^2 + 1587211*n_val + 944620)

def D_p27(n_val):
    return QQ(1) * (n_val+3)^4 * (n_val+4)^6 * (946*n_val^2 + 4515*n_val + 5399)

# P2.7: A_n u_{n+1} = B_n u_n - C_{n-1} u_{n-1} + D_{n-2} u_{n-2}
q_p27 = [QQ(0)] * N
q_p27[0] = QQ(-215040420000)
q_p27[1] = QQ(-167282265043404) / QQ(905)
q_p27[2] = QQ(-964185327658080) / QQ(6071)

for n_val in range(2, N-1):
    q_p27[n_val+1] = (B_p27(n_val) * q_p27[n_val] - C_p27(n_val-1) * q_p27[n_val-1] + D_p27(n_val-2) * q_p27[n_val-2]) / A_p27(n_val)

p_p27 = [QQ(0)] * N
p_p27[0] = QQ(-612218384750)
p_p27[1] = QQ(-9525021973931919) / QQ(18100)
p_p27[2] = QQ(-29561828382772029) / QQ(65380)

for n_val in range(2, N-1):
    p_p27[n_val+1] = (B_p27(n_val) * p_p27[n_val] - C_p27(n_val-1) * p_p27[n_val-1] + D_p27(n_val-2) * p_p27[n_val-2]) / A_p27(n_val)

# ============================================================
# Step 4: Compute errors using HIGH PRECISION
# ============================================================
print("\n=== Computing errors with mpmath ===")

# Convert to mpmath for high-precision error computation
L_val = mpmath.zeta(2) + mpmath.zeta(3)
print("L = ζ(2) + ζ(3) = %s" % mpmath.nstr(L_val, 30))

# Zudilin error: e_n^Z = p̃_n - L · a_n
e_zud = []
for n_val in range(min(40, N)):
    e = mpmath.mpf(p_zud[n_val]) - L_val * mpmath.mpf(a_list[n_val])
    e_zud.append(e)

print("\ne_n^Z = p̃_n - L·a_n:")
for n_val in range(min(15, len(e_zud))):
    print("  n=%2d: e^Z = %s" % (n_val, mpmath.nstr(e_zud[n_val], 20)))

# P2.7 error: e_n^{P27} = p_n - L · q_n
e_p27 = []
for n_val in range(min(40, N)):
    e = mpmath.mpf(p_p27[n_val]) - L_val * mpmath.mpf(q_p27[n_val])
    e_p27.append(e)

print("\ne_n^{P27} = p_n - L·q_n:")
for n_val in range(min(15, len(e_p27))):
    print("  n=%2d: e^P27 = %s" % (n_val, mpmath.nstr(e_p27[n_val], 20)))

# ============================================================
# Step 5: Error ratios
# ============================================================
print("\n=== Error ratio: e^{P27}_n / e^Z_n ===")
for n_val in range(min(25, len(e_zud))):
    if abs(e_zud[n_val]) > 1e-100 and abs(e_p27[n_val]) > 1e-100:
        ratio = e_p27[n_val] / e_zud[n_val]
        print("  n=%2d: e^P27/e^Z = %s" % (n_val, mpmath.nstr(ratio, 20)))

print("\n=== Scaled error ratio: e^{P27}_n / (64^{-n} · e^Z_n) ===")
for n_val in range(min(25, len(e_zud))):
    if abs(e_zud[n_val]) > 1e-100 and abs(e_p27[n_val]) > 1e-100:
        ratio = e_p27[n_val] / (mpmath.power(64, -n_val) * e_zud[n_val])
        print("  n=%2d: ratio = %s" % (n_val, mpmath.nstr(ratio, 20)))

# ============================================================
# Step 6: Error decay rates
# ============================================================
print("\n=== Error decay rates ===")
print("e^Z_{n+1} / e^Z_n (should → μ₁ ≈ 0.067):")
for n_val in range(min(20, len(e_zud)-1)):
    if abs(e_zud[n_val]) > 1e-100:
        ratio = e_zud[n_val+1] / e_zud[n_val]
        print("  n=%2d: %s" % (n_val, mpmath.nstr(ratio, 12)))

print("\ne^{P27}_{n+1} / e^{P27}_n (should → μ₁/64 ≈ 0.00105):")
for n_val in range(min(20, len(e_p27)-1)):
    if abs(e_p27[n_val]) > 1e-100:
        ratio = e_p27[n_val+1] / e_p27[n_val]
        print("  n=%2d: %s" % (n_val, mpmath.nstr(ratio, 12)))

# ============================================================
# Step 7: Gauge factor q_n^{P27} / (64^{-n} · a_n)
# ============================================================
print("\n=== Gauge: q^{P27}_n / (64^{-n} · a_n) = q^{P27}_n · 64^n / a_n ===")
gauge_vals = []
for n_val in range(min(25, N)):
    if a_list[n_val] != 0:
        g = mpmath.mpf(q_p27[n_val]) * mpmath.power(64, n_val) / mpmath.mpf(a_list[n_val])
        gauge_vals.append(g)
        print("  n=%2d: gauge = %s" % (n_val, mpmath.nstr(g, 15)))

print("\nGauge ratio g_{n+1}/g_n (should stabilize if g_n ~ n^α):")
for n_val in range(min(20, len(gauge_vals)-1)):
    if abs(gauge_vals[n_val]) > 1e-100:
        ratio = gauge_vals[n_val+1] / gauge_vals[n_val]
        print("  n=%2d: %s" % (n_val, mpmath.nstr(ratio, 12)))

# Log of |gauge| to find power law
import math
print("\nlog|gauge|/log(n) for n >= 2:")
for n_val in range(2, min(25, len(gauge_vals))):
    if abs(gauge_vals[n_val]) > 0:
        lg = float(mpmath.log(abs(gauge_vals[n_val])))
        ln = math.log(n_val)
        print("  n=%2d: log|g| = %.6f, log|g|/log(n) = %.6f" % (n_val, lg, lg/ln))
