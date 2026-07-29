#!/usr/bin/env sage
"""
Error ratio analysis with CORRECT P2.7 recurrence.

The recurrence is:
  u_{n+1} = (B_n/A_n) u_n - (C_{n-1}/A_{n-1}) u_{n-1} + (D_{n-2}/A_{n-2}) u_{n-2}

Key question: what is the relationship between e_n^{P27} and e_n^{Zudilin}?

Both errors are purely subdominant (c₀ = 0 for each).
Zudilin: |e^Z| ~ |μ₁|^n ≈ 0.067^n, σ₁ = -3/2
P2.7:    |e^P| ~ |μ₁/64|^n ≈ 0.00105^n, σ₁ ≈ -12

The 64-scaled ratio: e^P / (64^{-n} · e^Z) ~ n^{σ₁^P - σ₁^Z} = n^{-10.5}
"""
from ore_algebra import OreAlgebra, guess
import mpmath
mpmath.mp.dps = 500

N = 50

# ============================================================
# Compute a_n and rec_a
# ============================================================
a_list = []
for i in range(N):
    val = QQ(0)
    for k in range(i+1):
        val += binomial(i,k)^2 * binomial(i+k,i) * binomial(i+2*k,i)
    a_list.append(val)

Rn = PolynomialRing(QQ, 'nn')
nn_var = Rn.gen()
OS = OreAlgebra(Rn, 'Snn')
Snn = OS.gen()

rec_a_op = guess(a_list[:40], OS, order=3)
P_rec = [rec_a_op[j] for j in range(4)]

# ============================================================
# Zudilin companion p̃_n via rec_a
# ============================================================
p_zud = [QQ(0)] * N
p_zud[0] = QQ(0)
p_zud[1] = QQ(20)
p_zud[2] = QQ(7425) / QQ(16)

for n_val in range(N-3):
    val = sum(P_rec[j](nn=n_val) * p_zud[n_val + j] for j in range(3))
    p_zud[n_val + 3] = -val / P_rec[3](nn=n_val)

# ============================================================
# CORRECT P2.7 recurrence
# ============================================================
def A_p27(n):
    return QQ(1024) * (2*n+5)^4 * (2*n+7)^3 * (2*n+9)^3 * (946*n^2+6407*n+10860)

def B_p27(n):
    return QQ(128) * (2*n+7)^3 * (2*n+9)^3 * (104060*n^6 + 1745370*n^5 +
        12145238*n^4 + 44886481*n^3 + 92943995*n^2 + 102256019*n + 46709052)

def C_p27(n):
    return QQ(16) * (n+3)^4 * (2*n+9)^3 * (3784*n^5 + 57792*n^4 +
        351019*n^3 + 1059230*n^2 + 1587211*n + 944620)

def D_p27(n):
    return QQ(1) * (n+3)^4 * (n+4)^6 * (946*n^2 + 4515*n + 5399)

q_p27 = [QQ(0)] * N
q_p27[0] = QQ(-215040420000)
q_p27[1] = QQ(-167282265043404) / QQ(905)
q_p27[2] = QQ(-964185327658080) / QQ(6071)

p_p27 = [QQ(0)] * N
p_p27[0] = QQ(-612218384750)
p_p27[1] = QQ(-9525021973931919) / QQ(18100)
p_p27[2] = QQ(-29561828382772029) / QQ(65380)

for n in range(2, N-1):
    q_p27[n+1] = B_p27(n)/A_p27(n) * q_p27[n] - C_p27(n-1)/A_p27(n-1) * q_p27[n-1] + D_p27(n-2)/A_p27(n-2) * q_p27[n-2]
    p_p27[n+1] = B_p27(n)/A_p27(n) * p_p27[n] - C_p27(n-1)/A_p27(n-1) * p_p27[n-1] + D_p27(n-2)/A_p27(n-2) * p_p27[n-2]

# ============================================================
# Compute errors
# ============================================================
L = mpmath.zeta(2) + mpmath.zeta(3)

e_zud = []
e_p27 = []
for n_val in range(min(40, N)):
    ez = mpmath.mpf(p_zud[n_val]) - L * mpmath.mpf(a_list[n_val])
    e_zud.append(ez)

    ep = mpmath.mpf(p_p27[n_val].numerator()) / mpmath.mpf(p_p27[n_val].denominator()) - L * mpmath.mpf(q_p27[n_val].numerator()) / mpmath.mpf(q_p27[n_val].denominator())
    e_p27.append(ep)

# ============================================================
# Error decay rates
# ============================================================
print("=== Zudilin error decay |e^Z_{n+1}/e^Z_n| ===")
for n_val in range(min(25, len(e_zud)-1)):
    if abs(e_zud[n_val]) > 0:
        r = e_zud[n_val+1] / e_zud[n_val]
        print("n=%2d: ratio = %s, |ratio| = %s" % (n_val, mpmath.nstr(r, 12), mpmath.nstr(abs(r), 12)))

print("\n=== P2.7 error decay |e^P_{n+1}/e^P_n| ===")
for n_val in range(min(25, len(e_p27)-1)):
    if abs(e_p27[n_val]) > 0:
        r = e_p27[n_val+1] / e_p27[n_val]
        print("n=%2d: ratio = %s, |ratio| = %s" % (n_val, mpmath.nstr(r, 12), mpmath.nstr(abs(r), 12)))

# ============================================================
# Error ratios
# ============================================================
print("\n=== Raw ratio e^P_n / e^Z_n ===")
for n_val in range(min(25, len(e_zud))):
    if abs(e_zud[n_val]) > 0 and abs(e_p27[n_val]) > 0:
        r = e_p27[n_val] / e_zud[n_val]
        print("n=%2d: %s" % (n_val, mpmath.nstr(r, 15)))

print("\n=== 64-scaled ratio e^P_n / (64^{-n} · e^Z_n) ===")
for n_val in range(min(25, len(e_zud))):
    if abs(e_zud[n_val]) > 0 and abs(e_p27[n_val]) > 0:
        r = e_p27[n_val] / (mpmath.power(64, -n_val) * e_zud[n_val])
        print("n=%2d: %s" % (n_val, mpmath.nstr(r, 15)))

# ============================================================
# |e^P_n| · 64^n / |e^Z_n| — should behave like n^{-10.5}
# ============================================================
print("\n=== |e^P_n| · 64^n / |e^Z_n| ===")
ratios_scaled = []
for n_val in range(min(30, len(e_zud))):
    if abs(e_zud[n_val]) > 0 and abs(e_p27[n_val]) > 0:
        r = abs(e_p27[n_val]) * mpmath.power(64, n_val) / abs(e_zud[n_val])
        ratios_scaled.append((n_val, float(r)))
        print("n=%2d: %s" % (n_val, mpmath.nstr(r, 12)))

# Check power law: log ratio / log n
print("\nlog |ratio| / log n (should → ~-10.5 if ratio ~ n^{-10.5}):")
import math
for n_val, r in ratios_scaled:
    if n_val >= 2 and r > 0:
        print("n=%2d: log(r)/log(n) = %.4f" % (n_val, math.log(r)/math.log(n_val)))

# ============================================================
# Check successive log ratio — fits n^α with what α?
# ============================================================
print("\n=== Successive ratio (should → constant if power law) ===")
print("r_{n+1}/r_n · n/(n+1) for power law detection:")
for i in range(1, len(ratios_scaled)):
    n1, r1 = ratios_scaled[i-1]
    n2, r2 = ratios_scaled[i]
    if r1 > 0 and n1 > 0:
        actual_ratio = r2/r1
        print("n=%d→%d: r_{n+1}/r_n = %.6f" % (n1, n2, actual_ratio))

# ============================================================
# Also check the gauge q^P / (64^{-n} a_n) with correct recurrence
# ============================================================
print("\n=== Gauge: q^P_n · 64^n / a_n (CORRECT recurrence) ===")
for n_val in range(min(20, N)):
    if a_list[n_val] != 0:
        g = mpmath.mpf(q_p27[n_val]) * mpmath.power(64, n_val) / mpmath.mpf(a_list[n_val])
        print("n=%2d: %s" % (n_val, mpmath.nstr(g, 15)))

# ============================================================
# Try to identify the P2.7 recurrence in standard form
# ============================================================
# Multiply through to clear denominators:
# A_n A_{n-1} A_{n-2} u_{n+1} = B_n A_{n-1} A_{n-2} u_n
#   - C_{n-1} A_n A_{n-2} u_{n-1} + D_{n-2} A_n A_{n-1} u_{n-2}
#
# This is a 4-term recurrence with polynomial coefficients of degree ~ 36.
# Let me verify that guess() can find this recurrence.
print("\n=== Guessing recurrence from q_p27 ===")
q_list = [q_p27[i] for i in range(min(30, N))]
try:
    rec_q = guess(q_list, OS, order=3)
    print("Found order-3 recurrence!")
    for j in range(4):
        print("  P_%d: degree %d" % (j, rec_q[j].degree()))
    # Print leading and trailing coefficients
    for j in [0, 3]:
        poly = rec_q[j]
        print("  P_%d leading: %s * nn^%d" % (j, poly.leading_coefficient(), poly.degree()))
        print("  P_%d trailing: %s" % (j, poly(nn=0)))
except Exception as e:
    print("Guessing failed: %s" % e)
    # Try with more terms
    try:
        q_list2 = [q_p27[i] for i in range(min(40, N))]
        rec_q = guess(q_list2, OS, order=3)
        print("Found with 40 terms! Degrees:")
        for j in range(4):
            print("  P_%d: degree %d" % (j, rec_q[j].degree()))
    except Exception as e2:
        print("Still failed: %s" % e2)
