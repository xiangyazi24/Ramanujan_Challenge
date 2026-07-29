#!/usr/bin/env sage
"""
Search for nonlocal transforms from AESZ #209 (a_n) to P2.7 (q_n).
Try: binomial, Euler, Hadamard, weighted binomial.
"""
import mpmath
mpmath.mp.dps = 50

# ============================================================
# AESZ #209 sequence
# ============================================================
N = 30
a = [QQ(0)] * N
for n in range(N):
    val = QQ(0)
    for k in range(n+1):
        val += binomial(n,k)^2 * binomial(n+k,n) * binomial(n+2*k,n)
    a[n] = val

print("AESZ #209 first values:")
for i in range(10):
    print("  a_%d = %s" % (i, a[i]))

# ============================================================
# P2.7 sequences
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

q = [QQ(0)] * N
q[0] = QQ(-215040420000)
q[1] = QQ(-167282265043404) / QQ(905)
q[2] = QQ(-964185327658080) / QQ(6071)

for n in range(2, N-1):
    q[n+1] = B_p27(n)/A_p27(n)*q[n] - C_p27(n-1)/A_p27(n-1)*q[n-1] + D_p27(n-2)/A_p27(n-2)*q[n-2]

print("\nP2.7 first values:")
for i in range(8):
    print("  q_%d = %s" % (i, float(q[i])))

# ============================================================
# Pochhammer gauge
# ============================================================
def gamma_gauge(n):
    """γ_n = 2^{-20n} · (3)_n^4 · (4)_n^6 / [(5/2)_n^4 · (7/2)_n^3 · (9/2)_n^3]"""
    if n == 0:
        return QQ(1)
    num = QQ(2)^(-20*n)
    for j in range(n):
        num *= (3+j)^4 * (4+j)^6
    den = QQ(1)
    for j in range(n):
        den *= (QQ(5)/2+j)^4 * (QQ(7)/2+j)^3 * (QQ(9)/2+j)^3
    return num / den

print("\nGauged P2.7: v_n = q_n / γ_n")
v = [QQ(0)] * N
for i in range(min(15, N)):
    g = gamma_gauge(i)
    if g != 0:
        v[i] = q[i] / g
        print("  v_%d = %s = %e" % (i, v[i] if abs(v[i]) < 10^20 else "...", float(v[i])))

# Q₂₀₉ factor
def Q209(x):
    return 946*x^2 - 2623*x + 1830

print("\nQ209(n+83/22) values:")
for n in range(8):
    val = Q209(n + QQ(83)/22)
    print("  n=%d: Q209 = %s" % (n, val))

# Further gauge: w_n = v_n / Q209(n+83/22)
print("\nFully gauged: w_n = q_n / (γ_n · Q209(n+83/22))")
w = [QQ(0)] * N
for i in range(min(15, N)):
    g = gamma_gauge(i)
    Q_val = Q209(i + QQ(83)/22)
    if g != 0 and Q_val != 0:
        w[i] = q[i] / (g * Q_val)
        print("  w_%d = %e" % (i, float(w[i])))

# ============================================================
# Transform searches
# ============================================================

# 1. Binomial transform: b_n = Σ C(n,m) a_m
print("\n=== Binomial transform of a_n ===")
b_binom = []
for n in range(min(15, N)):
    val = sum(binomial(n, m) * a[m] for m in range(n+1))
    b_binom.append(val)
    print("  b_%d = %s" % (n, val))

# Compare with v_n
print("\nb_n / v_n:")
for n in range(min(10, len(b_binom))):
    if v[n] != 0 and b_binom[n] != 0:
        print("  n=%d: %e" % (n, float(b_binom[n] / v[n])))

# 2. Euler transform: b_n = Σ (-1)^{n-m} C(n,m) a_m
print("\n=== Euler transform of a_n ===")
b_euler = []
for n in range(min(15, N)):
    val = sum((-1)^(n-m) * binomial(n, m) * a[m] for m in range(n+1))
    b_euler.append(val)
    print("  b_%d = %s" % (n, val))

# 3. Double binomial: b_n = Σ C(n,m)^2 a_m
print("\n=== Double binomial transform ===")
b_double = []
for n in range(min(15, N)):
    val = sum(binomial(n, m)^2 * a[m] for m in range(n+1))
    b_double.append(val)
    print("  b_%d = %s" % (n, val))

# 4. Rising factorial weighted: b_n = Σ C(n,m) C(n+m,m) a_m
print("\n=== C(n,m)·C(n+m,m) transform ===")
b_rise = []
for n in range(min(15, N)):
    val = sum(binomial(n, m) * binomial(n+m, m) * a[m] for m in range(n+1))
    b_rise.append(val)
    print("  b_%d = %s" % (n, val))

# 5. Check if any transform matches v_n up to normalization
print("\n=== Ratio checks ===")
for name, b_list in [("binom", b_binom), ("euler", b_euler),
                      ("double", b_double), ("rise", b_rise)]:
    ratios = []
    for n in range(min(8, len(b_list))):
        if v[n] != 0 and b_list[n] != 0:
            ratios.append(float(b_list[n] / v[n]))
    print("%s: ratios = %s" % (name, ["%.6e" % r for r in ratios]))
    if len(ratios) >= 2:
        const = all(abs(r - ratios[0]) / max(abs(ratios[0]), 1e-30) < 0.01 for r in ratios[1:])
        print("  constant? %s" % const)

# 6. Try convolution: b_n = Σ a_m * a_{n-m}
print("\n=== Convolution a*a ===")
b_conv = []
for n in range(min(15, N)):
    val = sum(a[m] * a[n-m] for m in range(n+1))
    b_conv.append(val)
    print("  b_%d = %s" % (n, val))

# 7. Try Hadamard product: b_n = a_n * a_n
print("\n=== Hadamard square a_n^2 ===")
for n in range(8):
    print("  a_%d^2 = %s" % (n, a[n]^2))

# 8. Try to factor q_0
print("\n=== q_0 factored ===")
print("q_0 = %s" % q[0])
print("factored: %s" % factor(ZZ(q[0].numerator())))
print("\na_0 = %s, a_1 = %s, a_2 = %s" % (a[0], a[1], a[2]))
print("a_0·a_1·a_2·a_3 = %s" % (a[0]*a[1]*a[2]*a[3]))

# 9. Check: does q_0 divide into any simple expression?
q0 = ZZ(q[0])
# -215040420000 = -2^5 · 3^6 · 5^4 · 7^3 · 43
# Let's check: is q_0 related to factorials?
print("\nFactorial products:")
for p in range(1, 20):
    print("  %d! = %d" % (p, factorial(p)))

# 10. Try: does q_n / (some product of factorials) match a_n times something?
# The gauge γ_n involves (3)_n^4 (4)_n^6 / [(5/2)_n^4 (7/2)_n^3 (9/2)_n^3]
# At n=0, γ_0 = 1, so the normalization is all in q_0
# q_0 = -2^5 · 3^6 · 5^4 · 7^3 · 43
# The 43 is interesting — it appears in Q₂₀₉ factoring?
print("\n=== Q₂₀₉ at special values ===")
for x in range(20):
    val = 946*x^2 - 2623*x + 1830
    if val != 0 and val % 43 == 0:
        print("  Q209(%d) = %d = 43 * %d" % (x, val, val // 43))
print("  Q209(0) = %d" % Q209(0))
print("  Q209(83/22) = %s" % Q209(QQ(83)/22))

# 11. Explore: q_0 / Q209(83/22) · ...
q209_0 = Q209(QQ(83)/22)
print("\nQ209(83/22) = %s = %e" % (q209_0, float(q209_0)))
print("q_0 / Q209(83/22) = %s" % (q[0] / q209_0))

# 12. Try ore_algebra guess on gauged sequence
from ore_algebra import OreAlgebra, guess
Rn = PolynomialRing(QQ, 'nn')
OS = OreAlgebra(Rn, 'Snn')

print("\n=== Recurrence guess for gauged v_n ===")
v_list = [v[i] for i in range(min(20, N)) if v[i] != 0]
for order in [2, 3, 4]:
    try:
        rec = guess(v_list[:15], OS, order=order)
        print("Found order-%d recurrence for v_n!" % order)
        for j in range(order+1):
            print("  P_%d: degree %d" % (j, rec[j].degree()))
        break
    except:
        print("  order=%d: no relation" % order)
