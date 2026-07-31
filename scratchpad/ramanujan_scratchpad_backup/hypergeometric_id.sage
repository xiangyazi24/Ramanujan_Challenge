#!/usr/bin/env sage
"""
1. Verify the two-step error ratio e_n/e_{n-2} → -|r±|²
2. Compute the Pochhammer-gauged sequence v_n = q_n / h_n
3. Try to identify v_n as a hypergeometric series
"""
import mpmath
mpmath.mp.dps = 500

N = 60

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

p = [QQ(0)] * N
p[0] = QQ(-612218384750)
p[1] = QQ(-9525021973931919) / QQ(18100)
p[2] = QQ(-29561828382772029) / QQ(65380)

for n in range(2, N-1):
    q[n+1] = B_p27(n)/A_p27(n) * q[n] - C_p27(n-1)/A_p27(n-1) * q[n-1] + D_p27(n-2)/A_p27(n-2) * q[n-2]
    p[n+1] = B_p27(n)/A_p27(n) * p[n] - C_p27(n-1)/A_p27(n-1) * p[n-1] + D_p27(n-2)/A_p27(n-2) * p[n-2]

# Errors
L = mpmath.zeta(2) + mpmath.zeta(3)
e = []
for n_val in range(N):
    p_mp = mpmath.mpf(p[n_val].numerator()) / mpmath.mpf(p[n_val].denominator())
    q_mp = mpmath.mpf(q[n_val].numerator()) / mpmath.mpf(q[n_val].denominator())
    e.append(p_mp - L * q_mp)

# ============================================================
# 1. Two-step ratio e_n / e_{n-2}
# ============================================================
print("=== Two-step ratio e_n / e_{n-2} ===")
print("Expected: -|r±|² = -1/(4μ₀·64²) ≈ -1.110e-6")

# Compute -1/(4 μ₀ · 64²) exactly
# μ₀ is root of 4μ³ - 220μ² + 8μ - 1 = 0
# μ₀ ≈ 54.9637
# |r±|² = 1/(4μ₀ · 64²) = 1/(4 · 54.9637 · 4096) = 1/900698.06 ≈ 1.1103e-6
print("1/(4·54.9637·4096) = %e" % (1/(4*54.9637*4096)))

for n_val in range(2, min(40, len(e))):
    if abs(e[n_val-2]) > 0:
        r = e[n_val] / e[n_val-2]
        print("n=%2d: e_n/e_{n-2} = %s" % (n_val, mpmath.nstr(r, 15)))

# ============================================================
# 2. Pochhammer gauge
# ============================================================
print("\n=== Pochhammer gauge h_n ===")
def pochhammer(a, n):
    """Rising factorial (a)_n = a(a+1)...(a+n-1)"""
    return prod(a + j for j in range(n))

h = [QQ(0)] * N
for i in range(N):
    h[i] = QQ(2)^(-20*i) * pochhammer(QQ(3), i)^4 * pochhammer(QQ(4), i)^6 / \
           (pochhammer(QQ(5)/2, i)^4 * pochhammer(QQ(7)/2, i)^3 * pochhammer(QQ(9)/2, i)^3)

print("h_0 = %s" % h[0])
print("h_1 = %s" % h[1])
print("h_2 = %s" % h[2])

# Check D_n/A_n vs h_{n+1}/h_n
print("\nD_n/A_n vs h_{n+1}/h_n:")
for n_val in range(5):
    da = D_p27(n_val) / A_p27(n_val)
    hh = h[n_val+1] / h[n_val]
    print("  n=%d: D/A = %s, h_{n+1}/h_n = %s, match: %s" % (n_val, float(da), float(hh), da == hh))

# ============================================================
# 3. Gauged sequence v_n = q_n / h_n
# ============================================================
print("\n=== v_n = q_n / h_n ===")
v = [QQ(0)] * N
for i in range(N):
    if h[i] != 0:
        v[i] = q[i] / h[i]

for n_val in range(min(10, N)):
    print("v_%d = %s = %e" % (n_val, v[n_val], float(v[n_val])))

# Check ratios v_{n+1}/v_n
print("\nv_{n+1}/v_n:")
for n_val in range(min(15, N-1)):
    if v[n_val] != 0:
        r = float(v[n_val+1] / v[n_val])
        print("  n=%d: %.15f" % (n_val, r))

# ============================================================
# 4. Try to identify v_n as a hypergeometric sum
# ============================================================
# The Zudilin inner sum is a_n = Σ C(n,k)² C(n+k,k) C(n+2k,n)
# = ₄F₃(−n,−n,n+1,n+1; 1,1,1; 1) or similar.
# After Pochhammer gauge, v_n might be a similar ₄F₃.

# First, compute v_n directly from initial conditions and see structure
print("\n=== v_n factored ===")
for n_val in range(min(8, N)):
    if v[n_val] != 0:
        num = v[n_val].numerator()
        den = v[n_val].denominator()
        print("v_%d num = %s" % (n_val, factor(ZZ(num)) if abs(num) < 10^30 else "too large"))
        print("v_%d den = %s" % (n_val, factor(ZZ(den)) if abs(den) < 10^30 else "too large"))

# ============================================================
# 5. Compute D_n/A_n to verify Pochhammer gauge
# ============================================================
print("\n=== D_{n-2}/A_{n-2} pattern ===")
for n_val in range(5):
    da = D_p27(n_val) / A_p27(n_val)
    # Factor
    num_d = D_p27(n_val)
    den_a = A_p27(n_val)
    print("D_%d/A_%d = %s" % (n_val, n_val, da))

# ============================================================
# 6. Recurrence for v_n
# ============================================================
print("\n=== Guessing recurrence for v_n ===")
from ore_algebra import OreAlgebra, guess
Rn = PolynomialRing(QQ, 'nn')
nn_var = Rn.gen()
OS = OreAlgebra(Rn, 'Snn')
Snn = OS.gen()

v_list = v[:40]
for order in [3, 4, 5]:
    for deg in [12, 14, 16, 18, 20]:
        try:
            rec_v = guess(v_list, OS, order=order, degree=deg)
            print("Found for v: order=%d, degree=%d" % (order, deg))
            for j in range(order+1):
                print("  P_%d: degree %d" % (j, rec_v[j].degree()))
            print("  P_%d factored: %s" % (order, factor(rec_v[order])))
            print("  P_0 factored: %s" % factor(rec_v[0]))
            break
        except:
            pass
    else:
        continue
    break

# ============================================================
# 7. Compute the "well-poised" candidate
# ============================================================
# Try: w_n = Σ_{k=0}^n f(n,k) where f has Pochhammer structure
# related to the P2.7 parameters
print("\n=== Well-poised candidate ===")
# From Dauguet-Zudilin, the companion for ζ(2)+ζ(3) involves
# a ₅F₄ series with parameters related to the CY operator.
#
# The inner sum for AESZ #209 is:
# a_n = Σ C(n,k)² C(n+k,k) C(n+2k,k)
#     = ₄F₃(-n, -n, n+1, 1; 1, 1, 1; -1)  [not quite]
#
# The connection to ζ(2)+ζ(3) might come from a ε-derivative:
# d/dε [a_n(ε)]|_{ε=0} where a_n(ε) = Σ C(n,k)² C(n+k,k+ε) C(n+2k,n+ε)

# For now, let me check if q_n/h_n has a nice sum representation
# by examining the first few values.
print("v_0 = q_0/h_0 = q_0 = %s" % v[0])
print("v_1 = %s" % v[1])
print("v_2 = %s" % v[2])

# Check if v_n / (-215040420000) gives a clean sequence
c0 = v[0]
print("\nv_n / v_0:")
for n_val in range(min(8, N)):
    if c0 != 0:
        r = v[n_val] / c0
        print("  n=%d: %s" % (n_val, r))
