#!/usr/bin/env sage
"""
Shifted Barnes integral test — fixed version.
"""
import mpmath
mpmath.mp.dps = 150

def cotpi_mp(t):
    return mpmath.cos(mpmath.pi * t) / mpmath.sin(mpmath.pi * t)

def R_nu(nu, t):
    return (mpmath.gamma(nu+1-t)**3 * mpmath.gamma(t) /
            (mpmath.gamma(1-t)**3 * mpmath.gamma(nu+1)**2 * mpmath.gamma(t+nu+1)))

def integrand(nu, t):
    r = R_nu(nu, t)
    s = mpmath.pi / mpmath.sin(mpmath.pi * t)
    c = mpmath.pi * cotpi_mp(t) - 1
    return r * s**2 * c

def J_nu_contour(nu):
    def f(y):
        t = mpmath.mpf('0.5') + 1j * y
        return integrand(nu, t) / (2 * mpmath.pi * 1j)
    result = mpmath.quad(f, [-50, 50], error=True, maxdegree=8)
    return result

L = mpmath.zeta(2) + mpmath.zeta(3)

# ============================================================
# a_n and its recurrence
# ============================================================
a_list = []
for i in range(30):
    val = QQ(0)
    for k in range(i+1):
        val += binomial(i,k)^2 * binomial(i+k,i) * binomial(i+2*k,i)
    a_list.append(val)

from ore_algebra import OreAlgebra, guess
Rn = PolynomialRing(QQ, 'nn')
OS = OreAlgebra(Rn, 'Snn')

# Try guessing without specifying order
rec_a_op = guess(a_list[:25], OS)
print("rec_a: order %d" % rec_a_op.order())
P_rec = [rec_a_op[j] for j in range(rec_a_op.order()+1)]
for j in range(len(P_rec)):
    print("  P_%d: degree %d, value: %s" % (j, P_rec[j].degree(), P_rec[j]))

# Zudilin companion
# For order-2 recurrence: p̃_0 = 0, p̃_1 = ??
# The Zudilin companion for ζ(2)+ζ(3) with a_n is given by:
# p̃_n = Σ_{k=0}^n C(n,k)²C(n+k,k)C(n+2k,k) · H_n(k)
# where H_n(k) involves harmonic numbers.
#
# Let me compute p̃_n directly from the sum formula:
# p̃_n = Σ_{k=1}^n C(n,k)²C(n+k,k)C(n+2k,k) · [3·H_{n+k} - 3·H_k + (some digamma terms)]
#
# Actually, from the literature (Zudilin, Krattenthaler-Rivoal):
# The companion is more complex. Let me just use the known e_n values.
#
# For the SANITY CHECK, I'll compute J_n at integer n and verify it gives
# e_n = p̃_n - L·a_n for the Zudilin case. But I don't have the explicit
# p̃_n formula handy. Let me skip the sanity check and go directly to
# the shifted case — I can verify by checking if the J values satisfy
# any recurrence.

# ============================================================
# P2.7 errors
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

N_p27 = 12
q = [QQ(0)] * N_p27
q[0] = QQ(-215040420000)
q[1] = QQ(-167282265043404) / QQ(905)
q[2] = QQ(-964185327658080) / QQ(6071)

p = [QQ(0)] * N_p27
p[0] = QQ(-612218384750)
p[1] = QQ(-9525021973931919) / QQ(18100)
p[2] = QQ(-29561828382772029) / QQ(65380)

for n in range(2, N_p27-1):
    q[n+1] = B_p27(n)/A_p27(n) * q[n] - C_p27(n-1)/A_p27(n-1) * q[n-1] + D_p27(n-2)/A_p27(n-2) * q[n-2]
    p[n+1] = B_p27(n)/A_p27(n) * p[n] - C_p27(n-1)/A_p27(n-1) * p[n-1] + D_p27(n-2)/A_p27(n-2) * p[n-2]

print("\nP2.7 errors:")
e_p27 = []
for i in range(min(8, N_p27)):
    ep = mpmath.mpf(p[i].numerator()) / mpmath.mpf(p[i].denominator()) - \
         L * mpmath.mpf(q[i].numerator()) / mpmath.mpf(q[i].denominator())
    e_p27.append(ep)
    print("  e^P_%d = %s" % (i, mpmath.nstr(ep, 25)))

# ============================================================
# First: compute J_n for integer n = 0, 1, 2, 3 (Zudilin case)
# ============================================================
print("\n=== Zudilin J_n (integer ν) ===")
J_zud = []
for n_int in range(5):
    nu = mpmath.mpf(n_int)
    print("Computing J_%d..." % n_int)
    try:
        result = J_nu_contour(nu)
        j_val = result[0]
        J_zud.append(j_val)
        print("  J_%d = %s (err: %s)" % (n_int, mpmath.nstr(j_val, 25), mpmath.nstr(result[1], 5)))
    except Exception as ex:
        print("  FAILED: %s" % ex)
        J_zud.append(None)

# Check J_n / a_n ratio to extract the Zudilin error
print("\nJ_n / a_n (should give e^Z_n / a_n pattern):")
for i in range(min(len(J_zud), 5)):
    if J_zud[i] is not None and a_list[i] != 0:
        ratio = J_zud[i] / mpmath.mpf(a_list[i])
        print("  n=%d: J/a = %s" % (i, mpmath.nstr(ratio, 20)))

# ============================================================
# Shifted Barnes: J_{n+17/22}
# ============================================================
print("\n=== Shifted Barnes J_{n+17/22} ===")
J_shifted = []
for n_val in range(6):
    nu = n_val + mpmath.mpf(17)/22
    print("Computing J_{%d+17/22} (ν = %.6f)..." % (n_val, float(nu)))
    try:
        result = J_nu_contour(nu)
        j_val = result[0]
        J_shifted.append(j_val)
        print("  J = %s (err: %s)" % (mpmath.nstr(j_val, 25), mpmath.nstr(result[1], 5)))
    except Exception as ex:
        print("  FAILED: %s" % ex)
        J_shifted.append(None)

# Compare with P2.7 errors
print("\n=== Ratio J_{n+17/22} / e^P_n ===")
for i in range(min(len(J_shifted), len(e_p27))):
    if J_shifted[i] is not None and abs(e_p27[i]) > 1e-100:
        ratio = J_shifted[i] / e_p27[i]
        print("n=%d: %s" % (i, mpmath.nstr(ratio, 20)))

print("\n=== J_{n+17/22} / (64^{-n} · e^P_n) ===")
for i in range(min(len(J_shifted), len(e_p27))):
    if J_shifted[i] is not None and abs(e_p27[i]) > 1e-100:
        ratio = J_shifted[i] / (mpmath.power(64, -i) * e_p27[i])
        print("n=%d: %s" % (i, mpmath.nstr(ratio, 20)))

# Also try other fractional shifts
print("\n=== Try other shifts ===")
for shift_num in [5, 11, 17, 22, 39]:
    shift = mpmath.mpf(shift_num) / 22
    nu = 0 + shift
    try:
        result = J_nu_contour(nu)
        j_val = result[0]
        print("J_{%d/22} = %s" % (shift_num, mpmath.nstr(j_val, 20)))
    except Exception as ex:
        print("J_{%d/22} failed: %s" % (shift_num, ex))
