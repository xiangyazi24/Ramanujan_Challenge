#!/usr/bin/env sage
"""
Test shifted Barnes integral J_{n+17/22} against P2.7 errors.
"""
import mpmath
mpmath.mp.dps = 150

def cotpi(t):
    return mpmath.cos(mpmath.pi * t) / mpmath.sin(mpmath.pi * t)

def R_nu(nu, t):
    return (mpmath.gamma(nu+1-t)**3 * mpmath.gamma(t) /
            (mpmath.gamma(1-t)**3 * mpmath.gamma(nu+1)**2 * mpmath.gamma(t+nu+1)))

def integrand(nu, t):
    r = R_nu(nu, t)
    s = mpmath.pi / mpmath.sin(mpmath.pi * t)
    c = mpmath.pi * cotpi(t) - 1
    return r * s**2 * c

def J_nu_contour(nu):
    """Compute J_ν via vertical contour at Re(t) = 1/2."""
    def f(y):
        t = mpmath.mpf('0.5') + 1j * y
        val = integrand(nu, t)
        return val / (2 * mpmath.pi * 1j)

    result = mpmath.quad(f, [-50, 50], error=True, maxdegree=8)
    return result

# ============================================================
# Sanity check: Zudilin at integer ν
# ============================================================
print("=== Sanity check: Zudilin J_n for integer n ===")

L = mpmath.zeta(2) + mpmath.zeta(3)
print("ζ(2)+ζ(3) = %s" % mpmath.nstr(L, 30))

# Compute a_n
a_list = []
for i in range(30):
    val = QQ(0)
    for k in range(i+1):
        val += binomial(i,k)^2 * binomial(i+k,i) * binomial(i+2*k,i)
    a_list.append(val)

# Zudilin companion via explicit sum
# p̃_n = Σ_{k=0}^n C(n,k)²C(n+k,k)C(n+2k,k) · S_n(k) where S_n(k) involves harmonic sums
# Alternatively, use the recurrence with known initial conditions
# p̃_0 = 0, p̃_1 = 20, p̃_2 = 7425/16
#
# rec_a: P_0(n)a_n + P_1(n)a_{n+1} + P_2(n)a_{n+2} + P_3(n)a_{n+3} = 0
# a_n satisfies (n+1)⁴a_{n+1} = (2n+1)(3+34n²+51n³+17n⁴)a_n - n⁴a_{n-1}
# Wait, AESZ #209 is order 3, but the well-known recurrence for a_n is order 2!
# Let me check.

# Actually, AESZ #209 operator is order 4 as a differential operator (CY4).
# The recurrence for a_n is order 2 (Apéry-like):
# (n+1)⁴ a_{n+1} - (2n+1)(17n² + 17n + 5) a_n + n⁴ a_{n-1} = 0
# Wait, that's for a different sequence. Let me check AESZ #209 specifically.

# AESZ #209: a_n = Σ C(n,k)²C(n+k,n)C(n+2k,n)
# Actually this might have order 3 recurrence.

# Let me just compute using ore_algebra with more terms
from ore_algebra import OreAlgebra, guess
Rn = PolynomialRing(QQ, 'nn')
nn_var = Rn.gen()
OS = OreAlgebra(Rn, 'Snn')

rec_a_op = guess(a_list[:25], OS, order=3)
print("rec_a order: %d" % rec_a_op.order())
P_rec = [rec_a_op[j] for j in range(rec_a_op.order()+1)]
for j in range(len(P_rec)):
    print("  P_%d degree %d" % (j, P_rec[j].degree()))

# Compute Zudilin companion
N_zud = 10
p_zud = [QQ(0)] * N_zud
p_zud[0] = QQ(0)
p_zud[1] = QQ(20)
p_zud[2] = QQ(7425) / QQ(16)

for n_val in range(N_zud - rec_a_op.order()):
    val = sum(P_rec[j](nn=n_val) * p_zud[n_val + j] for j in range(rec_a_op.order()))
    p_zud[n_val + rec_a_op.order()] = -val / P_rec[rec_a_op.order()](nn=n_val)

print("\nZudilin errors:")
e_zud = []
for i in range(min(6, N_zud)):
    ez = mpmath.mpf(p_zud[i]) - L * mpmath.mpf(a_list[i])
    e_zud.append(ez)
    print("  e^Z_%d = %s" % (i, mpmath.nstr(ez, 25)))

# Now compute J_n for integer n=1,2,3 and compare
print("\n--- J_n via contour integral ---")
for n_int in [1, 2, 3]:
    try:
        result = J_nu_contour(mpmath.mpf(n_int))
        j_val = result[0]
        print("J_%d = %s (err: %s)" % (n_int, mpmath.nstr(j_val, 25), mpmath.nstr(result[1], 5)))
        if n_int < len(e_zud):
            ratio = j_val / e_zud[n_int]
            print("  J_%d / e^Z_%d = %s" % (n_int, n_int, mpmath.nstr(ratio, 15)))
    except Exception as ex:
        print("J_%d failed: %s" % (n_int, ex))

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

N_p27 = 10
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
# Shifted Barnes: J_{n+17/22}
# ============================================================
print("\n=== Shifted Barnes J_{n+17/22} ===")
J_shifted = []
for n_val in range(5):
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

# Compare
print("\n=== Ratio J_{n+17/22} / e^P_n ===")
for i in range(min(len(J_shifted), len(e_p27))):
    if J_shifted[i] is not None and abs(e_p27[i]) > 1e-100:
        ratio = J_shifted[i] / e_p27[i]
        print("n=%d: %s" % (i, mpmath.nstr(ratio, 20)))

# Also check if there's a simple gauge
print("\n=== J_{n+17/22} / (64^n · e^P_n) ===")
for i in range(min(len(J_shifted), len(e_p27))):
    if J_shifted[i] is not None and abs(e_p27[i]) > 1e-100:
        ratio = J_shifted[i] / (mpmath.power(64, i) * e_p27[i])
        print("n=%d: %s" % (i, mpmath.nstr(ratio, 20)))
