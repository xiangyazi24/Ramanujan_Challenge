#!/usr/bin/env sage
"""
Compute shifted Barnes integral via contour integration.
No ore_algebra needed.
"""
import mpmath
mpmath.mp.dps = 100

def cotpi_mp(t):
    return mpmath.cos(mpmath.pi * t) / mpmath.sin(mpmath.pi * t)

def R_nu(nu, t):
    """R_ν(t) = Γ(ν+1-t)³ · Γ(t) / [Γ(1-t)³ · Γ(ν+1)² · Γ(t+ν+1)]"""
    return (mpmath.gamma(nu+1-t)**3 * mpmath.gamma(t) /
            (mpmath.gamma(1-t)**3 * mpmath.gamma(nu+1)**2 * mpmath.gamma(t+nu+1)))

def integrand(nu, t):
    r = R_nu(nu, t)
    s = mpmath.pi / mpmath.sin(mpmath.pi * t)
    c = mpmath.pi * cotpi_mp(t) - 1
    return r * s**2 * c

def J_nu(nu, T=50):
    """J_ν via vertical contour at Re(t) = 1/2.
    (1/2πi) ∮ F(t) dt with dt = i·dy ⟹ integrand·i/(2πi) = integrand/(2π)
    """
    def f(y):
        t = mpmath.mpf('0.5') + 1j * y
        return integrand(nu, t) / (2 * mpmath.pi)
    return mpmath.quad(f, [-T, T], error=True, maxdegree=8)

L = mpmath.zeta(2) + mpmath.zeta(3)
print("ζ(2)+ζ(3) = %s" % mpmath.nstr(L, 30))

# ============================================================
# P2.7 errors (exact arithmetic)
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

N = 10
q = [QQ(0)] * N
q[0] = QQ(-215040420000)
q[1] = QQ(-167282265043404) / QQ(905)
q[2] = QQ(-964185327658080) / QQ(6071)

pp = [QQ(0)] * N
pp[0] = QQ(-612218384750)
pp[1] = QQ(-9525021973931919) / QQ(18100)
pp[2] = QQ(-29561828382772029) / QQ(65380)

for n in range(2, N-1):
    q[n+1] = B_p27(n)/A_p27(n)*q[n] - C_p27(n-1)/A_p27(n-1)*q[n-1] + D_p27(n-2)/A_p27(n-2)*q[n-2]
    pp[n+1] = B_p27(n)/A_p27(n)*pp[n] - C_p27(n-1)/A_p27(n-1)*pp[n-1] + D_p27(n-2)/A_p27(n-2)*pp[n-2]

e_p27 = []
print("\nP2.7 errors:")
for i in range(min(8, N)):
    ep = mpmath.mpf(pp[i].numerator())/mpmath.mpf(pp[i].denominator()) - \
         L * mpmath.mpf(q[i].numerator())/mpmath.mpf(q[i].denominator())
    e_p27.append(ep)
    print("  e_%d = %s" % (i, mpmath.nstr(ep, 20)))

# ============================================================
# Test 1: Zudilin J_n at integer n
# ============================================================
print("\n=== Zudilin J_n at integer ν ===")
for n_int in [0, 1, 2, 3]:
    nu = mpmath.mpf(n_int)
    print("J_%d: " % n_int, end="", flush=True)
    try:
        result = J_nu(nu)
        print("%s (err %s)" % (mpmath.nstr(result[0], 20), mpmath.nstr(result[1], 5)))
    except Exception as ex:
        print("FAILED: %s" % ex)

# ============================================================
# Test 2: Shifted J_{n+17/22}
# ============================================================
print("\n=== Shifted J_{n+17/22} ===")
J_sh = []
for n_val in range(6):
    nu = n_val + mpmath.mpf(17)/22
    print("J_{%d+17/22}: " % n_val, end="", flush=True)
    try:
        result = J_nu(nu)
        j_val = result[0]
        J_sh.append(j_val)
        print("%s (err %s)" % (mpmath.nstr(j_val, 20), mpmath.nstr(result[1], 5)))
    except Exception as ex:
        print("FAILED: %s" % ex)
        J_sh.append(None)

# ============================================================
# Compare
# ============================================================
print("\n=== J_{n+17/22} / e^P_n ===")
for i in range(min(len(J_sh), len(e_p27))):
    if J_sh[i] is not None and abs(e_p27[i]) > 1e-90:
        r = J_sh[i] / e_p27[i]
        print("n=%d: %s" % (i, mpmath.nstr(r, 15)))

print("\n=== J_{n+17/22} / (64^{-n} · e^P_n) ===")
for i in range(min(len(J_sh), len(e_p27))):
    if J_sh[i] is not None and abs(e_p27[i]) > 1e-90:
        r = J_sh[i] / (mpmath.power(64, -i) * e_p27[i])
        print("n=%d: %s" % (i, mpmath.nstr(r, 15)))

# ============================================================
# Also try shifts 5/22, 11/22, 1/22 — find which one fits
# ============================================================
print("\n=== J_{ν} for various ν ===")
for shift_num in [1, 5, 11, 17, 22, 39, 83]:
    shift = mpmath.mpf(shift_num) / 22
    print("J_{%d/22}: " % shift_num, end="", flush=True)
    try:
        result = J_nu(shift, T=30)
        print("%s" % mpmath.nstr(result[0], 15))
    except Exception as ex:
        print("FAILED: %s" % ex)
