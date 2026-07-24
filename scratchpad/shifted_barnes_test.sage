#!/usr/bin/env sage
"""
Test: Does the shifted Barnes integral J_{n+17/22} satisfy the P2.7 recurrence?

R_ν(t) = Γ(ν+1-t)³ · Γ(t) / [Γ(1-t)³ · Γ(ν+1)² · Γ(t+ν+1)]

J_ν = (1/2πi) ∮ R_ν(t) · [π/sin(πt)]² · [πcot(πt) - 1] dt

For integer ν=n, this gives the Zudilin companion: J_n = p̃_n - (ζ(2)+ζ(3))a_n

Hypothesis: J_{n+17/22} is proportional to e_n = p_n - (ζ(2)+ζ(3))q_n (P2.7 error)
"""
import mpmath
mpmath.mp.dps = 100

# ============================================================
# 1. Define the integrand for general ν
# ============================================================
def R_nu(nu, t):
    """R_ν(t) = Γ(ν+1-t)³ · Γ(t) / [Γ(1-t)³ · Γ(ν+1)² · Γ(t+ν+1)]"""
    return (mpmath.gamma(nu+1-t)**3 * mpmath.gamma(t) /
            (mpmath.gamma(1-t)**3 * mpmath.gamma(nu+1)**2 * mpmath.gamma(t+nu+1)))

def integrand(nu, t):
    """Full integrand: R_ν(t) · [π/sin(πt)]² · [πcot(πt) - 1]"""
    r = R_nu(nu, t)
    s = mpmath.pi / mpmath.sinpi(t)
    c = mpmath.pi * mpmath.cotpi(t) - 1
    return r * s**2 * c

# ============================================================
# 2. Compute via contour integration
# ============================================================
# The contour should separate left poles (t ≤ 0) from right poles.
# For ν non-integer, right poles are at t = ν+1+k (k=0,1,...) and integer t = 1,2,...
# For ν = 17/22 ≈ 0.773, ν+1 = 39/22 ≈ 1.773
# So contour at Re(t) = 1/2 separates t=0,-1,... from t=1,2,...,ν+1,...

def J_nu_contour(nu, N_poles=None):
    """Compute J_ν via vertical contour integral at Re(t) = 1/2.

    (1/2πi) ∫_{1/2-i∞}^{1/2+i∞} F(t) dt
    """
    def f(y):
        t = mpmath.mpf('0.5') + 1j * y
        return integrand(nu, t) / (2 * mpmath.pi * 1j)

    result = mpmath.quad(f, [-40, 40], error=True)
    return result

# ============================================================
# 3. Alternative: compute via residue sum
# ============================================================
def J_nu_residues(nu, max_poles=50):
    """Compute J_ν by summing residues at poles to the RIGHT of contour Re(t) = 1/2.

    Right poles come from:
    (a) Integer poles t = m (m = 1, 2, 3, ...) where R_ν and [π/sin πt]² interact
    (b) Gamma poles t = ν+1+k (k = 0, 1, 2, ...) from Γ(ν+1-t)³
    """
    total = mpmath.mpf(0)

    # Type (a): integer poles t = m ≥ 1
    # Near t = m: R_ν(t) is regular (for m ≤ floor(ν+1)) or has poles
    # [π/sin πt]² ~ 1/(t-m)², πcot πt ~ 1/(t-m) + O(t-m)
    # So [πcot πt - 1] ~ 1/(t-m) - 1 + O(t-m)
    # Full: R_ν · 1/(t-m)² · [1/(t-m) - 1] = R_ν/(t-m)³ - R_ν/(t-m)²
    #
    # Residue = (1/2) R_ν''(m) - R_ν'(m) + correction from cot expansion
    #
    # This is complex. Let me just use numerical differentiation.

    for m in range(1, max_poles + 1):
        # Compute residue of integrand at t = m numerically
        # The integrand has a pole of order ≤ 3 at integer points
        # Residue = (1/2πi) ∮ F(t) dt around t=m
        eps = mpmath.mpf('1e-15')

        def circle_integrand(theta):
            t = m + eps * mpmath.exp(1j * theta)
            return integrand(nu, t) * eps * 1j * mpmath.exp(1j * theta) / (2 * mpmath.pi * 1j)

        try:
            res = mpmath.quad(circle_integrand, [0, 2*mpmath.pi])
            if abs(res) > mpmath.mpf('1e-80'):
                total += res
        except:
            pass

    # Type (b): Gamma poles t = ν+1+k
    for k in range(max_poles):
        pole = nu + 1 + k
        if pole <= 0.5:  # skip if on wrong side of contour
            continue
        eps = mpmath.mpf('1e-15')

        def circle_integrand_b(theta):
            t = pole + eps * mpmath.exp(1j * theta)
            return integrand(nu, t) * eps * 1j * mpmath.exp(1j * theta) / (2 * mpmath.pi * 1j)

        try:
            res = mpmath.quad(circle_integrand_b, [0, 2*mpmath.pi])
            if abs(res) > mpmath.mpf('1e-80'):
                total += res
        except:
            pass

    return total

# ============================================================
# 4. First test: compute J_{17/22} and compare with e_0
# ============================================================
print("=== Testing shifted Barnes integral ===")
print("Computing J_{17/22} (ν = 17/22)...")

nu_shift = mpmath.mpf(17) / 22

# Try contour integral first
print("\n--- Contour integral method ---")
try:
    result = J_nu_contour(nu_shift)
    print("J_{17/22} = %s" % mpmath.nstr(result[0], 30))
    print("Error estimate: %s" % mpmath.nstr(result[1], 5))
except Exception as ex:
    print("Contour failed: %s" % ex)

# ============================================================
# 5. Compare with known values
# ============================================================
# Zudilin at integer ν: J_0 = 0 (trivially), J_1, J_2 are known
print("\n--- Zudilin at integer ν (sanity check) ---")
for n_int in [1, 2, 3]:
    try:
        result = J_nu_contour(mpmath.mpf(n_int))
        print("J_%d = %s" % (n_int, mpmath.nstr(result[0], 20)))
    except Exception as ex:
        print("J_%d failed: %s" % (n_int, ex))

# Compare with known a_n and p̃_n
a_list = []
for i in range(10):
    val = QQ(0)
    for k in range(i+1):
        val += binomial(i,k)^2 * binomial(i+k,i) * binomial(i+2*k,i)
    a_list.append(val)

print("\na_n values:")
for i in range(5):
    print("  a_%d = %s" % (i, a_list[i]))

L = mpmath.zeta(2) + mpmath.zeta(3)
print("\nζ(2)+ζ(3) = %s" % mpmath.nstr(L, 30))

# Zudilin companion p̃_n
from ore_algebra import OreAlgebra, guess
Rn = PolynomialRing(QQ, 'nn')
OS = OreAlgebra(Rn, 'Snn')
rec_a_op = guess(a_list[:10], OS, order=3)
P_rec = [rec_a_op[j] for j in range(4)]

p_zud = [QQ(0)] * 10
p_zud[0] = QQ(0)
p_zud[1] = QQ(20)
p_zud[2] = QQ(7425) / QQ(16)

for n_val in range(7):
    val = sum(P_rec[j](nn=n_val) * p_zud[n_val + j] for j in range(3))
    p_zud[n_val + 3] = -val / P_rec[3](nn=n_val)

print("\nZudilin errors e^Z_n = p̃_n - L·a_n:")
for i in range(5):
    ez = mpmath.mpf(p_zud[i]) - L * mpmath.mpf(a_list[i])
    print("  e^Z_%d = %s" % (i, mpmath.nstr(ez, 20)))

# ============================================================
# 6. P2.7 errors for comparison
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

p = [QQ(0)] * N
p[0] = QQ(-612218384750)
p[1] = QQ(-9525021973931919) / QQ(18100)
p[2] = QQ(-29561828382772029) / QQ(65380)

for n in range(2, N-1):
    q[n+1] = B_p27(n)/A_p27(n) * q[n] - C_p27(n-1)/A_p27(n-1) * q[n-1] + D_p27(n-2)/A_p27(n-2) * q[n-2]
    p[n+1] = B_p27(n)/A_p27(n) * p[n] - C_p27(n-1)/A_p27(n-1) * p[n-1] + D_p27(n-2)/A_p27(n-2) * p[n-2]

print("\nP2.7 errors e^P_n = p_n - L·q_n:")
e_p27 = []
for i in range(min(8, N)):
    ep = mpmath.mpf(p[i].numerator()) / mpmath.mpf(p[i].denominator()) - \
         L * mpmath.mpf(q[i].numerator()) / mpmath.mpf(q[i].denominator())
    e_p27.append(ep)
    print("  e^P_%d = %s" % (i, mpmath.nstr(ep, 20)))

# ============================================================
# 7. Compute J_{n+17/22} for n = 0, 1, 2, 3 and compare with e^P_n
# ============================================================
print("\n=== Shifted Barnes integral J_{n+17/22} ===")
J_shifted = []
for n_val in range(4):
    nu = n_val + mpmath.mpf(17)/22
    try:
        result = J_nu_contour(nu)
        J_shifted.append(result[0])
        print("J_{%d+17/22} = %s (err est: %s)" % (n_val, mpmath.nstr(result[0], 20), mpmath.nstr(result[1], 5)))
    except Exception as ex:
        print("J_{%d+17/22} failed: %s" % (n_val, ex))
        J_shifted.append(None)

# Compare ratios
print("\n=== Ratio J_{n+17/22} / e^P_n ===")
for i in range(min(len(J_shifted), len(e_p27))):
    if J_shifted[i] is not None and abs(e_p27[i]) > 0:
        ratio = J_shifted[i] / e_p27[i]
        print("n=%d: J/e = %s" % (i, mpmath.nstr(ratio, 15)))
