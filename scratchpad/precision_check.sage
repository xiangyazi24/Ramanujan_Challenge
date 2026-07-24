#!/usr/bin/env sage
"""
Ultra-high-precision check of whether p_n/q_n → ζ(2)+ζ(3) for P2.7.

The error_ratio.sage computation shows e_{n+1}/e_n → 0.859 (dominant root),
suggesting c₀(e) ≠ 0. But this could be a precision artifact.

With 2000 digits of L = ζ(2)+ζ(3), compute:
1. p_n/q_n - L at high precision
2. The convergence behavior
3. Whether e_n/q_n → 0 or → nonzero constant
"""
import mpmath
mpmath.mp.dps = 2000

# ============================================================
# Step 1: Compute P2.7 sequences in exact arithmetic
# ============================================================
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

for n in range(2, N-1):
    q[n+1] = (B_p27(n) * q[n] - C_p27(n-1) * q[n-1] + D_p27(n-2) * q[n-2]) / A_p27(n)

p = [QQ(0)] * N
p[0] = QQ(-612218384750)
p[1] = QQ(-9525021973931919) / QQ(18100)
p[2] = QQ(-29561828382772029) / QQ(65380)

for n in range(2, N-1):
    p[n+1] = (B_p27(n) * p[n] - C_p27(n-1) * p[n-1] + D_p27(n-2) * p[n-2]) / A_p27(n)

# ============================================================
# Step 2: Compute L = ζ(2)+ζ(3) at 2000 digits
# ============================================================
L = mpmath.zeta(2) + mpmath.zeta(3)
print("L = ζ(2)+ζ(3) (first 50 digits):")
print("  %s" % mpmath.nstr(L, 50))

# ============================================================
# Step 3: Compute p/q - L at ultra-high precision
# ============================================================
print("\n=== p_n/q_n - L at 2000-digit precision ===")
for n in range(min(40, N)):
    if q[n] != 0:
        # Convert to mpmath for high precision
        p_mp = mpmath.mpf(p[n].numerator()) / mpmath.mpf(p[n].denominator())
        q_mp = mpmath.mpf(q[n].numerator()) / mpmath.mpf(q[n].denominator())

        ratio = p_mp / q_mp
        diff = ratio - L
        if n <= 20:
            print("n=%2d: p/q - L = %s" % (n, mpmath.nstr(diff, 30)))

# ============================================================
# Step 4: Compute e_n = p_n - L·q_n and e_n/q_n
# ============================================================
print("\n=== e_n = p_n - L·q_n ===")
e_vals = []
for n in range(min(40, N)):
    p_mp = mpmath.mpf(p[n].numerator()) / mpmath.mpf(p[n].denominator())
    q_mp = mpmath.mpf(q[n].numerator()) / mpmath.mpf(q[n].denominator())
    e = p_mp - L * q_mp
    e_vals.append(e)
    if n <= 20:
        print("n=%2d: e = %s, e/q = %s" % (n, mpmath.nstr(e, 25), mpmath.nstr(e/q_mp, 25) if q_mp != 0 else "?"))

# ============================================================
# Step 5: Error decay rate
# ============================================================
print("\n=== e_{n+1}/e_n ===")
for n in range(min(35, len(e_vals)-1)):
    if abs(e_vals[n]) > 0:
        ratio = e_vals[n+1] / e_vals[n]
        print("n=%2d: %s" % (n, mpmath.nstr(ratio, 20)))

# ============================================================
# Step 6: q_{n+1}/q_n (Poincaré)
# ============================================================
print("\n=== q_{n+1}/q_n ===")
for n in range(min(35, N-1)):
    if q[n] != 0:
        r = float(q[n+1]/q[n])
        print("n=%2d: %.15f" % (n, r))

# ============================================================
# Step 7: Check absolute value pattern of e/q
# ============================================================
print("\n=== |e_n/q_n| pattern ===")
for n in range(min(35, len(e_vals))):
    if q[n] != 0:
        q_mp = mpmath.mpf(q[n].numerator()) / mpmath.mpf(q[n].denominator())
        ratio = abs(e_vals[n] / q_mp)
        print("n=%2d: |e/q| = %s" % (n, mpmath.nstr(ratio, 20)))
