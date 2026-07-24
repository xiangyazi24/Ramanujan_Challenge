#!/usr/bin/env sage
"""
FIXED precision check. The recurrence is:

  u_{n+1} = (B_n/A_n) u_n - (C_{n-1}/A_{n-1}) u_{n-1} + (D_{n-2}/A_{n-2}) u_{n-2}

NOT: u_{n+1} = (B_n u_n - C_{n-1} u_{n-1} + D_{n-2} u_{n-2}) / A_n

The denominators A_n, A_{n-1}, A_{n-2} are DIFFERENT for each term!
"""
import mpmath
mpmath.mp.dps = 2000

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

# CORRECT recurrence: different A denominators for each term
# u_{n+1} = (B_n/A_n) u_n - (C_{n-1}/A_{n-1}) u_{n-1} + (D_{n-2}/A_{n-2}) u_{n-2}
# For n >= 2
q = [QQ(0)] * N
q[0] = QQ(-215040420000)
q[1] = QQ(-167282265043404) / QQ(905)
q[2] = QQ(-964185327658080) / QQ(6071)

for n in range(2, N-1):
    q[n+1] = B_p27(n)/A_p27(n) * q[n] - C_p27(n-1)/A_p27(n-1) * q[n-1] + D_p27(n-2)/A_p27(n-2) * q[n-2]

p = [QQ(0)] * N
p[0] = QQ(-612218384750)
p[1] = QQ(-9525021973931919) / QQ(18100)
p[2] = QQ(-29561828382772029) / QQ(65380)

for n in range(2, N-1):
    p[n+1] = B_p27(n)/A_p27(n) * p[n] - C_p27(n-1)/A_p27(n-1) * p[n-1] + D_p27(n-2)/A_p27(n-2) * p[n-2]

# Check first few values match
print("=== First few values ===")
for n in range(5):
    print("q_%d = %s" % (n, q[n]))
    print("p_%d = %s" % (n, p[n]))

# Compare with OLD (wrong) recurrence
q_old = [QQ(0)] * 5
q_old[0] = QQ(-215040420000)
q_old[1] = QQ(-167282265043404) / QQ(905)
q_old[2] = QQ(-964185327658080) / QQ(6071)

for n in range(2, 4):
    q_old[n+1] = (B_p27(n) * q_old[n] - C_p27(n-1) * q_old[n-1] + D_p27(n-2) * q_old[n-2]) / A_p27(n)

print("\n=== Comparing old vs new at n=3 ===")
print("q_3 (correct): %s" % q[3])
print("q_3 (old/wrong): %s" % q_old[3])
print("Difference: %s" % (q[3] - q_old[3]))

# Compute L at 2000 digits
L = mpmath.zeta(2) + mpmath.zeta(3)
print("\nL = ζ(2)+ζ(3) (50 digits): %s" % mpmath.nstr(L, 50))

# p/q - L
print("\n=== p_n/q_n - L (FIXED recurrence) ===")
for n in range(min(30, N)):
    if q[n] != 0:
        p_mp = mpmath.mpf(p[n].numerator()) / mpmath.mpf(p[n].denominator())
        q_mp = mpmath.mpf(q[n].numerator()) / mpmath.mpf(q[n].denominator())
        diff = p_mp / q_mp - L
        if n <= 20 or n % 5 == 0:
            print("n=%2d: p/q - L = %s" % (n, mpmath.nstr(diff, 40)))

# Errors
print("\n=== e_n = p_n - L·q_n ===")
e_vals = []
for n in range(min(40, N)):
    p_mp = mpmath.mpf(p[n].numerator()) / mpmath.mpf(p[n].denominator())
    q_mp = mpmath.mpf(q[n].numerator()) / mpmath.mpf(q[n].denominator())
    e = p_mp - L * q_mp
    e_vals.append(e)
    if n <= 15:
        print("n=%2d: e = %s" % (n, mpmath.nstr(e, 30)))

# Error decay rate
print("\n=== e_{n+1}/e_n ===")
for n in range(min(30, len(e_vals)-1)):
    if abs(e_vals[n]) > 0:
        ratio = e_vals[n+1] / e_vals[n]
        print("n=%2d: %s" % (n, mpmath.nstr(ratio, 20)))

# q ratio (Poincaré)
print("\n=== q_{n+1}/q_n ===")
for n in range(min(20, N-1)):
    if q[n] != 0:
        r = float(q[n+1]/q[n])
        print("n=%2d: %.15f" % (n, r))

# |e/q| pattern
print("\n=== |e_n/q_n| ===")
for n in range(min(30, len(e_vals))):
    if q[n] != 0:
        q_mp = mpmath.mpf(q[n].numerator()) / mpmath.mpf(q[n].denominator())
        ratio = abs(e_vals[n] / q_mp)
        print("n=%2d: %s" % (n, mpmath.nstr(ratio, 25)))
