#!/usr/bin/env sage
"""
Search for the nonlocal transform T(n,k) that maps Zudilin's a_k to P2.7's q_n:
  q_n = Σ_k T(n,k) · a_k

Step 1: Compute many P2.7 q_n values using the recurrence
Step 2: See if q_n = linear combination of a_0,...,a_n with structured coefficients
Step 3: Identify the structure of T(n,k)
"""
from ore_algebra import OreAlgebra, guess
import mpmath
mpmath.mp.dps = 100

# ============================================================
# Step 1: Compute a_n (Zudilin inner sum)
# ============================================================
N = 40
a_vals = []
for i in range(N):
    s = QQ(0)
    for k in range(i+1):
        s += binomial(i,k)^2 * binomial(i+k,i) * binomial(i+2*k,i)
    a_vals.append(s)

# ============================================================
# Step 2: Compute P2.7 q_n values using the recurrence
# ============================================================
# Recurrence: A_n u_{n+1} = B_n u_n - C_{n-1} u_{n-1} + D_{n-2} u_{n-2}
def A(n):
    return QQ(1024) * (2*n+5)^4 * (2*n+7)^3 * (2*n+9)^3 * (946*n^2+6407*n+10860)

def B(n):
    return QQ(128) * (2*n+7)^3 * (2*n+9)^3 * (104060*n^6 + 1745370*n^5 +
        12145238*n^4 + 44886481*n^3 + 92943995*n^2 + 102256019*n + 46709052)

def C(n):
    return QQ(16) * (n+3)^4 * (2*n+9)^3 * (3784*n^5 + 57792*n^4 +
        351019*n^3 + 1059230*n^2 + 1587211*n + 944620)

def D(n):
    return QQ(1) * (n+3)^4 * (n+4)^6 * (946*n^2 + 4515*n + 5399)

# Initial conditions
q = [QQ(0)] * N
q[0] = QQ(-215040420000)
q[1] = QQ(-167282265043404) / QQ(905)
q[2] = QQ(-964185327658080) / QQ(6071)

# Forward recurrence
for n in range(2, N-1):
    q[n+1] = (B(n) * q[n] - C(n-1) * q[n-1] + D(n-2) * q[n-2]) / A(n)

# Also compute p_n
p = [QQ(0)] * N
p[0] = QQ(-612218384750)
p[1] = QQ(-9525021973931919) / QQ(18100)
p[2] = QQ(-29561828382772029) / QQ(65380)

for n in range(2, N-1):
    p[n+1] = (B(n) * p[n] - C(n-1) * p[n-1] + D(n-2) * p[n-2]) / A(n)

# Verify convergence to ζ(2)+ζ(3)
L = float(mpmath.zeta(2) + mpmath.zeta(3))
print("=== p_n/q_n convergence ===")
for n in range(min(15, N)):
    if q[n] != 0:
        r = float(p[n]/q[n])
        print("n=%2d: p/q = %.15f, diff = %.3e" % (n, r, abs(r - L)))

# ============================================================
# Step 3: Try to express q_n as Σ T(n,k) a_k
# ============================================================
# The simplest nonlocal transforms:
# (a) q_n = c_n · a_n (diagonal gauge) — already know this fails
# (b) q_n = Σ_{k=0}^n C(n,k) · r^k · a_k (binomial transform)
# (c) q_n = Σ_{k=0}^n C(2n,2k) · a_k (even binomial)

print("\n=== Search for nonlocal kernel ===")

# First, compute the ratio q_n / a_n to see the gauge factor
print("\nq_n / a_n:")
for n in range(min(10, N)):
    if a_vals[n] != 0:
        r = q[n] / a_vals[n]
        print("  n=%d: %.6e" % (n, float(r)))

# Check if q_n / (C(2n,n) * a_n) gives a pattern
print("\nq_n / (C(2n,n) * a_n) = q_n / A_n^{AESZ}:")
for n in range(min(10, N)):
    An = binomial(2*n,n) * a_vals[n]
    if An != 0:
        r = q[n] / An
        print("  n=%d: %s = %.6e" % (n, r, float(r)))

# Try binomial transform: does q_n = Σ C(n,k) · a_k?
print("\nBinomial transform Σ C(n,k) a_k:")
for n in range(min(10, N)):
    bt = sum(binomial(n,k) * a_vals[k] for k in range(n+1))
    print("  n=%d: bt = %s, q = %.6e, bt/q = %.6e" %
          (n, bt, float(q[n]), float(bt/q[n]) if q[n] != 0 else float('inf')))

# Try 64^n scaling: q_n / 64^n
print("\nq_n / 64^n:")
for n in range(min(10, N)):
    print("  n=%d: %.6e" % (n, float(q[n] / QQ(64)^n)))

# The ratio q_{n+1}/q_n (Poincaré)
print("\nq_{n+1}/q_n:")
for n in range(min(15, N-1)):
    if q[n] != 0:
        print("  n=%d: %.6f" % (n, float(q[n+1]/q[n])))

# ============================================================
# Step 4: Key insight from Q5048 — try the Zudilin v_n = 64^{-n} J_n^Z
# The scaled error ê = 64^{-n} e_n^Z should have rate ~0.00105
# The P2.7 error e_n^{P27} should also have rate ~0.00105
# If the RATES match but the recurrences differ, maybe
# e_n^{P27} = nonlocal transform of e_n^Z
# ============================================================

# Compute the Zudilin companion sequences using partial fractions
# R_n(t) = [(t-1)...(t-n)]^3 / [(n!)^2 * t*(t+1)*...*(t+n)]
print("\n=== Zudilin companion sequences ===")

def R_n(n_val, t):
    """R_n(t) = [(t-1)...(t-n)]^3 / [(n!)^2 * t*(t+1)*...*(t+n)]"""
    num = prod(t - j for j in range(1, n_val+1))^3
    den = factorial(n_val)^2 * prod(t + j for j in range(n_val+1))
    return num / den

# Compute b̃_n and b̃̃_n from the partial fraction decomposition
# R_n(t) has poles at t = 0, -1, ..., -n
# r_n(z) = Σ_ν z^ν R_n(ν) = a_n Li_1(z) - b_n
# r̃_n(z) = -Σ_ν z^ν R_n'(ν) = a_n Li_2(z) - b̃_n
# r̃̃_n(z) = ½ Σ_ν z^ν R_n''(ν) = a_n Li_3(z) - b̃̃_n

# The b̃_n can be computed from the partial fraction coefficients at the POLES
# At t = -m (m = 0,...,n):
# R_n(t) = Σ_m Σ_{j=1}^{order} A_{m,j} / (t+m)^j + polynomial

# Actually, it's easier to compute b̃_n from the RESIDUE formula:
# b̃_n = a_n * ζ(2) - r̃_n
# But r̃_n at z=1 diverges...

# Alternative: compute b̃_n from the explicit formula
# b̃_n = partial fraction constant from R_n'
# Use the relation: b̃_n = Σ_{m=0}^n Res_{t=-m} [R_n(t) π²/sin²(πt)]

# For computational purposes, b̃_n can be extracted from:
# a_n ζ(2) - b̃_n = lim_{z→1-} [Σ_ν z^ν R_n'(ν) + (Σ_ν z^ν/ν²)] ... too complicated

# Let me instead compute b̃_n + b̃̃_n DIRECTLY from the initial conditions
# Zudilin: p̃_n = b̃_n + b̃̃_n, with p̃_0 = 0, p̃_1 = 20, p̃_2 = 7425/16
# These satisfy rec_a (same as a_n).

# Verify: do p̃ values satisfy rec_a?
Rn = PolynomialRing(QQ, 'n')
n_var = Rn.gen()
OS = OreAlgebra(Rn, 'Sn')
Sn = OS.gen()
rec_a = guess(list(a_vals[:30]), OS, order=3)
P_rec = [rec_a[j] for j in range(4)]

# Compute p̃ using rec_a with given initial conditions
p_zud = [QQ(0)] * N
p_zud[0] = QQ(0)
p_zud[1] = QQ(20)
p_zud[2] = QQ(7425) / QQ(16)

for nn in range(N-3):
    val = sum(QQ(P_rec[j](n=nn)) * p_zud[nn+j] for j in range(3))
    p_zud[nn+3] = -val / QQ(P_rec[3](n=nn))

# Verify convergence of p̃_n/a_n to ζ(2)+ζ(3)
print("p̃_n / a_n → ζ(2)+ζ(3) ?")
for n in range(min(20, N)):
    if a_vals[n] != 0:
        r = float(p_zud[n] / a_vals[n])
        print("  n=%2d: p̃/a = %.15f, diff = %.3e" % (n, r, abs(r - L)))

# ============================================================
# Step 5: Core test — can q_n^{P27} be expressed in terms of a_k and p̃_k?
# ============================================================
# If q_n^{P27} = Σ T(n,k) a_k, then p_n^{P27} = Σ T(n,k) p̃_k
# And e_n^{P27} = p_n - L*q_n = Σ T(n,k) (p̃_k - L*a_k) = Σ T(n,k) e_k^Z → 0

# Try: q_n^{P27} = g(n) · Σ C(2n,2k) a_k
# Try: q_n^{P27} = g(n) · Σ C(n,k) 64^{-k} a_k
# Try: q_n^{P27} = g(n) · Σ h_k a_k  (Pochhammer-weighted)

# Most promising from Q5048: T(n,k) = C(αn+β, k) · c^k · Pochhammer ratio
# Let me search by solving for T(n,0), T(n,1), ... from q_n = Σ T(n,k) a_k

print("\n=== Solving for T(n,k) assuming q_n = Σ_{k=0}^n T(n,k) a_k ===")
# For each n, T(n,k) for k=0,...,n are unknowns. But we only have one equation per n.
# So T is underdetermined unless we impose structure.

# Instead: check if q_n / some_function = known transform of a_k
# Key: q_0 = -2^5 · 3^6 · 5^4 · 7^3 · 43 = -215040420000
# And a_0 = 1, so T(0,0) = q_0

# For n=1: q_1 = T(1,0)*a_0 + T(1,1)*a_1 = T(1,0) + 7*T(1,1)
# This has infinitely many solutions. We need structure.

# Try T(n,k) = g(n) · f(n,k) where f is a known combinatorial kernel
# and g(n) is a gauge.

# Try: q_n = q_0 · Σ C(n,k)^2 C(n+k,n) C(n+2k,n) · r^k / a_n · a_k
# i.e., q_n = q_0 · Σ (a_k F(n,k)/a_n) ... circular

# NEW IDEA: maybe q_n^{P27} = g(n) · C(2n,n) · a_n (i.e., proportional to A_n^{AESZ})
# with g(n) = some Pochhammer/polynomial gauge

# Check: q_n / A_n^{AESZ} = q_n / [C(2n,n) a_n]
print("\nq_n^{P27} / A_n^{AESZ} (= q_n / [C(2n,n) a_n]):")
gauge = []
for n in range(min(15, N)):
    An = binomial(2*n,n) * a_vals[n]
    g = q[n] / An
    gauge.append(g)
    print("  n=%d: %s" % (n, g))

# Check ratio of successive gauge values
print("\ng_{n+1}/g_n:")
for n in range(min(12, len(gauge)-1)):
    if gauge[n] != 0:
        r = gauge[n+1] / gauge[n]
        print("  n=%d: %s = %.6f" % (n, r, float(r)))

# Try to guess a recurrence for the gauge
print("\n=== Guessing recurrence for gauge ===")
try:
    rec_g = guess(gauge[:15], OS, order=2)
    print("Found recurrence for gauge! %s" % rec_g)
except:
    print("No order-2 recurrence found")
    try:
        rec_g = guess(gauge[:15], OS, order=3)
        print("Found order-3 recurrence for gauge! %s" % rec_g)
    except:
        print("No order-3 recurrence found either")

# Try to guess from the gauge values directly
print("\n=== Gauge values (exact) ===")
for n in range(min(8, len(gauge))):
    print("  g_%d = %s" % (n, gauge[n]))
    if n > 0:
        print("    = %s / %s" % (gauge[n].numerator(), gauge[n].denominator()))
