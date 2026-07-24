#!/usr/bin/env sage
"""
Ore algebra analysis — construct L_p27 directly, then compare with L_aesz.
"""
from ore_algebra import *
from ore_algebra import guess

R_poly.<n> = PolynomialRing(QQ)
R_ore.<Sn> = OreAlgebra(R_poly)

# ============================================================
# 1. Guess AESZ recurrence
# ============================================================
print("AESZ #209...")
aesz = []
for nn in range(60):
    aesz.append(binomial(2*nn,nn) * sum(binomial(nn,k)^2 * binomial(nn+k,nn) * binomial(nn+2*k,nn) for k in range(nn+1)))
L_aesz = guess(aesz, R_ore)
print(f"  Order: {L_aesz.order()}, Degree: {max(c.degree() for c in L_aesz.list())}")

# ============================================================
# 2. Construct L_p27 directly
# ============================================================
# The recurrence at index m (m >= 2):
#   A_m u_{m+1} = B_m u_m - (C_{m-1}/A_{m-1}) u_{m-1} + (D_{m-2}/A_{m-2}) u_{m-2}
# Shift m = n+2 to get operator on u_n, u_{n+1}, u_{n+2}, u_{n+3}:
#   A_{n+2} · A_n · A_{n+1} · u_{n+3}
#   - B_{n+2} · A_n · A_{n+1} · u_{n+2}
#   + C_{n+1} · A_n · u_{n+1}
#   - D_n · A_{n+1} · u_n = 0

print("\nConstructing L_p27...")

A_fn = R_poly(1024*(2*n+5)^4*(2*n+7)^3*(2*n+9)^3*(946*n^2+6407*n+10860))
B_fn = R_poly(128*(2*n+7)^3*(2*n+9)^3*(104060*n^6+1745370*n^5+12145238*n^4+44886481*n^3+92943995*n^2+102256019*n+46709052))
C_fn = R_poly(16*(n+3)^4*(2*n+9)^3*(3784*n^5+57792*n^4+351019*n^3+1059230*n^2+1587211*n+944620))
D_fn = R_poly((n+3)^4*(n+4)^6*(946*n^2+4515*n+5399))

# Shift: A(n+k) etc.
A0 = A_fn(n=n)        # A_n
A1 = A_fn(n=n+1)      # A_{n+1}
A2 = A_fn(n=n+2)      # A_{n+2}
B2 = B_fn(n=n+2)      # B_{n+2}
C1 = C_fn(n=n+1)      # C_{n+1}
D0 = D_fn(n=n)        # D_n

c3 = A0 * A1 * A2
c2 = -A0 * A1 * B2
c1 = A0 * A2 * C1
c0 = -A1 * A2 * D0

print(f"  Degree of c3: {c3.degree()}")
print(f"  Degree of c2: {c2.degree()}")
print(f"  Degree of c1: {c1.degree()}")
print(f"  Degree of c0: {c0.degree()}")

L_p27_raw = c3 * Sn^3 + c2 * Sn^2 + c1 * Sn + c0

# Verify against sequence
print("\n  Verifying L_p27 against q_n...")
q = [QQ(0)] * 100
q[0] = QQ(-215040420000)
q[1] = QQ(-167282265043404) / QQ(905)
q[2] = QQ(-964185327658080) / QQ(6071)

for nn in range(2, 99):
    An = QQ(A_fn(n=nn))
    Bn = QQ(B_fn(n=nn))
    Am1 = QQ(A_fn(n=nn-1))
    Cm1 = QQ(C_fn(n=nn-1))
    Am2 = QQ(A_fn(n=nn-2))
    Dm2 = QQ(D_fn(n=nn-2))
    q[nn+1] = Bn/An * q[nn] - Cm1/Am1 * q[nn-1] + Dm2/Am2 * q[nn-2]

# Check L_p27_raw annihilates q
print("  Checking annihilation...")
max_err = 0
for nn in range(90):
    val = QQ(c3(n=nn)) * q[nn+3] + QQ(c2(n=nn)) * q[nn+2] + QQ(c1(n=nn)) * q[nn+1] + QQ(c0(n=nn)) * q[nn]
    if val != 0:
        max_err += 1
print(f"  Non-zero residuals: {max_err} (should be 0)")

# Try to simplify: extract GCD of c0,c1,c2,c3
print("\n  Extracting common factor...")
g = gcd([c3, c2, c1, c0])
print(f"  GCD degree: {g.degree()}")
if g.degree() > 0:
    print(f"  GCD = {g.factor()}")
    c3s = c3 // g
    c2s = c2 // g
    c1s = c1 // g
    c0s = c0 // g
    L_p27 = c3s * Sn^3 + c2s * Sn^2 + c1s * Sn + c0s
    print(f"  Simplified degrees: {c3s.degree()}, {c2s.degree()}, {c1s.degree()}, {c0s.degree()}")
else:
    L_p27 = L_p27_raw
    print("  No common polynomial factor")

# ============================================================
# 3. GCRD of L_aesz and L_p27
# ============================================================
print("\n" + "="*60)
print("GCRD(L_aesz, L_p27)")
print("="*60)
try:
    G = L_aesz.gcrd(L_p27)
    print(f"  Order: {G.order()}")
    if G.order() > 0:
        print(f"  Degree: {max(c.degree() for c in G.list())}")
        print(f"  G = {G}")
    else:
        print("  Trivial GCRD — no common right factor")
except Exception as e:
    print(f"  Error: {e}")

# ============================================================
# 4. LCLM of L_aesz and L_p27
# ============================================================
print("\n" + "="*60)
print("LCLM(L_aesz, L_p27)")
print("="*60)
try:
    M = L_aesz.lclm(L_p27)
    print(f"  Order: {M.order()}")
    print(f"  Expected (if independent): {L_aesz.order() + L_p27.order()}")
    if M.order() < L_aesz.order() + L_p27.order():
        print("  *** ORDER DROP — shared solution space! ***")
        delta = L_aesz.order() + L_p27.order() - M.order()
        print(f"  Dimension of intersection: {delta}")
    else:
        print("  Independent solution spaces")
except Exception as e:
    print(f"  Error: {e}")

# ============================================================
# 5. Direct divisibility check: does L_aesz divide L_p27?
# ============================================================
print("\n" + "="*60)
print("Divisibility")
print("="*60)
try:
    Q, R = L_p27.quo_rem(L_aesz)
    print(f"  L_p27 = Q · L_aesz + R")
    print(f"  Q order: {Q.order()}, R order: {R.order()}")
    if R.is_zero():
        print("  *** L_aesz RIGHT-DIVIDES L_p27 ***")
        print(f"  Q = {Q}")
except Exception as e:
    print(f"  Error: {e}")

# ============================================================
# 6. Convert to differential operators
# ============================================================
print("\n" + "="*60)
print("Differential operator conversion")
print("="*60)

R_diff.<x> = PolynomialRing(QQ)
D_ore.<Dx> = OreAlgebra(R_diff)

print("Converting L_aesz to ODE...")
try:
    L_aesz_ode = L_aesz.to_D(D_ore)
    print(f"  Order: {L_aesz_ode.order()}")
    lc = L_aesz_ode.list()[-1]
    print(f"  Leading coeff degree: {lc.degree()}")
    print(f"  Leading coeff factored: {lc.factor()}")
except Exception as e:
    print(f"  Error: {e}")

# The ODE for L_p27 is much larger (degree 36+ polynomial coefficients).
# Skip for now and focus on the recurrence-level analysis.

# ============================================================
# 7. Factorization of L_p27
# ============================================================
print("\n" + "="*60)
print("Factorization of L_p27")
print("="*60)

print("Right factors of order 1...")
try:
    rf = L_p27.right_factors(1)
    print(f"  Found {len(rf)} right factors")
    for f in rf:
        print(f"  {f}")
except Exception as e:
    print(f"  {str(e)[:100]}")

print("\nRight factors of order 2...")
try:
    rf2 = L_p27.right_factors(2)
    print(f"  Found {len(rf2)} right factors")
except Exception as e:
    print(f"  {str(e)[:100]}")

# ============================================================
# 8. Check: is L_p27 a gauge (scalar multiplication) of L_aesz?
# ============================================================
print("\n" + "="*60)
print("Gauge search: L_p27 = gauge(L_aesz, r(n))?")
print("="*60)
# If u_n satisfies L_aesz, does v_n = r(n) · u_n satisfy L_p27?
# This means: L_p27(r(n) · u_n) = 0 for all u_n satisfying L_aesz.
# Equivalent to: L_p27 ∘ r(n) = Q · L_aesz for some Q.
# Where L_p27 ∘ r means: substitute u_n → r(n) u_n in L_p27.

# A hypergeometric gauge r(n+1)/r(n) = g(n) rational transforms:
# If v_n = r(n) u_n, then Sn · v_n = r(n+1) u_{n+1} = g(n) · r(n) · u_{n+1}
# The gauged operator is: L_{gauged} = r(n)^{-1} · L · r(n)

# For this, if r(n) = h_n (Pochhammer gauge), compute L_aesz_gauged:
print("  Computing Pochhammer-gauged L_aesz...")
# h_{n+1}/h_n = delta(n) = (n+3)^4(n+4)^6 / [1024(2n+5)^4(2n+7)^3(2n+9)^3]
delta_fn = R_poly((n+3)^4*(n+4)^6) / R_poly(1024*(2*n+5)^4*(2*n+7)^3*(2*n+9)^3)

# Gauged operator: if L = a3 Sn^3 + a2 Sn^2 + a1 Sn + a0
# then L_gauged = a3 δ(n)δ(n+1)δ(n+2) Sn^3 + a2 δ(n)δ(n+1) Sn^2 + a1 δ(n) Sn + a0
# Multiply through by denominators to get polynomial.

# Actually this is getting complicated. Let me check the NUMERICAL ratio instead.
# Compute v_n = h_n · a_n and check if L_p27 annihilates v_n.

from sage.all import rising_factorial
h = [QQ(0)] * 50
for nn in range(50):
    h[nn] = QQ(2)^(-20*nn)
    for j in range(nn):
        h[nn] *= QQ(3+j)^4 * QQ(4+j)^6
        h[nn] /= QQ(QQ(5)/2+j)^4 * QQ(QQ(7)/2+j)^3 * QQ(QQ(9)/2+j)^3

v = [h[nn] * aesz[nn] for nn in range(50)]
print(f"\n  v = h · a first values: {[float(v[i]) for i in range(5)]}")

# Check if L_p27 annihilates v
print("\n  Checking L_p27 · v = 0:")
annihilates = True
for nn in range(40):
    val = QQ(c3(n=nn)) * v[nn+3] + QQ(c2(n=nn)) * v[nn+2] + QQ(c1(n=nn)) * v[nn+1] + QQ(c0(n=nn)) * v[nn]
    if val != 0:
        if nn < 5:
            print(f"  n={nn}: val = {float(val):.6e} (NOT zero)")
        annihilates = False

if annihilates:
    print("  *** h_n · a_n satisfies L_p27! ***")
else:
    print("  h_n · a_n does NOT satisfy L_p27.")

# Try v_n = (64)^{-n} · a_n
v2 = [aesz[nn] / QQ(64)^nn for nn in range(50)]
print("\n  Checking L_p27 · (a_n/64^n) = 0:")
annihilates2 = True
for nn in range(40):
    val = QQ(c3(n=nn)) * v2[nn+3] + QQ(c2(n=nn)) * v2[nn+2] + QQ(c1(n=nn)) * v2[nn+1] + QQ(c0(n=nn)) * v2[nn]
    if val != 0:
        if nn < 3:
            print(f"  n={nn}: val = {float(val):.6e}")
        annihilates2 = False

if annihilates2:
    print("  *** a_n/64^n satisfies L_p27! ***")
else:
    print("  a_n/64^n does NOT satisfy L_p27.")

# Try v_n = h_n · a_n / 64^n
# Actually, let's just guess the operator for v = h·a
print("\n  Guessing operator for h_n · a_n...")
try:
    L_ha = guess([v[nn] for nn in range(45)], R_ore)
    print(f"  Order: {L_ha.order()}, Degree: {max(c.degree() for c in L_ha.list())}")

    # Check if this equals L_p27
    Q2, R2 = L_p27.quo_rem(L_ha)
    if R2.is_zero():
        print(f"  *** L_p27 = Q · L_ha! Q has order {Q2.order()} ***")
    else:
        print(f"  L_ha does not right-divide L_p27")

    G2 = L_ha.gcrd(L_p27)
    print(f"  GCRD(L_ha, L_p27) order: {G2.order()}")
except Exception as e:
    print(f"  {str(e)[:100]}")
