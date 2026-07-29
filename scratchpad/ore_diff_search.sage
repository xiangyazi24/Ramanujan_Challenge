#!/usr/bin/env sage
"""
Ore algebra search for differential intertwiner between AESZ #209 and P2.7.

Strategy:
1. Find the AESZ #209 recurrence operator in Q[n]<Sn>
2. Convert both AESZ and P2.7 to differential operators in Q(x)<Dx>
3. Compute LCLM and GCRD
4. Search for gauge/intertwiner
"""
from ore_algebra import *

# ============================================================
# 1. Set up the recurrence rings
# ============================================================
R_poly.<n> = PolynomialRing(QQ)
R_ore.<Sn> = OreAlgebra(R_poly)

# P2.7 recurrence: A(n) u_{n+1} = B(n) u_n - C(n-1)/A(n-1)·A(n)·u_{n-1} + D(n-2)/A(n-2)·A(n)·u_{n-2}
# In operator form: A(n)·S^3 - B(n)·S^2 + C(n)·S - D(n)   (shifted appropriately)
# Actually: the recurrence is
# u_{n+1} = (B_n/A_n) u_n - (C_{n-1}/A_{n-1}) u_{n-1} + (D_{n-2}/A_{n-2}) u_{n-2}
# Written as: A_n u_{n+1} - B_n u_n + C_{n-1}/A_{n-1}·A_n u_{n-1} - D_{n-2}/A_{n-2}·A_n u_{n-2} = 0
# To clear all denominators, multiply through to get polynomial coefficients.

# Actually, let's state the recurrence in standard form:
# P3(n) u_{n+3} + P2(n) u_{n+2} + P1(n) u_{n+1} + P0(n) u_n = 0
# where the shift is from the original recurrence indexed as:
# A_n u_{n+1} - B_n u_n + (C_{n-1}/A_{n-1}) A_n u_{n-1} - (D_{n-2}/A_{n-2}) A_n u_{n-2} = 0

# Let me use the compact form. First define the polynomial coefficients:
def A_poly(n):
    return 1024*(2*n+5)^4*(2*n+7)^3*(2*n+9)^3*(946*n^2+6407*n+10860)

def B_poly(n):
    return 128*(2*n+7)^3*(2*n+9)^3*(104060*n^6+1745370*n^5+12145238*n^4+44886481*n^3+92943995*n^2+102256019*n+46709052)

def C_poly(n):
    return 16*(n+3)^4*(2*n+9)^3*(3784*n^5+57792*n^4+351019*n^3+1059230*n^2+1587211*n+944620)

def D_poly(n):
    return (n+3)^4*(n+4)^6*(946*n^2+4515*n+5399)

# The recurrence u_{n+1} = (B/A)u_n - (C_{n-1}/A_{n-1})u_{n-1} + (D_{n-2}/A_{n-2})u_{n-2}
# Shift indices: let's write in terms of S_n acting on u_n.
# A_n·Sn·u_n - B_n·u_n + ... = 0
# In OreAlgebra, Sn shifts n by 1: Sn·f(n) = f(n+1)

# Clear denominators. The recurrence is:
# A(n) u(n+1) - B(n) u(n) + [C(n-1)/A(n-1)] · A(n) · u(n-1) - [D(n-2)/A(n-2)] · A(n) · u(n-2) = 0

# To get a polynomial operator, I need to multiply through by A(n-1)·A(n-2).
# Full clearing: A(n)·A(n-1)·A(n-2) u(n+1) - B(n)·A(n-1)·A(n-2) u(n)
#   + C(n-1)·A(n-2)·A(n) u(n-1) - D(n-2)·A(n-1)·A(n) u(n-2) = 0

# But this gives degree 12*3+12 = 48 polynomial coefficients, which is huge.
# Better: use the FACTORED form directly.

# Let me just use ore_algebra's guess function to find the operator from data.
print("Computing AESZ #209 sequence...")
aesz = []
for nn in range(50):
    aesz.append(sum(binomial(nn,k)^2 * binomial(nn+k,nn) * binomial(nn+2*k,nn) for k in range(nn+1)))

print(f"  a[0..4] = {aesz[:5]}")
print(f"  a[5..9] = {aesz[5:10]}")

# Guess the recurrence
print("\nGuessing AESZ #209 recurrence...")
from ore_algebra import guess
L_aesz = guess(aesz, R_ore)
print(f"  Order: {L_aesz.order()}")
print(f"  Degree: {max(c.degree() for c in L_aesz.list())}")

# Print the operator
print(f"\n  L_aesz = {L_aesz}")

# P2.7 sequence
print("\nComputing P2.7 sequence (exact rationals)...")
q = [QQ(0)] * 50
q[0] = QQ(-215040420000)
q[1] = QQ(-167282265043404) / QQ(905)
q[2] = QQ(-964185327658080) / QQ(6071)

A = lambda nn: R_poly(A_poly(nn))(nn)
B = lambda nn: R_poly(B_poly(nn))(nn)
C = lambda nn: R_poly(C_poly(nn))(nn)
D = lambda nn: R_poly(D_poly(nn))(nn)

for nn in range(2, 49):
    An = QQ(1024*(2*nn+5)^4*(2*nn+7)^3*(2*nn+9)^3*(946*nn^2+6407*nn+10860))
    Bn = QQ(128*(2*nn+7)^3*(2*nn+9)^3*(104060*nn^6+1745370*nn^5+12145238*nn^4+44886481*nn^3+92943995*nn^2+102256019*nn+46709052))
    An1 = QQ(1024*(2*(nn-1)+5)^4*(2*(nn-1)+7)^3*(2*(nn-1)+9)^3*(946*(nn-1)^2+6407*(nn-1)+10860))
    Cn1 = QQ(16*((nn-1)+3)^4*(2*(nn-1)+9)^3*(3784*(nn-1)^5+57792*(nn-1)^4+351019*(nn-1)^3+1059230*(nn-1)^2+1587211*(nn-1)+944620))
    An2 = QQ(1024*(2*(nn-2)+5)^4*(2*(nn-2)+7)^3*(2*(nn-2)+9)^3*(946*(nn-2)^2+6407*(nn-2)+10860))
    Dn2 = QQ(((nn-2)+3)^4*((nn-2)+4)^6*(946*(nn-2)^2+4515*(nn-2)+5399))
    q[nn+1] = Bn/An * q[nn] - Cn1/An1 * q[nn-1] + Dn2/An2 * q[nn-2]

print(f"  q[0..2] = {q[:3]}")
print(f"  q[3] = {q[3]}")

# Guess the P2.7 recurrence
print("\nGuessing P2.7 recurrence...")
L_p27 = guess(q[:40], R_ore)
print(f"  Order: {L_p27.order()}")
print(f"  Degree: {max(c.degree() for c in L_p27.list())}")

# ============================================================
# 2. Compare the two operators
# ============================================================
print("\n" + "="*60)
print("Comparing operators")
print("="*60)

# GCRD
print("\nComputing GCRD(L_aesz, L_p27)...")
try:
    G = L_aesz.gcrd(L_p27)
    print(f"  GCRD order: {G.order()}")
    if G.order() > 0:
        print(f"  GCRD: {G}")
    else:
        print("  GCRD is trivial (order 0) — operators share no common right factor")
except Exception as e:
    print(f"  Error: {e}")

# LCLM
print("\nComputing LCLM(L_aesz, L_p27)...")
try:
    M = L_aesz.lclm(L_p27)
    print(f"  LCLM order: {M.order()}")
    print(f"  Expected order (if independent): {L_aesz.order() + L_p27.order()}")
except Exception as e:
    print(f"  Error: {e}")

# ============================================================
# 3. Factor analysis
# ============================================================
print("\n" + "="*60)
print("Factor analysis")
print("="*60)

# Right factors of L_p27
print("\nRight factors of L_p27 (order 1):")
try:
    rf = L_p27.right_factors(1)
    if rf:
        for f in rf:
            print(f"  {f}")
    else:
        print("  None (irreducible at order 1)")
except Exception as e:
    print(f"  {e}")

# Left factors of L_p27
print("\nLeft factors of L_p27 (order 1):")
try:
    lf = L_p27.left_factors(1)
    if lf:
        for f in lf:
            print(f"  {f}")
    else:
        print("  None")
except Exception as e:
    print(f"  {e}")

# Check if L_p27 = Q · L_aesz · R for some operators Q, R
# This would mean every solution of L_aesz, after applying R, solves L_p27
# Equivalent to: L_aesz divides L_p27 from the right (after R gauge)
print("\nDoes L_aesz right-divide L_p27?")
try:
    quotient, remainder = L_p27.quo_rem(L_aesz)
    print(f"  Quotient order: {quotient.order()}")
    print(f"  Remainder order: {remainder.order()}")
    if remainder.is_zero():
        print("  YES — L_p27 = Q · L_aesz")
    else:
        print("  NO")
except Exception as e:
    print(f"  Error: {e}")

# ============================================================
# 4. Convert to differential operators
# ============================================================
print("\n" + "="*60)
print("Differential operator ring")
print("="*60)

R_diff.<x> = PolynomialRing(QQ)
D_ore.<Dx> = OreAlgebra(R_diff)

# Convert recurrence operator to differential operator
# The GF A(x) = Σ a_n x^n satisfies L_diff · A = 0
# where L_diff comes from the recurrence via n -> x·Dx (theta = x Dx)
print("\nConverting L_aesz to differential operator...")
try:
    L_aesz_diff = L_aesz.to_D(D_ore)
    print(f"  Order: {L_aesz_diff.order()}")
    print(f"  Max coeff degree: {max(c.degree() for c in L_aesz_diff.list())}")
except Exception as e:
    print(f"  Error: {e}")

print("\nConverting L_p27 to differential operator...")
try:
    L_p27_diff = L_p27.to_D(D_ore)
    print(f"  Order: {L_p27_diff.order()}")
    print(f"  Max coeff degree: {max(c.degree() for c in L_p27_diff.list())}")
except Exception as e:
    print(f"  Error: {e}")

# Singular points
print("\nSingular points of L_aesz (differential):")
try:
    lc = L_aesz_diff.list()[-1]
    print(f"  Leading coefficient: {lc.factor()}")
except Exception as e:
    print(f"  {e}")

print("\nSingular points of L_p27 (differential):")
try:
    lc = L_p27_diff.list()[-1]
    print(f"  Leading coefficient: {lc.factor()}")
except Exception as e:
    print(f"  {e}")

# ============================================================
# 5. Intertwiner search in differential ring
# ============================================================
print("\n" + "="*60)
print("Differential intertwiner search")
print("="*60)

# Search for T ∈ Q(x)<Dx> such that L_p27 · T ≡ 0 mod L_aesz
# i.e., T maps solutions of L_aesz to solutions of L_p27
# This means L_p27 divides L_aesz · T from the left
# Or equivalently: L_aesz · T ∈ L_p27 · Q(x)<Dx>

# A simpler check: does the LCLM have special structure?
print("\nDifferential LCLM order:")
try:
    M_diff = L_aesz_diff.lclm(L_p27_diff)
    print(f"  LCLM order: {M_diff.order()}")
    print(f"  Sum of orders: {L_aesz_diff.order() + L_p27_diff.order()}")
    if M_diff.order() < L_aesz_diff.order() + L_p27_diff.order():
        print("  *** LCLM order DROP — operators share solution space! ***")
    else:
        print("  No order drop — independent solution spaces")
except Exception as e:
    print(f"  Error: {e}")

print("\nDifferential GCRD:")
try:
    G_diff = L_aesz_diff.gcrd(L_p27_diff)
    print(f"  GCRD order: {G_diff.order()}")
except Exception as e:
    print(f"  Error: {e}")
