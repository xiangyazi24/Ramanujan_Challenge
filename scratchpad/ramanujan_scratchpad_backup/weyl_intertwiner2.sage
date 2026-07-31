#!/usr/bin/env sage
"""
ODE-level comparison of AESZ and P2.7 generating functions.
Guess ODEs from OGF coefficients, compute GCRD/LCLM at ODE level.
"""
from ore_algebra import *
from ore_algebra import guess

# ===== Sequences =====
N = 100

# AESZ inner sum
inner = []
for nn in range(N):
    inner.append(sum(binomial(nn,k)^2 * binomial(nn+k,nn) * binomial(nn+2*k,nn) for k in range(nn+1)))

# P2.7 sequence
def A_p(nn):
    nn = QQ(nn)
    return QQ(1024)*(2*nn+5)^4*(2*nn+7)^3*(2*nn+9)^3*(946*nn^2+6407*nn+10860)
def B_p(nn):
    nn = QQ(nn)
    return QQ(128)*(2*nn+7)^3*(2*nn+9)^3*(104060*nn^6+1745370*nn^5+12145238*nn^4+44886481*nn^3+92943995*nn^2+102256019*nn+46709052)
def C_p(nn):
    nn = QQ(nn)
    return QQ(16)*(nn+3)^4*(2*nn+9)^3*(3784*nn^5+57792*nn^4+351019*nn^3+1059230*nn^2+1587211*nn+944620)
def D_p(nn):
    nn = QQ(nn)
    return (nn+3)^4*(nn+4)^6*(946*nn^2+4515*nn+5399)

q27 = [QQ(-215040420000), QQ(-167282265043404)/905, QQ(-964185327658080)/6071]
for nn in range(2, N):
    q27.append(B_p(nn)/A_p(nn)*q27[nn] - C_p(nn-1)/A_p(nn-1)*q27[nn-1] + D_p(nn-2)/A_p(nn-2)*q27[nn-2])

# ===== Guess ODEs =====
R.<z> = PolynomialRing(QQ)
D.<Dz> = OreAlgebra(R)

print("Guessing ODE for AESZ inner OGF A(z) = Σ a_n z^n...")
# The ODE should have order 6 (= degree of recurrence) and degree 3 (= order of recurrence)
# Actually: order-3 recurrence with degree-6 coefficients -> ODE of order ≤ 6 with degree ≤ 3
for ord in range(3, 10):
    try:
        L_a = guess(inner[:N], D, order=ord)
        deg = max(c.degree() for c in L_a.list())
        print(f"  order={ord}: degree={deg}")
        break
    except:
        print(f"  order={ord}: no")

print(f"\nGuessing ODE for P2.7 OGF Q(z) = Σ q_n z^n...")
# Expected: order ≤ 18, degree ≤ 3
for ord in range(3, 22):
    try:
        L_q = guess([QQ(x) for x in q27[:N]], D, order=ord)
        deg = max(c.degree() for c in L_q.list())
        print(f"  order={ord}: degree={deg}")
        break
    except:
        if ord % 3 == 0:
            print(f"  order={ord}: no")

# ===== Also try: 64-gauged AESZ =====
v_seq = [inner[nn] / QQ(64)^nn for nn in range(N)]
print(f"\nGuessing ODE for v(z) = Σ (a_n/64^n) z^n...")
for ord in range(3, 10):
    try:
        L_v = guess(v_seq[:N], D, order=ord)
        deg = max(c.degree() for c in L_v.list())
        print(f"  order={ord}: degree={deg}")
        break
    except:
        print(f"  order={ord}: no")

# ===== GCRD at ODE level =====
print(f"\n{'='*60}")
print("GCRD and LCLM at ODE level")
print("="*60)

if 'L_a' in dir() and 'L_q' in dir():
    print(f"\nL_a: order={L_a.order()}, degree={max(c.degree() for c in L_a.list())}")
    print(f"L_q: order={L_q.order()}, degree={max(c.degree() for c in L_q.list())}")

    try:
        G = L_a.gcrd(L_q)
        print(f"\nGCRD(L_a, L_q): order={G.order()}")
        if G.order() > 0:
            print("*** NONTRIVIAL ODE-LEVEL GCRD! ***")
            print(f"  degree={max(c.degree() for c in G.list())}")
    except Exception as e:
        print(f"GCRD error: {e}")

    try:
        M = L_a.lclm(L_q)
        print(f"\nLCLM(L_a, L_q): order={M.order()}")
        expected = L_a.order() + L_q.order()
        if M.order() < expected:
            print(f"*** ORDER DROP by {expected - M.order()}! ***")
    except Exception as e:
        print(f"LCLM error: {e}")

if 'L_v' in dir() and 'L_q' in dir():
    print(f"\n--- 64-gauged AESZ vs P2.7 ---")
    print(f"L_v: order={L_v.order()}, degree={max(c.degree() for c in L_v.list())}")

    try:
        G2 = L_v.gcrd(L_q)
        print(f"\nGCRD(L_v, L_q): order={G2.order()}")
        if G2.order() > 0:
            print("*** NONTRIVIAL! ***")
    except Exception as e:
        print(f"GCRD error: {e}")

    try:
        M2 = L_v.lclm(L_q)
        print(f"LCLM(L_v, L_q): order={M2.order()}")
        expected2 = L_v.order() + L_q.order()
        if M2.order() < expected2:
            drop = expected2 - M2.order()
            print(f"*** ORDER DROP by {drop}! Shared solution space dim = {drop} ***")
    except Exception as e:
        print(f"LCLM error: {e}")

# ===== Singular loci =====
print(f"\n{'='*60}")
print("Singular loci comparison")
print("="*60)

if 'L_a' in dir():
    lc_a = L_a.list()[-1]
    tc_a = L_a.list()[0]
    print(f"L_a leading: {lc_a.factor()}")
    print(f"L_a trailing: {tc_a.factor()}")

if 'L_q' in dir():
    lc_q = L_q.list()[-1]
    tc_q = L_q.list()[0]
    print(f"\nL_q leading: {lc_q.factor()}")
    print(f"L_q trailing: {tc_q.factor()}")

if 'L_v' in dir():
    lc_v = L_v.list()[-1]
    tc_v = L_v.list()[0]
    print(f"\nL_v leading: {lc_v.factor()}")
    print(f"L_v trailing: {tc_v.factor()}")
