#!/usr/bin/env sage
"""
Find the recurrence satisfied by the Cauchy-twisted Zudilin basis,
then search for Ore intertwiner with P2.7.
"""
from ore_algebra import *
from ore_algebra import guess

R.<n> = PolynomialRing(QQ)
A.<Sn> = OreAlgebra(R)

N = 100

# AESZ inner sum
inner = []
for nn in range(N):
    inner.append(sum(binomial(nn,k)^2 * binomial(nn+k,nn) * binomial(nn+2*k,nn) for k in range(nn+1)))

# Zudilin b2 and b3 companions
# Use Zudilin recurrence to compute them
def z3(nn): return ZZ(2*(946*nn^2-731*nn+153)*(2*nn+1)*(nn+1)^3)
def z2_c(nn): return ZZ(-2*(104060*nn^6+127710*nn^5+12788*nn^4-34525*nn^3-8482*nn^2+3298*nn+1071))
def z1(nn): return ZZ(2*nn*(3784*nn^5-1032*nn^4-1925*nn^3+853*nn^2+328*nn-184))
def z0(nn): return ZZ(-(946*nn^2+1161*nn+368)*nn*(nn-1)^3)

def zudilin_rec(init):
    u = list(init)
    for nn in range(2, N):
        u.append((-z2_c(nn)*u[nn] - z1(nn)*u[nn-1] - z0(nn)*u[nn-2]) / z3(nn))
    return u

a = zudilin_rec([QQ(1), QQ(7), QQ(163)])
b2 = zudilin_rec([QQ(0), QQ(23)/2, QQ(2145)/8])
b3 = zudilin_rec([QQ(0), QQ(17)/2, QQ(3135)/16])

# Verify a matches inner sum
for nn in range(min(N, 30)):
    assert a[nn] == inner[nn], f"Mismatch at n={nn}"
print("Zudilin a matches inner sum. ✓")

# c_j(1) for the Cauchy twist
c = [QQ(1)]
for nn in range(N+5):
    cm1 = c[nn-1] if nn >= 1 else QQ(0)
    cm2 = c[nn-2] if nn >= 2 else QQ(0)
    c.append(((440*nn+660)*c[nn] + (160-16*nn)*cm1 + (2*nn-43)*cm2) / (8*(nn+1)))

print(f"c[0:5] = {c[:5]}")
print(f"c[5] = {c[5]}")

# Cauchy twist: (C₁ f)_n = (1/64^n) Σ_{m=0}^n c_{n-m} f_m
def cauchy_twist(seq):
    return [sum(c[nn-m]*seq[m] for m in range(nn+1)) / QQ(64)^nn
            for nn in range(N)]

ta = cauchy_twist(a)
tb2 = cauchy_twist(b2)
tb3 = cauchy_twist(b3)

print(f"\nTwisted a: ta[0:5] = {[float(x) for x in ta[:5]]}")
print(f"Growth ta[20]/ta[19] = {float(ta[20]/ta[19]):.10f}")

# Guess the recurrence for the twisted sequence ta
print(f"\nGuessing recurrence for Cauchy-twisted a_n...")
L_twist = guess(ta, A)
print(f"L_twist: order={L_twist.order()}, degree={max(c.degree() for c in L_twist.list())}")

# P2.7 operator
A_fn = R(1024*(2*n+5)^4*(2*n+7)^3*(2*n+9)^3*(946*n^2+6407*n+10860))
B_fn = R(128*(2*n+7)^3*(2*n+9)^3*(104060*n^6+1745370*n^5+12145238*n^4+44886481*n^3+92943995*n^2+102256019*n+46709052))
C_fn = R(16*(n+3)^4*(2*n+9)^3*(3784*n^5+57792*n^4+351019*n^3+1059230*n^2+1587211*n+944620))
D_fn = R((n+3)^4*(n+4)^6*(946*n^2+4515*n+5399))

A0, A1, A2 = A_fn(n=n), A_fn(n=n+1), A_fn(n=n+2)
B2 = B_fn(n=n+2)
C1 = C_fn(n=n+1)
D0 = D_fn(n=n)

c3_p = A0 * A1 * A2
c2_p = -A0 * A1 * B2
c1_p = A0 * A2 * C1
c0_p = -A1 * A2 * D0
g = gcd([c3_p, c2_p, c1_p, c0_p])
L_p27 = (c3_p // g) * Sn^3 + (c2_p // g) * Sn^2 + (c1_p // g) * Sn + (c0_p // g)
print(f"L_p27: order={L_p27.order()}, degree={max(c.degree() for c in L_p27.list())}")

# GCRD of twisted recurrence with P2.7
print(f"\n{'='*60}")
print("GCRD(L_twist, L_p27)")
print("="*60)
try:
    G = L_twist.gcrd(L_p27)
    print(f"  Order: {G.order()}")
    if G.order() > 0:
        print(f"  *** NONTRIVIAL GCRD! ***")
        print(f"  Degree: {max(c.degree() for c in G.list())}")
    else:
        print("  Trivial (order 0)")
except Exception as e:
    print(f"  Error: {e}")

# LCLM
print(f"\nLCLM(L_twist, L_p27)")
try:
    M_lclm = L_twist.lclm(L_p27)
    print(f"  Order: {M_lclm.order()}")
    print(f"  Expected max: {L_twist.order() + L_p27.order()}")
    if M_lclm.order() < L_twist.order() + L_p27.order():
        drop = L_twist.order() + L_p27.order() - M_lclm.order()
        print(f"  *** ORDER DROP by {drop}! ***")
except Exception as e:
    print(f"  Error: {e}")

# Also guess recurrence for tb2 and tb3
print(f"\nGuessing recurrence for twisted b2...")
try:
    L_tb2 = guess(tb2, A)
    print(f"  order={L_tb2.order()}, degree={max(c.degree() for c in L_tb2.list())}")
    if L_tb2 == L_twist:
        print("  *** Same as L_twist! ***")
except Exception as e:
    print(f"  Error: {e}")

print(f"\nGuessing recurrence for twisted b3...")
try:
    L_tb3 = guess(tb3, A)
    print(f"  order={L_tb3.order()}, degree={max(c.degree() for c in L_tb3.list())}")
    if L_tb3 == L_twist:
        print("  *** Same as L_twist! ***")
except Exception as e:
    print(f"  Error: {e}")

# Poincaré analysis of L_twist
print(f"\nPoincaré analysis of L_twist:")
leads = [c.leading_coefficient() for c in L_twist.list()]
print(f"  Leading coefficients: {leads}")
poincare = sum(leads[i] * var('mu')^i for i in range(len(leads)))
print(f"  Poincaré polynomial: {poincare}")
from sage.all import solve as solve_eq
mu = var('mu')
roots = solve_eq(poincare, mu)
for r in roots:
    try:
        val = complex(r.rhs().n())
        print(f"    Root: {val} (modulus {abs(val):.10f})")
    except:
        print(f"    Root: {r}")

# Birkhoff formal powers
print(f"\nFormal powers of L_twist vs L_p27:")
for op, name in [(L_twist, 'L_twist'), (L_p27, 'L_p27')]:
    coeffs = op.list()
    ord = op.order()
    leading = coeffs[ord]
    trailing = coeffs[0]
    # Formal power at dominant root: σ = -deg(leading) + ...
    # Actually: σ is determined by the next-order expansion
    # For polynomial-coefficient recurrences, the formal power is related to
    # the factored form of leading/trailing
    print(f"  {name}:")
    print(f"    Leading coeff degree: {leading.degree()}")
    print(f"    Trailing coeff degree: {trailing.degree()}")
    try:
        print(f"    Leading factored: {leading.factor()}")
    except:
        pass
