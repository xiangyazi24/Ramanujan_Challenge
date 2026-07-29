#!/usr/bin/env sage
"""
GCRD of 64-gauged inner operator with L_p27.
Key idea: b_n/64^n has same Poincaré roots as q_n.
"""
from ore_algebra import *
from ore_algebra import guess

R.<n> = PolynomialRing(QQ)
A.<Sn> = OreAlgebra(R)

# Inner sum (AESZ #209 without C(2n,n))
inner = []
for nn in range(80):
    inner.append(sum(binomial(nn,k)^2 * binomial(nn+k,nn) * binomial(nn+2*k,nn) for k in range(nn+1)))

L_inner = guess(inner, A)
print(f"L_inner: order={L_inner.order()}, degree={max(c.degree() for c in L_inner.list())}")

# 64-gauged sequence: v_n = b_n / 64^n
v_seq = [inner[nn] / QQ(64)^nn for nn in range(80)]
print(f"\nv_n = b_n/64^n: v[0..4] = {[float(v_seq[i]) for i in range(5)]}")
print(f"  Growth v[20]/v[19] = {float(v_seq[20]/v_seq[19]):.10f}  (target: ~0.859)")

L_gauged = guess(v_seq, A)
print(f"\nL_gauged: order={L_gauged.order()}, degree={max(c.degree() for c in L_gauged.list())}")

# Also compute: 64-gauged via operator algebra
# If L = Σ c_i(n) Sn^i annihilates b_n
# Then b_n/64^n is annihilated by L' with L'(n) = Σ c_i(n) · 64^(-i) · Sn^i
# i.e., just rescale the Sn coefficients
coeffs = L_inner.list()
L_gauged_alg = sum(coeffs[i] * QQ(64)^(len(coeffs)-1-i) * Sn^i for i in range(len(coeffs)))
print(f"\nL_gauged (algebraic): {L_gauged_alg}")
print(f"  Matches guessed? {L_gauged == L_gauged_alg or L_gauged_alg.is_zero()}")

# Verify L_gauged annihilates v_n
print("\nVerification:")
ok = True
for nn in range(70):
    val = sum(QQ(L_gauged.list()[i](n=nn)) * v_seq[nn+i] for i in range(L_gauged.order()+1))
    if val != 0:
        ok = False
        if nn < 3:
            print(f"  n={nn}: val = {float(val):.6e}")
print(f"  L_gauged annihilates v_n: {ok}")

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
print(f"\nL_p27: order={L_p27.order()}, degree={max(c.degree() for c in L_p27.list())}")

# GCRD computation
print("\n" + "="*60)
print("GCRD(L_gauged, L_p27)")
print("="*60)
try:
    G = L_gauged.gcrd(L_p27)
    print(f"  Order: {G.order()}")
    if G.order() > 0:
        print(f"  *** NONTRIVIAL GCRD! ***")
        print(f"  Degree: {max(c.degree() for c in G.list())}")
        print(f"  G = {G}")
    else:
        print("  Trivial (order 0)")
except Exception as e:
    print(f"  Error: {e}")

# LCLM
print("\nLCLM(L_gauged, L_p27)")
try:
    M = L_gauged.lclm(L_p27)
    print(f"  Order: {M.order()}")
    print(f"  Expected: {L_gauged.order() + L_p27.order()}")
    if M.order() < L_gauged.order() + L_p27.order():
        drop = L_gauged.order() + L_p27.order() - M.order()
        print(f"  *** ORDER DROP by {drop}! Shared solution space dim {drop} ***")
except Exception as e:
    print(f"  Error: {e}")

# Also try: GCRD of L_inner directly with L_p27
print("\n" + "="*60)
print("GCRD(L_inner, L_p27)")
print("="*60)
try:
    G2 = L_inner.gcrd(L_p27)
    print(f"  Order: {G2.order()}")
except Exception as e:
    print(f"  Error: {e}")

# Right divisibility: L_p27 = Q · L_gauged + R?
print("\nRight division L_p27 / L_gauged:")
try:
    Q_div, R_div = L_p27.quo_rem(L_gauged)
    print(f"  Quotient order: {Q_div.order()}")
    print(f"  Remainder order: {R_div.order()}")
    if R_div.is_zero():
        print("  *** L_gauged RIGHT-DIVIDES L_p27! ***")
except Exception as e:
    print(f"  Error: {e}")

# Left divisibility: L_p27 = L_gauged · Q + R?
print("\nLeft division (experimental):")
try:
    Q_ldiv, R_ldiv = L_p27.quo_rem(L_gauged)
    print(f"  Already computed above.")
except:
    pass

# Exterior square of L_gauged
print("\n" + "="*60)
print("Exterior square Λ²(L_gauged)")
print("="*60)

# The exterior square of an order-3 operator L is another order-3 operator
# whose solutions are y_i(n)·y_j(n+1) - y_j(n)·y_i(n+1) (Casoratians)
# In ore_algebra, this is the symmetric product with respect to the
# difference algebra structure.

# Compute the Casoratian sequences
# For the inner sum, compute b_n, b'_n, b''_n (three independent solutions)
# Actually, we only have ONE explicit solution (the inner sum itself).
# For the exterior square, we can compute it algebraically.

# The exterior square operator can be computed via the formula:
# If L = Σ a_i Sn^i (order 3), the exterior square Λ²L has order 3
# and its Poincaré roots are μ_i·μ_j (products of pairs).

# For our Poincaré roots x₀, x₊, x₋:
# Λ² roots: x₀·x₊, x₀·x₋, x₊·x₋
# Since x₊x₋ = 1/(4x₀) (from Vieta), |x₊x₋| = 1/(4·55) ≈ 0.00455
# And x₀x₊, x₀x₋ have modulus x₀·|x₊| ≈ 55·0.067 ≈ 3.71

print("  Poincaré roots of Λ²:")
print(f"    x₀·x₊ and x₀·x₋ have modulus ≈ {55*0.0674:.4f}")
print(f"    x₊·x₋ = 1/(4·x₀) ≈ {1/(4*54.96):.6f}")

# Try to compute the exterior square using ore_algebra
# ore_algebra may have a symmetric_product method
try:
    print("\n  Computing symmetric_product(L_gauged, L_gauged)...")
    # The exterior square is related to the symmetric product
    # Λ²(L) solutions: {y₁y₂' - y₂y₁'} where y₁,y₂ are solutions of L
    # This equals the 2nd exterior power of the difference module.
    # In ore_algebra, one approach: compute the annihilator of det(y_i(n+j))_{i∈{a,b}, j∈{0,1}}
    # for independent solutions y_a, y_b.

    # Simplest: compute Casoratian c_n = y₁(n)y₂(n+1) - y₁(n+1)y₂(n) numerically
    # and guess its recurrence.

    # We have v_n = b_n/64^n as one solution. We need another.
    # Use the recurrence to generate a second independent solution from different ICs.

    # Second solution: start with v'_0=0, v'_1=0, v'_2=1
    v2 = [QQ(0)] * 80
    v2[0] = QQ(0)
    v2[1] = QQ(0)
    v2[2] = QQ(1)
    coeffs_g = L_gauged.list()
    for nn in range(2, 77):
        # L_gauged: c3 v_{n+3} + c2 v_{n+2} + c1 v_{n+1} + c0 v_n = 0
        c3_val = QQ(coeffs_g[3](n=nn))
        c2_val = QQ(coeffs_g[2](n=nn))
        c1_val = QQ(coeffs_g[1](n=nn))
        c0_val = QQ(coeffs_g[0](n=nn))
        if c3_val == 0:
            print(f"    c3({nn})=0, skipping")
            continue
        v2[nn+3] = -(c2_val * v2[nn+2] + c1_val * v2[nn+1] + c0_val * v2[nn]) / c3_val

    # Casoratian: cas_n = v_n · v2_{n+1} - v_{n+1} · v2_n
    cas = [v_seq[nn] * v2[nn+1] - v_seq[nn+1] * v2[nn] for nn in range(75)]

    print("  Guessing exterior square operator...")
    L_ext = guess(cas, A)
    print(f"  Λ² order: {L_ext.order()}, degree: {max(c.degree() for c in L_ext.list())}")

    # Check if Λ² has right factors of order 1 (hypergeometric solutions)
    print("  Right factors of Λ² (order 1):")
    rf_ext = L_ext.right_factors(1)
    print(f"    Found {len(rf_ext)}")
    for f in rf_ext:
        print(f"    {f}")
        # Extract the hypergeometric ratio
        coeffs_f = f.list()
        if len(coeffs_f) == 2:
            ratio = -coeffs_f[0] / coeffs_f[1]
            print(f"    Hypergeometric ratio u_{'{n+1}'}/u_n = {ratio}")

except Exception as e:
    import traceback
    traceback.print_exc()
