#!/usr/bin/env python3
"""
P2.5: Compute the Birkhoff dominant left eigenvector V+ at n=0.
Then check if p·V+ / (q·V+) = G (Catalan's constant).

Method:
1. Compute M_H(n) for large N
2. Forward iterate both q-row and p-row through N steps
3. The dominant projection is extracted by:
   α+ = lim (q·M_prod) / (e₁·M_prod) (component-wise, any component)
   β+ = lim (p·M_prod) / (e₁·M_prod)
   Then L = β+ / α+
4. Numerically verify L = G to 200+ digits
5. Try PSLQ on the individual projections to find closed forms
"""
from mpmath import mp, mpf, matrix as mpmat, catalan, pi, log, sqrt, euler
from mpmath import identify, pslq

mp.dps = 300

def M_entries(n):
    m11 = (-2*n-5)*(n+3)**2*(136*n**4+1424*n**3+5548*n**2+9551*n+6141)
    m12 = 384*n**6+6384*n**5+44168*n**4+162698*n**3+336377*n**2+369933*n+169011
    m13 = -(480*n**4+4980*n**3+19210*n**2+32690*n+20730)
    m21 = (n+2)**2*(n+3)**2*(4*n+10)*(48*n**3+386*n**2+1017*n+879)
    m22 = (n+2)**2*(-272*n**5-3848*n**4-21732*n**3-61184*n**2-85761*n-47808)
    m23 = (n+2)**2*(320*n**3+2540*n**2+6610*n+5640)
    m31 = (-4*n-10)*(n+2)**2*(n+3)**2*(32*n**4+302*n**3+1037*n**2+1530*n+813)
    m32 = (n+2)**2*(192*n**6+2984*n**5+19116*n**4+64452*n**3+120256*n**2+117279*n+46476)
    m33 = (n+2)**2*(-16*n**5-408*n**4-2912*n**3-8884*n**2-12254*n-6240)
    return [[m11,m12,m13],[m21,m22,m23],[m31,m32,m33]]

def delta_H(n):
    return -2*(n+2)**2*(n+3)**2*(2*n+5)*(2*n+7)**2

G = catalan

print("=== Computing Birkhoff eigenvector for P2.5 ===")
print(f"G = {mp.nstr(G, 50)}")

NMAX = 200

# Forward iterate q-row and p-row through the NORMALIZED CMF M_H(n)
q = [mpf(33750), mpf(-36000), mpf(9000)]
p = [mpf(30921), mpf(-32972), mpf(8240)]

# We need to track the evolution of both rows under M_H(n)
# q_N = q · M_H(0) · M_H(1) · ... · M_H(N-1)
# p_N = p · M_H(0) · M_H(1) · ... · M_H(N-1)

# For computing L = lim P/Q, we just need the ratio at each step
# P_{N,j}/Q_{N,j} = (p · M_prod · e_j) / (q · M_prod · e_j)

# But for the Birkhoff eigenvector, we need to project onto the dominant mode.
# The dominant mode grows like ξ_+^N where ξ_+ = 17+12√2.
# The ratio P/Q converges because both numerator and denominator are dominated
# by the same eigenvector V_+.

# Direct approach: just compute P_{N,1}/Q_{N,1} to 300 digits
q_row = list(q)
p_row = list(p)

print(f"\nComputing {NMAX} CMF steps at {mp.dps}-digit precision...")

ratios = []
for N in range(NMAX):
    M = M_entries(mpf(N))
    d = mpf(delta_H(N))
    MH = [[M[i][j]/d for j in range(3)] for i in range(3)]

    q_new = [sum(q_row[i]*MH[i][j] for i in range(3)) for j in range(3)]
    p_new = [sum(p_row[i]*MH[i][j] for i in range(3)) for j in range(3)]

    q_row = q_new
    p_row = p_new

    # Ratio for column 1 (index 0)
    if q_row[0] != 0:
        ratio = p_row[0] / q_row[0]
        ratios.append(ratio)
        if N % 20 == 19:
            err = ratio - G
            if err != 0:
                digits = -mp.log10(abs(err))
            else:
                digits = mp.dps
            print(f"  N={N+1}: P/Q - G has {mp.nstr(digits, 6)} correct digits")

print(f"\n=== Final ratio analysis ===")
final_ratio = ratios[-1]
err = final_ratio - G
digits_match = -mp.log10(abs(err))
print(f"At N={NMAX}: P/Q matches G to {mp.nstr(digits_match, 8)} digits")

# Now try to understand the INDIVIDUAL dominant projections
# Instead of forward iteration (which mixes all modes),
# use the RATIO of consecutive values to extract the dominant eigenvalue
print(f"\n=== Poincaré analysis of normalized sequence ===")
xi_plus = 17 + 12*sqrt(2)
xi_minus = 17 - 12*sqrt(2)
rho = xi_minus
print(f"ξ₊ = {mp.nstr(xi_plus, 15)} = 17 + 12√2")
print(f"ξ₋ = {mp.nstr(xi_minus, 15)} = 17 - 12√2 = ρ")

# The normalized Q̂_N = q · Π M_H(m) · e₁
# Let's track this and check the Poincaré ratio Q̂_{N+1}/Q̂_N
# This should → ξ_+ as N → ∞ (dominant root)
print(f"\nRecomputing Q̂_N from scratch for Poincaré analysis...")
q_row2 = [mpf(33750), mpf(-36000), mpf(9000)]
qhat_vals = [q_row2[0]]  # Q̂_0 = first component of q

for N in range(NMAX):
    M = M_entries(mpf(N))
    d = mpf(delta_H(N))
    MH = [[M[i][j]/d for j in range(3)] for i in range(3)]
    q_new = [sum(q_row2[i]*MH[i][j] for i in range(3)) for j in range(3)]
    q_row2 = q_new
    qhat_vals.append(q_row2[0])

# Poincaré ratios
print("Poincaré ratios Q̂_{N+1}/Q̂_N:")
for N in [10, 20, 50, 100, 150, 190]:
    if N < len(qhat_vals) - 1:
        r = qhat_vals[N+1] / qhat_vals[N]
        err = r - xi_plus
        if err != 0:
            print(f"  N={N}: ratio = {mp.nstr(r, 20)}, error from ξ₊: {mp.nstr(err, 10)}")

# Now the key computation: extract the dominant LEFT eigenvector
# Using the matrix power method on M_H^T
# Actually, a simpler approach: the ratio P_{N,j}/Q_{N,j} already gives L.
# But to get individual projections, we need different initial vectors.

# Let's try: compute the 3 solutions of the scalar recurrence starting
# from different initial conditions, and decompose Q̂ and P̂ in terms of them.
print(f"\n=== Scalar recurrence decomposition ===")
print("Computing 3 fundamental solutions of the scalar recurrence...")

# We already have the recurrence coefficients from the earlier script.
# But rather than recompute them, use a different approach:
# Start from 3 different initial vectors in the CMF and track column 1.

# Solution 1: initial vector e₁ = (1,0,0)
e1_row = [mpf(1), mpf(0), mpf(0)]
# Solution 2: initial vector e₂ = (0,1,0)
e2_row = [mpf(0), mpf(1), mpf(0)]
# Solution 3: initial vector e₃ = (0,0,1)
e3_row = [mpf(0), mpf(0), mpf(1)]

u1, u2, u3 = [e1_row[0]], [e2_row[0]], [e3_row[0]]

for N in range(NMAX):
    M = M_entries(mpf(N))
    d = mpf(delta_H(N))
    MH = [[M[i][j]/d for j in range(3)] for i in range(3)]

    e1_new = [sum(e1_row[i]*MH[i][j] for i in range(3)) for j in range(3)]
    e2_new = [sum(e2_row[i]*MH[i][j] for i in range(3)) for j in range(3)]
    e3_new = [sum(e3_row[i]*MH[i][j] for i in range(3)) for j in range(3)]

    e1_row, e2_row, e3_row = e1_new, e2_new, e3_new
    u1.append(e1_row[0])
    u2.append(e2_row[0])
    u3.append(e3_row[0])

# q = (33750, -36000, 9000), so Q̂_N = 33750·u1_N - 36000·u2_N + 9000·u3_N
# p = (30921, -32972, 8240), so P̂_N = 30921·u1_N - 32972·u2_N + 8240·u3_N

# Verify:
for N in [0, 1, 10, 50]:
    qhat_check = 33750*u1[N] - 36000*u2[N] + 9000*u3[N]
    err = abs(qhat_check - qhat_vals[N])
    if err > mpf('1e-280'):
        print(f"  MISMATCH at N={N}: err = {mp.nstr(err, 5)}")
    else:
        print(f"  N={N}: decomposition verified")

# The dominant eigenvalue ξ₊ ≈ 33.97. For large N, u_i(N) ~ A_i · ξ₊^N · N^τ.
# The ratios u_i(N+1)/u_i(N) → ξ₊ for all i (unless i is in a subdominant subspace).
# But the ratios u2/u1 and u3/u1 converge to constants as N → ∞.

print(f"\n=== Asymptotic ratios of fundamental solutions ===")
for N in [50, 100, 150, 190]:
    r21 = u2[N] / u1[N]
    r31 = u3[N] / u1[N]
    print(f"  N={N}: u2/u1 = {mp.nstr(r21, 30)}, u3/u1 = {mp.nstr(r31, 30)}")

# These limiting ratios give the components of the dominant RIGHT eigenvector
# v₊ = (1, lim u2/u1, lim u3/u1)
# Then Q̂_N ~ (q·v₊) · w₊(0)·e₁ · ξ₊^N · ...
# And L = (p·v₊) / (q·v₊)

r21_inf = u2[NMAX] / u1[NMAX]
r31_inf = u3[NMAX] / u1[NMAX]
print(f"\nLimiting ratios (from N={NMAX}):")
print(f"  u2/u1 → {mp.nstr(r21_inf, 50)}")
print(f"  u3/u1 → {mp.nstr(r31_inf, 50)}")

# Dominant eigenvector v₊ ≈ (1, r21_inf, r31_inf)
# q·v₊ = 33750 + (-36000)*r21_inf + 9000*r31_inf
# p·v₊ = 30921 + (-32972)*r21_inf + 8240*r31_inf

q_dot_v = 33750 + (-36000)*r21_inf + 9000*r31_inf
p_dot_v = 30921 + (-32972)*r21_inf + 8240*r31_inf

L_computed = p_dot_v / q_dot_v
err_L = L_computed - G
digits_L = -mp.log10(abs(err_L)) if err_L != 0 else mp.dps

print(f"\nq·v₊ = {mp.nstr(q_dot_v, 50)}")
print(f"p·v₊ = {mp.nstr(p_dot_v, 50)}")
print(f"L = p·v₊ / q·v₊ = {mp.nstr(L_computed, 50)}")
print(f"G = {mp.nstr(G, 50)}")
print(f"L - G: {mp.nstr(err_L, 10)}")
print(f"L matches G to {mp.nstr(digits_L, 6)} digits")

# PSLQ on the individual projections
print(f"\n=== PSLQ analysis of q·v₊ and p·v₊ ===")

# Try to express q_dot_v and p_dot_v in terms of known constants
# Known: G, π, log(2), √2, π²/6 = ζ(2), etc.
constants = {
    'G': G,
    'pi': pi,
    'log2': log(2),
    'sqrt2': sqrt(2),
    'pi_sq': pi**2,
    'euler': euler,
}

for name, val in [("q·v₊", q_dot_v), ("p·v₊", p_dot_v)]:
    print(f"\n  {name} = {mp.nstr(val, 30)}")

    # Try: val = a + b*G for integer a,b
    result = pslq([val, G, 1])
    if result:
        print(f"    PSLQ {name} vs (G, 1): {result}")
        a, b, c = result
        if a != 0:
            print(f"    → {name} = ({-b}*G + {-c}) / {a}")

    # Try: val = a + b*G + c*pi
    result = pslq([val, G, pi, 1])
    if result:
        print(f"    PSLQ {name} vs (G, π, 1): {result}")

    # Try: val = a + b*G + c*√2
    result = pslq([val, G, sqrt(2), 1])
    if result:
        print(f"    PSLQ {name} vs (G, √2, 1): {result}")

# Also try: is r21_inf or r31_inf a known algebraic number?
print(f"\n=== PSLQ on eigenvector components ===")
for name, val in [("r21=u2/u1", r21_inf), ("r31=u3/u1", r31_inf)]:
    print(f"\n  {name} = {mp.nstr(val, 30)}")

    # Try minimal polynomial of degree 2
    result = pslq([val**2, val, 1])
    if result:
        print(f"    PSLQ degree-2: {result[0]}x² + {result[1]}x + {result[2]} = 0")

    # Try minimal polynomial of degree 3
    result = pslq([val**3, val**2, val, 1])
    if result:
        print(f"    PSLQ degree-3: {result}")

    # Try: val = a + b*√2 for rational a,b
    result = pslq([val, sqrt(2), 1])
    if result:
        print(f"    PSLQ vs (√2, 1): {result}")

# Also try: is q_dot_v / p_dot_v exactly 1/G?
print(f"\n  q·v₊ / p·v₊ = {mp.nstr(q_dot_v / p_dot_v, 30)}")
print(f"  1/G = {mp.nstr(1/G, 30)}")

# Key insight: if the eigenvector components are algebraic (involving √2),
# then q·v₊ and p·v₊ would be linear combinations of 1 and √2 times the
# initial row entries, and their ratio would need to equal G.
# This can only happen if G is algebraic — contradiction since G is transcendental.
# Therefore, the eigenvector components CANNOT be purely algebraic.
# They must involve G (or some other period) themselves.

# Try: r21 = a + b*G + c*G*√2 + d*√2
result = pslq([r21_inf, G, G*sqrt(2), sqrt(2), 1])
if result:
    print(f"\n  PSLQ r21 vs (G, G√2, √2, 1): {result}")

result = pslq([r31_inf, G, G*sqrt(2), sqrt(2), 1])
if result:
    print(f"  PSLQ r31 vs (G, G√2, √2, 1): {result}")

print("\n=== Done ===")
